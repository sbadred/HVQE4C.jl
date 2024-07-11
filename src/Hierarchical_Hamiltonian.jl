using HMatrices
using StaticArrays
using SparseArrays
using ITensors
"""
      partition_H(level,Op::SecondQuantizationOperators,n;tol::Float64=1e-12)

Build a compressed hierarchical Hamiltonian based a certain level k

"""
function partition_H(level, number_leaves, n, h, v; tol::Float64=1e-12)
  @assert n / (2^level) <= number_leaves
  H = ITensors.OpSum()
  tree = generate_binary_tree(number_leaves, n)
  leaves = Vector{SVector{1,Int64}}[]
  collect_leaves_at_level!(tree, 1, level, leaves)
  H_1 = compute_normal(leaves, h, v)


  leaves = Vector{SVector{1,Int64}}[]
  collect_leaves_at_level!(tree, 1, 2, leaves)
  H = compute_interactions_without_compression(leaves[1], leaves[2], h, v, atol=tol)

  for leaf = 2:(level-1) #2^(level-1)
    leaves = Vector{SVector{1,Int64}}[]
    collect_leaves_at_level!(tree, 1, leaf + 1, leaves)
    H += compute_interactions_without_compression(leaves[1], leaves[2], h, v, atol=tol)
    for iter = 3:2:length(leaves)
      @show iter
      H += compute_interactions_without_compression(leaves[iter], leaves[iter+1], h, v, atol=tol)
    end
  end
  H += H_1
  H = electron_to_fermion(H)
  return H
end


"""
        generate_binary_tree(level,n)

Generate binray tree from 1D interval of size 'n' at depth  'level' a
"""
function generate_binary_tree(level,n)
    points = SVector.(range(1, n; step=1))
    splitter = HMatrices.GeometricSplitter(; nmax = level)
    splitter = HMatrices.CardinalitySplitter(level)
    tree = HMatrices.ClusterTree(points, splitter)
    Base.summary(tree)
    return tree
end

"""
         collect_leaves_at_level!(tree::ClusterTree{1, Int64}, current_level::Int, target_level::Int, leaves)

Generate a set of leaves at the current level of  a binary tree of depth n
"""
function collect_leaves_at_level!(tree::ClusterTree{1, Int64}, current_level::Int, target_level::Int, leaves)
    # Base case: If the node is nothing, return
    if tree === nothing
        return
    end

    # Check if the current node is a leaf and at the target level
    if current_level == target_level #&&  HMatrices.isleaf(tree)
        elements=tree._elements[tree.index_range]
        push!(leaves, elements)
    else
        # Recursively traverse the left and right subtrees
        collect_leaves_at_level!(tree.children[1], current_level + 1, target_level, leaves)
        collect_leaves_at_level!(tree.children[2], current_level + 1, target_level, leaves)
    end

end

"""
         compute_normal(tree::ClusterTree{1, Int64},level,h, v; atol::Float64=1e-14)

Compute the Hamiltonian in the spin basis without block interaction at specific levels from the tree 
"""
function compute_normal(leaves,h, v; atol::Float64=1e-14)
   
    hamiltonian = _molecular_orbital_hamiltonian(leaves,h, v,atol=atol)
    return hamiltonian
end




function _molecular_orbital_hamiltonian(intervals,h, v; atol::Float64=1e-14)
    # Representation of the second quantized quantum chemistry Hamiltonian.
    v *= 0.5; #permutedims(v , (3, 2, 1, 4))
    hamiltonian = ITensors.OpSum()
    for id_ in eachindex(intervals)
        for id1 in intervals[id_] 
          for id2 in intervals[id_]
            i=id1[1]; j=id2[1]
            @show i,j

            if norm(h[i, j]) > atol
              ITensors.add!(hamiltonian, h[i, j], "c†↑", i, "c↑", j)
              ITensors.add!(hamiltonian, h[i, j], "c†↓", i, "c↓", j)
            end
          end
        end
    
          for id1 in intervals[id_], id2 in intervals[id_], id3 in intervals[id_], id4 in intervals[id_]
            i=id1[1]; j=id2[1]; k=id3[1]; l=id4[1]
            if norm(v[i, j, k, l]) > atol
              #if (i ≠ j) && (k ≠ l) # Otherwise the terms are exactly zero
                ITensors.add!(hamiltonian, v[i, j, k, l], "c†↑", i, "c†↑", k, "c↑", l, "c↑", j)
                ITensors.add!(hamiltonian, v[i, j, k, l], "c†↓", i, "c†↓", k, "c↓", l, "c↓", j)
              #end
              ITensors.add!(hamiltonian, v[i, j, k, l], "c†↑", i, "c†↓", k, "c↓", l, "c↑", j)
              ITensors.add!(hamiltonian, v[i, j, k, l], "c†↓", i, "c†↑", k, "c↑", l, "c↓", j)
            end
          end
    end
    
    return hamiltonian
end

"""
    electron_to_fermion(hamiltonian::OpSum)
Map an OpSum from spinfull to spinless fermions.
"""
function electron_to_fermion(hamiltonian::ITensors.OpSum)
  fermion_hamiltonian = ITensors.OpSum()
  # loop over MPOTerms
  for k in 1:length(hamiltonian)
    h = hamiltonian[k]
    c = ITensors.coefficient(h)
    sites = first.(ITensors.sites.(ITensors.terms(h)))
    O = ITensors.name.(ITensors.terms(h))

    if O == ["Id"]
      fermion_hamiltonian += c, "Id", 1
    else
      ops_and_sites = []
      # loop over each single-site operator
      for (j, o) in enumerate(O)
        # the fermion is placed at twice the site number + 1 if ↓
        fermionsite = 2 * sites[j] + Int(o[end] == '↓') - 1
        ops_and_sites = vcat(ops_and_sites, (String(strip(o, o[end])), fermionsite)...)
      end
      fermion_hamiltonian += (c, ops_and_sites...)
    end
  end
  return fermion_hamiltonian
end


function compute_interactions_without_compression(intervals_left, intervals_right, h, v; atol::Float64=1e-14)
  v *= 0.5
  hamiltonian_1 = ITensors.OpSum()
  hamiltonian_2 = ITensors.OpSum()

  for id_left in eachindex(intervals_left)
    i = intervals_left[id_left][1]
    for s in 1:2
      H_21, H_21T = ITensors.OpSum(), ITensors.OpSum()
      op = s == 1 ? "↑" : "↓"
      ITensors.add!(H_21T, 2, "c†$op", i)
      ITensors.add!(H_21, 2, "c$op", i)
      R, RT = evaluate_R(i, intervals_right, h, v, s)
      hamiltonian_1 += Ops.expand(H_21T * R) + Ops.expand(RT * H_21)
    end

    for id_left2 in eachindex(intervals_left)
      j = intervals_left[id_left2][1]
      P = evaluate_P1(i, j, v, intervals_right)
      for s in 1:2
        op = s == 1 ? "↑" : "↓"
        hamiltonian_2 += (2, "c†$op", i, "c$op", j) * P
        for s_ in 1:2
          P2, PT, P_ = evaluate_P2(i, j, v, s, s_, intervals_right)
          op_ = s_ == 1 ? "↑" : "↓"
          hamiltonian_2 += (1, "c†$op", i, "c†$op_", j) * P2
          hamiltonian_2 += PT * (1, "c$op_", i, "c$op", j)
          hamiltonian_2 += (-2, "c†$op_", i, "c$op", j) * P_
        end
      end
    end
  end
  hamiltonian_1 += hamiltonian_2




  for id_right in eachindex(intervals_right)
    i = intervals_right[id_right][1]
    for s in 1:2
      H_21, H_21T = ITensors.OpSum(), ITensors.OpSum()
      op = s == 1 ? "↑" : "↓"
      ITensors.add!(H_21T, 2, "c†$op", i)
      ITensors.add!(H_21, 2, "c$op", i)
      R, RT = evaluate_R(i, intervals_left, h, v, s)
      hamiltonian_1 += Ops.expand(H_21T * R) + Ops.expand(RT * H_21)
    end
  end

  return hamiltonian_1
end



function evaluate_R(id_left, intervals_right, h, v, s; atol::Float64=1e-14)
  R = ITensors.OpSum()
  RT = ITensors.OpSum()
  i = id_left[1]

  for j_ in eachindex(intervals_right)
    j = intervals_right[j_][1]
    if norm(h[i, j]) > atol
      op = s == 1 ? "↑" : "↓"
      ITensors.add!(R, h[i, j] / 4, "c$op", j)
      @show i,j,s
    end

    for k_ in eachindex(intervals_right), l_ in eachindex(intervals_right)
      k, l = intervals_right[k_][1], intervals_right[l_][1]
      if norm(v[i, j, k, l]) > atol
        ITensors.add!(R, v[i, j, k, l], "c†↑", k, "c↑", l, "c$op", j)
        ITensors.add!(R, v[i, j, k, l], "c†↓", k, "c↓", l, "c$op", j)
      end
    end
  end

  for j_ in eachindex(intervals_right)
    j = intervals_right[j_][1]
    if norm(h[i, j]) > atol
      op = s == 1 ? "↑" : "↓"
      ITensors.add!(RT, h[i, j] / 4, "c†$op", j)
    end

    for k_ in eachindex(intervals_right), l_ in eachindex(intervals_right)
      k, l = intervals_right[k_][1], intervals_right[l_][1]
      if norm(v[i, j, k, l]) > atol
        ITensors.add!(RT, v[j, i, l, k], "c†$op", j, "c†↑", l, "c↑", k)
        ITensors.add!(RT, v[j, i, l, k], "c†$op", j, "c†↓", l, "c↓", k)
      end
    end
  end

  return R, RT
end




function evaluate_P1(i, j, v, intervals_right)
  P = ITensors.OpSum()

  for  k_ in eachindex(intervals_right), l_ in eachindex(intervals_right), s in 1:2
    k, l = intervals_right[k_][1], intervals_right[l_][1]
    op = s == 1 ? "↑" : "↓"
    ITensors.add!(P, v[i, j, k, l], "c†$op", k, "c$op", l)
  end

  return P
end




function evaluate_P2(i, j, v, s, s_, intervals_right)
  P = ITensors.OpSum()
  PT = ITensors.OpSum()
  P2 = ITensors.OpSum()

  for k_ in eachindex(intervals_right), l_ in eachindex(intervals_right)
    k, l = intervals_right[k_][1], intervals_right[l_][1]
    op = s == 1 ? "↑" : "↓"
    op_ = s_ == 1 ? "↑" : "↓"

    if (s == 1 && s_ == 2) || (s == 2 && s_ == 1)
      ITensors.add!(P, v[i, k, j, l], "c$op_", l, "c$op", k)
      ITensors.add!(PT, v[l, j, k, i], "c†$op", l, "c†$op_", k)
      ITensors.add!(P2, v[k, j, i, l], "c†$op", k, "c$op_", l)
    else
      ITensors.add!(P, v[i, k, j, l], "c$op", l, "c$op", k)
      ITensors.add!(PT, v[l, j, k, i], "c†$op", l, "c†$op", k)
      ITensors.add!(P2, v[k, j, i, l], "c†$op", k, "c$op", l)
    end
  end
  return P, PT, P2
end



