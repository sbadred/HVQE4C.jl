using HMatrices
using StaticArrays
using SparseArrays
using ITensors


"""

Define a function to handle interaction computation based on the compression flag

"""
function compute_interactions(leaves,iter, h, v; use_compression::Bool=false, tol::Float64=1e-12, δ_compress::Float64=1e-12)
  if use_compression
      return compute_interactions_compression(leaves[iter], leaves[iter+1], h, v, atol=tol,δ_compress=δ_compress)
  else
      return compute_interactions_without_compression(leaves[iter], leaves[iter+1], h, v, atol=tol)
  end
end

"""
      partition_H(level,Op::SecondQuantizationOperators,n;tol::Float64=1e-12)

Build a compressed hierarchical Hamiltonian based a certain level k

"""
function partition_H(level, number_leaves, n, h, v; use_compression::Bool=false, tol::Float64=1e-12, δ_compress::Float64=1e-12)
  @assert n / (2^level) <= number_leaves
  H = ITensors.OpSum()
  tree = generate_binary_tree(number_leaves, n)
  leaves = Vector{SVector{1,Int64}}[]
  collect_leaves_at_level!(tree, 1, level, leaves)
  H_1 = compute_normal(leaves, h, v)


  leaves = Vector{SVector{1,Int64}}[]
  collect_leaves_at_level!(tree, 1, 2, leaves)
  # Initialize H based on the flag
  H =compute_interactions(leaves,1, h, v, use_compression=use_compression, tol=tol, δ_compress=δ_compress)

  for leaf = 2:(level-1)
    leaves = Vector{SVector{1,Int64}}[]
    collect_leaves_at_level!(tree, 1, leaf + 1, leaves)
    H += compute_interactions(leaves,1, h, v, use_compression=use_compression, tol=tol, δ_compress=δ_compress)
    for iter = 3:2:length(leaves)
      H +=compute_interactions(leaves,iter, h, v, use_compression=use_compression, tol=tol, δ_compress=δ_compress)
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
  #hamiltonian_1 += hamiltonian_2




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
    end

    for k_ in eachindex(intervals_right), l_ in eachindex(intervals_right)
      k, l = intervals_right[k_][1], intervals_right[l_][1]
      if norm(v[i, j, k, l]) > atol
        op = s == 1 ? "↑" : "↓"
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
        op = s == 1 ? "↑" : "↓"
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



"""
Compress Hamiltonian through the compression of integrals
"""
function compute_interactions_compression(intervals_left, intervals_right, h, v; atol::Float64=1e-14,  δ_compress::Float64=1e-12)
  v_=v*0.5  
  hamiltonian = ITensors.OpSum()
 
  hamiltonian=reshape_and_evaluate(h, v_, intervals_left, intervals_right,atol=atol,δ_compress=δ_compress);

  hamiltonian+=reshape_and_evaluate(h, v_, intervals_right, intervals_left,atol=atol,δ_compress=δ_compress);
  hamiltonian+=evaluate_B_compressed(intervals_left, intervals_right,v_; atol=atol)

  return hamiltonian
end



function evaluate_R_compressed(intervals_left, intervals_right,u_1,v_1, u_2, v_2; atol::Float64=1e-14)
  hamiltonian=ITensors.OpSum()
  low_rank=size(u_1,2)
  for α=1:low_rank
    for s in 1:2
        Op1= ITensors.OpSum(); Op1_T= ITensors.OpSum();Op2= ITensors.OpSum();Op2_T = ITensors.OpSum()
        op = s == 1 ? "↑" : "↓"
        add_ops_1e!(Op1_T, Op1, Op2_T, Op2, u_1, v_1, α, intervals_left, intervals_right, op)
        hamiltonian+= Ops.expand(Op1_T * Op2) + Ops.expand(Op2_T * Op1)
    end
  end

  low_rank=size(u_2,3)
  for α = 1:low_rank
      Op1, Op1_T = add_ops_2e!(u_2, α, intervals_left, intervals_right)
      Op2, Op2_T = add_ops_2e!(permutedims(v_2, [2, 3, 1]), α, intervals_right, intervals_right)
      hamiltonian +=  2*(Ops.expand(Op2*Op1) +  Ops.expand(Op1_T*Op2_T ))

  end
  return hamiltonian
end

function evaluate_B_compressed(intervals_left, intervals_right,v; atol::Float64=1e-14)
  hamiltonian=ITensors.OpSum()
  u_2, v_2 = compress_integrals(v, [intervals_left, intervals_left, intervals_right, intervals_right], δ_compress=δ_compress)

  L_left_start,L_left_end = intervals_left[1][1],intervals_left[end][1]
  L_right_start, L_right_end=intervals_right[1][1], intervals_right[end][1]
  
  u_2 = reshape(u_2, L_left_end-L_left_start+1, L_left_end-L_left_start+1, :); 
  v_2 = reshape(v_2, :, L_right_end-L_right_start+1, L_right_end-L_right_start+1)
 
  low_rank=size(u_2,3)
  for α = 1:low_rank
      Op1, ~ = add_ops_2e!(u_2, α, intervals_left, intervals_left)
      Op2, ~ = add_ops_2e!(permutedims(v_2, [2, 3, 1]), α, intervals_right, intervals_right)
      hamiltonian +=  2*Ops.expand(Op1*Op2)
  end

  u_2, v_2 = compress_integrals(v, [intervals_left, intervals_right, intervals_right, intervals_left], δ_compress=δ_compress)

  u_2 = reshape(u_2, L_left_end-L_left_start+1, L_right_end-L_right_start+1, :); 
  v_2 = reshape(v_2, :, L_right_end-L_right_start+1, L_left_end-L_left_start+1)
 
  low_rank=size(u_2,3)
  for α = 1:low_rank
      Op_1, ~ = add_ops_2e!(u_2, α, intervals_left, intervals_right)
      Op_2, ~ = add_ops_2e!(permutedims(v_2, [2, 3, 1]), α, intervals_right, intervals_left)
      hamiltonian +=  -Ops.expand(Op_1*Op_2)
  end

  u_2, v_2 = compress_integrals(v, [intervals_right,intervals_left, intervals_left,intervals_right], δ_compress=δ_compress)

  u_2 = reshape(u_2, L_right_end-L_right_start+1,L_left_end-L_left_start+1,  :); 
  v_2 = reshape(v_2, :, L_left_end-L_left_start+1,L_right_end-L_right_start+1)
 
  low_rank=size(u_2,3)
  for α = 1:low_rank
      Op1, ~ = add_ops_2e!(u_2, α, intervals_right, intervals_left)
      Op_2, ~ = add_ops_2e!(permutedims(v_2, [2, 3, 1]), α, intervals_right, intervals_left)
      hamiltonian +=  Ops.expand(Op1*Op_2)

      Op2, ~ = add_ops_2e!(permutedims(v_2, [2, 3, 1]), α, intervals_left, intervals_right)
      hamiltonian +=  -Ops.expand(Op1*Op2)
      
  end

  Operator=ITensors.OpSum()
  for id_left in eachindex(intervals_left)
    for id_right1 in eachindex(intervals_right)
      for id_right2 in eachindex(intervals_right)
        j = intervals_left[id_left][1]
        i = intervals_right[id_right1][1]
        l = intervals_right[id_right2][1]
        ITensors.add!(Operator, v[i, j, j, l], "c†↑", i, "c↑", l)
        ITensors.add!(Operator, v[i, j, j, l], "c†↓", i, "c↓", l)
      end
    end
  end

  for id_right in eachindex(intervals_right)
    for id_left1 in eachindex(intervals_left)
      for id_left2 in eachindex(intervals_left)
        j = intervals_left[id_right][1]
        i = intervals_right[id_left1][1]
        k = intervals_right[id_left2][1]
        ITensors.add!(Operator, v[i, j, k, j], "c†↑", i, "c↑", k)
        ITensors.add!(Operator, v[i, j, k, j], "c†↓", i, "c↓", k)
      end
    end
  end
  hamiltonian+=Operator
  return hamiltonian
end
"""
Compress integrals according to a given threshold
"""
function compress_integrals(ints,list; δ_compress::Float64=1e-12)
  dim=length(list)
  N=size(ints)[1]
  if dim==2
    M=ints[list[1][1][1]:list[1][end][1],list[2][1][1]:list[2][end][1]]
  else
    M=ints[list[1][1][1]:list[1][end][1],list[2][1][1]:list[2][end][1],list[3][1][1]:list[3][end][1],list[4][1][1]:list[4][end][1]]
    dim_tensor=size(M)
    M=reshape(M,prod(dim_tensor[1:2]),:)
  end
  u,s,v = LinearAlgebra.svd(M)
  #truncated_rank
  s_trunc = sv_trunc(s,δ_compress)
  tail = diagm(s_trunc) * v[:, 1:length(s_trunc)]'
  return u[:,1:length(s_trunc)], tail
end


function sv_trunc(s::Array{Float64}, tol; degenerate=true, degenerate_eps=1e-14)
  if tol == 0.0
      return s
  else
      d = length(s)
      i = 0
      weight = 0.0
      norm2 = dot(s, s)
      while (i < d) && weight < tol^2#*norm2
          weight += s[d-i]^2
          i += 1
      end

      if degenerate && (d - i + 1) != d
          k = i - 1
          while k > 0
              if abs(s[d-k+1] - s[d-i+1]) / abs(s[d-k+1]) > degenerate_eps
                  break
              else
                  i -= 1
                  k -= 1
              end
          end
      end
      return s[1:(d-i+1)]
  end
end


function add_ops_1e!(Si_T, Si, Sj_T, Sj, u, v, α, intervals_left, intervals_right, op)
  idx1=1
  idx2=1
  for id_left in eachindex(intervals_left)
      i = intervals_left[id_left][1]
      ITensors.add!(Si_T, 0.5 * u[idx1, α], "c†$op", i)
      ITensors.add!(Si, 0.5 * u[idx1, α], "c$op", i)
      idx1+=1
  end
  for id_right in eachindex(intervals_right)
      j = intervals_right[id_right][1]
      ITensors.add!(Sj_T, v[α, idx2], "c†$op", j)
      ITensors.add!(Sj, v[α, idx2], "c$op", j)
      idx2+=1
  end
end


function add_ops_2e!(v, α, intervals_left, intervals_right)
  Operator=ITensors.OpSum()
  Operator_T=ITensors.OpSum()
  
  for s in 1:2
    idx1=1
    idx2 = 1
    op = s == 1 ? "↑" : "↓"
    for id_left in eachindex(intervals_left)
      for id_right in eachindex(intervals_right)
        i = intervals_left[id_left][1]
        j = intervals_right[id_right][1]
        ITensors.add!(Operator, v[idx1, idx2, α], "c†$op", i, "c$op", j)
        ITensors.add!(Operator_T, v[idx1, idx2, α], "c†$op", j, "c$op", i)
        idx2 += 1
      end
      idx2 = 1
      idx1 += 1
    end
  end
  return Operator,Operator_T
end

function reshape_and_evaluate(h, v, intervals1, intervals2;atol::Float64=1e-14,δ_compress::Float64=1e-14)
  u_1, v_1 = compress_integrals(h, [intervals1, intervals2], δ_compress=δ_compress)
  u_2, v_2 = compress_integrals(v, [intervals1, intervals2, intervals2, intervals2], δ_compress=δ_compress)
  @show size(u_2)
  L_left_start,L_left_end = intervals1[1][1],intervals1[end][1]
  L_right_start, L_right_end=intervals2[1][1], intervals2[end][1]
  
  u_2 = reshape(u_2, L_left_end-L_left_start+1, L_right_end-L_right_start+1, :); 
  v_2 = reshape(v_2, :, L_right_end-L_right_start+1, L_right_end-L_right_start+1)
  hamiltonian=evaluate_R_compressed(intervals1, intervals2, u_1,v_1,u_2, v_2; atol=atol)
  return hamiltonian
end
