using Revise
using HVQE4C


@testset "Compute_normal H" begin 
    #Generate one and two-electron integrals
    MOLECULE = "h2"
    δ = 1e-12
    number_leaves= 1
    level=2

    data = "data/fcidump_files/FCIDUMP." * MOLECULE
    Vnn, sites,N, h, v = read_electron_integral_tensors(data)
    n=sites
    @assert n/(2^level) <= number_leaves

    tree =generate_binary_tree(number_leaves,n)
    leaves=Vector{SVector{1, Int64}}[]
    collect_leaves_at_level!(tree,1,level,leaves)
    H_1=compute_normal(leaves,h, v,atol=δ)
end

@testset "Compute Interaction block" begin 
    tol = 1e-12
    n=2
    number_leaves= 1
    level=3
    h=rand(n,n)
    h=h+h'
    v=rand(n^2,n^2); 
    v=reshape(v+v',n,n,n,n)

    tree =generate_binary_tree(number_leaves,n)
    leaves=Vector{SVector{1, Int64}}[]
    collect_leaves_at_level!(tree,1,level,leaves)
    H_2=compute_interactions_without_compression(leaves,h, v, atol=tol)
end

@testset "Test the partition of H at level k" begin 
    MOLECULE = "h4"
    δ = 1e-12
    number_leaves= 1
    level=2

    data = "data/fcidump_files/FCIDUMP." * MOLECULE
    Vnn, sites,N, h, v = read_electron_integral_tensors(data)
    n=sites
    K=2*n
    @assert n/(2^level) <= number_leaves

  
    H=partition_H(level,number_leaves,n,h,v,use_compression=false,tol=δ)
    sitetype = "Fermion"
    sites= ITensors.siteinds(sitetype, K)
    HMPO=ITensors.MPO(H, sites,cutoff=0.0)

    H_partition=ITensors.contract(HMPO).tensor
    new_indices = vcat(2*K-1:-2:1, 2*K:-2:2)

    H_partition=reshape(permutedims(ITensors.array(H_partition),new_indices),(2^K, 2^K))
    
    #evaluate H_exact
    H_exact= molecular_hamiltonian_matrix(h, v)

    @assert norm(H_partition-H_exact) <= 1e-12

end


@testset "Test the partition of compression H at level k" begin 
    MOLECULE = "h4"
    δ = 1e-12
    number_leaves= 1
    level=2

    data = "data/fcidump_files/FCIDUMP." * MOLECULE
    Vnn, sites,N, h, v = read_electron_integral_tensors(data)
    n=sites
    K=2*n
    @assert n/(2^level) <= number_leaves

  
    H=partition_H(level,number_leaves,n,h,v,use_compression=true,tol=δ)
    sitetype = "Fermion"
    sites= ITensors.siteinds(sitetype, K)
    HMPO=ITensors.MPO(H, sites,cutoff=0.0)

    H_partition=ITensors.contract(HMPO).tensor
    new_indices = vcat(2*K-1:-2:1, 2*K:-2:2)

    H_partition=reshape(permutedims(ITensors.array(H_partition),new_indices),(2^K, 2^K))
    
    #evaluate H_exact
    H_exact= molecular_hamiltonian_matrix(h, v)

    @assert norm(H_partition-H_exact) <= 1e-12

end


@testset "test first block of interaction" begin 
    δ = 1e-12
    number_leaves = 1
    level = 2
    n=2
    h=rand(n,n)
    h=h+h'
    v=rand(n^2,n^2); 
    v=reshape(v+v',n,n,n,n)
    h[1,1]=0;h[2,2]=0;
    v[1,1,1,1]=0;v[2,2,2,2]=0

    K = 2 * n
    @assert n / (2^level) <= number_leaves
    tree = generate_binary_tree(number_leaves, n)
    leaves = Vector{SVector{1,Int64}}[]
    collect_leaves_at_level!(tree, 1, level, leaves)
    L = leaves[1] 
    R = leaves[2]

    operators = SecondQuantizationOperators(n)
    cup = operators.cup
    cdagup = operators.cdagup
    cdn = operators.cdn
    cdagdn = operators.cdagdn

    MPO_test = spzeros(4^n, 4^n)
    for (left, right) in Set([ (L, L), (R, R)])
        for id1 in left
            for id2 in right
                i = id1[1]
                j = id2[1]

                for s in [1, 2]
                    if s == 1
                        MPO_test += h[i, j] * cdagup[i] * cup[j]
                    else
                        MPO_test += h[i, j] * cdagdn[i] * cdn[j]
                    end
                end
            end
        end
    end
    v_ = 0.5 * v
    for (idx1, idx2, idx3, idx4) in Set([(L, L, L, L), (R, R, R, R)])
        for id1 in idx1, id2 in idx2, id3 in idx3, id4 in idx4
            i = id1[1]
            j = id2[1]
            k = id3[1]
            l = id4[1]
            for s in [1, 2], s_ in [1, 2]
                if s == 1 && s_ == 1
                    MPO_test += v_[i, j, k, l] * (cdagup[i] * cdagup[k] * cup[l] * cup[j])
                elseif s == 1 && s_ == 2
                    MPO_test += v_[i, j, k, l] * (cdagup[i] * cdagdn[k] * cdn[l] * cup[j])
                elseif s == 2 && s_ == 1
                    MPO_test += v_[i, j, k, l] * (cdagdn[i] * cdagup[k] * cup[l] * cdn[j])
                else
                    MPO_test += v_[i, j, k, l] * (cdagdn[i] * cdagdn[k] * cdn[l] * cdn[j])
                end
            end
        end
    end



    for (left, right) in Set([(L, R), (R, L)])#, (L, L), (R, R)])
        for id1 in left
            for id2 in right
                i = id1[1]
                j = id2[1]

                for s in [1, 2]
                    if s == 1
                        MPO_test += h[i, j] * cdagup[i] * cup[j]
                    else
                        MPO_test += h[i, j] * cdagdn[i] * cdn[j]
                    end
                end
            end
        end
    end


    v_ = 0.5 * v
    for (idx1, idx2, idx3, idx4) in Set([(L, L, L, R), (L, L, R, L), (L, R, L, L), (R, L, L, L),(R, R, R, L), (R, R, L, R), (R, L, R, R), (L, R, R, R)])#])#)#])#, (L, L, L, L), (R, R, R, R)])
        for id1 in idx1, id2 in idx2, id3 in idx3, id4 in idx4
            i = id1[1]
            j = id2[1]
            k = id3[1]
            l = id4[1]
            for s in [1, 2], s_ in [1, 2]
                if s == 1 && s_ == 1
                    MPO_test += v_[i, j, k, l] * (cdagup[i] * cdagup[k] * cup[l] * cup[j])
                elseif s == 1 && s_ == 2
                    MPO_test += v_[i, j, k, l] * (cdagup[i] * cdagdn[k] * cdn[l] * cup[j])
                elseif s == 2 && s_ == 1
                    MPO_test += v_[i, j, k, l] * (cdagdn[i] * cdagup[k] * cup[l] * cdn[j])
                else
                    MPO_test += v_[i, j, k, l] * (cdagdn[i] * cdagdn[k] * cdn[l] * cdn[j])
                end
            end
        end
    end

    v_ = 0.5 * v
    for (idx1, idx2, idx3, idx4) in Set([(L, R, L, R), (L, R, R, L), (L, L, R, R), (R, R, L, L), (R, R, L, L), (R, L, L, R), (R, L, R, L)])
        for id1 in idx1, id2 in idx2, id3 in idx3, id4 in idx4
            i = id1[1]
            j = id2[1]
            k = id3[1]
            l = id4[1]
            @show i,j,k,l
            for s in [1, 2], s_ in [1, 2]
                if s == 1 && s_ == 1
                    MPO_test += v_[i, j, k, l] * (cdagup[i] * cdagup[k] * cup[l] * cup[j])
                elseif s == 1 && s_ == 2
                    MPO_test += v_[i, j, k, l] * (cdagup[i] * cdagdn[k] * cdn[l] * cup[j])
                elseif s == 2 && s_ == 1
                    MPO_test += v_[i, j, k, l] * (cdagdn[i] * cdagup[k] * cup[l] * cdn[j])
                else
                    MPO_test += v_[i, j, k, l] * (cdagdn[i] * cdagdn[k] * cdn[l] * cdn[j])
                end
            end
        end
    end


    #Evaluate All
    H=partition_H(level,number_leaves,n,h,v,tol=δ)
    sitetype = "Fermion"
    sites = ITensors.siteinds(sitetype, K)
    HMPO = ITensors.MPO(H, sites, cutoff=0.0)

    H_partition = ITensors.contract(HMPO).tensor
    new_indices = vcat(2*K-1:-2:1, 2*K:-2:2)

    H_partition = reshape(permutedims(ITensors.array(H_partition), new_indices), (2^K, 2^K))
    @assert norm(H_partition-MPO_test) <= 1e-13


end


@testset "test integrals compression with only screening" begin 
    MOLECULE = "h4"
    δ = 1e-12
    data = "data/fcidump_files/FCIDUMP." * MOLECULE
    Vnn, sites,N, h, v = read_electron_integral_tensors(data)
    number_leaves= 1
    level=2
    n=sites
    @assert n / (2^level) <= number_leaves

    tree = generate_binary_tree(number_leaves, n)
    leaves = Vector{SVector{1,Int64}}[]
    collect_leaves_at_level!(tree, 1, level, leaves)
    L=leaves[1]
    R=leaves[2]
    
    tols=[1e-2, 1e-3,1e-4,1e-6,1e-8,1e-10,1e-12,1e-14]
    for tol in tols
       #1electron integrals
       list=[L,R]
       U_1,V_1=compress_integrals(h,list,δ_compress=tol)
       println("Frobenuis Norm 1: ", norm(h[1:list[1][end][1],1:list[2][end][1]]-U_1*V_1))

       list=[L,R,R,R]
       U_2,V_2=compress_integrals(v,[L,R,R,R],δ_compress=tol)
       M=v[1:list[1][end][1],1:list[2][end][1],1:list[3][end][1],1:list[4][end][1]]
       dim_tensor=size(M)
       M=reshape(M,prod(dim_tensor[1:2]),:)
       println("Frobenuis Norm 2: ", norm(M-U_2*V_2))

    end

end












@testset "Test the partition of compression H at level k with compression of intergrals" begin 
    MOLECULE = "h4"
    δ = 1e-12
   
    number_leaves= 1
    level=3

    data = "data/fcidump_files/FCIDUMP." * MOLECULE
    Vnn, sites,N, h, v = read_electron_integral_tensors(data)
    n=sites
    K=2*n
    @assert n/(2^level) <= number_leaves

  
    #without compression
    
    H=partition_H(level,number_leaves,n,h,v,use_compression=false,tol=δ)
    sitetype = "Fermion"
    sites= ITensors.siteinds(sitetype, K)
    HMPO=ITensors.MPO(H, sites,cutoff=0.0)

    H_partition=ITensors.contract(HMPO).tensor
    new_indices = vcat(2*K-1:-2:1, 2*K:-2:2)

    H_partition=reshape(permutedims(ITensors.array(H_partition),new_indices),(2^K, 2^K))
    
    #with compression
    list=[]
    tols = [1e-2, 1e-3, 1e-4, 1e-6, 1e-8, 1e-10, 1e-12, 1e-14]
    for δ_compress in tols
        @show δ_compress
        H_compressed = partition_H(level, number_leaves, n, h, v, use_compression=true, tol=δ, δ_compress=δ_compress)
        sitetype = "Fermion"
        sites = ITensors.siteinds(sitetype, K)
        HMPO = ITensors.MPO(H_compressed, sites, cutoff=0.0)

        H_compressed_partition = ITensors.contract(HMPO).tensor
        new_indices = vcat(2*K-1:-2:1, 2*K:-2:2)

        H_compressed_partition = reshape(permutedims(ITensors.array(H_compressed_partition), new_indices), (2^K, 2^K))
        push!(list,norm(H_partition - H_compressed_partition)./norm(H_partition))
        #@assert norm(H_partition - H_compressed_partition) <= 1e-12
    end

end