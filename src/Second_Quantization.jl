"""
    SecondQuantizationOperators
This struct contains the complete set of second quantization operator
matrices for an arbitrary number of orbitals (sites). Each operator is
stored as a sparse matrix of size (4^site_count, 4^site_count).
The operator matrices are stored in a set of arrays where each element
of an array corresponds to the operator acting on the orbital.
For example `cup` contains `site_count` annihilation operators, where
`cup[1]` annihilates an up electron on the first orbital.
All the operator matrices correctly constructed using the Jordan-Wigner
transformation. There are 6 operator arrays in the struct:
    id: identity operators
    f: Jordan-Wigner operators
    cup: up-electron annihilation operators
    cdagup: up-electron creation operators
    cdn: down-electron annihilation operators
    cdagdn: down-electron creation operators
"""

struct SecondQuantizationOperators
    id::Array
    f::Array
    cup::Array
    cdagup::Array
    cdn::Array
    cdagdn::Array
    Number_operator_up::Array
    Number_operator_down::Array
    function SecondQuantizationOperators(site_count::Int)

        id_one_site = sparse(I(4))
        f_one_site = sparse([1 0 0 0; 0 -1 0 0; 0 0 -1 0; 0 0 0 1])
        cdagup_one_site = sparse([0 0 0 0; 0 0 0 0; 1 0 0 0; 0 1 0 0])
        cup_one_site = cdagup_one_site'
        cdagdn_one_site = sparse([0 0 0 0; 1 0 0 0; 0 0 0 0; 0 0 -1 0])
        cdn_one_site = cdagdn_one_site'
        Number_one_site_up = sparse([0 0 0 0; 0 1 0 0; 0 0 0 0; 0 0 0 1])
        Number_one_site_down = sparse([0 0 0 0; 0 0 0 0; 0 0 1 0; 0 0 0 1])

        id = [id_one_site]
        f = [f_one_site]

        for i in 2:site_count-1
            push!(id, kron(id_one_site, id[i-1]))
            push!(f, kron(f_one_site, f[i-1]))
        end

        cup = [kron(cup_one_site, id[end])]
        cdagup = [kron(cdagup_one_site, id[end])]
        cdn = [kron(cdn_one_site, id[end])]
        cdagdn = [kron(cdagdn_one_site, id[end])]
        Number_operator_up = [kron(Number_one_site_up, id[end])]
        Number_operator_down = [kron(Number_one_site_down, id[end])]
        for i in 1:site_count-2
            push!(cup, kron(f[i], cup_one_site, id[end-i]))
            push!(cdagup, kron(f[i], cdagup_one_site, id[end-i]))
            push!(cdn, kron(f[i], cdn_one_site, id[end-i]))
            push!(cdagdn, kron(f[i], cdagdn_one_site, id[end-i]))
            push!(Number_operator_up, kron(id[i], Number_one_site_up, id[end-i]))
            push!(Number_operator_down, kron(id[i], Number_one_site_down, id[end-i]))
        end

        push!(cup, kron(f[end], cup_one_site))
        push!(cdagup, kron(f[end], cdagup_one_site))
        push!(cdn, kron(f[end], cdn_one_site))
        push!(cdagdn, kron(f[end], cdagdn_one_site))
        push!(Number_operator_up, kron(id[end], Number_one_site_up))
        push!(Number_operator_down, kron(id[end], Number_one_site_down))

        number = cdagup_one_site * cup_one_site + cdagdn_one_site * cdn_one_site
        @assert number == [0 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 2]

        @assert norm(cup[1] * cup[1] + cup[1] * cup[1]) == 0
        @assert norm(cup[2] * cup[2] + cup[2] * cup[2]) == 0
        @assert norm(cup[1] * cup[2] + cup[2] * cup[1]) == 0

        @assert cup[1] * cdagup[1] + cdagup[1] * cup[1] == I(4^site_count)
        @assert cup[2] * cdagup[2] + cdagup[2] * cup[2] == I(4^site_count)
        @assert norm(cup[1] * cdagup[2] + cdagup[2] * cup[1]) == 0
        @assert norm(cup[2] * cdagup[1] + cdagup[1] * cup[2]) == 0

        new(id, f, cup, cdagup, cdn, cdagdn, Number_operator_up, Number_operator_down)

    end

end


"""
    one_electron_hamiltonian_matrix(one_electron_tensor::Array{Float64, 2},
                                    operators=undef)
Construct the one-electron part of the molecular hamiltonian matrix. This
operator encompasses the electron kinetic energy and nuclear attraction
part of the molecular hamiltonian.
"""
function one_electron_hamiltonian_matrix(one_electron_tensor::Array{Float64,2},
    operators=undef)

    site_count = size(one_electron_tensor)[1]

    if operators == undef
        operators = SecondQuantizationOperators(site_count)
    end

    cup = operators.cup
    cdagup = operators.cdagup
    cdn = operators.cdn
    cdagdn = operators.cdagdn

    matrix = spzeros(4^site_count, 4^site_count)

    for c in CartesianIndices(size(one_electron_tensor))
        t = one_electron_tensor[c]
        matrix += t * (cdagup[c[1]] * cup[c[2]])
        matrix += t * (cdagdn[c[1]] * cdn[c[2]])
    end

    return matrix

end

"""
    two_electron_hamiltonian_matrix(two_electron_tensor::Array{Float64, 4},
                                    operators=undef)
Construct the two-electron part of the molecular hamiltonian matrix. This
operator encompasses the electron replulstion part of the molecular hamiltonian.
This function is optionally multithreaded. To exploit this, simply launch Julia
with several threads. For example,
    julia --threads 4
"""
function two_electron_hamiltonian_matrix(two_electron_tensor::Array{Float64,4},
    operators=undef)

    site_count = size(two_electron_tensor)[1]
    if operators == undef
        operators = SecondQuantizationOperators(site_count)
    end
    cup = operators.cup
    cdagup = operators.cdagup
    cdn = operators.cdn
    cdagdn = operators.cdagdn

    matrix = spzeros(4^site_count, 4^site_count)
    #mutex = SpinLock()
    for c in CartesianIndices(size(two_electron_tensor))
        thread_matrix = spzeros(4^site_count, 4^site_count)
        u = 1 / 2 * two_electron_tensor[c]
        if u == 0.0
            continue
        end

        #thread_matrix += u * (cdagup[c[1]] * cdagup[c[3]] * cup[c[4]] * cup[c[2]])
        thread_matrix += u * (cdagup[c[1]] * cdagup[c[3]] * cup[c[4]] * cup[c[2]])
        thread_matrix += u * (cdagup[c[1]] * cdagdn[c[3]] * cdn[c[4]] * cup[c[2]])
        thread_matrix += u * (cdagdn[c[1]] * cdagup[c[3]] * cup[c[4]] * cdn[c[2]])
        thread_matrix += u * (cdagdn[c[1]] * cdagdn[c[3]] * cdn[c[4]] * cdn[c[2]])
        #lock(mutex) do
        matrix += thread_matrix
        #end


    end
    return matrix
end


"""
    molecular_hamiltonian_matrix(one_electron_tensor::Array{Float64, 2},
                                 two_electron_tensor::Array{Float64, 4})
Construct the molecular hamiltonian matrix. This is simply a sum of the
matrices generated by one_electron_hamiltonian_matrix() and
two_electron_hamiltonian_matrix().
The two-electron part of this function is optionally multithreaded.
To exploit this, simply launch Julia with several threads. For example,
    julia --threads 4
"""
function molecular_hamiltonian_matrix(one_electron_tensor::Array{Float64,2},
    two_electron_tensor::Array{Float64,4})

    site_count = size(one_electron_tensor)[1]
    operators = SecondQuantizationOperators(site_count)
    one_electron_matrix = one_electron_hamiltonian_matrix(one_electron_tensor, operators)
    two_electron_matrix = two_electron_hamiltonian_matrix(two_electron_tensor, operators)

    return one_electron_matrix + two_electron_matrix

end

function molecular_hamiltonian_matrix(one_electron_tensor::Array{Float64,2})

    site_count = size(one_electron_tensor)[1]
    operators = SecondQuantizationOperators(site_count)
    one_electron_matrix = one_electron_hamiltonian_matrix(one_electron_tensor, operators)

    return one_electron_matrix

end

"""
    split_indices(all_indices_original)
Split array of indices into an array of arrays where each array contains
the set of indices that should be processed by a specific thread.
This information is used to distribute work amongst several threads.
"""
function split_indices(all_indices_original)

    @assert length(all_indices_original) >= Threads.nthreads()

    all_indices = copy(all_indices_original)

    thread_indices = []
    chunk_length = convert(Int64, ceil(length(all_indices) / Threads.nthreads()))

    for i in 1:Threads.nthreads()-1
        push!(thread_indices, splice!(all_indices, 1:chunk_length))
    end

    push!(thread_indices, all_indices)
    return thread_indices

end