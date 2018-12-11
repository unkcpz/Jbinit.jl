include("spglibwrapper.jl")

"""
Kpoints describe bloch wave vector of electronic states
"""
struct Kpoints
    N::Int64
    mesh::Tuple{Int64, Int64, Int64}
    k::Array{Float64, 2}    # 不可约的k
    wk::Vector{Float64}   # 各个k的权重
    RecVecs::Array{Float64, 2} # copy of reciprocal vectors
end

function Kpoints(cell::Cell, mesh::Array{Int64, 1}, is_shift::Array{Int64, 1})
    N, grid, mapping = get_ir_reciprocal_mesh(cell, mesh, is_shift)

    Ntotal = length(mapping)
    @assert Ntotal == prod(mesh)

    counter = zeros(Int64, Ntotal, 1)
    wka = zeros(Int64, Ntotal, 1)

    for i in 1:Ntotal
        xi = mapping[i]
        counter[xi] += 1
    end
    wka = counter .* inv(Ntotal)

    wk = Vector{Float64}(undef, N)
    kgrid = Array{Int64, 2}(undef, N, 3)
    idx = 1
    for i in 1:Ntotal
        if counter[i] != 0
            wk[idx] = wka[i]
            kgrid[idx, :] = grid[i, :]
            idx += 1
        end
    end
    k = kgrid ./ mesh

    @assert size(wk, 1) == N
    @assert size(k, 1) == N

    tuple_mesh = (mesh[1], mesh[2], mesh[3])
    RecVecs = 2π*inv(cell.latt')

    return Kpoints(N, tuple_mesh, k, wk, RecVecs)
end # Kpoints function
