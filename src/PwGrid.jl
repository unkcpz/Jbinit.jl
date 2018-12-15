include("kpoints.jl")

using FFTW
"""
The type for set of G-vectors for describing density and potentials
"""
struct GVectors
    Ng::Int64
    G::Array{Float64, 2}
    G_length::Vector{Float64}
    # idx_g2r::Array{Int64, 1}
end # GVectors struct

"""
The type for set of G-vectors for describing wave function
"""
struct GVectorsWF
    N_max::Int64             # maximum(Ngk)
    Ngw::Vector{Int64}    # number of GvectorsWF for each K-points
    # idx_gw2g::Array{Array{Int64, 1}, 1}
    # idx_gw2r::Array{Array{Int64, 1}, 1}
    kpoints::Kpoints
end # GVectorsWF struct

"""
PWGrid describe plane wave basis set for a given periodic unit cell
"""
struct PWGrid
    ecutwfc::Float64    # energy cut of wavefunction
    ecutrho::Float64    #
    Ns::Tuple{Int64, Int64, Int64}
    lattice::Array{Float64, 2}
    recLatt::Array{Float64, 2}
    CellVolume::Float64
    r::Array{Float64, 2}
    gvec::GVectors
    gvecwf::GVectorsWF
    planfw  # return by FFTW.plan_fft()
    planbw  # return by FFTW.plan_ifft()
end # PWGrid struct

function PWGrid(ecutwfc::Float64, lattice::Array{Float64, 2}; kpoints=nothing)

    ecutrho = 4.0*ecutwfc
    recLatt = 2π * inv(transpose(lattice))

    lattLen = Vector{Float64}(undef, 3)
    for i in 1:3
        lattLen[i] = norm(lattice[i, :])
    end

    Ns = map(x -> 2*round(Int64, √(ecutrho/2)x/π) + 1, lattLen)

    r = init_grid_real(lattice, Ns)

    gvec = init_Gvectors(ns, recLatt, ecutrho)

    if kpoints == nothing
        kpoints = KPoints(1, (1,1,1), zeros(1,3), [1.0], recLattice)
    end

    gvecw = init_gvecw(ecutwfc, gvec, kpoints)

    planfw = plan_fft(zeros(Ns))
    planbw = plan_ifft(zeros(Ns))

    return PWGrid(ecutwfc, ecutrho, (Ns[1], Ns[2], Ns[3]),
                  lattice, recLatt, CellVolume, r,
                  gvec, gvecw,
                  planfw, planbw)

end # PwGrid module

"""
实空间中的点阵坐标,以一个3D矩阵表示,每个矩阵元表示该位置实际的分数坐标
"""
function init_grid_real(lattice, Ns)
    r = Array{Vector{Float64}, 3}(undef, Ns[1], Ns[2], Ns[3])
    # b is the smallest grid unit
    b = lattice ./ hcat(Ns, Ns, Ns)
    for I in CartesianIndices(r)
        scale = collect(Tuple(I)) .- 1
        r[I] = transpose(b) * scale
    end

    return r
end # init_grid_real function

"""
倒空间中满足|g|^2 < 2ecutrho的点, 及其到Γ的距离
"""
function init_GVectors(recLatt::Array{Float64, 2}, Ns::Vector{Int64}, ecutrho::Float64)
    iter_points = Iterators.product(
        ceil(Int64, -Ns[1]/2):ceil(Int64, Ns[1]/2)-1,
        ceil(Int64, -Ns[2]/2):ceil(Int64, Ns[2]/2)-1,
        ceil(Int64, -Ns[3]/2):ceil(Int64, Ns[3]/2)-1
    )

    G_tmp = []
    G_length = []
    for (idx, i) in enumerate(iter_points)
        scale = collect(i)
        g = transpose(scale) * recLatt
        gl = norm(g)
        if gl < √(2ecutrho)
            push!(G_tmp, g)
            push!(G_length, gl)
        end
    end

    Ng = length(G_tmp)
    G = Array{Float64, 2}(undef, Ng, 3)
    for i = 1:Ng
        G[i,:] = transpose(G_tmp[i])
    end

    G_length = convert(Vector{Float64}, G_length)

    return GVectors(Ng, G, G_length)
end

"""
对于每一个k点,给出该k点需要的平面波数量
"""
function init_GVectorsWF(ecutwfc::Float64, gvec::GVectors, kpoints::Kpoints)
    G = gvec.G
    Ng = gvec.Ng

    kpts = kpoints.k
    Nkpt = kpoints.N

    Ngw = Vector{Int64}(undef, Nkpt)
    Gk_l = Vector{Float64}(undef, Ng)
    for ik = 1:Nkpt
        for ig = 1:Ng
            Gk = G[ig, :] .+ kpts[ik, :]
            Gk_l[ig] = norm(Gk)
        end
        Ngw[ik] = length(findall(Gk_l .< ecutwfc))
    end

    N_max = maximum(Ngw)

    return GVectorsWF(N_max, Ngw, kpoints)
end
