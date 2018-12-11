include("kpoints.jl")

using FFTW
"""
The type for set of G-vectors for describing density and potentials
"""
struct GVectors
    Ng::Int64
    G::Array{Float64, 2}
    G2::Array{Float64, 1}
    idx_g2r::Array{Int64, 1}
end # GVectors struct

"""
The type for set of G-vectors for describing wave function
"""
struct GVectorsWF
    Ngwx::Int64             # maximum(Ngk)
    Ngw::Array{Int64, 1}    # number of GvectorsWF for each K-points
    idx_gw2g::Array{Array{Int64, 1}, 1}
    idx_gw2r::Array{Array{Int64, 1}, 1}
    kpoints::Kpoints
end # GVectorsWF struct

"""
PWGrid describe plane wave basis set for a given periodic unit cell
"""
struct PWGrid
    ecutwfc::Float64    # energy cut of wavefunction
    ecutrho::Float64    #
    Ns::Tuple{Int64, Int64, Int64}
    LatVecs::Array{Float64, 2}
    RecVecs::Array{Float64, 2}
    CellVolume::Float64
    r::Array{Float64, 2}
    gvec::GVectors
    gvecwf::GVectorsWF
    planfw  # return by FFTW.plan_fft()
    planbw  # return by FFTW.plan_ifft()
end # PWGrid struct

# function PWGrid(ecutwfc::Float64, LatVecs::Array{Float64, 2}; kpoints=nothing)
#
#     ecutrho = 4.0*ecutwfc
#
#
#     return PWGrid(ecutwfc, ecutrho, Ns, LatVecs, RecVecs, CellVolume, r, gvec, gvecw,
#                   planfw, planbw)
#
# end # PwGrid module

"""
实空间中的点阵坐标
"""
function init_grid_R(latt, ns)
    R = Array{Float64, 2}(undef, prod(ns), 3)
    # b is the smallest grid unit
    b = latt .* (1 ./ hcat(ns, ns, ns))
    idx = 1
    for i in Iterators.product(0:ns[1]-1, 0:ns[2]-1, 0:ns[3]-1)
        scale = collect(i)
        R[idx, :] = scale' * b
        idx += 1
    end

    return R
end # init_grid_R function
