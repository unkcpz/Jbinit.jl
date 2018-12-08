module PwGrid

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
    kpoints::KPoints
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


end # PwGrid module
