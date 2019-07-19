using FFTW

struct PlaneWaveBasis{T <: Real}
    # Lattice and reciprocal lattice vectors in columns
    lattice::Mat3{T}
    recip_lattice::Mat3{T}
    unit_cell_volume::T
    recip_cell_volume::T

    # Ecut::T

    grid_size::Vec3{Int}
    idx_DC::Int  # Index of the DC component in the rectangular grid
    # kpoints::Vector{Vec3{T}}
    # basis_wf::Vector{Vector{Vec3{Int}}}

    # kweights::Vector{T}

    FFT::FFTW.cFFTWPlan{Complex{T},-1,true,3}
    iFFT::AbstractFFTs.ScaledPlan{Complex{Float64},FFTW.cFFTWPlan{Complex{T},1,true,3},T}
end

function PlaneWaveBasis(lattice::AbstractMatrix{T}, grid_size) where {T <: Real}
    lattice = Mat3{T}(lattice)
    recip_lattice = 2π * inv(lattice')

    @assert(mod.(grid_size, 2) == ones(3),
            "grid_size needs to be a 3D Array with all entries odd, such that the " *
            "symmetry of a real quantity can be properly represented.")
    idx_DC = LinearIndices((Int.(grid_size)..., ))[ceil.(Int, grid_size ./ 2)...]

    # Plan a FFT, spending some time on finding an optimal algorithm
    # for the machine on which the computation runs
    fft_size = nextpow.(2, grid_size)  # Optimise FFT grid
    tmp = Array{Complex{T}}(undef, fft_size...)
    fft_plan = plan_fft!(tmp, flags=FFTW.MEASURE)
    ifft_plan = plan_ifft!(tmp, flags=FFTW.MEASURE)

    pw = PlaneWaveBasis{T}(lattice, recip_lattice, det(lattice), det(recip_lattice),
                           grid_size, idx_DC, fft_plan, ifft_plan)
end
"""
Reset the kpoints of an existing Plane-wave basis and change the basis accordingly.
"""
function set_kpoints!(pw::PlaneWaveBasis{T}, kpoints, kweights; Ecut=pw.Ecut) where T
    @assert(length(kpoints) == length(kweights),
            "Lengths of kpoints and length of kweights need to agree")
    @assert sum(kweights) ≈ 1 "kweights are assumed to be normalized."
    max_E = sum(abs2, pw.recip_lattice * floor.(Int, pw.grid_size ./ 2)) / 2
    @assert(Ecut ≤ max_E, "Ecut should be less than the maximal kinetic energy " *
            "the grid supports (== $max_E)")

    resize!(pw.kpoints, length(kpoints)) .= kpoints
    resize!(pw.kweights, length(kweights)) .= kweights
    resize!(pw.basis_wf, length(kpoints))

    # Update basis_wf: For each k-Point select those G coords,
    # satisfying the energy cutoff
    for (ik, k) in enumerate(kpoints)
        energy(q) = sum(abs2, pw.recip_lattice * q) / 2
        p = [G for G in basis_ρ(pw) if energy(k + G) ≤ Ecut]
        pw.basis_wf[ik] = [G for G in basis_ρ(pw) if energy(k + G) ≤ Ecut]
    end
    pw
end

function G_to_r!(pw::PlaneWaveBasis, f_fourier, f_real; gcoords=basis_ρ(pw))
    @assert length(f_fourier) == length(gcoords)
    @assert(size(f_real) == size(pw.iFFT),
            "Size mismatch between f_real(==$(size(f_real)) and " *
            "FFT size(==$(size(pw.iFFT))")
    fft_size = [size(pw.FFT)...]

    # Pad the input data, then perform an FFT on the rectangular cube
    f_real .= 0
    for (ig, G) in enumerate(gcoords)
        idx_fft = 1 .+ mod.(G, fft_size)
        f_real[idx_fft...] = f_fourier[ig]
    end
    mul!(f_real, pw.iFFT, f_real)

    # IFFT has a normalization factor of 1/length(ψ),
    # but the normalisation convention used in this code is
    # e_G(x) = e^iGx / sqrt(|Γ|), so we need to use the factor
    # below in order to match both conventions.
    f_real .*= length(pw.iFFT)
end

function r_to_G!(pw::PlaneWaveBasis, f_real, f_fourier; gcoords=basis_ρ(pw))
    @assert length(f_fourier) == length(gcoords)
    @assert(size(f_real) == size(pw.FFT),
            "Size mismatch between f_real(==$(size(f_real)) and " *
            "FFT size(==$(size(pw.FFT))")
    fft_size = [size(pw.FFT)...]

    # Do FFT on the full FFT plan, but truncate the resulting frequency
    # range to the part defined by the idx_to_fft array
    f_fourier_extended = pw.FFT * f_real  # This destroys data in f_real
    f_fourier .= 0
    for (ig, G) in enumerate(gcoords)
        idx_fft = 1 .+ mod.(G, fft_size)
        f_fourier[ig] = f_fourier_extended[idx_fft...]
    end
    # Again adjust normalisation as in G_to_r
    f_fourier .*= 1 / length(pw.FFT)
end
