using FFTW

struct PlaneWaveBasis{T <: Real}
    # Lattice and reciprocal lattice vectors in columns
    lattice::Mat3{T}
    recip_lattice::Mat3{T}
    unit_cell_volume::T
    recip_cell_volume::T

    grid_size::Vec3{Int}
    idx_DC::Int  # Index of the DC component in the rectangular grid

    FFT::FFTW.cFFTWPlan{Complex{T},-1,true,3}
    iFFT::AbstractFFTs.ScaledPlan{Complex{Float64},FFTW.cFFTWPlan{Complex{T},1,true,3},T}
end

function PlaneWaveBasis(lattice::AbstractMatrix{T}, grid_size) where {T <: Real}
    lattice = SMatrix{3, 3, T, 9}(lattice)
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

# [15,15,15] -> [-7:7,-7:7,-7:7] generator
function basis_ρ(grid_size)
    start = -ceil.(Int, (grid_size .- 1) ./ 2)
    stop  = floor.(Int, (grid_size .- 1) ./ 2)
    (Vector{Int}([i, j, k]) for i=start[1]:stop[1], j=start[2]:stop[2], k=start[3]:stop[3])
end

function basis_ρ(pw::PlaneWaveBasis) = basis_ρ(pw.grid_size)

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

function compute_hartree(pwb::PlaneWaveBasis, ρ)
    # Solve the Poisson equation ΔV = -4π ρ in Fourier space,
    # i.e. Multiply elementwise by 4π / |G|^2.
    values_Y = [4π * ρ[ig] / sum(abs2, pwb.recip_lattice * G)
                for (ig, G) in enumerate(basis_ρ(pwb))]

    # Zero the DC component (i.e. assume a compensating charge background)
    values_Y[pwb.idx_DC] = 0

    return values_Y
end

function compute_hartree!(precomp, pwb::PlaneWaveBasis, ρ)
    # Solve the Poisson equation in Fourier space,
    values_Y = compute_hartree(pwb, ρ)

    # Fourier-transform and store in values_real
    T = eltype(pwb.lattice)
    values_real = similar(precomp, Complex{T})
    G_to_r!(pwb, values_Y, values_real)

    if maximum(imag(values_real)) > 100 * eps(T)
        throw(ArgumentError("Expected potential on the real-space grid B_ρ to be entirely" *
                            " real-valued, but the present density gives rise to a " *
                            "maximal imaginary entry of $(maximum(imag(values_real)))."))
    end
    recomp .= real(values_real)
end
