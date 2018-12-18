include("PwGrid.jl")
include("wfinit.jl")

function calc_rhoe(
    Nelec::Float64,
    pw::PwGrid, socc::Array{Float64, 2},
    bwf::BlochWF, Nspin::Int64;
    renormalize=true
)
    CellVolume = pw.CellVolume
    Ns = pw.Ns
    kpoints = pw.gvecwf.kpoints
    Nkpt = kpoints.N
    wk = kpoints.wk
    Npoints = prod(Ns)
    Nstates = bwf.Nstates
    psiks = bwf.psiks

    gvecwf = pw.gvecwf

    pfft = plan_ifft(zeros(ComplexF64, Ns[1], Ns[2], Ns[3], Nstates))
    rhoe = zeros(Float64, Ns[1], Ns[2], Ns[3], Nspin)
    psiG = zeros(ComplexF64, Ns[1], Ns[2], Ns[3], Nstates)
    psiR = zeros(ComplexF64, Ns[1], Ns[2], Ns[3], Nstates)

    for ispin = 1:Nspin
        for ik = 1:Nkpt
            psiG[:, :, :, :] .= 0.0 + im*0.0
            iks = ik + (ispin -1)*Nkpt
            psi = psiks[iks]
            Is = gvecwf.kgw_p[iks]
            for i in eachindex(Is)
                psiG[Is[i], :] = psi[i, :]
            end

            psiR = pfft * psiG
            ortho_sqrt!(psiR)
            psiR = √(Npoints/CellVolume) * psiR
            for in = 1:Nstates
                psi = psiR[:, :, :, in]
                rhoe += real(conj(psi).*psi)
            end
        end
    end
    @show rhoe
    return rhoe
end
