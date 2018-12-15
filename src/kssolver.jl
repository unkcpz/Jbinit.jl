include("hamiltonian.jl")
include("wfinit.jl")

"""
Solves Kohn-Sham problem using trditional self-consistent field(SCF)
iterations with density mixing.
"""
function scf_solve!(Ham::Hamiltonian;
                initwfc=nothing, savewfc=false,
                initrhoe=:gaussian,
                is_band=true, is_energy=true, is_gap=true,
                check_rhoe=false,
                use_smearing=false, kT=1e-3,
                update_psi="LOBPCG", cheby_degree=8,
                mix_method="simple", MIXDIM=5,
                ETOT_CONV_THR=1-6)
    pw = Ham.pw

    Ngw = pw.gvecw.Ngw
    kpoints = pw.gvecw.kpoints
    Nkpt = kpoints.Nkpt
    wk = kpoints.wk

    # 实空间分割和实空间点密度
    Ns = pw.Ns
    Npoints = prod(Ns)
    CellVolume = pw.CellVolume
    dVol = CellVolume/Npoints

    electrons = Ham.electrons
    Nelec = elctrons.Nelec
    socc = elctrons.socc
    Nstates = electrons.Nstates
    Nstates_occ = electrons.Nstates_occ
    Nspin = electrons.Nspin

    Nkspin = Nkpt*Nspin

    # Random guess of wave function
    if startingwfc = nothing
        psiks = rand_bwf(pw, electrons)
    else
        psiks = startingwfc
    end

    # Calculated electron density from wf psik and update hamiltonian
    rhoe = zeros(Float64, Npoints, Nspin) # 2维矩阵是实空间中格点数量x自旋数量

    # TODO: Nspin == 2
    rhoe[:, :] = calc_rhoe(Nelec, pw, socc, psiks, Nspin)
    update!(Ham, rhoe)

end
