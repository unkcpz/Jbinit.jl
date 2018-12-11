include("hamiltonian.jl")

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
    wk = Ham.pw.gvecw.kpoints.wk


end
