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

    cpsi = zeros(ComplexF64, Npoints, Nstates)
    psiR = zeros(ComplexF64, Npoints, Nstates)
    rhoe = zeros(Float64, Npoints, Nspin)
