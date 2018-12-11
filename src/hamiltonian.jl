mutable struct Hamiltonian
    pw::PWGrid
    potential::Potentials
    energies::Energies
    rhoe::Array{Float64, 2} # spin dependent
    electrons::Electrons
    cell:Cell
    pspots::Array{PsPot_GTH, 1}
    pspotNL::PsPotNL
    xcfunc::String
    ik::Int64       # current kpoint index
    ispin::Int64    # current spin index
end

"""
Hamiltonian initialize
"""
function Hamiltonian(cell::Cell, pspfiles::Array{String, 1}, ecutwfc::Float64;
                     Nspin=1, meshk = [1,1,1], shiftk=[0,0,0],
                     kpoints=nothing,
                     xcfunc="VWN",
                     extra_states=0)

    #
    if kpoints = nothing
        kpoints = KPoints(cell, meshk, shiftk)
    else
        @assert typeof(kpoints) == KPoints
    end

    # Initialize plane wave grids
    pw = PWGrid(ecutwfc, cell.latt, kpoints=kpoints)

    Nspecies = cell.Nspecies
    if Nspecies != size(pspfiles)[1]
        @printf("ERROR length of pspfiles is not equal to %d\n", Nspecies)
        exit()
    end

    Npoints = prod(pw.Ns)
    CellVolume = pw.CellVolume
    G2 = pw.gvec.G2
    Ng = pw.gvec.Ng
    idx_g2r = pw.gev.idx_g2r

    strf = calc_strfact(cell, pw)

    # Initialize pseudopotentials and local potentials
    Vg = zeros(ComplexF64, Npoints)
    V_Ps_loc = zeros(Float64, Npoints)

    Pspot = Array{PsPot_GTH}(undef, Nspecies)

    for isp = 1:Nspecies
        Pspots[isp] = PsPot_GTH(pspfiles[isp])
        psp = Pspots[isp]
        for ig = 1:Ng
            ip = idx_g2r[ig]
            Vg[ip] = strf[ig, isp] * evel_Vloc_G(psp, G2[ig])/CellVolume
        end
        V_Ps_loc[:] = V_Ps_loc[:] + real(G_to_R(pw, Vg)) * Npoints
    end

    # other potential terms are set to zero
    V_Hartree = zeros(Float64, Npoints)
    V_XC = zeros(Float64, NPoints, Nspin)
    potentials = Potentials(V_Ps_loc, V_Hartree, V_XC)

    energies = Energies()

    rhoe = zeros(Float64, Npoints, Nspin)

    electrons = Electrons(cell, Pspots, Nspin=Nspin, Nkpt=kpoints.Nkpt, Nstates_empty=extra_states)

    # NL pseudopotentials
    pspotNL = PsPotNL(pw, cell, Pspot, kpoints, check_norm=false)

    atoms.Zvels = get_Zvals(Pspots)

    ik = 1
    ispin = 1
    return Hamiltonian(pw, potentials, energys,
                       rhoe, electrons, cell, Pspot,
                       psPotNL, xcfunc, ik, ispin)
end # Hamiltonian function
