module Hamiltonian

mutable struct Hamiltonian
    pw::PWGrid
    potential::Potentials
    energies::Energies
    rhoe::Array{Float64, 2} # spin dependent
    electrons::Electrons
    atoms::Atoms
    pspots::Array{PsPot_GTH, 1}
    pspotNL::PsPotNL
    xcfunc::String
    ik::Int64       # current kpoint index
    ispin::Int64    # current spin index
end

end # Hamiltonian
