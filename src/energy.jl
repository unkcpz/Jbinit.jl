"""
The type collecte all energy terms:
kinetic, local and nonlocal pseudopotential,
hartree, xc,
nuclei-nuclei terms
end electronic entropy term
"""
mutable struct Energies
    kinetic::Float64
    ps_loc::Float64
    ps_nloc::Float64
    hartree::Float64
    xc::Float64
    psp_core::Float64
    mTS::Float64
    cc::Float64 # energy of core-core interaction
end

"""
sum get total energy by summing all of its components
"""
function elec_energy(energies::Energies)
    return energies.kinetic + energies.ps_loc + energies.ps_nloc +
           energies.hartree + energies.xc + energies.psp_core + energies.mTS
end

function total_energy(energies::Energies)
    return elec_energy(energies) + energies.cc
end
