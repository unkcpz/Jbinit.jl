module Hamiltonian

mutable struct Hamiltonian
    pw::PWGrid
    potential::Potentials
    energies::Energies
end

end # Hamiltonian
