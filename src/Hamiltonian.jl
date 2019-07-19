struct Hamiltonian
    basis         # Plane-wave basis object
    kinetic       # Kinetic energy operator
    pot_local     # Local potential term
    pot_nonlocal  # Non-local potential term
    pot_hartree   # Hartree potential term
    pot_xc        # Exchange-correlation term
end


"""
    Hamiltonian(basis::PlaneWaveBasis; kinetic=Kinetic(basis), pot_local=nothing,
                pot_nonlocal=nothing, pot_hartree=nothing, pot_xc=nothing)

Hamiltonian discretized in a plane-wave `basis`. Terms which are not specified
(left as `nothing`) are ignored during application. If only a basis object is specified,
a free-electron Hamiltonian is constructed. The `kinetic` `pot_local` and `pot_nonlocal`
terms shall not contain a non-linearity in the density.
"""
function Hamiltonian(basis; kinetic=Kinetic(basis), pot_local=nothing, pot_nonlocal=nothing,
            pot_hartree=nothing, pot_xc=nothing)
    Hamiltonian(basis, kinetic, pot_local, pot_nonlocal, pot_hartree, pot_xc)
end
