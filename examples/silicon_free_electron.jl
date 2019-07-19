using PyCall
using JPDFT
mg = pyimport("pymatgen")
symmetry = pyimport("pymatgen.symmetry")
elec_structure = pyimport("pymatgen.electronic_structure")
plotter = pyimport("pymatgen.electronic_structure.plotter")

#
# Calculation parameters
#
kgrid = [3, 3, 3]
Ecut = 5  # Hartree
n_bands = 10
kline_density = 20


#
# Setup silicon structure in pymatgen
#
a = 5.431020504 * mg.units.ang_to_bohr
A = mg.ArrayWithUnit(a / 2 .* [[0 1 1.];
                               [1 0 1.];
                               [1 1 0.]], "bohr")
lattice = mg.lattice.Lattice(A)
recip_lattice = lattice.reciprocal_lattice
structure = mg.Structure(lattice, ["Si", "Si"], [ones(3)/8, -ones(3)/8])

# Get k-Point mesh for Brillouin-zone integration
spgana = symmetry.analyzer.SpacegroupAnalyzer(structure)
bzmesh = spgana.get_ir_reciprocal_mesh(kgrid)
# unkcpz: TODO kpoints as a julia type
kpoints = [mp[1] for mp in bzmesh]
kweigths = [mp[2] for mp in bzmesh]
kweigths = kweigths / sum(kweigths)

#
# Basis and Hamiltonian in DFTK
#
# Construct basis: transpose is required, since pymatgen uses rows for the
# lattice vectors and DFTK uses columns
grid_size = JPDFT.determine_grid_size(A', Ecut, kpoints=kpoints) * ones(Int, 3)
basis = PlaneWaveBasis(A', grid_size)

# Construct a free-electron Hamiltonian
ham = Hamiltonian(basis)

#
# Band structure calculation in DFTK
#
# Get the kpoints at which the band structure should be computed
symm_kpath = symmetry.bandstructure.HighSymmKpath(structure)
kpoints, klabels = symm_kpath.get_kpoints(kline_density, coords_are_cartesian=true)
println("Computing bands along kpath:\n     $(join(symm_kpath.kpath["path"][1], " -> "))")

# TODO Maybe think about some better mechanism here:
#      This kind of feels implicit, since it also replaces the kpoints
#      from potential other references to the ham or PlaneWaveBasis object.
kweigths = ones(length(kpoints)) ./ length(kpoints)

# Compute bands:
band_data = lobpcg(ham, n_bands, prec=PreconditionerKinetic(ham, α=0.5))
if ! band_data.converged
    println("WARNING: Not all k-points converged.")
end
