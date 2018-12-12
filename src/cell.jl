struct Cell
    lattice::Array{Float64, 2}
    positions::Array{Float64, 2}
    symbols::Vector{String}
    Natoms::Int64
    atoms::Vector{Int}
    compositions::Dict{String, Int64}
end

# function Cell(latt::Array{Float64, 2}, positions::Array{Float64, 2}, atoms::Vector{Int})
#     Natoms = length(atoms)
#     return Cell(latt, positions, atoms, Natoms)
# end
function Cell(lattice::Array{Float64, 2}, positions::Array{Float64, 2},
              symbols::Vector{String})
    #
    Natoms = length(symbols)

    atoms = Vector{Int64}(undef, Natoms)
    compositions = Dict{String, Int64}()
    for i in 1:Natoms
        s_key = symbols[i]
        atoms[i] = ZATOMS[s_key]
        if haskey(compositions, s_key)
            compositions[s_key] += 1
        else
            compositions[s_key] = 1
        end
    end
    @assert length(atoms) == Natoms

    return Cell(lattice, positions, symbols, Natoms, atoms, compositions)
end

ZATOMS = Dict(
"H"  => 1,
"He" => 2,
"Li" => 3,
"Be" => 4,
"B"  => 5,
"C"  => 6,
"N"  => 7,
"O"  => 8,
"F"  => 9,
"Ne" => 10,
"Na" => 11,
"Mg" => 12,
"Al" => 13,
"Si" => 14,
"P"  => 15,
"S"  => 16,
"Cl" => 17,
"Ar" => 18,
"K"  => 19,
"Ca" => 20,
"Sc" => 21,
"Ti" => 22,
"V"  => 23,
"Cr" => 24,
"Mn" => 25,
"Fe" => 26,
"Co" => 27,
"Ni" => 28,
"Cu" => 29,
"Zn" => 30,
"Ga" => 31,
"Ge" => 32,
"As" => 33,
"Se" => 34,
"Br" => 35,
"Kr" => 36,
"Rb" => 37,
"Sr" => 38,
"Y"  => 39,
"Zr" => 40,
"Nb" => 41,
"Mo" => 42,
"Tc" => 43,
"Ru" => 44,
"Rh" => 45,
"Pd" => 46,
"Ag" => 47,
"Cd" => 48,
"In" => 49,
"Sn" => 50,
"Sb" => 51,
"Te" => 52,
"I"  => 53,
"Xe" => 54,
"Cs" => 55,
"Ba" => 56,
"La" => 57,
"Ce" => 58,
"Pr" => 59,
"Nd" => 60,
"Pm" => 61,
"Sm" => 62,
"Eu" => 63,
"Gd" => 64,
"Tb" => 65,
"Dy" => 66,
"Ho" => 67,
"Er" => 68,
"Tm" => 69,
"Yb" => 70,
"Lu" => 71,
"Hf" => 72,
"Ta" => 73,
"W"  => 74,
"Re" => 75,
"Os" => 76,
"Ir" => 77,
"Pt" => 78,
"Au" => 79,
"Hg" => 80,
"Tl" => 81,
"Pb" => 82,
"Bi" => 83,
"Po" => 84,
"At" => 85,
"Rn" => 86,)

# ZATOMS = Dict(
# "H"  => 1,"He" => 2,
# "Li" => 3,"Be" => 4,"B" => 5,"C"  => 6,"N"  => 7,"O"  => 8,"F"  => 9,"Ne" => 10,
# "Na" => 11,"Mg" => 12,"Al" => 13,"Si" => 14,"P"  => 15,"S"  => 16,"Cl" => 17,"Ar" => 18,
# "K"  => 19,"Ca" => 20,"Sc" => 21,"Ti" => 22,"V"  => 23,"Cr" => 24,"Mn" => 25,"Fe" => 26,
# "Co" => 27,"Ni" => 28,"Cu" => 29,"Zn" => 30,"Ga" => 31,"Ge" => 32,"As" => 33,"Se" => 34,"Br" => 35,"Kr" => 36,
# "Rb" => 37,"Sr" => 38,"Y"  => 39,"Zr" => 40,"Nb" => 41,"Mo" => 42,"Tc" => 43,"Ru" => 44,
# "Rh" => 45,"Pd" => 46,"Ag" => 47,"Cd" => 48,"In" => 49,"Sn" => 50,"Sb" => 51,"Te" => 52,"I"  => 53,"Xe" => 54,
# "Cs" => 55,"Ba" => 56,"La" => 57,"Ce" => 58,"Pr" => 59,"Nd" => 60,"Pm" => 61,"Sm" => 62,
# "Eu" => 63,"Gd" => 64,"Tb" => 65,"Dy" => 66,"Ho" => 67,"Er" => 68,"Tm" => 69,"Yb" => 70,"Lu" => 71,"Hf" => 72,
# "Ta" => 73,"W"  => 74,"Re" => 75,"Os" => 76,"Ir" => 77,"Pt" => 78,"Au" => 79,"Hg" => 80,
# "Tl" => 81,"Pb" => 82,"Bi" => 83,"Po" => 84,"At" => 85,"Rn" => 86,
# )
