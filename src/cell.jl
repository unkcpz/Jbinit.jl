struct Cell
    latt::Array{Float64, 2}
    positions::Array{Float64, 2}
    atoms::Vector{Int}
    Natoms::Int64
end

function Cell(latt::Array{Float64, 2}, positions::Array{Float64, 2}, atoms::Vector{Int})
    Natoms = length(atoms)
    return Cell(latt, positions, atoms, Natoms)
end
