include("../src/cell.jl")

latt = Array{Float64, 2}([0.5 0.0 14.0; 0.0 0.5 14.0; 0.0 0.0 14.0])*5.4
positions = Array{Float64, 2}([0.0 0.0 0.0; 0.5 0.5 0.5])
atoms = Vector{Int}([14, 14])
cell = Cell(latt, positions, atoms)

using Test

@testset "cell testset" begin
    @test cell.Natoms == 2
end
