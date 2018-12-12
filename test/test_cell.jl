include("../src/cell.jl")

using Test

@testset "cell testset" begin
    latt = Array{Float64, 2}([0.5 0.0 14.0; 0.0 0.5 14.0; 0.0 0.0 14.0])*5.4
    positions = Array{Float64, 2}([0.0 0.0 0.0; 0.1 0.1 0.1; 0.5 0.5 0.5])
    symbols = Vector{String}(["Si", "C", "Si"])
    cell = Cell(latt, positions, symbols)
    @test cell.Natoms == 3
    @test cell.atoms == [14, 6, 14]
    @test cell.compositions == Dict("Si" => 2, "C" => 1)
end
