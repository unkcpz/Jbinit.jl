include("../src/spglibwrapper.jl")

using Test

latt = Array{Float64, 2}([0 0.5 0.5; 0.5 0.0 0.5; 0.5 0.5 0])*5.4
positions = Array{Float64, 2}([0.875 0.875 0.875; 0.125 0.125 0.125])
atoms = Vector{Int}([1; 1])
cell = Cell(latt, positions, atoms)
@testset "spglibwrapper testset" begin
    Nirk, grid, mapping = get_ir_reciprocal_mesh([8;8;8], cell, [0;0;0])
    @test Nirk == 59
end
