include("../src/spglibwrapper.jl")

using Test

latt = Array{Float64, 2}([0.5 0.0 14.0; 0.0 0.5 14.0; 0.0 0.0 14.0])*5.4
positions = Array{Float64, 2}([0.0 0.0 0.0])
atoms = Vector{Int}([1])
cell = Cell(latt, positions, atoms)
@testset "spglibwrapper testset" begin
    Nirk, = get_ir_reciprocal_mesh([3;3;1], cell, [0;0;0])
    @test Nirk == 3
    Nirk, = get_ir_reciprocal_mesh([2;2;1], cell, [0;0;0])
    @test Nirk == 3
    Nirk, = get_ir_reciprocal_mesh([7;7;1], cell, [0;0;0])
    @test Nirk == 10
end
