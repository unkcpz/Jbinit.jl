include("../src/spglibwrapper.jl")

using Test

latt = Array{Float64, 2}([0.5 0.0 14.0; 0.0 0.5 14.0; 0.0 0.0 14.0])*5.4
positions = Array{Float64, 2}([0.0 0.0 0.0])
symb = Vector{String}(["H"])
cell = Cell(latt, positions, symb)
@testset "spglibwrapper testset" begin
    Nirk, = get_ir_reciprocal_mesh(cell, [3,3,1], [0,0,0])
    @test Nirk == 3
    Nirk, = get_ir_reciprocal_mesh(cell, [2,2,1], [0,0,0])
    @test Nirk == 3
    Nirk, grid, mapping = get_ir_reciprocal_mesh(cell, [7,7,1], [0,0,0])
    @test typeof(grid) == Array{Int64, 2}
    @test typeof(mapping) == Array{Int64, 1}
    @test Nirk == 10
end
