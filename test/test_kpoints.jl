include("../src/kpoints.jl")

using Test
using LinearAlgebra

latt = Array{Float64, 2}([0.5 0.0 14.0; 0.0 0.5 14.0; 0.0 0.0 14.0])*5.4
positions = Array{Float64, 2}([0.0 0.0 0.0])
atoms = Vector{Int}([1])
cell = Cell(latt, positions, atoms)
kp = Kpoints(cell, [2,2,1], [0,0,0])
@testset "kpoints testset" begin
    @test kp.N == 3
    @test kp.mesh == (2, 2, 1)
    @test norm(kp.k - Array{Float64, 2}([0 0 0;0.5 0 0; 1 1 0])) < 1e-6
    @test norm(kp.wk - Vector{Float64}([0.25,0.5,0.25])) < 1e-6
    @test norm(kp.RecVecs - Array{Float64, 2}([2.32711 0.0 0.0;
                                               0.0 2.32711 0.0;
                                               -2.32711 -2.32711 0.08311])) < 1e-5
end
