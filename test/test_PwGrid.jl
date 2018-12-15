include("../src/PwGrid.jl")

using Test
using LinearAlgebra

@testset "init_grid_R test" begin
    latt = Array{Float64, 2}([10 10 10; 5 0 5; 5 5 0])
    ns = Vector{Int64}([10, 8, 6])
    r = init_grid_real(latt, ns)

    @test size(r) == (10, 8, 6)
    @test r[1, 1, 1] == [0.0, 0.0, 0.0]
    @test r[2, 1, 1] == [1.0, 1.0, 1.0]
end

"""
"""
latt = [5.0 0 0; 0 5.0 0; 0 0 5.0]
recLatt = 2π * inv(transpose(latt))
lattLen = Vector{Float64}(undef, 3)
for i in 1:3
    lattLen[i] = norm(latt[i, :])
end
ecutrho = 120.0
Ns = map(x -> 2*round(Int64, √(ecutrho/2)x/π) + 1, lattLen)
gvec = init_GVectors(recLatt, Ns, ecutrho)
@testset "init_GVectors testset" begin
    @test size(gvec.G, 1) == size(gvec.G_length, 1) == gvec.Ng
end

latt = [5.0 0 0; 0 5.0 0; 0 0 5.0]
positions = Array{Float64, 2}([0.0 0.0 0.0])
symb = Vector{String}(["H"])
cell = Cell(latt, positions, symb)
kpoints = Kpoints(cell, [2,2,1], [0,0,0])
ecutwfc = 30.0
gvecwf = init_GVectorsWF(ecutwfc, gvec, kpoints)
@testset "init_GVectorsWF testset" begin
    @test length(gvecwf.Ngw) == gvecwf.kpoints.N
    @test gvecwf.N_max == 7809 # ??? useless ut
end
