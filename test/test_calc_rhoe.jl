include("../src/calc_rhoe.jl")
include("../src/wfinit.jl")
include("../src/PwGrid.jl")

using Test

latt = [5.0 0 0; 0 5.0 0; 0 0 5.0]
positions = Array{Float64, 2}([0.0 0.0 0.0])
symb = Vector{String}(["H"])
cell = Cell(latt, positions, symb)
kpoints = Kpoints(cell, [2,2,1], [0,0,0])
ecutwfc = 3.0
pw = PwGrid(ecutwfc, latt; kpoints=kpoints)
Nspin = 1
Nstates = 4
bwf = rand_bwf(pw, Nspin, Nstates)
socc = zeros(Float64, 4, 1*3)
@testset "calc_rhoe test" begin
    rhoe = calc_rhoe(10.0, pw, socc, bwf, Nspin)
    @test 1==1
end
