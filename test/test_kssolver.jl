include("../src/kssolver.jl")

using Test

@testset "solver KS testset" begin
    latt = [2.0 2.0 0.0;
            0.0 2.0 2.0;
            2.0 0.0 2.0]
    pos = [0.0 0.0 0.0]
    symb = ["Si"]
    c = Cell(latt, pos, symb)

    pspfiles = ["../pps/pade_gth/Si-q4.gth"]
    ecutwfc_Ry = 30.0
    ham = Hamiltonian(c, pspfiles, ecutwfc_Ry*0.5, meshk=[3, 3, 3])
    scf_solve!(ham)
    @test 1 == 1
end
