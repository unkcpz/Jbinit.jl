include("../src/pspot.jl")

using Test
using LinearAlgebra

@testset "pspot testset" begin
    ce_psp_string =
    """
    Ce GTH-PBE-q30
        6   12   11    1
    0.20000000  2
    148.23700910  3.91725495
        4
         0.24060899  2
         2.21770923  4.21351645
                 -5.43962635
         0.16613857  3
         -4.53697472  11.78843522  0.01475880
         -13.92369350  -0.03492570
                 0.02482188
         0.17161109  1
         -15.00199037
         0.26526997  1
         -4.75946878
    """
    ce_psp = Pspot(ce_psp_string)
    @test ce_psp.symb == "Ce"
    @test ce_psp.zval == 30
    @test ce_psp.rlocal == 0.2
    @test norm(ce_psp.rc - [0.240609, 0.166136, 0.171611, 0.26527]) < 1e-5
    @test norm(ce_psp.c - [148.23700910, 3.91725495, 0, 0]) < 1e-5
    @test norm(ce_psp.h[1,2,2] - 11.78843522) < 1e-5
    @test ce_psp.lmax == 4
    @test ce_psp.Nproj_l == [2,3,1,1]
end
