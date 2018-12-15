include("../src/wfinit.jl")

using Test

@testset "wfinit testset" begin
    bwf = rand_bwf([1,2,2,6], 2, 20)
    @test bwf.Nstates == 20
    @test length(bwf.psiks) == 4*2
    @test size(bwf.psiks[1]) == (1, 20)
    @test size(bwf.psiks[4]) == (6, 20)
end
