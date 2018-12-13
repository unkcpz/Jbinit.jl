include("../src/wfinit.jl")

using Test

@testset "wfinit testset" begin
    psiks = rand_bwf([1,2,2,6], 2, 20)
    @test length(psiks) == 4*2
    @test size(psiks[1]) == (1, 20)
    @test size(psiks[4]) == (6, 20)
end
