include("../src/energy.jl")

using Test

@testset "energy testset" begin
    s = Energies(0, 0, 0, 0, 0, 0, 0, 1)
    @test elec_energy(s) == 0.0
    @test total_energy(s) == 1.0
end
