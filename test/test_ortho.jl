include("../src/ortho.jl")

using Test
using LinearAlgebra
using Random

import .ortho:ortho_sqrt, ortho_sqrt!

atol = 1e-6
seed = reinterpret(Int32, Random.GLOBAL_RNG.seed)

p = ortho_sqrt(rand(ComplexF64, 2, 10))
@testset "ortho_sqrt test with seed = $seed" begin
    @test norm(p*p' - Matrix{ComplexF64}(I, 2, 2)) ≤ atol
    @test norm(p[1,:]'*p[1,:] - (1.0 + 0.0im)) ≤ atol
    @test norm(p[1,:]'*p[2,:]) ≤ atol
end

p = rand(ComplexF64, 2, 10)
ortho_sqrt!(p)
@testset "ortho_sqrt! test with seed = $seed" begin
    @test norm(p*p' - Matrix{ComplexF64}(I, 2, 2)) ≤ atol
    @test norm(p[1,:]'*p[1,:] - (1.0 + 0.0im)) ≤ atol
    @test norm(p[1,:]'*p[2,:]) ≤ atol
end
