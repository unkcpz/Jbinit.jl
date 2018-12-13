include("../src/PwGrid.jl")

using Test
using LinearAlgebra

@testset "init_grid_R test" begin
    latt = Array{Float64, 2}([10 10 10; 5 0 5; 5 5 0])
    ns = Vector{Int64}([10, 8, 6])
    r = init_grid_real(latt, ns)

    function has_row(v, m)
        return any(Bool[v == m[i, :]' for i=1:size(m, 1)])
    end

    @test size(r) == (480, 3)
    @test has_row([0 0 0], r)
    @test has_row([3 3 3], r)
end

@testset "init_Gvectors testset" begin
    latt = [5.0 0 0; 0 5.0 0; 0 0 5.0]

    recLatt = 2π * inv(transpose(latt))
    lattLen = Vector{Float64}(undef, 3)
    for i in 1:3
        lattLen[i] = norm(latt[i, :])
    end
    ecutrho = 120.0
    Ns = map(x -> 2*round(Int64, √(ecutrho/2)x/π) + 1, lattLen)

    gvec = init_Gvectors(recLatt, Ns, ecutrho)
    @test size(gvec.G, 1) == size(gvec.G_length, 1) == gvec.Ng
end
