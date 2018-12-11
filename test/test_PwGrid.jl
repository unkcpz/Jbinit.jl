include("../src/PwGrid.jl")

using Test

@testset "init_grid_R test" begin
    latt = Array{Float64, 2}([10 10 10; 5 0 5; 5 5 0])
    ns = Vector{Int64}([10, 8, 6])
    r = init_grid_R(latt, ns)
    
    function has_row(v, m)
        return any(Bool[v == m[i, :]' for i=1:size(m, 1)])
    end

    @test size(r) == (480, 3)
    @test has_row([0 0 0], r)
    @test has_row([3 3 3], r)
end
