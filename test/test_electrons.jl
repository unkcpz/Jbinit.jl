include("../src/electrons.jl")
include("../src/cell.jl")
include("../src/pspot.jl")

using Test

o_gth_str =
"""
O GTH-PBE-q6
2  4  0  0
0.24455430  2
-16.66721480  2.48731132
    2
     0.22095592  1
     18.33745811
     0.21133247  0
"""

c_gth_str =
"""
C GTH-PBE-q4
2  2  0  0
0.33847124  2
-8.80367398  1.33921085
    2
     0.30257575  1
     9.62248665
     0.29150694  0
"""

@testset "electrons testset" begin
    latt = [5.0 0 0; 0 5.0 0; 0 0 5.0]
    pos = [0 0 0; 0.25 0.25 0.25; 0.5 0.5 0.5]
    symb = ["C", "O", "C"]
    c = Cell(latt, pos, symb)

    o_gth = Pspot(o_gth_str)
    c_gth = Pspot(c_gth_str)
    p_dict = Dict{String, Pspot}(["C" => c_gth, "O" => o_gth])
    @test get_Nelec(c, p_dict) == 14

end
