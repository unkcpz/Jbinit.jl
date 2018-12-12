struct Pspot
    filename::String
    symb::String
    zval::Int64
    rlocal::Float64
    rc::Vector{Float64}
    c::Vector{Float64}
    h::Array{Float64, 3}
    lmax::Int64
    Nproj_l::Vector{Int64}
    rcut_nlo::Vector{Int64}
end
