include("pspot.jl")
include("cell.jl")

mutable struct Electrons
    Nelec::Float64
    Nstates::Int64
    Nstates_occ::Int64
    socc::Array{Float64, 2}
    ebands::Array{Float64, 2}
    Nspin::Int64
end

"""
Creates an instance of `Electrons`
"""
function Electrons(cell::Cell, pspots_dict::Dict{String, Pspot};
                   Nspin=1, Nkpt=1, Nstates=nothing, Nstates_empty=0)
    #
    Nelec = get_Nelec(cell, pspots_dict)

    # Nstates manually from Nelec
    Nstates_occ = round(Int64, Nelec/2)
    if 2Nstates_occ < Nelec
        Nstates_occ += 1
    end
    Nstates = Nstates_occ + Nstates_empty

    socc = zeros(Float64, Nstates, Nkpt*Nspin)
    ebands = zeros(Float64, Nstates, Nkpt*Nspin)

    is_odd = round(Int64, Nelec) % 2 == 1
    if Nspin == 1
        for i = 1:Nstates_occ-1
            socc[i, :] .= 2.0
        end
        if is_odd
            socc[Nstates_occ, :] .= 1.0
        else
            socc[Nstates_occ, :] .= 2.0
        end
    elseif Nspin == 2
        for i = 1:Nstates_occ-1
            socc[i, :] .= 1.0
        end
        if is_odd
            # only spin up
            socc[Nstates_occ, 1:Nkpt] .= 1.0
        else
            socc[Nstates_occ, :] .= 1.0
        end
    end

    sum_socc = sum(socc)/Nkpt
    if abs(sum_socc - Nelec) > eps()
        error("sum(socc) = $sum_socc and Nelec = $Nelec is different\n")
    end

    return Electrons(Nelec, Nstates, Nstates_occ, socc, ebands, Nspin)
end

"""
return number of electrons for a given `cell::Cell` and `pspots::Array{pspot_GTH, 1}`
"""
function get_Nelec(cell::Cell, pspots_dict::Dict{String, Pspot})
    Nelec = 0.0
    com = cell.compositions
    for (k, v) in com
        Nelec += pspots_dict[k].zval * v
    end
    return Nelec
end
