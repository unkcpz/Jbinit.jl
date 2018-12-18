include("PwGrid.jl")
include("electrons.jl")
include("ortho.jl")

mutable struct BlochWF
    Nstates::Int64
    psiks::Vector{Array{ComplexF64, 2}}
end

"""
返回一个Vector每一个元素是一个二维Nbasis*Nstates的矩阵,
每一行代表一个正交基
"""
function rand_bwf(Ngw::Vector{Int64}, Nspin::Int64, Nstates::Int64)
    # Ngw 的每个元素为表示该k点平面波的数量
    # 所有的本征态的数量
    Nkpt = length(Ngw)
    Nkspin = Nspin * Nkpt

    # 初始化每一个态
    psiks = Vector{Array{ComplexF64, 2}}(undef, Nkspin)

    for ispin = 1:Nspin
        for ik = 1:Nkpt
            ikspin = ik + (ispin-1) * Nkpt
            psiks[ikspin] = rand_WF(Ngw[ik], Nstates)
        end
    end

    return BlochWF(Nstates, psiks)
end
"""
rand_bwf function generate a random bloch wave function
"""
function rand_bwf(pw::PwGrid, Nspin::Int64, Nstates::Int64)
    Ngw = pw.gvecwf.kgw_n
    return rand_bwf(Ngw, Nspin, Nstates)
end # rand_bwf function

function rand_bwf(pw::PwGrid, electrons::Electrons)
    Nspin = electrons.Nspin
    Nstates = electrons.Nstates
    return rand_bwf(pw, Nspin, Nstates)
end # rand_bwf function

function rand_WF(Nbasis, Nstates)
    return ortho_sqrt(rand(ComplexF64, Nbasis, Nstates))
end # rand_WF function
