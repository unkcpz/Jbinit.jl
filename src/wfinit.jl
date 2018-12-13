const BlochWaveFunc = Array{Array{ComplexF64, 2}, 1}
"""
rand_bwf function generate a random bloch wave function
"""
function rand_bwf(pw::PWGrid, Nstates::Int64, Nspin::Int64)
    Nkpt = pw.gvecw.kpoints.Nkpt
    Ngw = pw.gvecw.Ngw

    # 所有的本征态的数量
    Nkspin = Nspin * Nkpt

    # 初始化每一个态
    psiks = BlochWaveFunc(undef, Nkspin)

    for ispin = 1:Nspin
        for ik = 1:Nkpt
            ikspin = ik + (ispin-1) * Nkpt
            psiks[ikspin] = rand_WF(Ngw[ik], Nstates)
        end
    end

    return psiks
end # rand_bwf function

function rand_bwf(pw::PWGrid, electrons::Electrons)
    Nspin = electrons.Nspin
    Nstates = electrons.Nstates
    return rand_bwf(pw, Nstates, Nspin)
end # rand_bwf function

function rand_WF(Nbasis, Nstates)
    return ortho_sqrt(rand(ComplexF64, Nbasis, Nstates))
end # rand_WF function
