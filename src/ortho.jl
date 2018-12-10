module ortho

"""
ortho_sqrt function return a orthoganal psi φ,
φ><φ = I
"""
function ortho_sqrt(psi::Array{ComplexF64, 2})
    udagger = inv(sqrt(psi'*psi))
    return psi*udagger
end # ortho_sqrt function

function ortho_sqrt!(psi::Array{ComplexF64, 2})
    udagger = inv(sqrt(psi'*psi))
    psi[:,:] = psi*udagger
end # ortho_sqrt! function

end # ortho module
