"""
ortho_sqrt function return a orthoganal psi φ,
φ><φ = I
对给定的一组波函数进行幺正变换,使得输出的每一行是一个正交基
"""
function ortho_sqrt(psi::Array{ComplexF64, 2})
    udagger = inv(sqrt(psi'*psi))
    return psi*udagger
end # ortho_sqrt function

function ortho_sqrt!(psi::Array{ComplexF64, 2})
    udagger = inv(sqrt(psi'*psi))
    psi[:,:] = psi*udagger
end # ortho_sqrt! function
