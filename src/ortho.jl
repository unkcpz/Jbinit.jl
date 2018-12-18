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

function ortho_sqrt(psi::Array{ComplexF64, 4})
    s = size(psi)
    Np = s[1]*s[2]*s[3]
    psi_tmp = reshape(psi, Np, s[4])
    udagger = inv(sqrt(psi_tmp'*psi_tmp))
    return reshape(psi_tmp*udagger, s[1], s[2], s[3], s[4])
end

function ortho_sqrt!(psi::Array{ComplexF64, 4})
    s = size(psi)
    Np = s[1]*s[2]*s[3]
    psi_tmp = reshape(psi, Np, s[4])
    udagger = inv(sqrt(psi_tmp'*psi_tmp))
    psi = reshape(psi_tmp*udagger, s[1], s[2], s[3], s[4])
end
