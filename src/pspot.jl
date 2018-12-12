struct Pspot
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

"""
Initialize Pspot from filename
"""
function Pspot(content::String)
    c = zeros(Float64, 4)
    rc = zeros(Float64, 4)
    h = zeros(Float64, 3, 3, 4)
    Nproj_l = zeros(Int64, 4)
    rcut_nlo = zeros(Float64, 4)

    lines = split(content, "\n")

    idx = 0
    while idx <= length(lines)
        if idx > 5
            @assert typeof(lmax) == Int64
            for i = 1:lmax
                s = split(lines[idx])
                rc[i] = parse(Float64, s[1])
                n = Nproj_l[i] = parse(Int64, s[2])
                for j in 1:n
                    s = split(lines[idx+j])
                    jk = parse.(Float64, s)
                    for k = 1:length(jk)
                        h[j,k,i] = jk[k]
                    end
                end
                idx += n + 1
            end
        elseif idx == 1
            global symb = split(lines[idx])[1]
        elseif idx == 2
            s = split(lines[idx])
            n_s = parse(Int64, s[1])
            n_p = parse(Int64, s[2])
            n_d = parse(Int64, s[3])
            n_f = parse(Int64, s[4])
            global zval = n_s + n_p + n_d + n_f
        elseif idx == 3
            s = split(lines[idx])
            global rlocal = parse(Float64, s[1])
            global n_local = parse(Int64, s[2])
        elseif idx == 4
            s = split(lines[idx])
            for i = 1:n_local
                c[i] = parse(Float64, s[i])
            end
        elseif idx == 5
            global lmax = parse(Int64, split(lines[idx])[1])
        end # end if
        idx += 1
    end

    return Pspot(symb, zval, rlocal, rc, c, h, lmax, Nproj_l, rcut_nlo)
end
