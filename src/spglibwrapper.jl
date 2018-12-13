using Libdl

macro checked_lib(libname, path)
    if Libdl.dlopen_e(path) == C_NULL
        error("Unable to load \n\n$libname ($path)\n\n")
    end
    quote
        const $(esc(libname)) = $path
    end
end


@checked_lib LIBSYMSPG "/usr/lib/libsymspg.so"

include("cell.jl")

function get_ir_reciprocal_mesh(cell::Cell,
                                mesh::Vector{Int64},
                                is_shift::Vector{Int64})
    is_time_reversal = 1
    symprec = 1e-5

    ## have to ??? why??
    cmesh = Base.cconvert(Vector{Cint}, mesh)
    # cmesh = mesh

    Nkpts = prod(mesh)
    grid = Array{Cint}(undef, 3, Nkpts)
    mapping = Vector{Cint}(undef, Nkpts)
    Nirk = ccall((:spg_get_ir_reciprocal_mesh, LIBSYMSPG),
        Cint,
        (Ref{Cint}, Ref{Cint}, Ptr{Cint}, Ptr{Cint},
        Cint, Ref{Cdouble}, Ref{Cdouble}, Ptr{Cint}, Cint, Cdouble),
        grid, mapping, cmesh, is_shift,
        is_time_reversal, cell.lattice, cell.positions, cell.atoms, cell.Natoms, symprec)

    Nirk = convert(Int64, Nirk)
    grid = transpose(grid)
    grid = convert(Array{Int64, 2}, grid)
    mapping = convert(Array{Int64, 1}, mapping)
    mapping = mapping .+ 1

    return Nirk, grid, mapping
end # spg_get_ir_reciprocal_mesh function
