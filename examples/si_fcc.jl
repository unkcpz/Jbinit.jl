function main()
    latt = [2.0 2.0 0.0;
            0.0 2.0 2.0;
            2.0 0.0 2.0]
    pos = [0.0 0.0 0.0]
    atoms = [14]
    c = Cell(latt, pos, atoms)

    pspfiles = ["../pps/pade_gth/Si-q4.gth"]
    ecutwfc_Ry = 30.0
    Ham = Hamiltonian(c, pspfiles, ecutwfc_Ry*0.5, meshk=[3, 3, 3])
    println(Ham)

    KS_solve_SCF!(Ham)
end # main function
