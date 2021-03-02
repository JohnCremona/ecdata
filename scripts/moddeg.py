def get_modular_degree(E, label):
    degphi_magma = 0
    degphi_sympow = 0
    #return E.modular_degree(algorithm='sympow')
    try:
        degphi_magma = E.modular_degree(algorithm='magma')
    except RuntimeError:
        print("{}: degphi via magma failed".format(label))
        try:
            degphi_sympow = E.modular_degree(algorithm='sympow')
        except RuntimeError:
            print("{}: degphi via sympow failed".format(label))
    if degphi_magma:
        if degphi_sympow:
            if degphi_magma == degphi_sympow:
                return degphi_magma
            else:
                print("{}: degphi = {} from magma but {} from sympow!".format(
                    label, degphi_magma, degphi_sympow))
                return degphi_magma
        else:
            return degphi_magma
    else:
        if degphi_sympow:
            return degphi_sympow
        else:
            print("{}: no success in computing degphi via magma or sympow".format(label))
            return 0
