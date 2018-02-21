import sys



def calculate_mua_mus(Aup,kup,Adown,kdown,g1,g2,StO2,C_Hb,C_Collagen,thickness):
    #   the absorption of hemoglobain used is from the research of S.A. Prahl
    #   "Tabulated molar extinction coefficient for hemoglobin in water"
    #   http://omlc.ogi.edu/spectra/hemoglobin/summary.html

    # USE DICTIONARY STRUCTURE TO EASILY FIND WAVELENGTH THAT IS WANTED
    mua_oxy = {}
    mua_deoxy = {}
    mua_up = {}
    mua_collagen = {}
    ###############################################################################
    # READ HB ABSORPTION FROM FILE

    with open('input/mua_Hb.txt', 'r') as mua:
        row = mua.read().split('\n')
        del row[len(row) - 1]
        for element in row:
            element = element.split('\t')
            # mua_oxy.append([element[0],element[1]])
            # mua_deoxy.append([element[0],element[2]])
            mua_oxy[element[0]] = element[1]
            mua_deoxy[element[0]] = element[2]

    # READ UPPER ABSORPTION FROM FILE
    with open('input/mua_up.txt', 'r') as mua:
        row = mua.read().split('\n')
        del row[len(row) - 1]
        for element in row:
            element = element.split('\t')
            mua_up[element[0]] = element[1]

    # READ UPPER ABSORPTION FROM FILE
    with open('input/mua_collagen.txt', 'r') as mua:
        row = mua.read().split('\n')
        del row[len(row) - 1]
        for element in row:
            element = element.split('\t')
            mua_collagen[element[0]] = element[1]

    # READ WAVELENGTH FROM FILE
    with open('input/wavelength.txt', 'r') as wave:
        wavelength = wave.read().split('\n')
        del wavelength[len(wavelength) - 1]
        wavelength = [365] + wavelength
        wavelength = [int(wavelength[i]) for i in range(len(wavelength))]








    #   load hemoglobin mua that is in specific wavelength
    oxy=[]
    deoxy=[]
    upper=[]
    collagen=[]
    for i in wavelength:
        oxy.append(float(mua_oxy[str(i)]))
        deoxy.append(float(mua_deoxy[str(i)]))
        upper.append(float(mua_up[str(i)]))
        collagen.append(float(mua_collagen[str(i)]))
    #   compute mus,mua in specific wavelength
    muaUP = [upper[i] for i in range(len(wavelength))]
    musUP = [Aup*pow(i,-kup)/(1-g1) for i in wavelength]

    muaDOWN = [2.303/64500*C_Hb*(StO2*oxy[i]\
                     +(1-StO2)*deoxy[i]) \
                     +C_Collagen*collagen[i]\
    for i in range(len(wavelength))]
    musDOWN = [Adown*pow(i,-kdown)/(1-g2) for i in wavelength]


    # unit: [1/cm]




    return wavelength,(muaUP,musUP,muaDOWN,musDOWN,thickness)



