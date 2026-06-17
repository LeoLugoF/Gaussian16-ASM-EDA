import re

breverse = False
filename = ""

def NAtoms(lcoordinates):
    na = 1
    kn = 1
    for coors in lcoordinates:
        knp = int(coors.split()[0])
        if (knp - kn) < 0:
            na = kn
            break
        kn = knp
    return na

def OrderCartesianCoorAndEnergies(totalcoors,scfenergies,NAtoms,CriticalX):
    # At CriticalX starts the new forward/reverse calculation.
    FirstIRCpath = []
    Firstscfenergies = []
    SecondIRCpath = []
    Secondscfenergies = []
    for i in range(0,int(CriticalX)):
        FirstIRCpath.append(totalcoors[i])
        Firstscfenergies.append(scfenergies[i])
    for i in range(int(CriticalX),len(scfenergies)):
        SecondIRCpath.append(totalcoors[i])
        Secondscfenergies.append(scfenergies[i])
    SecondIRCpath.reverse()
    Secondscfenergies.reverse()

    TotalIRCpath =  SecondIRCpath + FirstIRCpath 
    Totalscfenergies = Secondscfenergies + Firstscfenergies
    
    global breverse
    global filename

    if breverse:
        TotalIRCpath.reverse()
        Totalscfenergies.reverse()
        savefilename = filename.split('.')[0]+"_reversed" + ".xyz"
    else:
        savefilename = filename.split('.')[0] + ".xyz"
 
    file = open(savefilename,"w")
    for i in range(0,len(scfenergies)):
        file.write(str(NAtoms)+'\n')
        file.write(Totalscfenergies[i]+'\n')
        file.write(TotalIRCpath[i])

    file.close()


def ReadCartesianCoorAndEnergies(data):
    energylinefromchk = re.compile(r'Energy From Chk'+r'.*',re.IGNORECASE).findall(data)
    scfenergies = []
    if (len(energylinefromchk) > 0):
        scfenergies = [energylinefromchk[0].split()[4]]

    TSCoorRCFC =     re.compile( r'^\s+\d+\s+[-+]?\d+\.\d+\s+[-+]?\d+\.\d+\s+[+-]?\d+\.\d+',re.MULTILINE).findall(data)
    if (len(TSCoorRCFC) == 0):
         TSCoorRCFC = re.compile( r'^\s+\w\,\d\,[+-]?\d+\.\d+\,[+-]?\d+\.\d+\,[+-]?\d+\.\d+',re.MULTILINE).findall(data)
         for i in range(0,len(TSCoorRCFC)):
             TSCoorRCFC[i] = TSCoorRCFC[i].replace(',',' ')
             TSCoorRCFC[i] = TSCoorRCFC[i].replace(' 0 ',' ')
    raw_totalcoors = re.compile( r'^\s+\d+\s+\d+\s+\d+\s+[-+]?\d+\.\d+\s+[-+]?\d+\.\d+\s+[+-]?\d+\.\d+',re.MULTILINE).findall(data)
    raw_scfenergies = re.findall(r'SCF Done:'+r'.*?'+ r'A.U. after',data)
    raw_steps = re.compile(r'Path Number:'+r'\s+\d+',re.MULTILINE).findall(data)
    #search at what step the forward or reverse IRC calculation finishes
    steps = [int(x.split()[-1]) for x in raw_steps]
    CriticalX = 0
    for i in range(0,len(steps)):
        if (steps[i] == 2):
            CriticalX = i 
            break

    
    #search number of atoms
    na = NAtoms(raw_totalcoors)
    
    for energies in raw_scfenergies:
        scfenergies.append(energies.split()[4])
    
    #Process Cartesian Coordinates into list
    xo = 0
    xf = na
    totalcoors = []
    while xo < len(raw_totalcoors):
        ln = raw_totalcoors[xo:xf]
        RCcoor = []
        if (xo == 0) and (len(energylinefromchk) > 0):
            for cooro in TSCoorRCFC:
                xyzline = cooro.split()
                xyzline = '     '.join(xyzline)
                RCcoor.append(xyzline)
            totalcoors.append('\n'.join(RCcoor))
            totalcoors[-1] = totalcoors[-1] + "\n"
            RCcoor= []
        for coor in ln:
            xyzline = coor.split()
            xyzline.pop(0)
            xyzline.pop(1)
            xyzline = '     '.join(xyzline)
            RCcoor.append(xyzline)
        totalcoors.append('\n'.join(RCcoor))
        totalcoors[-1] = totalcoors[-1] + "\n"
        xf+=na
        xo+=na

    OrderCartesianCoorAndEnergies(totalcoors,scfenergies,na,CriticalX)



def organizeIRC(File):
    global filename
    global breverse
    breverse = True
    filename = File
    file = open(File,"r")
    data = file.read()
    file.close()
    ReadCartesianCoorAndEnergies(data)

    breverse = False
    filename = File
    file = open(File,"r")
    data = file.read()
    file.close()
    ReadCartesianCoorAndEnergies(data)


