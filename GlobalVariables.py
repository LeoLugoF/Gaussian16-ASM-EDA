import os,re

def PrintConf() -> bool:
    if(os.path.exists('ASMG_conf.txt')):
        return True
    else:
        #Generates ASMG_conf.txt file
        f = open('ASMG_conf.txt','w')
        f.write("FILE USED TO GENERATE INPUTS\OUTPUTS\n\n")

        f.write("--Gaussian parameters--\n")
        f.write("%nprocshared=8\n")
        f.write("%mem=8gb\n")
        f.write("Method-Functional=pbe1pbe\n")
        f.write("Basis-set=6-311g(d)\n")
        f.write("Extra-command=\n")
        f.write("IOp(3/74)=\n")
        f.write("-----------------------\n\n")

        f.write("--File Paths---------\n")
        f.write("Path of XYZ file of IRC=/home/user/path.xyz\n")
        f.write("Path of Multiwfn=\n")
        f.write("---------------------\n\n")

        f.write("--X axis parameters--\n")
        f.write("Bond1=1,2\n")
        f.write("Bond2=\n")
        f.write("Angle1=1,2,3\n")
        f.write("Angle2=1,2,3\n")
        f.write("Pyramidalization1=1,2,3,3,2,4,4,2,1\n")
        f.write("Pyramidalization2=\n")
        f.write("Centroids1=1,2,3*1,2\n")
        f.write("Centroids2=\n")
        f.write("Dihedral1=1,2,3,4\n")
        f.write("Hirshfeld1=2\n")
        f.write("Hirshfeld2=3\n")
        f.write("--------------------------------------------\n\n")

        f.write("Distance between IRC points (if 1, all IRC points are considered)=5\n")
        f.write("Charge and multiplicity of molecular complex=0,1\n")
        f.write("More fragments can be added; only add the following lines: i.e. --Fragment3--\n")
        f.write("--Fragment1--\n")
        f.write("Charge and multiplicity of fragment 1=0,1\n")
        f.write("Atom indexes of fragment 1=1-6,10\n")
        f.write("SCF energy of optimized fragment 1=0.000000\n")
        f.write("--Fragment2--\n")
        f.write("Charge and multiplicity of fragment 2=0,1\n")
        f.write("Atom indexes of fragment 2=12-15,18\n")
        f.write("SCF energy of optimized fragment 2=0.000000\n")
        f.close()
    return False

class ASMG_parameters:
    #Parameters to perform a quantum calculation
    g09nprocshared=""
    g09mem=""
    g09functional=""
    g09basisset=""
    gaussianversion="g16"

    #Files needed
    xyzfilepath=""
    multiwfnpath=""

    #Charge and multiplicity of Molecular Complex, Fragment 1, ..., Fragment N : { 0:[0,1],1:[0,1],... }
    dCM={}

    #Space between the newly selected IRC points for ASA calculations
    IRCgap=1

    #Index of atoms of M.C, F1, ..., FN.
    dAIndex={}

    #X axis parameters: dBonds = {1:[1,2],...}, 
    #                   dPyramidalization = {1:[1,2,3,3,2,4,4,2,1]}, dangle = {1:[1,2,3]}
    #                   dcentroids = {0 : { 0 : [1,2,3] , 1: [2,3] }, .. }
    dbonds={}
    dpyramidalization={}
    dangle={}
    ddihedral={}
    dhirshfeld={}
    # dcentroid[n1][n2]
    # {0: ['11', '5'], 1: ['13', '3']}
    # n1 = atoms that conform the centroid-centroid distance
    # n2 = 0 : centroid_sub1, 1 : centroid_sub2 (given in array)
    dcentroids={}
    # dprintxaxis[distance_n][n2]
    # { distance_n : [b1,b2,...,bn], ...  
    # distance_n = Proyected coordinate (distance, angle, ...). Starts at 0. (angle_n, pyramidalization_n, hirshfeld_n)
    # n2 = X point position
    dprintxaxis={}

    #Functional Parameters for SCF decomposition
    dfunctionalps={'tpssh':"-35",'pbe1pbe':"-13",'b3lyp':"-5",'hf':'100'}
    #functionalparameter
    functionalpar=""

    #Calculated Hirshfeld charges for each molecular complex
    # dcalchirshfeld[n1][n2][n3]
    # { RC_1 : [[hirshfeldcharge1,spinpopulation1,...],...], ..., nRC }
    # n1 = RC number (only the calculated one)
    # n2 = Atom Index
    # n3 = 0 : Hirshfeld cahrge, 1 : spin population, ...,
    dcalchirshfeld={}
    typeofcharge = 0
    #Coordinates of molecular complex : 
    # [  [ [Iatom1 x11 y11 z11],[Iatom2 x21 y21 z21],... ]  , ...]
    # Molecular complex 1                        Molecular complex 2
    # lcoors[n1][n2]
    # CONTAINS ALL THE COORDINATES OF THE IRC.
    # n1 = IRC point (starts at 0)
    # n2 = Atom Index : Gives in one line the atomic number and cartesian coordinates 
    lMCcoors=[]

    #Number of atoms
    natoms=0
    # SCF energy of fragment 1, 2, ..., n
    # fragment 1 corresponds to index 0
    # { 0 :  -1000.00, ..., n: -2000.000 }
    dSCFe={}

    #IRC point of TS
    TSXpoint=0

    #Command needed for each fragment
    extracommand=""
    fragmentcommand="em=gd3bj ExtraLinks=L608 nosymm"
    promolcommand="em=gd3bj ExtraLinks=L608 nosymm scf=maxcyc=-1 guess=read iop(4/6=222)"
    mcomplexcommand="em=gd3bj ExtraLinks=L608 nosymm guess=read pop=hirshfeld"

    #Decomposed energies
    # { 0 : { 0 : [Eels,Ex,Ec,Edisp,Etotal] } }
    # dDecomposedE[n1][n2][n3]
    # n1 = 0 : Promol, 1 : fragment1, 2 : fragment 2.
    # n2 = number of file; according to its point in the IRC.
    # n3 = Type of energy
    dDecomposedE = {}
    #Energy of MComplex in its final stage and frozen state
    # { 0 : [Escffinalstate,Escffrozenstate], }
    # dEscf[n1][n2]
    # n1 = number of file; according to its point in the IRC.
    # n2 = 1) Escf of final state and 2) of frozen state
    dEscf = {}

    # EDA enegies
    #{ 0 : [Eels,Ex,Ec,Edisp,Eorb,Epauli,Eint,Etotal] }
    # dEDA [n1][n2]
    # n1 = number of file; according to its point in the IRC.
    # n2 = Type of energy
    dEDA = {}

    # Fragment strain
    # { 0 : { 0 : 10.1 }, ...}
    # dDecomposedE[n1][n2]
    # n1 = 0 : total strain, 1 : fragment1, 2 : fragment 2 ... 
    # n2 = number of file; according to its point in the IRC.
    dFragmentStrain = {}

    lCorrectseqoffragments = []

    def __init__(self):
        #Reads ASMG_Parameters.txt file

        file = open('ASMG_conf.txt','r')
        output = file.read()
        file.close()

        self.g09nprocshared=re.compile(r'^\%nprocshared=\d+',re.MULTILINE).findall(output)[0].split('=')[-1]
        self.g09mem=re.compile(r'^\%mem=\d+.*',re.MULTILINE).findall(output)[0].split('=')[-1]
        self.g09functional=re.compile(r'^\b(Method-Functional=.*)',re.MULTILINE).findall(output)[0].split('=')[-1]
        self.g09basisset=re.compile(r'^\b(Basis-set=.*)',re.MULTILINE).findall(output)[0].split('=')[-1]
        self.extracommand=re.compile(r'^\b(Extra-command=.*)',re.MULTILINE).findall(output)[0].split('=')[-1]

        lbonds = re.compile(r'^\b(Bond.*)',re.MULTILINE).findall(output)
        i = 0
        for bond in lbonds:
            rawbond = bond.split('=')[-1]
            _rl = rawbond.split(',')
            if len(_rl) == 0:
                break
            if _rl[0] == '':
                break
            rl = [int(number) for number in _rl]
            self.dbonds[i]=rl
            
            i=i+1

        lang = re.compile(r'^\b(Angle.*)',re.MULTILINE).findall(output)
        i = 0
        for ang in lang:
            rawang = ang.split('=')[-1]
            _rl = rawang.split(',')
            if len(_rl) == 0:
                break
            if _rl[0] == '':
                break
            rl = [int(number) for number in _rl]
            self.dangle[i]=rl
            i=i+1

        lang = re.compile(r'^\b(Dihedral.*)',re.MULTILINE).findall(output)
        i = 0
        for ang in lang:
            rawang = ang.split('=')[-1]
            _rl = rawang.split(',')
            if len(_rl) == 0:
                break
            if _rl[0] == '':
                break
            rl = [int(number) for number in _rl]
            self.ddihedral[i]=rl
            i=i+1

        lpyr = re.compile(r'^\b(Pyramidalization.*)',re.MULTILINE).findall(output)
        i = 0
        for pyr in lpyr:
            rawpyr = pyr.split('=')[-1]
            _rl = rawpyr.split(',')
            if len(_rl) == 0:
                break
            if _rl[0] == '':
                break
            rl = [int(number) for number in _rl]
            self.dpyramidalization[i]=rl
            i=i+1

        _lcentroid = re.compile(r'\b(Centroids.*)',re.MULTILINE).findall(output)
        i = 0
        for _centroids in _lcentroid:
            rawcentroid = _centroids.split('=')[-1]
            rc = rawcentroid.split('*')
            if len(rc) == 2:
                self.dcentroids[i] = {}
                _rc2 = rc[0].split(',')
                if len(_rc2) == 0:
                    break
                if _rc2[0] == '':
                    break
                rc2 = [int(number) for number in _rc2]
                self.dcentroids[i][0] = rc2
                _rc3 = rc[1].split(',')
                if len(_rc3) == 0:
                    break
                if _rc3[0] == '':
                    break
                rc3 = [int(number) for number in _rc3]
                self.dcentroids[i][1] = rc3
            i=i+1

        lhirsh = re.compile(r'^\b(Hirshfeld.*)',re.MULTILINE).findall(output)
        i = 0
        for hirsh in lhirsh:
            FAtoms = []
            RawList=hirsh.split("=")[-1].split(",")
            for a in RawList:
                b = a.strip().split('-')[0]
                c = a.strip().split('-')[-1]
                if b != c:
                    FAtoms = FAtoms + list(range(int(b),int(c)+1))
                else:
                    FAtoms.append(int(b))

            if len(FAtoms) == 0:
                break
            self.dhirshfeld[i]=FAtoms
            i=i+1

        self.xyzfilepath=re.compile(r'\b(XYZ file.*)',re.MULTILINE).findall(output)[0].split('=')[-1]
        self.multiwfnpath=re.compile(r'\b(Multiwfn.*)',re.MULTILINE).findall(output)[0].split('=')[-1]
        self.IRCgap=int(re.compile(r'\b(IRC points.*)',re.MULTILINE).findall(output)[0].split('=')[-1])
        
        _lCM = re.compile(r'\b(Charge and mul.*)',re.MULTILINE).findall(output)
        i = 0
        for _cm in _lCM:
            _cm : str
            cm = _cm.split('=')[-1].split(',')
            self.dCM[i] = cm
            i=i+1

        _SCFe = re.compile(r'\b(SCF energy of.*)',re.MULTILINE).findall(output)
        i = 0
        for _scf in _SCFe:
            _scf : str
            scf = _scf.split('=')[-1]
            self.dSCFe[i] = float(scf)
            i=i+1

        _lai = re.compile(r'\b(Atom indexes of.*)',re.MULTILINE).findall(output)
        i = 0
        for _ai in _lai:
            FAtoms = []
            RawList=_ai.split("=")[-1].split(",")
            for a in RawList:
                b = a.strip().split('-')[0]
                c = a.strip().split('-')[-1]
                if b != c:
                    FAtoms = FAtoms + list(range(int(b),int(c)+1))
                else:
                    FAtoms.append(int(b))

            if len(FAtoms) == 0:
                break
            self.dAIndex[i]=FAtoms
            i=i+1

        _IOP = re.compile(r'\b(IOp.*)',re.MULTILINE).findall(output)[0].split('=')[-1].strip()
        if _IOP == '':
            for functional in self.dfunctionalps.keys():
                if self.g09functional.strip() in functional:
                    self.functionalpar = self.dfunctionalps[functional]
        else:
            self.functionalpar = _IOP

        return


    def debug(self):
        print("g09nprocshared=",self.g09nprocshared)
        print("g09mem=",self.g09mem)
        print("g09functional=",self.g09functional)
        print("g09basisset=",self.g09basisset)
        print("Extra-command=",self.extracommand)
        print("xyzfilepath=",self.xyzfilepath)
        print("multiwfnpath=",self.multiwfnpath)
        print("dCM=",self.dCM)
        print("IRCgap=",self.IRCgap)
        print("dAIndex=",self.dAIndex)
        print("dSCFe=",self.dSCFe)
        print("dbonds=",self.dbonds)
        print("dpyramidalization=",self.dpyramidalization)
        print("dcentroids=",self.dcentroids)
        print("dangle=",self.dangle)
        print("ddihedral=",self.ddihedral)
        print("dhirshfeld=",self.dhirshfeld)
        print("functionalpar=",self.functionalpar)
        return
