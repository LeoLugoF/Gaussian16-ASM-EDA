import re, os
import sys, glob
import GlobalVariables
import subprocess, math
from itertools import permutations

class ReadFileProp:

    def __init__(self):

        return

    def GetlistofcalculatedfilesClass(self,GVariables:GlobalVariables.ASMG_parameters) -> list:
        TSXpoint = GVariables.TSXpoint
        spacepoints = GVariables.IRCgap
        lMCcoor = GVariables.lMCcoors

        nRC1 = list(range(int(TSXpoint),-1,-int(spacepoints)))
        nRC1.reverse()
        nRC2 = list(range(int(TSXpoint)+int(spacepoints),len(lMCcoor),int(spacepoints)))
        nRC = nRC1 + nRC2
        nRC.sort()
        return nRC

    def Getlistofcalculatedfiles(self,typeoffile:str,extension:str) -> list:
        #Returns a list with the indexes of each fragment calculation#
        #extension=*.log
        #filename=fragment1
        nRC = [ int(x.split('_')[-1].split('.')[0]) for x in glob.glob(extension) if typeoffile in x]
        nRC = list(dict.fromkeys(nRC))
        nRC.sort()
        return nRC

    def GenerateFragmentGJF(self,GVariables:GlobalVariables.ASMG_parameters):
        nRC = self.GetlistofcalculatedfilesClass(GVariables)
        for ircpoint in nRC:
            for fragment in GVariables.dAIndex.keys():
                sfnxyz = self.definefragmentsXYZ(GVariables.lMCcoors[ircpoint],GVariables.dAIndex[fragment])
                self.WriteGJFfile(fragment+1,"fragment"+str(fragment+1)+"_"+str(ircpoint),
                                  sfnxyz,GVariables.fragmentcommand,
                                  GVariables)
        return

    def GeneratePromolGJF(self,GVariables:GlobalVariables.ASMG_parameters):
        nRC = self.Getlistofcalculatedfiles('fragment1',"*.log")
        for ircpoint in nRC:
            sfnxyz = self.definefragmentsXYZ(GVariables.lMCcoors[ircpoint],[])
            self.WriteGJFfile(0,"promol_"+str(ircpoint),
                                sfnxyz,GVariables.promolcommand,
                                GVariables)
        return

    def GenerateMComplexGJF(self,GVariables:GlobalVariables.ASMG_parameters):
        nRC = self.Getlistofcalculatedfiles('promol',"*.log")
        for ircpoint in nRC:
            sfnxyz = self.definefragmentsXYZ(GVariables.lMCcoors[ircpoint],[])
            self.WriteGJFfile(0,"mcomplex_"+str(ircpoint),
                                sfnxyz,GVariables.mcomplexcommand,
                                GVariables)
        return

    def GenerateformartedChk(self,GVariables:GlobalVariables.ASMG_parameters):
        nRC = self.Getlistofcalculatedfiles('fragment1',"*.log")
        for ircpoint in nRC:
            for fragment in GVariables.dAIndex.keys():
                filenamechk = "fragment"+str(fragment+1)+"_"+str(ircpoint)+".chk"
                os.system("formchk "+filenamechk + " > /dev/null")
        return

    def GenerateJobFile(self,typeoffile,GVariables:GlobalVariables.ASMG_parameters):
        #Type of file=
        # 1) fragment, 2) promol or 3) mcomplex
        nRC = self.Getlistofcalculatedfiles(typeoffile,"*.gjf")
        File = open(typeoffile+"_job.sh",'w')
        File.write("#!/bin/bash\n")
        if 'fragment' in typeoffile:
            for fragment in GVariables.dAIndex.keys():
                snfragment = str(fragment + 1)
                for RC in nRC:
                    filename='fragment'+snfragment+"_"+str(RC)
                    File.write(GVariables.gaussianversion+" < "+
                               filename+".gjf > "+filename+".log\n")

        else:
            for RC in nRC:
                filename=typeoffile+"_"+str(RC)
                if 'promol' in typeoffile:
                    self.GeneratePromolchkwithMultiwfn(GlobalVariables.ASMG_parameters,RC)
                result = os.system(GVariables.gaussianversion+" < "+
                            filename+".gjf > "+filename+".log\n")
                if result != 0:
                    print("ERROR! Input Density of "+filename+".chk file is not normalized by Multiwfn.")
                    print("Continuing the next calculations...")

        File.close()
        return

    def RunJobFile(self,typeoffile):
        #Type of file=
        # 1) fragment, 2) promol or 3) mcomplex
        bashinput=""
        if "fragment" in typeoffile:
            bashinput = "fragment_job.sh"
        if "promol" in typeoffile:
            bashinput = "promol_job.sh"
        if "mcomplex" in typeoffile:
            bashinput = "mcomplex_job.sh"

        os.system("bash "+bashinput)
        #File = open("ASMG_details.log","w")
        #File.write("All fragments were calculated!")
        #File.close()

        return

    def WriteGJFfile(self,nCM:int,filenamewithoutext:str,xyz:str,commands:str,GVariables:GlobalVariables.ASMG_parameters):
        File = open(filenamewithoutext+".gjf","w")
        File.write("%nprocshared="+GVariables.g09nprocshared+"\n")
        File.write("%mem="+GVariables.g09mem+"\n")
        if 'fragment' in filenamewithoutext:
            File.write("%chk="+filenamewithoutext+".chk"+"\n")
            File.write("#p "+GVariables.g09functional+'/'+GVariables.g09basisset+" "+commands+" "+GVariables.extracommand+"\n\n")
            File.write(filenamewithoutext+"\n\n")
            File.write(GVariables.dCM[nCM][0]+" "+GVariables.dCM[nCM][1]+"\n")
            File.write(xyz+"\n")
            File.write(GVariables.functionalpar+" "+"5"+"\n\n")
        else:
            File.write("%chk="+"mcomplex_"+filenamewithoutext.split('_')[-1]+".chk"+"\n")
            File.write("%oldchk="+"promol_"+filenamewithoutext.split('_')[-1]+".chk"+"\n")
            File.write("#p "+GVariables.g09functional+'/'+GVariables.g09basisset+" "+commands+" "+GVariables.extracommand+" geom=check\n\n")
            File.write(filenamewithoutext+"\n\n")
            File.write(GVariables.dCM[nCM][0]+" "+GVariables.dCM[nCM][1]+"\n")
            File.write("\n")
            File.write(GVariables.functionalpar+" "+"5"+"\n\n")
        File.close()


        return

    def readIRC(self,xyzFilePath):
        File = open(xyzFilePath,"r")
        output = File.read()
        File.close()
        raw_totalcoors = re.compile(r'^\w+\s+[-+]?\d+\.\d+\s+[-+]?\d+\.\d+\s+[-+]?\d+\.\d+',re.MULTILINE).findall(output)
        if len(raw_totalcoors) == 0:
            raw_totalcoors = re.compile(r'\w+\s+[-+]?\d+\.\d+\s+[-+]?\d+\.\d+\s+[-+]?\d+\.\d+',re.MULTILINE).findall(output)    
        natoms = int(output.split('\n')[0])
        #Process Cartesian Coordinates into list
        xo = 0
        xf = natoms
        #cartesian coordinates of Molecular complexes 
        totalcoors = []

        #Generate cartesian coordinates of molecular complexes
        while xo < len(raw_totalcoors):
            ln = raw_totalcoors[xo:xf]
            RCcoor = []
            for coor in ln:
                xyzline = coor.split()
                xyzline = '     '.join(xyzline)
                RCcoor.append(xyzline)
            totalcoors.append(RCcoor)
            #totalcoors.append('\n'.join(RCcoor))
            #totalcoors[-1] = totalcoors[-1] + "\n\n\n\n\n"
            xf+=natoms
            xo+=natoms

            # Get TS X point
        energies = []
        lines = output.split('\n')
        for line in lines:
            raw_line = line.split('.')
            if len(raw_line) == 2:
                energies.append(float(line))

        TSXindex = int(energies.index(max(energies)))

        return totalcoors, TSXindex

    def definefragmentsXYZ(self,lxyzcoors,lAtomIndexes) -> str:
        newcoors=[]
        i=1
        for coors in lxyzcoors:
            lcoors = coors.split()
            if lAtomIndexes == []:
                lcoors[0] = lcoors[0]+"   "
                newcoors.append('    '.join(lcoors))
                continue

            if i in lAtomIndexes:
                lcoors[0] = lcoors[0]+"   "
            #else:
            #    lcoors[0] = lcoors[0]+"-Bq"
                newcoors.append('    '.join(lcoors))
            i = i + 1
        stringnewcoors='\n'.join(newcoors)+"\n"
        return stringnewcoors

    def get_all_permutations(self,input_list):
        # Use permutations to generate all possible permutations
        all_permutations = permutations(input_list)

        # Convert the permutations to a list and return
        # Returns : [(1,2,3,...,n),...]
        return list(all_permutations)

    def testforcorrectpromolchk(self,GVariables:GlobalVariables.ASMG_parameters, IRCpoint):
        testpoint = IRCpoint
        multiwfnpath="Multiwfn"
        if GVariables.multiwfnpath != "":
            multiwfnpath=GVariables.multiwfnpath
        ltotalfragments=list(GVariables.dAIndex.keys())
        permutations = self.get_all_permutations(ltotalfragments)
        for permutation in permutations:
            sfragments = ""
            for nf in permutation:
                if nf == permutation[0]:
                    continue
                sfragments = sfragments + "fragment" + str(nf+1) + "_" + str(testpoint) + ".fchk\n"
            multiwfnargs = 'fragment' + str(permutation[0]+1) + "_" + str(testpoint) + ".fchk\n"
            multiwfnargs = multiwfnargs + '100\n19\n3\n' + str(len(ltotalfragments)) + "\n" + sfragments + '0\n\q\n'
            p = subprocess.Popen(multiwfnpath, universal_newlines=True,stdin=subprocess.PIPE,stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)
            stdO = p.communicate(multiwfnargs)
            os.system("mv combine.fch "+" promol_"+str(testpoint)+".fch")
            os.system("unfchk "+ " promol_"+str(testpoint)+".fch" + " promol_"+str(testpoint)+".chk > /dev/null")
            result = os.system(GVariables.gaussianversion + " < promol_"+str(testpoint)+".gjf > " +  " promol_"+str(testpoint)+".log > /dev/null 2>&1" )
            #/dev/null 2>&1
            

            
            if result == 0:
                GVariables.lCorrectseqoffragments = list(permutation).copy()
                break
        return

    def GeneratePromolchkwithMultiwfn(self,GVariables:GlobalVariables.ASMG_parameters,RC):
        nRC = self.Getlistofcalculatedfiles('fragment1',"*.log")
        #multiwfnpath = GVariables.multiwfnpath
        multiwfnpath="Multiwfn"
        if GVariables.multiwfnpath != "":
            multiwfnpath=GVariables.multiwfnpath

        if(RC == nRC[0]):
            self.testforcorrectpromolchk(GlobalVariables.ASMG_parameters,nRC[0])

        sfragments=""
        n = len(GVariables.lCorrectseqoffragments)
        if (n == 0):
            if(RC != nRC[0]):
                self.testforcorrectpromolchk(GlobalVariables.ASMG_parameters,RC)
            return 

        for i in GVariables.lCorrectseqoffragments:
            if i == GVariables.lCorrectseqoffragments[0]:
                continue
            sfragments = sfragments + "fragment" + str(i+1) + "_" + str(RC) + ".fchk\n"
        multiwfnargs= 'fragment'+ str(GVariables.lCorrectseqoffragments[0]+1) + '_' + str(RC) + '.fchk\n' + '100\n19\n3\n' + str(n) + "\n" + sfragments + '0\n\q\n' #add \n-
        p = subprocess.Popen(multiwfnpath, universal_newlines=True,stdin=subprocess.PIPE,stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)
        stdO = p.communicate(multiwfnargs)
        os.system("mv combine.fch "+" promol_"+str(RC)+".fch")
        os.system("unfchk "+ " promol_"+str(RC)+".fch" + " promol_"+str(RC)+".chk > /dev/null")
        return 

    def ReadEnergiesFromOutputs(self,GVariables:GlobalVariables.ASMG_parameters):
        for nfragment in GVariables.dAIndex.keys():
            i = nfragment+1
            GVariables.dDecomposedE[i]={}
            nRC = self.Getlistofcalculatedfiles('fragment','*.log')
            for nfile in nRC:
                filename='fragment'+str(i)+"_"+str(nfile)+".log"
                GVariables.dDecomposedE[i][nfile]=self.DecompInteractionEnergy(filename)
        GVariables.dDecomposedE[0]={}
        nRC = self.Getlistofcalculatedfiles('promol','*.log')
        for nfile in nRC:
            filename='promol_'+str(nfile)+".log"
            GVariables.dDecomposedE[0][nfile]=self.DecompInteractionEnergy(filename)
        return

    def read_energy(self,File) -> float:
        #Reads SCF energy from outputfile
        inp = open(File,'r')
        gout = inp.read()
        inp.close()
        energyline = re.compile(r'SCF Done:'+r'.*?'+ r'A.U. after', re.IGNORECASE|re.DOTALL).findall(gout) 

        if (len(energyline) == 0):
            print("Cannot read SCF energy of file: " + File)
            sys.exit()

        energy = energyline[0].split()
        return float(energy[4])

    def DecompInteractionEnergy(self,i : int,EInt : float) -> list:
        filename = "Complex_"+str(i)+".log"
        ofile = open(filename,'r')
        output = ofile.read()
        ofile.close()
        firstscf = re.compile(r'^\s\w\=\s\-\d+\.\d+\s+',re.MULTILINE).findall(output)[0].strip()
        firstscf = float(firstscf.split()[-1])
        scf = self.read_energy(filename)
        Eorb = round((scf - firstscf)*627.51,2)
        Epaupele = round((EInt - Eorb),2)
        return [Eorb,Epaupele]

    def GetDecomposedEnergies(self,filename) -> list: 
        File = open(filename,'r')
        output = File.read()
        File.close()
        
        regexd = re.compile(r'(N-N=.*)',re.MULTILINE).findall(output)
        if (len(regexd) == 0):
            return [0,0,0,0,0]
        rawenergies=regexd[-1]
        #Get ENTJ
        ENVT = re.findall(r'=\s*([+-]?\d+\.\d+)[DdEe]([+-]?\d+)', rawenergies)
        ENVT = [float(f"{base}e{exp}") for base, exp in ENVT]
        #Get EJ
        regexdEJ=re.findall(r'EJ=\s*([-+]?\d+\.\d+)', output)
        EJ=float(regexdEJ[-1])
        #Get ENTVJ
        ENTVJ=ENVT[0]+ENVT[1]+ENVT[2]+EJ
        #Get Ex and Ec
        matches = re.search(r'Ex=\s*([-+]?\d+\.\d+).*?Ec=\s*([-+]?\d+\.\d+)', output)
        Ex=float(matches.group(1))
        Ec=float(matches.group(2))
        #Insert Energies
        ldecomposedE = []
        ldecomposedE.append(ENTVJ)
        ldecomposedE.append(Ex)
        ldecomposedE.append(Ec)
        Edisp=float(re.compile(r'(Dispersion energy.*)',re.MULTILINE).findall(output)[0].split()[-2])
        Escf =float(re.compile(r'SCF Done:'+r'.*?'+ r'A.U. after', re.IGNORECASE|re.DOTALL).findall(output)[0].split()[4])
        ldecomposedE.append(Edisp)
        ldecomposedE.append(Escf)
        return ldecomposedE

    def GetFinalandFrozenstateSCF(self,nfilename:int) -> list:
        File_MComplex = open('mcomplex_'+str(nfilename)+'.log','r')
        out_MComplex = File_MComplex.read()
        File_MComplex.close()
        regexd = re.compile(r'^\s\w\=\s\-\d+\.\d+\s+',re.MULTILINE).findall(out_MComplex)
        if len(regexd) > 0:
            firstscf = regexd[0].strip()
            firstscf = float(firstscf.split()[-1])
            scf = self.read_energy('mcomplex_'+str(nfilename)+'.log')
        else:
            firstscf=0
            scf=0

        return [scf,firstscf]

    def ReadEnergiesFromAllFiles(self,GVariables:GlobalVariables.ASMG_parameters):
        #Read energies of fragment and promol files
        #Promol
        nRC = self.Getlistofcalculatedfiles('fragment1',"*.log")
        for i in range(0,len(GVariables.dAIndex.keys())+1):
            GVariables.dDecomposedE[i] = {}
        #Fragments, final and frozen state
        for RC in nRC:
            promol_filename = 'promol_'+str(RC)+'.log'
            GVariables.dDecomposedE[0][RC] = self.GetDecomposedEnergies(promol_filename)
            GVariables.dEscf[RC] = self.GetFinalandFrozenstateSCF(RC)
            for i in range(0,len(GVariables.dAIndex.keys())):
                n = i+1
                fragment_filename = 'fragment' + str(n) + '_'+str(RC)+'.log'
                GVariables.dDecomposedE[n][RC] = self.GetDecomposedEnergies(fragment_filename)
        return

    def CalculateEDA(self,GVariables:GlobalVariables.ASMG_parameters):
        nRC = self.GetlistofcalculatedfilesClass(GVariables)
        for RC in nRC:
            GVariables.dEDA[RC] = []
            for typeofenergy in range(0,4):
                # 0,1,2,3 = Eels+disp,Ex,Ec,Edisp
                sumfragments_energies = 0.0
                for i in GVariables.dAIndex.keys():
                    nfrag = i + 1
                    sumfragments_energies = sumfragments_energies + GVariables.dDecomposedE[nfrag][RC][typeofenergy]
                En = round((GVariables.dDecomposedE[0][RC][typeofenergy] - sumfragments_energies)*627.51,2)
                GVariables.dEDA[RC].append(En)

        for RC in nRC:
            Eorb = GVariables.dEscf[RC][0]- GVariables.dEscf[RC][1]
            Eorb = round(Eorb*627.51,2)
            Epauli = GVariables.dEscf[RC][1] - GVariables.dDecomposedE[0][RC][4]
            Epauli = round(Epauli*627.51,2)
            # 4,5 = Eorb,Epauli
            GVariables.dEDA[RC].append(Eorb)
            GVariables.dEDA[RC].append(Epauli)

            #Correct electrostatic energy by resting to it the dispersion energy
            GVariables.dEDA[RC][0] = round(GVariables.dEDA[RC][0] - GVariables.dEDA[RC][3],2)

        #Calculate interaction energy
        for RC in nRC:
            sumfragments_energies = 0.0
            for i in GVariables.dAIndex.keys():
                nfrag = i + 1
                sumfragments_energies = sumfragments_energies + GVariables.dDecomposedE[nfrag][RC][4]
            Eint = round(( GVariables.dEscf[RC][0] - sumfragments_energies)*627.51,2)
            GVariables.dEDA[RC].append(Eint)
        return

    def CalculateStrain(self,GVariables:GlobalVariables.ASMG_parameters):
        nRC = self.GetlistofcalculatedfilesClass(GVariables)

        GVariables.dFragmentStrain[0] = {}
        for i in GVariables.dAIndex.keys():
            nfrag = i + 1
            GVariables.dFragmentStrain[nfrag] = {}

        for RC in nRC:
            TotalStrain = 0.0
            TotalSCF_fragments = 0.0
            for i in GVariables.dAIndex.keys():
                nfrag = i + 1
                fragment_strain = round((GVariables.dDecomposedE[nfrag][RC][4] - GVariables.dSCFe[i])*627.51,2)
                GVariables.dFragmentStrain[nfrag][RC] = fragment_strain
                TotalStrain = TotalStrain + fragment_strain
                TotalSCF_fragments = TotalSCF_fragments + GVariables.dSCFe[i]
            GVariables.dFragmentStrain[0][RC] = round(TotalStrain,2)
            #Calculate Total Energy
            Etotal = round((GVariables.dEscf[RC][0] - TotalSCF_fragments)*627.51,2)
            GVariables.dEDA[RC].append(Etotal)

        return

    def CalculateEDAandStrain(self,GVariables:GlobalVariables.ASMG_parameters):
        self.CalculateEDA(GVariables)
        self.CalculateStrain(GVariables)
        return

    def PrintEnergies(self,GVariables:GlobalVariables.ASMG_parameters):
        File = open("ASMG_energies.txt",'w')
        nRC = self.GetlistofcalculatedfilesClass(GVariables) #TOTAL SPACES 8
        headline = "mcomplex_n  Eels    Ex      Ec      Edisp   Eorb    Epauli  Eint  Etotal  TStrain "
        fragments_line=""
        for i in GVariables.dAIndex.keys():
            nfrags = i + 1
            fragments_line = fragments_line + "Strain" + str(nfrags)+" "
        headline = headline + fragments_line + "\n" 
        File.write(headline)

        for RC in nRC:
            el = str(RC)
            for energy in GVariables.dEDA[RC]:
                el = el + "       " + str(energy)
            for strain in GVariables.dFragmentStrain.keys():
                el = el + "       " + str( GVariables.dFragmentStrain[strain][RC])
            File.write(el+"\n")

        File.close()

        os.system("column -t ASMG_energies.txt > ASMG-Energies.txt")
        os.system("rm ASMG_energies.txt")

        return

    def GetXYZCoorsInlist(self,IRCpoint:int,lAIndexes:list,GVariables:GlobalVariables.ASMG_parameters):
        # OutCoors = [ [x1,y1,z1], [x2,y2,z2], ... , [xn,yn,zn] ]
        # where n = len(lAIndexes) 
        OutMatrix = []

        if lAIndexes == []:
            numberofatoms = len(GVariables.lMCcoors[IRCpoint])
            for Ai in range(0,numberofatoms):
                a=GVariables.lMCcoors[IRCpoint][int(Ai)].split()
                a.pop(0)
                b=[]
                b.append(float(a[0]))
                b.append(float(a[1]))
                b.append(float(a[2]))
                OutMatrix.append(b)
            return OutMatrix
        for Ai in lAIndexes:
            a=GVariables.lMCcoors[IRCpoint][int(Ai)-1].split()
            a.pop(0)
            b=[]
            b.append(float(a[0]))
            b.append(float(a[1]))
            b.append(float(a[2]))
            OutMatrix.append(b)
        return OutMatrix

    def GetDistanceXYZ(self,Coor1 : list, Coor2 : list):
        lAi1 = Coor1
        lAi2 = Coor2
        Sum=pow(float(lAi1[0])-float(lAi2[0]),2)+pow(lAi1[1]-lAi2[1],2)+pow(lAi1[2]-lAi2[2],2)
        distance=round(pow(Sum,0.5),4)
        return distance

    def ReadDistances(self,GVariables:GlobalVariables.ASMG_parameters):
        nX = 1
        nRC = self.GetlistofcalculatedfilesClass(GVariables)
        for lbonds in GVariables.dbonds.values():
            GVariables.dprintxaxis["Bond"+str(nX)] = []
            for RC in nRC:
                lAi = self.GetXYZCoorsInlist(RC,lbonds,GVariables)
                distance=self.GetDistanceXYZ(lAi[0],lAi[1])
                GVariables.dprintxaxis["Bond"+str(nX)].append(distance)
            nX=nX+1
        return

    def GetAngleXYZ(self,Coor1 : list, Coor2 : list, Coor3 : list):
        c1=Coor1
        c2=Coor2
        c3=Coor3
        v1=[float(c1[0])-float(c2[0]),
            float(c1[1])-float(c2[1]),
            float(c1[2])-float(c2[2])]
        v2=[float(c3[0])-float(c2[0]),
            float(c3[1])-float(c2[1]),
            float(c3[2])-float(c2[2])]
        v1mag=0
        for i in v1:
            v1mag=v1mag+math.pow(i,2)
        v1mag=math.sqrt(v1mag)
        v2mag=0
        for j in v2:
            v2mag=v2mag+math.pow(j,2)
        v2mag=math.sqrt(v2mag)
    
        v1norm=[v1[0]/v1mag,v1[1]/v1mag,v1[2]/v1mag] 
        v2norm=[v2[0]/v2mag,v2[1]/v2mag,v2[2]/v2mag]

        dotProduct=0
        for i in [0,1,2]:
            dotProduct=dotProduct+v1norm[i]*v2norm[i]

        angle=round((math.acos(dotProduct)*180)/math.pi,4)
    
        return angle

    def ReadAngles(self,GVariables:GlobalVariables.ASMG_parameters):
        nX = 1
        nRC = self.GetlistofcalculatedfilesClass(GVariables)
        for langle in GVariables.dangle.values():
            GVariables.dprintxaxis["Angle"+str(nX)] = []
            for RC in nRC:
                lAi = self.GetXYZCoorsInlist(RC,langle,GVariables)
                angle = self.GetAngleXYZ(lAi[0],lAi[1],lAi[2])
                GVariables.dprintxaxis["Angle"+str(nX)].append(angle)
            nX = nX + 1
        return

    def ReadDihedral(self,GVariables:GlobalVariables.ASMG_parameters):
        nX = 1
        nRC = self.GetlistofcalculatedfilesClass(GVariables)
        for ldihedral in GVariables.ddihedral.values():
            GVariables.dprintxaxis["Dihedral"+str(nX)] = []
            for RC in nRC:
                lAi = self.GetXYZCoorsInlist(RC,ldihedral,GVariables)
                dangle = self.GetDihedralAngle(lAi[0],lAi[1],lAi[2],lAi[3])
                GVariables.dprintxaxis["Dihedral"+str(nX)].append(round(dangle,4))
            nX = nX + 1
        return

    def ReadPyramidalizationAngle(self,GVariables:GlobalVariables.ASMG_parameters):
        nX = 1
        nRC = self.GetlistofcalculatedfilesClass(GVariables)
        for lpyr in GVariables.dpyramidalization.values():
            GVariables.dprintxaxis["Pyramidalization"+str(nX)] = []
            for RC in nRC:
                lAi = self.GetXYZCoorsInlist(RC,lpyr,GVariables)
                theta12=self.GetAngleXYZ(lAi[0],lAi[1],lAi[2])
                theta23=self.GetAngleXYZ(lAi[3],lAi[4],lAi[5])
                theta31=self.GetAngleXYZ(lAi[6],lAi[7],lAi[8])

                cos12=math.cos(theta12*math.pi/180)
                cos23=math.cos(theta23*math.pi/180)
                cos31=math.cos(theta31*math.pi/180)
                sin12=math.sin(theta12*math.pi/180)

                x1=-1
                x2=-cos12
                y2=sin12
                x3=-cos31
                y3=(cos23-cos12*cos31)/sin12
                x33=math.sqrt(1-math.pow(x3,2)-math.pow(y3,2))
                x1y2z3=x1*y2*x33
                D1=math.pow(y2*x33,2)
                D2=math.pow((x2-x1)*x33,2)
                D3=math.pow((x2-x1)*y3-y2*(x3-x1),2)
                cos_sp=x1y2z3/math.sqrt(D1+D2+D3)
                theta_sp=math.acos(cos_sp)*180/math.pi
                Pyramidalization_tetha=theta_sp-90
                Pyramidalization_tetha=round(Pyramidalization_tetha,4)
                GVariables.dprintxaxis["Pyramidalization"+str(nX)].append(Pyramidalization_tetha)
            nX = nX + 1
        return

    def ReadCentroidsDistance(self,GVariables:GlobalVariables.ASMG_parameters):
        nX = 1
        nRC = self.GetlistofcalculatedfilesClass(GVariables)
        for dcent in GVariables.dcentroids.values():
            GVariables.dprintxaxis["Centroids"+str(nX)] = []
            for RC in nRC:
                lAi1 = self.GetXYZCoorsInlist(RC,dcent[0],GVariables)
                lAi2 = self.GetXYZCoorsInlist(RC,dcent[1],GVariables)

                centroid1 = [0.0,0.0,0.0]
                for Ai1 in lAi1:
                    for i in range(0,3):
                        centroid1[i] = centroid1[i]+(Ai1[i]/len(lAi1))
                centroid2 = [0.0,0.0,0.0]
                for Ai2 in lAi2:
                    for i in range(0,3):
                        centroid2[i] = centroid2[i]+(Ai2[i]/len(lAi2))
                centroid_distance = self.GetDistanceXYZ(centroid1,centroid2)
                GVariables.dprintxaxis["Centroids"+str(nX)].append(centroid_distance)
            nX = nX + 1
        return

    def euclidean_distance(self,point_a, point_b):
        return math.sqrt(sum((a - b) ** 2 for a, b in zip(point_a, point_b)))

    def rmsd(self,matrix_a, matrix_b):

        n = len(matrix_a)
        squared_distances = [self.euclidean_distance(matrix_a[i], matrix_b[i]) ** 2 for i in range(n)]
        mean_squared_distance = sum(squared_distances) / n
        rmsd_value = math.sqrt(mean_squared_distance)

        return rmsd_value

    def ReadRmsd(self,GVariables:GlobalVariables.ASMG_parameters):
        nRC = self.GetlistofcalculatedfilesClass(GVariables)
        for fragment in GVariables.dAIndex.keys():
            lA1 = []
            GVariables.dprintxaxis["RMSD_F"+str(fragment+1)] = []
            for RC in nRC:
                lAi = self.GetXYZCoorsInlist(RC,GVariables.dAIndex[fragment],GVariables) 
                if RC == nRC[0]:
                    lA1 = lAi.copy()
                rmsd = self.rmsd(lAi,lA1)
                GVariables.dprintxaxis["RMSD_F"+str(fragment+1)].append(rmsd)
        lA1 = []
        GVariables.dprintxaxis["RMSD_MC"] = []
        for RC in nRC:
            lAi = self.GetXYZCoorsInlist(RC,[],GVariables) 
            if RC == nRC[0]:
                lA1 = lAi.copy()
            rmsd = self.rmsd(lAi,lA1)
            GVariables.dprintxaxis["RMSD_MC"].append(rmsd)
        return

    def calculate_cross_product(self,v1,v2):
        return [v1[1]*v2[2] - v1[2]*v2[1], v1[2]*v2[0] - v1[0]*v2[2], v1[0]*v2[1] - v1[1]*v2[0]]

    def calculate_dot_product(self,v1,v2):
        return sum(x*y for x, y in zip(v1, v2))

    def GetDihedralAngle(self,atom1,atom2,atom3,atom4):

        b1 = [x - y for x, y in zip(atom1, atom2)]
        b2 = [x - y for x, y in zip(atom3, atom2)]
        b3 = [x - y for x, y in zip(atom4, atom3)]

        # Normalize the vectors
        b1_magnitude = sum(x**2 for x in b1)**0.5
        b2_magnitude = sum(x**2 for x in b2)**0.5
        b3_magnitude = sum(x**2 for x in b3)**0.5

        b1 = [x / b1_magnitude for x in b1]
        b2 = [x / b2_magnitude for x in b2]
        b3 = [x / b3_magnitude for x in b3]

        # Calculate the normals to the planes
        n1 = self.calculate_cross_product(b1, b2)
        n2 = self.calculate_cross_product(b2, b3)

        # Calculate the dihedral angle
        m1 = self.calculate_cross_product(n1, b2)
        x = self.calculate_dot_product(m1, n2)
        y = self.calculate_dot_product(n1, n2)

        return 180.0 / 3.14159265358979323846 * math.atan2(x, y)



    def GetHirshfeldCharges(self,GVariables:GlobalVariables.ASMG_parameters):
        nRC = self.GetlistofcalculatedfilesClass(GVariables)
        for RC in nRC:
            MCfile = open("mcomplex_"+str(RC)+".log","r")
            output = MCfile.read()
            MCfile.close()

            lrawhirshfeld = re.compile(r'^\s+\d+\s+\w+\s+[-+]?\d+.\d+\s+\s+[-+]?\d+.\d+\s+[+-]?\d+.\d+\s+[-+]?\d+.\d+\s+[-+]?\d+.\d+\s+[-+]?\d+.\d+',re.MULTILINE).findall(output)
            if(len(lrawhirshfeld) == 0):
                continue
            totalatoms = int(lrawhirshfeld[-1].split()[0])
            lhirshfeld = lrawhirshfeld[-totalatoms:] 
            lhirshfelf_patitioned = []
            for hirsh in lhirshfeld:
                a = [ float(x) for x in hirsh.split() if len(x.split(".")) == 2 ] 
                lhirshfelf_patitioned.append(a)
            GVariables.dcalchirshfeld[RC] = []
            GVariables.dcalchirshfeld[RC] = lhirshfelf_patitioned.copy()
        return

    def SumHirshfeldCharges(self,lAHirshfeld:list,RC:int,GVariables:GlobalVariables.ASMG_parameters):
        sumHirshfeld = 0.0
        for rawatomIn in lAHirshfeld:
            if RC not in GVariables.dcalchirshfeld:
                return 'ERROR'
            Hirshfeld = GVariables.dcalchirshfeld[RC][rawatomIn-1][GVariables.typeofcharge]
            sumHirshfeld = round(sumHirshfeld + Hirshfeld,4)
        return sumHirshfeld


    def ReadHirshfeld(self,GVariables:GlobalVariables.ASMG_parameters):
        self.GetHirshfeldCharges(GVariables)
        nRC = self.GetlistofcalculatedfilesClass(GVariables)
        #Sum Hirshfeld Charges
        i = 1
        for lAiHirshfeld in GVariables.dhirshfeld.values():
            GVariables.dprintxaxis["Hirshfeld"+str(i)] = []
            for RC in nRC:
                GVariables.dprintxaxis["Hirshfeld"+str(i)].append(self.SumHirshfeldCharges(lAiHirshfeld,RC,GVariables))

            i = i + 1
        return

    def PrintXFile(self,GVariables:GlobalVariables.ASMG_parameters):

        height = len(list(GVariables.dprintxaxis.values())[-1])

        XFile = open("ASMG_Xaxis.txt","w")
        header = "mcomplex_n    "
        for key in GVariables.dprintxaxis.keys():
            header = header + key + "    "
        XFile.write(header+"\n")
        nRC = self.GetlistofcalculatedfilesClass(GVariables)
        for i in range(0,height):
            line_xpar = str(nRC[i]) + "    "
            for key in GVariables.dprintxaxis.keys():
                line_xpar = line_xpar + str(GVariables.dprintxaxis[key][i]) + "    "
            XFile.write(line_xpar+"\n")
        XFile.close()

        os.system("column -t ASMG_Xaxis.txt > ASMG-Xaxis.txt")
        os.system("rm ASMG_Xaxis.txt")

        return

    def PrintXvalues(self,GVariables:GlobalVariables.ASMG_parameters):
        self.ReadDistances(GVariables)
        self.ReadAngles(GVariables)
        self.ReadDihedral(GVariables)
        self.ReadPyramidalizationAngle(GVariables)
        self.ReadCentroidsDistance(GVariables)
        self.ReadHirshfeld(GVariables)
        self.ReadRmsd(GVariables)
        self.PrintXFile(GVariables)
        return
