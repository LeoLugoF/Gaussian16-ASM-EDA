import GlobalVariables as GV
from GenReadFileProp import ReadFileProp

class calcASMGfiles:
    #Generate instance of the GV class
    GlobalVars: GV.ASMG_parameters
    readfileprop: ReadFileProp
    #This is the first module
    def __init__(self):
        bConfExists:bool = GV.PrintConf()
        self.readfileprop = ReadFileProp()
        # If ASMG_conf.txt exists, then read and store variables 
        # from configuration file.
        if (bConfExists == True):
            self.GlobalVars = GV.ASMG_parameters()
            #self.GlobalVars.debug()
            # Generate and calculate fragment*.gjf from IRC xyz file.
            self.GlobalVars.lMCcoors,self.GlobalVars.TSXpoint = self.readfileprop.readIRC(self.GlobalVars.xyzfilepath)
        return

    def M2_frags(self,run=True):
        self.readfileprop.GenerateFragmentGJF(self.GlobalVars)
        self.readfileprop.GenerateJobFile('fragment',self.GlobalVars)
        if run:
            self.readfileprop.RunJobFile('fragment')
        return

    def M3_promol(self):
        self.readfileprop.GenerateformartedChk(self.GlobalVars)
        self.readfileprop.GeneratePromolGJF(self.GlobalVars)
        self.readfileprop.GenerateJobFile('promol',self.GlobalVars)
        self.readfileprop.RunJobFile('promol')
        return

    def M4_Mcomplex(self):
        self.readfileprop.GenerateMComplexGJF(self.GlobalVars)
        self.readfileprop.GenerateJobFile('mcomplex',self.GlobalVars)
        self.readfileprop.RunJobFile('mcomplex')
        return

    def M5_energies(self):
        self.readfileprop.ReadEnergiesFromAllFiles(self.GlobalVars)
        self.readfileprop.CalculateEDAandStrain(self.GlobalVars)
        self.readfileprop.PrintEnergies(self.GlobalVars)
        self.readfileprop.PrintXvalues(self.GlobalVars)

        return
