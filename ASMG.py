import ASMGmodules as ASMG
import argparse, sys, os
from IRC import organizeIRC
from GlobalVariables import PrintConf

#Parser
parser = argparse.ArgumentParser(description='Activation Strain Model for Gaussian : Author: Leonardo Israel Lugo Fuentes')
group = parser.add_mutually_exclusive_group()
group.add_argument('--IRC',action='store_true',help="Orders the IRC *.log (calculated with LQA) from reactants to products and transforms it into a *.xyz file")
group.add_argument('--r1', action='store_true',help="Generate inputs of fragments.")
group.add_argument('--r2',action='store_true',help="Calculate promol and mcomplex files.")
group.add_argument('--r3',action='store_true',help="Calculate mcomplex files.")
group.add_argument('--rdata',action='store_true',help="Generate energetic and x-axis data.")
group.add_argument('--run',action='store_true',help="Calculate fragments, promol, mcomplex and generate energetic/x-axis data).")

args = parser.parse_args()

if not os.path.exists('ASMG_conf.txt'):
    PrintConf()

if args.r1:
    ASMG_module = ASMG.calcASMGfiles()
    ASMG_module.M2_frags(False)

if args.r2:
    ASMG_module = ASMG.calcASMGfiles()
    ASMG_module.M3_promol()
    ASMG_module.M4_Mcomplex()

if args.r3:
    ASMG_module = ASMG.calcASMGfiles()
    ASMG_module.M4_Mcomplex()

if args.rdata:
    ASMG_module = ASMG.calcASMGfiles()
    ASMG_module.M5_energies()

if args.run:
    ASMG_module = ASMG.calcASMGfiles()
    ASMG_module.M2_frags()
    ASMG_module.M3_promol()
    ASMG_module.M4_Mcomplex()
    ASMG_module.M5_energies()

if args.IRC:
    FilePath=input("File Path of IRC *.log:")
    organizeIRC(FilePath)
    sys.exit()

#If not options given, display help menu
args_list = vars(args)
default_args = vars(parser.parse_args([]))
if args_list == default_args:
    parser.print_help()




