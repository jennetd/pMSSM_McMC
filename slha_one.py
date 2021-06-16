#this file should contain the mcmc decision function and the point generator
import os
import random
from math import sqrt,log
import argparse
import os
from ROOT import *
import numpy as np
import utils
import likelihood
import subprocess
from datetime import datetime
import pyslha

# set up the parameter ranges 
# positive definite only: signs dealt with below
parameter_values ={}

parameter_values["tb"] = 59.6803
parameter_values["Mh3"] = 10512.8
parameter_values["mu"] = 2733.62
parameter_values["M1"] = 5612.28
parameter_values["M2"] = -5704.96
parameter_values["M3"] = 3819.43
parameter_values["Mr1"] = 5745.33
parameter_values["Ml1"] = 18991.6
parameter_values["Mr3"] = 19199.8
parameter_values["Ml3"] = 1506.45
parameter_values["Mq1"] = 1750.39
parameter_values["Mq3"] = 5885.49
parameter_values["Mu1"] = 2717.19
parameter_values["Mu3"] = 2566.18
parameter_values["Md1"] = 23291.2
parameter_values["Md3"] = 5523.24
parameter_values["Al"] = 725.269
parameter_values["Ab"] = 2402.21
parameter_values["At"] = 3420.42

# paths and executables
homedir = "/pMSSM_McMC"
packagedir = homedir+"/packages/"
spnexe = packagedir+"SPheno-4.0.4/bin/SPheno"
fhexe = packagedir+"FeynHiggs-2.16.1/x86_64-Linux/bin/FeynHiggs"
sisoexe = packagedir+"superiso_v4.0/slha.x"
#sisochi2exe = packagedir+"superiso_v4.0/slha_chi2.x" #use all of the non-controversial low-energy results in superiso chi2 calculation. takes approximately 20s/call
sisochi2exe = packagedir+"superiso_v4.0/slha_chi2_reduced.x"#use only branching ratios in superiso chi2. takes approximately 8s/call
hbexe = "./packages/higgsbounds/build/HiggsBounds"
hbchi2exe = "./packages/higgsbounds/build/example_programs/HBwithLHClikelihood_SLHA"
hbslhaexe="./packages/higgsbounds/build/example_programs/HBSLHAinputblocksfromFH"
hsexe = "./packages/higgssignals/build/HiggsSignals"
mmgsexe = packagedir+"micromegas_5.2.4/MSSM/main"
gm2exe = packagedir+"GM2Calc-1.7.3/build/bin/gm2calc.x"

def generate_point():
    output_point = {}

    for parameter in parameter_values.keys():
        output_point[parameter] = parameter_values[parameter]

    output_point["mtop"] = 173.1
    output_point["mbottom"] = 4.18
    output_point["alpha_s"] = 0.1181

    output_point["scale"] = sqrt(output_point["Mq3"]*output_point["Mu3"])

    print(output_point)

    return output_point

# SPheno
def run_spheno(inpath,devnull):
    cmd = " ".join([spnexe,inpath,devnull])
    os.system(cmd)
    error = open("Messages.out","r").read()
#    print(error)
    return len(error) == 0

# FeynHiggs
def run_feynhiggs(devnull):
    """
    Replace the SPheno Higgs sector with the FeynHiggs one
    fhin = FeynHiggs output file
    spnin = SPheno output file
    """

    fhin = "SPheno.spc.fh-001"#again, writing a lot of files to disk
    spnin = "SPheno.spc"
    cmd = " ".join([fhexe,spnin,devnull])
    
    t1 = datetime.now()
    os.system(cmd)# run FeynHiggs                                           
    print("FeynHiggs",str(datetime.now()-t1))

    if not os.path.exists(fhin):#check to see whether FeynHiggs produced an output
        print "could not find Feynhiggs output, skipping point"
        return False

    # Replace Higgs parameters from SPheno with FeynHiggs
    dspnin = pyslha.read(spnin)
    dfhin = pyslha.read(fhin)

    dspnin.blocks["DMASS"] = dfhin.blocks["DMASS"]
    dspnin.blocks["ALPHA"] = dfhin.blocks["ALPHA"]

    # If the Higgs mass is too terrible
#    if dspnin.blocks["MASS"][25] < 123.6 or dspnin.blocks["MASS"][25] > 126.6:
#        print "awful Higgs mass, skipping point"
#        return False

    particles_to_replace = [24,25,35,36,37]
    for particle in particles_to_replace:
        dspnin.blocks["MASS"][particle] = dfhin.blocks["MASS"][particle]
        dspnin.decays[particle] = dfhin.decays[particle]

    pyslha.writeSLHAFile("SPheno.spc",dspnin)

    return True

# helper function for reading values from superiso output
def read_superiso_out(search_str,siso_out):
    variable_str = siso_out[siso_out.find(search_str)+len(search_str):]
    variable = float(variable_str[:variable_str.find("\n")].strip())
    return variable

# superiso
def run_superiso(slhapath):

    t1 = datetime.now()
    siso_call = subprocess.Popen([sisoexe,str(slhapath)], stdout=subprocess.PIPE)
    siso_out = siso_call.stdout.read()
    print("superiso",str(datetime.now()-t1))

    if len(siso_out)<10:
        print "something went wrong with siso call!"
#        print siso_out
    returndict = {"superiso_stdout":{"value":siso_out,"special_case":""}}

    # get the individual observables from stdout
    try:
        returndict["Delta0_B_to_K0star_gamma"] = {"value":read_superiso_out("delta0(B->K* gamma)",siso_out)}
        returndict["BR_B0_K0star_gamma"] = {"value":read_superiso_out("BR(B0->K* gamma)",siso_out)}
        returndict["BR_Bs_to_mu_mu"] = {"value":read_superiso_out("BR(Bs->mu mu)",siso_out)}
        returndict["BR_Bd_to_mu_mu"] = {"value":read_superiso_out("BR(Bd->mu mu)",siso_out)}
#        returndict["BR_b_to_s_mu_mu"] = {"value":read_superiso_out()}
#        returndict["BR_b_to_s_e_e"] = {"value":read_superiso_out()}
        returndict["BR_b_to_s_gamma"] = {"value":read_superiso_out("BR(b->s gamma)",siso_out)}

    except:
        print "something went wrong with siso call, printing output"
        print siso_out
        print "rejecting candidate point"
        return -1

    return returndict

def run_superiso_chi2(slhapath):

    t1 = datetime.now()
    siso_chi2_call = subprocess.Popen([sisochi2exe,str(slhapath)], stdout=subprocess.PIPE)
    siso_chi2_out = siso_chi2_call.stdout.read()
    print("superiso chi2",str(datetime.now()-t1))

    if len(siso_chi2_out)<10:
        print "something went wrong with siso chi2 call!"
        print siso_chi2_out
    #special_case key in dictionary tells the likelihood function that this observable should not be handled in a standard way
    returndict = {"superiso_chi2_stdout":{"value":siso_chi2_out,"special_case":""}}
    chi2 = siso_chi2_out[siso_chi2_out.find("chi2"):]
    chi2 = float(chi2[chi2.find("=")+1:chi2.find("\n")].strip())
    returndict["siso_chi2"]={"value":chi2,"special_case":""}
    ndf = siso_chi2_out[siso_chi2_out.find("n_obs"):]
    ndf = int(ndf[ndf.find("=")+1:ndf.find("\n")].strip())
    returndict["siso_chi2_ndf"]={"value":ndf,"special_case":""}
#    print siso_chi2_out
    return returndict

# HiggsSignals
def run_higgssignals(slhapath):

    t1 = datetime.now()
    hb_slha_call = subprocess.Popen(hbslhaexe + " " +slhapath, stdout=subprocess.PIPE, shell=True)
    hb_slha_out = hb_slha_call.stdout.read()
    os.system("cp "+slhapath+".fh "+slhapath)
    
    hs_call = subprocess.Popen(hsexe+" latestresults 1 SLHA 3 1 "+slhapath, stdout=subprocess.PIPE, shell=True)
    hs_out = hs_call.stdout.read()
    print("HiggsSignals",str(datetime.now()-t1))

    returndict = {"hb_slha_stdout":{"value":hb_slha_out,"special_case":""},
                  "hs_stdout":{"value":hs_out,"special_case":""}}
    return returndict
    
# HiggsBounds
def run_higgsbounds(slhapath):

    t1 = datetime.now()
    hb_call = subprocess.Popen(hbexe+" LandH SLHA 3 1 "+slhapath, stdout=subprocess.PIPE,shell=True)
    hb_out = hb_call.stdout.read()
    print(hb_out)
    print("HiggsBounds",str(datetime.now()-t1))

    returndict = {"hb_stdout":{"value":hb_out,"special_case":""}}
    return returndict

def run_higgsbounds_chi2(slhapath):

    os.system("cp "+slhapath+" "+slhapath+".1")
    t1 = datetime.now()
    hb_call = subprocess.Popen(hbchi2exe+" 1 "+slhapath, stdout=subprocess.PIPE,shell=True)
    hb_out = hb_call.stdout.read()
    print("HiggsBounds chi2",str(datetime.now()-t1))

    returndict = {"hb_chi2_stdout":{"value":hb_out,"special_case":""}}
    
    os.system("cp "+slhapath+".1 "+slhapath)

    with open("Mh125_HBwithLHClikelihood.dat", "r") as hb_outfile:
         content = hb_outfile.read().split()

    try:
        returndict["llh_CMS8"]        = {"value":float(content[1]),"special_case":""}
        returndict["llh_CMS13"]       = {"value":float(content[3]),"special_case":""}
        returndict["llh_ATLAS20"]     = {"value":float(content[5]),"special_case":""}
    except:
        print "something went wrong with higgsbounds chi2 call, printing output"
        print hb_out
        print "rejecting candidate point"
        return -1

#    print(returndict["llh_CMS8"]["value"],returndict["llh_CMS13"]["value"],returndict["llh_ATLAS20"]["value"])
    
    return returndict

# MicroMegas
def run_micromegas(slhapath):

    t1 = datetime.now()
    print "calling micromegas"
    micromegas_call = subprocess.Popen(mmgsexe+" "+str(slhapath), stdout=subprocess.PIPE,shell=True)
    print "processing micromegas output"
    micromegas_out = micromegas_call.stdout.read()
    print "I got the output! yay! This is it:"
    print("Micromegas",str(datetime.now()-t1))

    print micromegas_out

    #if any of these quantities are used to binarily reject candidate points, micromegas should be run first and terminate if, and as soon as, a rejection criterion is fulfilled
    returndict = {"micromegas_stdout":{"value":micromegas_out,"special_case":""}}
#    ztoinv_excluded = micromegas_out.find("Excluded by Z->invisible") != -1
#    returndict["ztoinv_excluded"] = {"value":ztoinv_excluded,"special_case":""}
#    lep_excluded = micromegas_out.find("Excluded by LEP  by e+,e- -> DM q qbar. Cross section =")!=-1
#    returndict["lep_excluded"] = {"value":lep_excluded,"special_case":""}
#    masslim = micromegas_out.find("MassLimits OK")!=-1
#    returndict["masslim"] = {"value":masslim,"special_case":""}#if true, mass limits are not ok
#    omegah2 = micromegas_out[micromegas_out.find("Omega=")+len("Omega="):]
#    omegah2 = float(omegah2.split("\n")[0].strip())
#    returndict["omegah2"]={"value":omegah2,"special_case":""}
#    omegaxf = float(micromegas_out[micromegas_out.find("Xf=")+len("Xf="):micromegas_out.find("Omega=")].strip())
#    returndict["omegaxf"] = {"value":omegaxf,"special_case":""}
    return returndict

# GM2Calc                                                                                  
def run_gm2calc(slhapath):

    t1 = datetime.now()

    d = pyslha.read(slhapath)
    d.blocks["GM2CalcConfig"] = pyslha.Block("GM2CalcConfig")
    d.blocks["GM2CalcConfig"][0] = 4 # output format
    d.blocks["GM2CalcConfig"][1] = 2 # loop order (0, 1 or 2)
    d.blocks["GM2CalcConfig"][2] = 1 # disable/enable tan(beta) resummation (0 or 1)
    d.blocks["GM2CalcConfig"][3] = 0 # force output (0 or 1)
    d.blocks["GM2CalcConfig"][4] = 0 # verbose output (0 or 1)
    d.blocks["GM2CalcConfig"][5] = 1 # calculate uncertainty

    pyslha.writeSLHAFile(slhapath,d)
    
    gm2_call = subprocess.Popen(gm2exe+" --slha-input-file="+str(slhapath), stdout=subprocess.PIPE,shell=True)
    gm2_out = gm2_call.stdout.read()
    print("GM2Calc",str(datetime.now()-t1))

    returndict = {"gm2calc_stdout":{"value":gm2_out,"special_case":""}}

    try:
        blocks = gm2_out.split("Block")
        gm2_str = blocks[-1].split()[2]
        gm2_unc_str = blocks[-1].split()[6]
        returndict["Delta_a_mu_x1E11"] = {"value":float(gm2_str)*pow(10,11),"uncertainty":float(gm2_unc_str)*pow(10,11)}
    
    except:
        print "something went wrong with siso call, printing output"
        print gm2_out
        print "rejecting candidate point"
        return -1
    
    return returndict
    
def setup_tree(outtree):
    for branch in tree_branches.keys():
        if tree_branches[branch]["dtype"] == "TString":#the container in this case is just a placeholder, since a new TString is created for each new tree entry. I am not sure if this can be done differently
            outtree.Branch(branch,tree_branches[branch]["container"])
        else:
            outtree.Branch(branch,tree_branches[branch]["container"],branch+"/"+tree_branches[branch]["dtype"])
        

def prepare_fill(point,outtree):
    point_info = {}
    tvals = {}
    for key,val in point.items():
        if type(val) != str:
            point_info[key] = val
        else:#strings are handled differently
            tvals[key] = TString(val)
            outtree.SetBranchAddress(key,tvals[key])
            #            point_info[key]=val

    d = pyslha.read("SPheno.spc")
    if "EXTPAR" in d.blocks:
        point_info["M1"] = d.blocks["EXTPAR"][1]
        point_info["M2"] = d.blocks["EXTPAR"][2]
        point_info["M3"] = d.blocks["EXTPAR"][3]
        point_info["At"] = d.blocks["EXTPAR"][11]
        point_info["Ab"] = d.blocks["EXTPAR"][12]
        point_info["Al"] = d.blocks["EXTPAR"][13]
        point_info["mu"] = d.blocks["EXTPAR"][23]
        point_info["tb"] = d.blocks["EXTPAR"][25]
        point_info["Mh3"] = d.blocks["EXTPAR"][26]
        point_info["Ml1"] = d.blocks["EXTPAR"][31]
        point_info["Ml2"] = d.blocks["EXTPAR"][32]
        point_info["Ml3"] = d.blocks["EXTPAR"][33]
        point_info["Mr1"] = d.blocks["EXTPAR"][34]
        point_info["Mr2"] = d.blocks["EXTPAR"][35]
        point_info["Mr3"] = d.blocks["EXTPAR"][36]
        point_info["Mq1"] = d.blocks["EXTPAR"][41]
        point_info["Mq2"] = d.blocks["EXTPAR"][42]
        point_info["Mq3"] = d.blocks["EXTPAR"][43]
        point_info["Mu1"] = d.blocks["EXTPAR"][44]
        point_info["Mu2"] = d.blocks["EXTPAR"][45]
        point_info["Mu3"] = d.blocks["EXTPAR"][46]
        point_info["Md1"] = d.blocks["EXTPAR"][47]
        point_info["Md2"] = d.blocks["EXTPAR"][48]
        point_info["Md3"] = d.blocks["EXTPAR"][49]

        for parameter in ["mu","M1","M2","Al","Ab","At"]:
            point_info[parameter+"_sign"] = parameter_signs[parameter]

    with open("SPheno.spc","r") as spnin:
        slhacont = spnin.read()
        slhafile = slhacont
    slha_file = TString(slhafile)
    outtree.SetBranchAddress("slha_file",slha_file)

    for key in tree_branches.keys():#exclude strings
        if tree_branches[key]["dtype"]=="TString":continue
        tree_branches[key]["container"][0]=point_info[key]

    return point_info

def run(arguments):
    outdir = arguments.output
    devnull = '> /dev/null'

    utils.clean()

    candidate = generate_point()#generate a point from the last point                       
    spnin = utils.write_spheno_input(candidate)#write the input for spheno                            
    spnerr = run_spheno(spnin,devnull) #run spheno, check if viable point                             
    fherr = run_feynhiggs(devnull)

    gm2_obs = run_gm2calc(slhapath="SPheno.spc")
    hs_obs = run_higgssignals(slhapath="SPheno.spc")
    hb_obs = run_higgsbounds(slhapath="SPheno.spc")
    hb_obs = run_higgsbounds_chi2(slhapath="SPheno.spc")

    # save slha
    outfile = args.output+"/pMSSM_MCMC_one.shla"
    os.system("cp SPheno.spc "+outfile)        

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i","--input",help="If runmode resume, specify input root file",default = None)
    parser.add_argument("-o","--output",help="Specify an output directory",required=True)
    parser.add_argument("-s","--save_interval",default = 300,help = "How many points to generate before starting a new root file",type=int)
    args=parser.parse_args()
    if not os.path.exists(args.output):
        print "Output directory "+args.output+" does not exist. Please specify an existing output directory"
        exit()
    run(args)

    

