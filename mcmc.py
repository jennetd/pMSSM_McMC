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
parameter_ranges ={}

parameter_ranges["tb"] = (1,60)
parameter_ranges["Mh3"] = (100,25000)
parameter_ranges["mu"] = (80,25000) # can be <0
parameter_ranges["M1"] = (1,25000) # can be <0
parameter_ranges["M2"] = (70,25000) # can be <0
parameter_ranges["M3"] = (200,50000)
for parameter in ["Ml1","Mr1","Ml3","Mr3"]:
    parameter_ranges[parameter] = (90,25000)
for parameter in ["Mq1","Mu1","Md1"]:
    parameter_ranges[parameter] = (200,50000)
for parameter in ["Mq3","Mu3","Md3"]:
    parameter_ranges[parameter] = (100,50000)
for parameter in ["Ab","Al"]: # can be <0
    parameter_ranges[parameter] = (1,7000)
parameter_ranges["At"] = (1,7000) # can be <0

# Malte's ranges
#for parameter in ["mu","M1","M2"]:
#    parameter_ranges[parameter] = (0,4000)
#for parameter in ["M3","Mq1","Mq3","Mu1","Mu3","Md1","Md3"]:
#    parameter_ranges[parameter] = (0,10000)
#for parameter in ["Ml1","Mr1","Ml3","Mr3","Mh3"]:
#    parameter_ranges[parameter] = (0,4000)
#for parameter in ["At","Ab","Al"]:
#    parameter_ranges[parameter] = (0,7000)
#parameter_ranges["tb"] = (2,60)

# signs
parameter_signs = {}
for parameter in ["tb","Mh3","M3","Ml1","Mr1","Ml3","Mr3","Mq1","Mu1","Md1","Mq3","Mu3","Md3"]:
    parameter_signs[parameter] = 1
for parameter in ["mu","M1","M2","Al","Ab","At"]:
    parameter_signs[parameter] = np.random.randint(2, size=1)[0]*2 - 1
#    parameter_signs[parameter] = 1
    
# width coefficient of the gaussian for the mcmc step. 
# This coefficient is multiplied by the parameter range to give the width of the gaussian.
width_coefficient = .01
base = np.e

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

#containers for the tree branches. Numpy arrays are used as an interface between python types and the root branches
tree_branches = {}
tree_branches["slha_file"]={"container":TString(),"dtype":"TString"}
tree_branches["likelihood"]= {"container":np.zeros(1,dtype = float),"dtype":"D"}
tree_branches["iteration_index"] = {"container":np.zeros(1,dtype = int),"dtype":"I"}
tree_branches["accepted_index"] = {"container":np.zeros(1,dtype = int),"dtype":"I"}
tree_branches["chain_index"] = {"container":np.zeros(1,dtype = int),"dtype":"I"}

tree_branches["mtop"] = {"container":np.zeros(1,dtype = float),"dtype":"D"}
tree_branches["mbottom"] = {"container":np.zeros(1,dtype = float),"dtype":"D"}
tree_branches["alpha_s"] = {"container":np.zeros(1,dtype = float),"dtype":"D"}
tree_branches["mhiggs"] = {"container":np.zeros(1,dtype = float),"dtype":"D"}
tree_branches["mW"] = {"container":np.zeros(1,dtype = float),"dtype":"D"}
for param in parameter_ranges.keys():
    tree_branches[param] = {"container":np.zeros(1,dtype=float),"dtype":"D"}
for param in ["mu","M1","M2","Al","Ab","At"]:
    tree_branches[param+"_sign"] = {"container":np.zeros(1,dtype=int),"dtype":"I"}

tree_branches["superiso_chi2_stdout"]={"container":TString(),"dtype":"TString"}
tree_branches["superiso_stdout"]={"container":TString(),"dtype":"TString"}
tree_branches["Delta0_B_to_K0star_gamma"] = {"container":np.zeros(1,dtype=float),"dtype":"D"}
tree_branches["BR_B0_K0star_gamma"] = {"container":np.zeros(1,dtype=float),"dtype":"D"}
tree_branches["BR_Bs_to_mu_mu"] = {"container":np.zeros(1,dtype=float),"dtype":"D"}
tree_branches["BR_Bd_to_mu_mu"] = {"container":np.zeros(1,dtype=float),"dtype":"D"}
#tree_branches["BR_b_to_s_mu_mu"] = {"container":np.zeros(1,dtype=float),"dtype":"D"}
#tree_branches["BR_b_to_s_e_e"] = {"container":np.zeros(1,dtype=float),"dtype":"D"}
tree_branches["BR_b_to_s_gamma"] = {"container":np.zeros(1,dtype=float),"dtype":"D"}
tree_branches["siso_chi2"]={"container":np.zeros(1,dtype=float),"dtype":"D"}
tree_branches["siso_chi2_ndf"]={"container":np.zeros(1,dtype=int),"dtype":"I"}

tree_branches["hb_stdout"]={"container":TString(),"dtype":"TString"}
tree_branches["hb_ch1"]={"container":np.zeros(1,dtype=int),"dtype":"I"}
tree_branches["hb_ch1_exclusion"]={"container":np.zeros(1,dtype=int),"dtype":"I"}
tree_branches["hb_ch2"]={"container":np.zeros(1,dtype=int),"dtype":"I"}
tree_branches["hb_ch2_exclusion"]={"container":np.zeros(1,dtype=int),"dtype":"I"}
tree_branches["hb_ch3"]={"container":np.zeros(1,dtype=int),"dtype":"I"}
tree_branches["hb_ch3_exclusion"]={"container":np.zeros(1,dtype=int),"dtype":"I"}
tree_branches["hb_ch4"]={"container":np.zeros(1,dtype=int),"dtype":"I"}
tree_branches["hb_ch4_exclusion"]={"container":np.zeros(1,dtype=int),"dtype":"I"}

tree_branches["hb_chi2_stdout"]={"container":TString(),"dtype":"TString"}
tree_branches["llh_CMS8"]={"container":np.zeros(1,dtype=float),"dtype":"D"}
tree_branches["llh_CMS13"]={"container":np.zeros(1,dtype=float),"dtype":"D"}
tree_branches["llh_ATLAS20"]={"container":np.zeros(1,dtype=float),"dtype":"D"}

tree_branches["hb_slha_stdout"] = {"container":TString(),"dtype":"TString"}

tree_branches["hs_stdout"]={"container":TString(),"dtype":"TString"}
tree_branches["hs_chi2"]={"container":np.zeros(1,dtype=float),"dtype":"D"}
tree_branches["hs_chi2_ndf"]={"container":np.zeros(1,dtype=int),"dtype":"I"}

tree_branches["gm2calc_stdout"]={"container":TString(),"dtype":"TString"}
tree_branches["Delta_a_mu_x1E11"]={"container":np.zeros(1,dtype = float),"dtype":"D"}

#tree_branches["omegah2"] = {"container":np.zeros(1,dtype=float),"dtype":"D"}
#tree_branches["micromegas_stdout"]={"container":TString(),"dtype":"TString"}

def generate_point(input_point = {}):
    """
    @param input_point: if using mode 2, give a dictionary keying the pMSSM parameter values of the previous point in the chain
    """
    # mode 1: generate point from flat prior
    from_scratch = len(input_point) == 0
    output_point = {}

    if from_scratch:
        for parameter,parameterrange in parameter_ranges.items():
            parametervalue = parameter_signs[parameter]*random.uniform(parameterrange[0],parameterrange[1])
            output_point[parameter] = parametervalue
        output_point["mtop"] = 173.1
        output_point["mbottom"] = 4.18
        output_point["alpha_s"] = 0.1181

    else: #mode 2: generate point from previous point

        parameter_ranges["At"] = (1,sqrt(input_point["Mq3"]*input_point["Mu3"]))
    
        for parameter,parameterrange in parameter_ranges.items():
            in_range = False

            if True: # constant width, log mean
                width = width_coefficient*log(parameterrange[1]-parameterrange[0],base)
                #log(width_coefficient*(parameterrange[1]-parameterrange[0]),base)
                mean = log(abs(input_point[parameter]),base)
            
                while not in_range:
                    parametervalue = pow(base,random.gauss(mean,width))                
                    in_range = parametervalue > parameterrange[0] and parametervalue < parameterrange[1]
                
            if False: # constant width, lin mean
                width = width_coefficient*(parameterrange[1]-parameterrange[0])
                mean = input_point[parameter]
                
                while not in_range:
                    parametervalue = random.gauss(mean,width)
                    in_range = parametervalue > parameterrange[0] and parametervalue < parameterrange[1]

            output_point[parameter] = np.sign(input_point[parameter])*parametervalue
            
        for param in ["mu","M1","M2","Al","Ab","At"]:
            output_point[param+"_sign"] = parameter_signs[param]

        output_point["mtop"]= random.gauss(input_point["mtop"],0.9*width_coefficient)
        output_point["mbottom"]= random.gauss(input_point["mbottom"],((0.18+0.04)/2)*width_coefficient)
        output_point["alpha_s"] = random.gauss(input_point["alpha_s"],0.0011*width_coefficient)

    output_point["scale"] = sqrt(output_point["Mq3"]*output_point["Mu3"])

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

def get_observables(slhapath):
    """
    get the observables from the slha file
    """
    returndict = {}
    
    d = pyslha.read(slhapath)

    if "MASS" in d.blocks:
        returndict["mtop"] = {"value":float(d.blocks["MASS"][6])}
        returndict["mhiggs"] = {"value":float(d.blocks["MASS"][25])}
        returndict["mW"] = {"value":float(d.blocks["MASS"][24])}
    else:
        print("Error: no block MASS in SLHA file")

    if "SMINPUTS" in d.blocks:
        returndict["mbottom"] = {"value":float(d.blocks["SMINPUTS"][5])}
        returndict["alpha_s"] = {"value":float(d.blocks["SMINPUTS"][3])}
    else:
        print("Error: no block SMINPUTS in SLHA file")

    if "HIGGSSIGNALSRESULTS" in d.blocks:
        try:
            returndict["hs_chi2"] = {"value":float(d.blocks["HIGGSSIGNALSRESULTS"][17]),"special_case":""}
            returndict["hs_chi2_ndf"] = {"value":float(d.blocks["HIGGSSIGNALSRESULTS"][8]),"special_case":""}
        except:
            print "something went wrong with HiggsSignals"
            print "rejecting candidate point"
            print(d)
            return -1
                                            
    else:
        print("Error: no block HIGGSSIGNALSRESULTS in SLHA file")
        
    if "HIGGSBOUNDSRESULTS" in d.blocks:
            returndict["hb_ch1"] = {"value":int(d.blocks["HIGGSBOUNDSRESULTS"][0,1]),"special_case":""}
            returndict["hb_ch1_exclusion"] = {"value":int(d.blocks["HIGGSBOUNDSRESULTS"][0,2]),"special_case":""}
            returndict["hb_ch2"] = {"value":int(d.blocks["HIGGSBOUNDSRESULTS"][1,1]),"special_case":""}
            returndict["hb_ch2_exclusion"] = {"value":int(d.blocks["HIGGSBOUNDSRESULTS"][1,2]),"special_case":""}
            returndict["hb_ch3"] = {"value":int(d.blocks["HIGGSBOUNDSRESULTS"][2,1]),"special_case":""}
            returndict["hb_ch3_exclusion"] = {"value":int(d.blocks["HIGGSBOUNDSRESULTS"][2,2]),"special_case":""}
            returndict["hb_ch4"] = {"value":int(d.blocks["HIGGSBOUNDSRESULTS"][3,1]),"special_case":""}
            returndict["hb_ch4_exclusion"] = {"value":int(d.blocks["HIGGSBOUNDSRESULTS"][3,2]),"special_case":""}
    else:
        print("Error: no block HIGGSBOUNDSRESULTS in SLHA file")
                
    return returndict

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
    mode = arguments.mode
    chainix = arguments.chain_index
    tend = arguments.npoints
    start = 0
    move_every = arguments.save_interval
    outdir = arguments.output
    inpath = arguments.input
    devnull = '> /dev/null'

    lastaccepted = {}
    notaccepted = {}
    #select mode, use argparse for this
    move = 1
    #run the thing
    #mode 1: start from new point
    if mode == "new":
        outname = "pMSSM_MCMC_"+str(chainix)+"_"+str(start)+"to"+str(min(start+move_every,start+tend))+".root"
        outroot = TFile(outname,"recreate")
        outtree = TTree("mcmc","mcmc")
        setup_tree(outtree)
        tree_branches["chain_index"][0] = chainix
        finite_lh = False
        while not finite_lh:
            utils.clean()
            spnerr = False
            n_spnerr = 0
            while not spnerr:#find a viable point                                                        
                n_spnerr += 1
                utils.clean()
                candidate = generate_point(lastaccepted)#generate a point from the last point                       
                spnin = utils.write_spheno_input(candidate)#write the input for spheno                            
                spnerr = run_spheno(spnin,devnull) #run spheno, check if viable point                             

            print("SPheno errs:",n_spnerr)
            if not run_feynhiggs(devnull):#run feynhiggs, replace higgs sector                                
                continue

            gm2_obs = run_gm2calc(slhapath="SPheno.spc")
            hs_obs = run_higgssignals(slhapath="SPheno.spc")
            hb_obs = run_higgsbounds(slhapath="Spheno.spc")
            hb_obs = run_higgsbounds_chi2(slhapath="SPheno.spc")

#            os.system("cp SPheno.spc mmgsin.slha")
#            mmgs_obs = run_micromegas(slhapath="mmgsin.slha")#micromegas seems to consume the input file?!?!
#            os.system("mv mmgsin.slha SPheno.spc")

            print "getting the stuff from the slha file"
            observables = get_observables(slhapath = "SPheno.spc") #get observables for the likelihood
            siso_obs = run_superiso("SPheno.spc")
            if siso_obs == -1:
                continue
            for obs in siso_obs:
                observables[obs] = siso_obs[obs]
            siso_chi2_obs = run_superiso_chi2("SPheno.spc")
            for obs in siso_chi2_obs:
                observables[obs] = siso_chi2_obs[obs]
            for obs in gm2_obs:
                observables[obs] = gm2_obs[obs]
            for obs in hb_obs:
                observables[obs] = hb_obs[obs]
#            for obs in mmgs_obs:
#                observables[obs] = mmgs_obs[obs]

            _l = likelihood.get_likelihood(observables)#get likelihood
            finite_lh = _l != 0
        lastaccepted["likelihood"]=_l
        lastaccepted["iteration_index"] = 1
        lastaccepted["accepted_index"] = 1
        lastaccepted["chain_index"] = chainix
        
        # write point to root, start loop
        for obs in observables.keys():
            lastaccepted[obs] = observables[obs]["value"]

#        print(type(observables["Delta_a_mu_x1E11"]["value"]))
        lastaccepted = prepare_fill(lastaccepted,outtree)#add the rest of the point info, fill the tree branches
        outtree.Fill()
        # save slha
        outfile = args.output+"/pMSSM_MCMC_"+str(lastaccepted["accepted_index"])+".shla"
        os.system("cp SPheno.spc "+outfile)        

    #mode 2: continue from previous point/root file?
    elif mode == "resume":
        lastaccepted = utils.get_point_from_rootfile(inpath,chainix)
    start = lastaccepted["iteration_index"]+1
    if mode == "resume":
        outname = "pMSSM_MCMC_"+str(chainix)+"_"+str(start)+"to"+str(min(start+move_every,start+tend))+".root"
        print "Creating file "+outname
        outroot = TFile(outname,"recreate")
        outtree = TTree("mcmc","mcmc")
        setup_tree(outtree)
        tree_branches["chain_index"][0] = chainix

    #run
    print "reached run loop"
    for iter_ix in range(start,start+tend+1):
        print iter_ix
        if move == move_every-1 and iter_ix < start+tend-1:
            outtree.BuildIndex("chain_index","iteration_index")
            outtree.Write()
            outroot.Close()
            print "Made "+str(move_every)+" iterations, moving "+outname+" to storage"
            os.system(" ".join(["mv",outname,outdir]))
            outname = "pMSSM_MCMC_"+str(chainix)+"_"+str(iter_ix)+"to"+str(min(iter_ix+move_every,start+tend+1))+".root"
            print "Creating file "+outname
            outroot = TFile(outname,"recreate")
            outtree = TTree("mcmc","mcmc")
            setup_tree(outtree)
            move = -1
        finite_lh = False
        while not finite_lh:
            utils.clean()
            spnerr = False
            n_spnerr = 0
            while not spnerr:#find a viable point
                n_spnerr += 1
                utils.clean()
                candidate = generate_point(lastaccepted)#generate a point from the last point
                spnin = utils.write_spheno_input(candidate)#write the input for spheno
                spnerr = run_spheno(spnin,devnull) #run spheno, check if viable point

            print("SPheno errs:",n_spnerr)
            if not run_feynhiggs(devnull):#run feynhiggs, replace higgs sector
                continue

            gm2_obs = run_gm2calc(slhapath="SPheno.spc")
            hs_obs = run_higgssignals(slhapath="SPheno.spc")
            hb_obs = run_higgsbounds(slhapath="SPheno.spc")
            hb_obs = run_higgsbounds_chi2(slhapath="SPheno.spc")

#            os.system("cp SPheno.spc mmgsin.slha")
#            mmgs_obs = run_micromegas(slhapath="mmgsin.slha")
#            os.system("mv mmgsin.slha SPheno.spc")

            observables = get_observables(slhapath = "SPheno.spc") #get observables for the likelihood
            siso_obs = run_superiso("SPheno.spc")
            if siso_obs == -1:
                continue
            for obs in siso_obs:
                observables[obs] = siso_obs[obs]
            siso_chi2_obs = run_superiso_chi2("SPheno.spc")
            for obs in siso_chi2_obs:
                observables[obs] = siso_chi2_obs[obs]
            for obs in gm2_obs:
                observables[obs] = gm2_obs[obs]
            for obs in hb_obs:
                observables[obs] = hb_obs[obs]
                
#            for obs in mmgs_obs:
#                observables[obs] = mmgs_obs[obs]

            _l = likelihood.make_decision(observables,lastaccepted["likelihood"])#get likelihood
            finite_lh = _l != 0

        if _l<0: # point is not accepted
            move +=1
            if iter_ix == start+tend:
                print "Made all "+str(tend)+" iterations, moving "+outname+" to storage"
                outtree.BuildIndex("chain_index","iteration_index")
                outtree.Write()
                outroot.Close()
                os.system(" ".join(["mv",outname,outdir]))
            
            notaccepted["likelihood"] = _l
            notaccepted["iteration_index"] = iter_ix
            notaccepted["accepted_index"] = -1
            notaccepted["chain_index"] = chainix

            for obs in observables.keys():
                notaccepted[obs] = observables[obs]["value"]
                
            notaccepted = prepare_fill(notaccepted,outtree)
            outtree.Fill() # save it even though not accepted
            continue

        lastaccepted["likelihood"]=_l
        lastaccepted["iteration_index"] = iter_ix
        lastaccepted["accepted_index"] =lastaccepted["accepted_index"]+1
        lastaccepted["chain_index"] = chainix

        for obs in observables.keys():
            lastaccepted[obs] = observables[obs]["value"]

        lastaccepted = prepare_fill(lastaccepted,outtree)#add the rest of the point info, fill the tree branches
        outtree.Fill()
        # save slha
        outfile = args.output+"/pMSSM_MCMC_"+str(lastaccepted["accepted_index"])+".shla"
        os.system("cp SPheno.spc "+outfile)

        if iter_ix == start+tend:
            print "Made all "+str(tend)+" iterations, moving "+outname+" to storage"
            outtree.BuildIndex("chain_index","iteration_index")
            outtree.Write()
            outroot.Close()
            os.system(" ".join(["mv",outname,outdir]))
        move+=1

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-m","--mode",choices=["new","resume"],help="Specify run mode, either new or resume",required=True)
    parser.add_argument("-i","--input",help="If runmode resume, specify input root file",default = None)
    parser.add_argument("-o","--output",help="Specify an output directory",required=True)
    parser.add_argument("-n","--npoints",help = "How many points to run the MCMC for",default=1000,type=int)
    parser.add_argument("-c","--chain_index",help = "chain index for the chain",type = int,default = 1)
    parser.add_argument("-s","--save_interval",default = 300,help = "How many points to generate before starting a new root file",type=int)
    args=parser.parse_args()
    if not os.path.exists(args.output):
        print "Output directory "+args.output+" does not exist. Please specify an existing output directory"
        exit()
    if args.mode =="resume":
        if args.input == None:
            parser.error("Need to specify an input file in order to resume")
        if args.chain_index!=None:
            print "Ignoring specified chain index, taking from input file"            
    if args.mode == "new":
        if args.input!=None:
            print "Ignoring input file, as we are starting a new chain"
    run(args)

    

