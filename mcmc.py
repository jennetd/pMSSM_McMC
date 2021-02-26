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

# signs
parameter_signs = {}
for parameter in ["tb","Mh3","M3","Ml1","Mr1","Ml3","Mr3","Mq1","Mu1","Md1","Mq3","Mu3","Md3"]:
    parameter_signs[parameter] = 1
for parameter in ["mu","M1","M2","Al","Ab","At"]:
    parameter_signs[parameter] = np.random.randint(2, size=1)[0]*2 - 1

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
hbexe = packagedir+"higgsbounds/build/HiggsBounds"
hbchi2exe = packagedir+"higgsbounds/build/example_programs/HBwithLHClikelihood_SLHA"
hsexe =packagedir+"higgssignals/build/HiggsSignals"
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

#tree_branches["hb_stdout"]={"container":TString(),"dtype":"TString"}
#tree_branches["hb_ch1"]={"container":np.zeros(1,dtype=int),"dtype":"I"}
#tree_branches["hb_ch1_exclusion"]={"container":np.zeros(1,dtype=int),"dtype":"I"}
#tree_branches["hb_ch2"]={"container":np.zeros(1,dtype=int),"dtype":"I"}
#tree_branches["hb_ch2_exclusion"]={"container":np.zeros(1,dtype=int),"dtype":"I"}
#tree_branches["hb_ch3"]={"container":np.zeros(1,dtype=int),"dtype":"I"}
#tree_branches["hb_ch3_exclusion"]={"container":np.zeros(1,dtype=int),"dtype":"I"}
#tree_branches["hb_ch4"]={"container":np.zeros(1,dtype=int),"dtype":"I"}
#tree_branches["hb_ch4_exclusion"]={"container":np.zeros(1,dtype=int),"dtype":"I"}

#tree_branches["hb_chi2_stdout"]={"container":TString(),"dtype":"TString"}
#tree_branches["llh_CMS8"]={"container":np.zeros(1,dtype=float),"dtype":"D"}
#tree_branches["llh_CMS13"]={"container":np.zeros(1,dtype=float),"dtype":"D"}
#tree_branches["llh_ATLAS20"]={"container":np.zeros(1,dtype=float),"dtype":"D"}

#tree_branches["hs_stdout"]={"container":TString(),"dtype":"TString"}
#tree_branches["hs_chi2"]={"container":np.zeros(1,dtype=float),"dtype":"D"}
#tree_branches["hs_chi2_ndf"]={"container":np.zeros(1,dtype=int),"dtype":"I"}

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

            # constant width, log mean
            width = log(width_coefficient*(parameterrange[1]-parameterrange[0]),base)
            mean = log(abs(input_point[parameter]),base)
            
            while not in_range:
                parametervalue = pow(base,random.gauss(mean,width))                
                in_range = parametervalue > parameterrange[0] and parametervalue < parameterrange[1]
                
            output_point[parameter] = np.sign(input_point[parameter])*parametervalue
            
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
    os.system(cmd)# run FeynHiggs
    if not os.path.exists(fhin):#check to see whether FeynHiggs produced an output
        print "could not find Feynhiggs output, skipping point"
        return False
    #fhin part: Get the new Higgs parameters
    masses = {"Mh0":"","MHH":"","MA0":"","MHp":"","MW":""}#dictionary to collect the FeynHiggs Higgs masses
    dmassblock= "Block"# No replacing needs to be done here, can just append the whole block at the end
    alpha = -1# for the FeynHiggs alpha values
    decaytabs = {"25":"","35":"","36":"","37":"","24":""}#dictionary to collect the FeynHiggs decay tables
    with open(fhin,"r") as feynin:
        fhtab = feynin.read()
        blocks = fhtab[:fhtab.find("DECAY")].split("BLOCK")#put all the block in a list
        decays = fhtab[fhtab.find("DECAY"):].split("DECAY")#put the decay tables in a list
        for block in blocks:# loop over the slha blocks
            blocklines = block.split("\n")
            blockname =blocklines[0].strip()
            if blockname == "MASS":#we are only interested in the MASS block
                for bline in blocklines[1:-1]: #iterate over the lines with the particle masses
                    ptc = bline[bline.find("#")+1:].strip()#particle name
                    val = bline.split("   ")[-2].strip()#particle mass
                    if ptc in masses.keys():#only interested in the Higgs bosons
                        masses[ptc] = val.lower()#write value to dictionay
            if blockname == "DMASS":#just append the block
                dmassblock+=block
            if blockname == "ALPHA":#This block only contains the one value
                alpha = blocklines[1]
                alpha = alpha[:alpha.find("#")].strip()
        for decay in decays:#loop over decay blocks for the different particles
            pdgid = decay.split("\n")[0]
            pdgid = pdgid[:pdgid.find(".")-1].strip()#this is the pdgid to which the decay block belongs
            if pdgid in decaytabs:#only interested in Higgs boson decay tables
                decaytabs[pdgid] = "DECAY"+decay#save the respective decay blocks
    #spnin part: find the obsolete Higgs-related parts in the SPheno generated SLHA file
    translate = {"h0":"Mh0","H0":"MHH","A0":"MA0","H)+":"MHp","W+":"MW"}# dictionary to translate from SPheno naming convention to FeynHiggs naming convention
    with open(spnin,"r") as tmpin:
        spntab = tmpin.read()#original file content
        newcont = spntab#copy of the original file content. This will be modified with the FeynHiggs masses and decays, then replace the original content
        blocks = spntab[:spntab.find("DECAY")].split("Block")#SPheno SLHA blocks
        decays = spntab[spntab.find("DECAY"):].split("DECAY")#SPheno decay tables
        replaceblock = "\n".join([dmassblock,"# Higgs mixing\n"])#Puts the DMASS block below the MASS block (and above the Higgs mixing part of the SLHA)
        newcont = newcont.replace("# Higgs mixing\n",replaceblock)#insert the DMASS block
        #beware: the order of replacement in this string might matter
        for block in blocks:#loop over SLHA blocks
            blocklines = block.split("\n")
            blockname =blocklines[0][:blocklines[0].find("#")].strip()
            if blockname == "MASS":#only interested in MASS block
                for bline in blocklines[1:-2]:#crop the mass block to relevant part, then loop over it
                    ptc = bline[bline.find("#")+1:].strip()#Particle name of the mass block, beware SPheno naming convention differs from FeynHiggs
                    val = bline.split("  ")[-2].strip()# particle mass
                    if ptc in translate.keys():#Only interested in Higgs boson masses
                        replaceline = bline.replace(val,masses[translate[ptc]])#replace the line containing the mass in SPheno with the FeynHiggs line
                        newcont = newcont.replace(bline,replaceline)#replace the SPheno line
            if blockname == "alpha":#replace the alpha value
                spnalpha = blocklines[1]
                spnalpha = spnalpha[:spnalpha.find("#")].strip()
                replaceline = blocklines[1].replace(spnalpha,alpha)
                newcont = newcont.replace(blocklines[1],replaceline)
        #now replace the decay tables
        for decay in decays:
            pdgid = decay.split("\n")[0]
            pdgid = pdgid[:pdgid.find(".")-1].strip()
            if pdgid in decaytabs:#only interested in Higgs boson tables
                newcont = newcont.replace("DECAY"+decay,decaytabs[pdgid])#replace the decay table for the Higgs bosons
        #write out new file
    with open(spnin,"w") as outfile:
        outfile.write(newcont)#overwrite SPheno file
#    os.system("rm "+fhin)# clean up
    return True

def get_observables(slhapath):
    """
    get the observables from the slha file
    """
    returndict = {}
    with open(slhapath,"r") as slhain:
        slhacont = slhain.read()
        blocks = slhacont.split("Block")[1:]
        #get the masses
        mblock = blocks[14]
        masses = mblock.split("\n")[2:-2]

        returndict["mtop"] = {"value":float(" ".join(masses[0].split()).split()[1])}
        returndict["mhiggs"] = {"value":float(" ".join(masses[4].split()).split()[1])}
        returndict["mW"] = {"value":float(" ".join(masses[2].split()).split()[1])}

        #get the higgs mass uncertainty
        dmblock = blocks[15]
        dmasses = dmblock.split("\n")[1:-1]

        try:
            returndict["mhiggs"]["uncertainty"]=float(" ".join(dmasses[0].split()).split()[1])
        except:
            print " ".join(dmasses[0].split()).split()
            exit()
        #get alpha_s
        smblock = blocks[5]
        sminputs = smblock.split("\n")[1:-1]
        returndict["alpha_s"] = {"value":float(" ".join(sminputs[2].split()).split()[1])}
        returndict["mbottom"] = {"value":float(" ".join(sminputs[4].split()).split()[1])}

        # get higgs bounds exclusions
#        hbblock = blocks[30]
#        hbblock = hbblock.split("(rank = 0: global result)")[-1]
#        hbinputs = hbblock.split("\n")[1:-1]
        
#        try:
#            returndict["hb_ch1"] = {"value":int(hbinputs[0].split()[2]),"special_case":""}
#            returndict["hb_ch1_exclusion"] = {"value":int(hbinputs[1].split()[2]),"special_case":""}
#            returndict["hb_ch2"] = {"value":int(hbinputs[5].split()[2]),"special_case":""}
#            returndict["hb_ch2_exclusion"] = {"value":int(hbinputs[6].split()[2]),"special_case":""}
#            returndict["hb_ch3"] = {"value":int(hbinputs[10].split()[2]),"special_case":""}
#            returndict["hb_ch3_exclusion"] = {"value":int(hbinputs[11].split()[2]),"special_case":""}
#            returndict["hb_ch4"] = {"value":int(hbinputs[15].split()[2]),"special_case":""}
#            returndict["hb_ch4_exclusion"] = {"value":int(hbinputs[16].split()[2]),"special_case":""}

#        except:
#            print("Problem getting exclusions from HiggsBounds")
#            returndict["hb_ch1"] = {"value":-1,"special_case":""}
#            returndict["hb_ch1_exclusion"] = {"value":-1,"special_case":""}
#            returndict["hb_ch2"] = {"value":-1,"special_case":""}
#            returndict["hb_ch2_exclusion"] = {"value":-1,"special_case":""}
#            returndict["hb_ch3"] = {"value":-1,"special_case":""}
#            returndict["hb_ch3_exclusion"] = {"value":-1,"special_case":""}
#            returndict["hb_ch4"] = {"value":-1,"special_case":""}
#            returndict["hb_ch4_exclusion"] = {"value":-1,"special_case":""}

        # get higgs signals chi2
#        hsblock = blocks[30]
#        hsinputs = hsblock.split("\n")[1:-1]
#        returndict["hs_chi2"] = {"value":float(" ".join(hsinputs[48].split()).split()[1]),"special_case":""}
#        returndict["hs_chi2_ndf"] = {"value":float(" ".join(hsinputs[39].split()).split()[1]),"special_case":""}
        
    return returndict

# helper function for reading values from superiso output
def read_superiso_out(search_str,siso_out):
    variable_str = siso_out[siso_out.find(search_str)+len(search_str):]
    variable = float(variable_str[:variable_str.find("\n")].strip())
    return variable

# superiso
def run_superiso(slhapath):
    siso_call = subprocess.Popen([sisoexe,str(slhapath)], stdout=subprocess.PIPE)
    siso_out = siso_call.stdout.read()
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
    siso_chi2_call = subprocess.Popen([sisochi2exe,str(slhapath)], stdout=subprocess.PIPE)
    siso_chi2_out = siso_chi2_call.stdout.read()
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
    hs_call = subprocess.Popen(hsexe+" latestresults 1 SLHA 3 1 "+slhapath, stdout=subprocess.PIPE,shell=True)
    hs_out = hs_call.stdout.read()

    returndict = {"hs_stdout":{"value":hs_out,"special_case":""}}
    return returndict
    
# HiggsBounds
def run_higgsbounds(slhapath):

    hb_call = subprocess.Popen(hbexe+" LandH SLHA 3 1 "+slhapath, stdout=subprocess.PIPE,shell=True)
    hb_out = hb_call.stdout.read()
    
    returndict = {"hb_stdout":{"value":hb_out,"special_case":""}}
    return returndict

def run_higgsbounds_chi2(slhapath):

    os.system("ln -s "+slhapath+" "+slhapath+".1")
    hb_call = subprocess.Popen(hbchi2exe+" 1 "+slhapath, stdout=subprocess.PIPE,shell=True)
    hb_out = hb_call.stdout.read()
    
    returndict = {"hb_chi2_stdout":{"value":hb_out,"special_case":""}}

    with open("Mh125_HBwithLHClikelihood.dat", "r") as hb_outfile:
         content = hb_outfile.read().split()

    try:
        returndict["llh_CMS8"]        = {"value":float(content[12]),"special_case":""}
        returndict["llh_CMS13"]       = {"value":float(content[14]),"special_case":""}
        returndict["llh_ATLAS20"]     = {"value":float(content[18]),"special_case":""}
    except:
        print "something went wrong with higgsbounds chi2 call, printing output"
        print hb_out
        print "rejecting candidate point"
        return -1

#    print(returndict["llh_CMS8"]["value"],returndict["llh_CMS13"]["value"],returndict["llh_ATLAS20"]["value"])
    
    return returndict

# MicroMegas
def run_micromegas(slhapath):

    print "calling micromegas"
    micromegas_call = subprocess.Popen(mmgsexe+" "+str(slhapath), stdout=subprocess.PIPE,shell=True)
    print "processing micromegas output"
    micromegas_out = micromegas_call.stdout.read()
    print "I got the output! yay! This is it:"

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

    with open(slhapath, "a+") as file_object:
        file_object.write("\n")
        file_object.write("Block GM2CalcConfig\n")
        file_object.write("     0     4     # output format \n")
        file_object.write("     1     2     # loop order (0, 1 or 2) \n")
        file_object.write("     2     1     # disable/enable tan(beta) resummation (0 or 1) \n")
        file_object.write("     3     0     # force output (0 or 1) \n")
        file_object.write("     4     0     # verbose output (0 or 1) \n")
        file_object.write("     5     1     # calculate uncertainty \n")

    gm2_call = subprocess.Popen(gm2exe+" --slha-input-file="+str(slhapath), stdout=subprocess.PIPE,shell=True)
    gm2_out = gm2_call.stdout.read()

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

    with open("SPheno.spc","r") as spnin:
        slhacont = spnin.read()
        slhafile = slhacont
        blocks = slhacont.split("Block")[1:]
        params = blocks[4].split("\n")[2:-1]
        point_info["M1"] = float(" ".join(params[0].split()).split()[1])
        point_info["M2"] = float(" ".join(params[1].split()).split()[1])
        point_info["M3"] = float(" ".join(params[2].split()).split()[1])
        point_info["At"] = float(" ".join(params[3].split()).split()[1])
        point_info["Ab"] = float(" ".join(params[4].split()).split()[1])
        point_info["Al"] = float(" ".join(params[5].split()).split()[1])
        point_info["mu"] = float(" ".join(params[6].split()).split()[1])
        point_info["tb"] = float(" ".join(params[7].split()).split()[1])
        point_info["Mh3"] = float(" ".join(params[8].split()).split()[1])
        point_info["Ml1"] = float(" ".join(params[9].split()).split()[1])
        point_info["Ml2"] = float(" ".join(params[10].split()).split()[1])
        point_info["Ml3"] = float(" ".join(params[11].split()).split()[1])
        point_info["Mr1"] = float(" ".join(params[12].split()).split()[1])
        point_info["Mr2"] = float(" ".join(params[13].split()).split()[1])
        point_info["Mr3"] = float(" ".join(params[14].split()).split()[1])
        point_info["Mq1"] = float(" ".join(params[15].split()).split()[1])
        point_info["Mq2"] = float(" ".join(params[16].split()).split()[1])
        point_info["Mq3"] = float(" ".join(params[17].split()).split()[1])
        point_info["Mu1"] = float(" ".join(params[18].split()).split()[1])
        point_info["Mu2"] = float(" ".join(params[19].split()).split()[1])
        point_info["Mu3"] = float(" ".join(params[20].split()).split()[1])
        point_info["Md1"] = float(" ".join(params[21].split()).split()[1])
        point_info["Md2"] = float(" ".join(params[22].split()).split()[1])
        point_info["Md3"] = float(" ".join(params[23].split()).split()[1])
        sms =  blocks[5].split("\n")[1:-1]
        point_info["alpha_s"] = float(" ".join(sms[2].split()).split()[1])
        point_info["mbottom"] = float(" ".join(sms[4].split()).split()[1])
        point_info["mtop"] = float(" ".join(sms[5].split()).split()[1])
        masses = blocks[14].split("\n")[2:-2]
        point_info["mhiggs"] = float(" ".join(masses[4].split()).split()[1])
        point_info["mW"] = float(" ".join(masses[2].split()).split()[1])
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
            while not spnerr:#find a viable point
                utils.clean()
                candidate = generate_point()#generate a point from flat prior
                spnin = utils.write_spheno_input(candidate)#write the input for spheno
                spnerr = run_spheno(spnin,devnull) #run spheno, check if viable point

            if not run_feynhiggs(devnull):#run feynhiggs, replace higgs sector
                continue

            gm2_obs = run_gm2calc(slhapath="SPheno.spc")

#            hb_obs = run_higgsbounds(slhapath="SPheno.spc")
#            hb_obs = run_higgsbounds_chi2(slhapath="SPheno.spc")
#            hs_obs = run_higgssignals(slhapath="SPheno.spc")

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
#            for obs in hb_obs:
#                observables[obs] = hb_obs[obs]
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
            while not spnerr:#find a viable point
                utils.clean()
                candidate = generate_point(lastaccepted)#generate a point from the last point
                spnin = utils.write_spheno_input(candidate)#write the input for spheno
                spnerr = run_spheno(spnin,devnull) #run spheno, check if viable point

            if not run_feynhiggs(devnull):#run feynhiggs, replace higgs sector
                continue

            gm2_obs = run_gm2calc(slhapath="SPheno.spc")

#            hb_obs = run_higgsbounds(slhapath="SPheno.spc")
#            hb_obs = run_higgsbounds_chi2(slhapath="SPheno.spc")
#            hs_obs = run_higgssignals(slhapath="SPheno.spc")

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
#            for obs in hb_obs:
#                observables[obs] = hb_obs[obs]
                
#            for obs in mmgs_obs:
#                observables[obs] = mmgs_obs[obs]

            _l = likelihood.make_decision(observables,lastaccepted["likelihood"])#get likelihood
            finite_lh = _l != 0
        if _l<0:
            move +=1
            if iter_ix == start+tend:
                print "Made all "+str(tend)+" iterations, moving "+outname+" to storage"
                outtree.BuildIndex("chain_index","iteration_index")
                outtree.Write()
                outroot.Close()
                os.system(" ".join(["mv",outname,outdir]))
            continue #point was not accepted

        lastaccepted["likelihood"]=_l
        lastaccepted["iteration_index"] = iter_ix
        lastaccepted["accepted_index"] =lastaccepted["accepted_index"]+1
        lastaccepted["chain_index"] = chainix

        for obs in observables.keys():
            lastaccepted[obs] = observables[obs]["value"]

        lastaccepted = prepare_fill(lastaccepted,outtree)#add the rest of the point info, fill the tree branches
        outtree.Fill()
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

    

