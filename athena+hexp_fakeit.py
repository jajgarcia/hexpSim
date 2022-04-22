#!/usr/local/bin/python3
#
# fakeit.py
#
# Simulate data using the fakit command
# Uses the NuSTAR response files, and the relxill model
# 
# This version uses the HEX-P responses and the parameters
# similar to GX 339-4 in the hard state
#
# Requires: xspec
#
import sys 
from xspec import *
from optparse import OptionParser
import os,os.path
import glob
import numpy as np
#
# ------------------------------------------------------------------------------
#
# MAIN PROGRAM
#
#
#
version='0.1a'
date='--Mon Aug 23 11:21:57 MDT 2021'
author='Javier Garcia'
#
ul=[]
ul.append("usage: %prog [options] PREFIX")
ul.append("")
ul.append("Get total counts in different bands for a given observation")
ul.append("PREFIX can be a single PHA file or a group (e.g. *.pha)")
usage=""
for u in ul: usage+=u+'\n'

parser=OptionParser(usage=usage)
parser.add_option("-v","--version",action="store_true",dest="version",default=False,help="show version number")
parser.add_option("-i","--input",dest="inpfile",default="",help="specify spectrum file (default: meg_-1.pha.gz)")
parser.add_option("-o","--output",dest="outfile",default="",help="specify alternative output file (default: counts.out)")

(options,args)=parser.parse_args()

#-----

# Reponse files
rmf1 = 'XIFU_CC_BASELINECONF_IMPROVED_SPECTRAL_RESOLUTION_2018_10_10.rmf'
arf1 = 'XIFU_CC_BASELINECONF_IMPROVED_SPECTRAL_RESOLUTION_2018_10_10.arf'
bgd1 = ''

rmf2 = 'HEXP_v001.rmf'
arf2 = 'HEXP_Ni_EffArea_D2F20_220keV.arf'
#bgd2 = 'HEXP_bkg_220keV.fits'
bgd2 = ''

#specfile
specbase1='sim-xifu'
specbase2='sim-hexp'

# Source flux
f1mC = 3.e-11  # Flux for 1 mCrab
sflux = 1.e-1*f1mC

# No chatter
Xset.chatter = 0

# Fit statistics
Fit.statMethod = "cstat"

# Query
Fit.query = 'yes'

# Set abundances and cross sections
Xset.abund = "wilm"
Xset.xsect = "vern"

# Parameter of interest
spin_vals=['0.0','0.2','0.4','0.6','0.8','0.9','0.95','0.998']

# Numer of simulations per case
maxiter = 1

# Source and Background exposure times (sec)
stime=1.e4   # 10 ks 
btime=0.

# Load local models
AllModels.lmod("relxill")

# Output
outfile='spfit-spin.error'
f = open(outfile, 'w')

# Output
f.write('# Flux (2-10)='+str(sflux)+' Exposure (s)='+str(stime)+'\n')

for spin in spin_vals:

  specfile1=specbase1+spin+'.fak'
  grpspecfile1=specbase1+spin+'_grp.fak'

  specfile2=specbase2+spin+'.fak'
  grpspecfile2=specbase2+spin+'_grp.fak'

  for it in range(maxiter):

    print ('*** spin',spin,'iter',it)

    # Define the Model
    m1 = Model("tbabs*relxillCp")

    m1(1).values = "0.5 -1"      # Tbabs   Nh
    m1(2).values = "5. -1"       # relxill Index1
    m1(3).values = "3. -1"       # relxill Index2
    m1(4).values = "15. -1"      # relxill Rbr
    m1(5).values = spin+" 1"     # relxill a
    m1(6).values = "45. 1"       # relxill Incl
    m1(7).values = "-1. -1"      # relxill Rin
    m1(8).values = "400. -1"     # relxill Rout
    m1(9).values = "1.e-2 -1"    # relxill z
    m1(10).values = "2.0 1"      # relxill gamma
    m1(11).values = "2. 1"       # relxill logxi
    m1(12).values = "3. -1"      # relxill Afe
    m1(13).values = "100. 1"     # relxill kTe
    m1(14).values = "1.0 -1"     # relxill refl_frac
    m1(15).values = "1. 1"       # relxill norm

    # Calculate 2-10 keV flux
    AllModels.calcFlux("2. 10.")
    flux=AllModels(1).flux[0]

    # New normalization
    norm = sflux/flux
    m1(15).values = str(norm)      # relxill norm

    # Reset variables
    cval=[]
    wval=[]

    # First simulate Athena XIFU
    #response, arf, background, exposure, correction, backExposure, fileName
    fs1 = FakeitSettings(rmf1,arf1,bgd1,stime,1.,btime,specfile1)
    AllData.fakeit(1, fs1)

    # Unload data
    AllData.clear()

    # Then simulate HEX-P
    #response, arf, background, exposure, correction, backExposure, fileName
    fs2 = FakeitSettings(rmf2,arf2,bgd2,stime,1.,btime,specfile2)
    AllData.fakeit(1, fs2)

    # Unload data
    AllData.clear()

    # Group the data
    bashcommand='ftgrouppha infile='+specfile1+' outfile='+grpspecfile1+' grouptype=opt respfile='+rmf1
    os.system(bashcommand)

    bashcommand='ftgrouppha infile='+specfile2+' outfile='+grpspecfile2+' grouptype=opt respfile='+rmf2
    os.system(bashcommand)

    # Unload the model
    AllModels.clear()

    #
    # Load data
    AllData('1:1 '+grpspecfile1+' 2:2 '+grpspecfile2)

    s1 = AllData(1)
    s2 = AllData(2)

    # Ignore/notice data
    s1.ignore("0.-0.4,10.-**")
    s2.ignore("0.-2.,220.-**")

    # Define the Model
    AllModels += "tbabs*relxillCp"
    m1 = AllModels(1)
    m2 = AllModels(2)

    m1(1).values = "0.5 -1"      # Tbabs   Nh
    m1(2).values = "5. -1"       # relxill Index1
    m1(3).values = "3. -1"       # relxill Index2
    m1(4).values = "15. -1"      # relxill Rbr
    m1(5).values = spin+" 1"     # relxill a
    m1(6).values = "45. 1"      # relxill Incl
    m1(7).values = "-1. -1"      # relxill Rin
    m1(8).values = "400. -1"     # relxill Rout
    m1(9).values = "1.e-2 -1"    # relxill z
    m1(10).values = "2.0 1"     # relxill gamma
    m1(11).values = "2. -1"      # relxill logxi
    m1(12).values = "3. -1"      # relxill Afe
    m1(13).values = "100. 1"    # relxill kTe
    m1(14).values = "1.0 -1"     # relxill refl_frac
    m1(15).values = "1. 1"      # relxill norm

    m2(15).values = "=15"

    # Fit
    Fit.renorm()
    Fit.perform()

    # Calculate Errors
    Fit.error("maximum 3. 2.706 5")

    # Save values
    parval=m1(5).values[0]
    uplim=m1(5).error[1]
    lolim=m1(5).error[0]

    #If pegged at the upper boundary
    if (uplim==0.):
       uplim=0.998
       parval=lolim

    #If pegged at the lower boundary
    if (lolim==0.):
       lolim=-0.998
       parval=uplim

    #cval.append(parval)
    #wval.append(1./((uplim-lolim)/2.)**2.) # There should be a better way!

    # Unload data
    AllData.clear()

    # Unload the model
    AllModels.clear()

  # Calculate the weighted average and error
  #avgVal = np.average(np.array(cval),weights=np.array(wval),returned=True)

  # Output
  #f.write(str(spin)+' '+str(mod*avgVal[0])+' '+str(np.sqrt(1./avgVal[1]))+'\n')

  # Output
  f.write(str(spin)+' '+str(parval)+' '+str(lolim)+' '+str(uplim)+'\n')

f.close()
#
#
sys.exit()
# ------------------------------------------------------------------------------
