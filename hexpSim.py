#!/usr/local/bin/python3
#
# hexpSim.py
#
# Simulate data using the 'fakit' command
# Uses the HEX-P response files, and the relxill model
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
version='0.3a'
date='-- Thu Sep 15 00:37:06 MDT 2022'
author='Javier Garcia'
#
ul=[]
ul.append("usage: %prog [options]")
ul.append("")
ul.append("Run a series of HEX-P simulations for BH spin measurements with RELXILL")
#ul.append("PREFIX can be a single PHA file or a group (e.g. *.pha)")
usage=""
for u in ul: usage+=u+'\n'

parser=OptionParser(usage=usage)
parser.add_option("-v","--version",action="store_true",dest="version",default=False,help="show version number")
parser.add_option("-e","--exposure",dest="exposure",default="10",help="source exposure in ks (default: 10)")
parser.add_option("-s","--scale",dest="scale",default="100",help="scale factor for background exposure w.r.t. source (default: 100)")
parser.add_option("-f","--flux",dest="flux",default="10",help="source flux (2-10 keV) in mCrab (default: 10)")
parser.add_option("-i","--ier",dest="iter",default="10",help="Number of simulations (default: 10)")

(options,args)=parser.parse_args()

if options.version:
  print ('hexpSim.py version:',version,date)
  print ('Author:',author)
  sys.exit()

exposure=float(options.exposure)
scale=float(options.scale)
flux=float(options.flux)
maxiter=int(options.iter)


#-----

# Reponse files
rmf1 = "HEXP_LET_v02.rmf"
arf1 = "HEXP_LET_PSFcor_v02.arf"
bgd1 = "HEXP_LET_LEO_bkg_5arcsec_v02.fits"

rmf2 = "HEXP_HET_v01.rmf"
arf2 = "HEXP_HET_2x_PSFcor_v01.arf"
bgd2 = "HEXP_HET_2x_bkg_5arcsec_v01.fits"

#specfile
specbase1='sim-HEX-P_LET'
specbase2='sim-HEX-P_HET'

# Source flux
f1mC = 3.e-11  # Flux for 1 mCrab
sflux = flux*f1mC

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
#spin_vals=['0.0','0.5','0.9','0.95','0.998']
#spin_vals=['0.95']

# Numer of simulations per case
#maxiter = 100

# Source and Background exposure times (sec)
stime=exposure*1.e3
btime=scale*stime

# Load local models
AllModels.lmod("relxill")

# Output
outfile='spfit-spin.error'
f = open(outfile, 'w')

# Output
f.write('# Flux (2-10 keV, erg/cm^2/s) = '+str(sflux)+'    Exposure (ks) = '+str(exposure)+'\n')

for spin in spin_vals:

  specfile1=specbase1+spin+'.fak'
  grpspecfile1=specbase1+spin+'_grp.fak'

  specfile2=specbase2+spin+'.fak'
  grpspecfile2=specbase2+spin+'_grp.fak'

  for it in range(maxiter):

    print ('*** spin',spin,'iter',it+1)

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

    # Simulate LEO spectrum
    #response, arf, background, exposure, correction, backExposure, fileName
    fs1 = FakeitSettings(rmf1,arf1,bgd1,stime,1.,btime,specfile1)
    AllData.fakeit(1, fs1)

    # Unload data
    AllData.clear()

    # Simulate HEO spectrum
    #response, arf, background, exposure, correction, backExposure, fileName
    fs2 = FakeitSettings(rmf2,arf2,bgd2,stime,1.,btime,specfile2)
    AllData.fakeit(1, fs2)

    # Unload data
    AllData.clear()

    # Grop the data
    bashcommand='ftgrouppha infile='+specfile1+' outfile='+grpspecfile1+' grouptype=opt respfile='+rmf1
    os.system(bashcommand)

    bashcommand='ftgrouppha infile='+specfile2+' outfile='+grpspecfile2+' grouptype=opt respfile='+rmf2
    os.system(bashcommand)

    # Unload the model
    AllModels.clear()

    # Unload data
    AllData.clear()

    #
    # Load data
    AllData('1:1 '+grpspecfile1+' 2:2 '+grpspecfile2)

    s1 = AllData(1)
    s2 = AllData(2)

    # Ignore/notice data
    s1.ignore("0.-0.3,10.5-**")
    s2.ignore("0.-2.0,195.-**")

    # Define the Model
    AllModels += "const*tbabs*relxillCp"
    m1 = AllModels(1)
    m2 = AllModels(2)

    m1(1).values = "1. -1"       # Constant
    m1(2).values = "0.5 -1"      # Tbabs   Nh
    m1(3).values = "5. -1"       # relxill Index1
    m1(4).values = "3. -1"       # relxill Index2
    m1(5).values = "15. -1"      # relxill Rbr
    m1(6).values = spin+" 1"     # relxill a
    m1(7).values = "45. 1"       # relxill Incl
    m1(8).values = "-1. -1"      # relxill Rin
    m1(9).values = "400. -1"     # relxill Rout
    m1(10).values = "1.e-2 -1"   # relxill z
    m1(11).values = "2.0 1"      # relxill gamma
    m1(12).values = "2. -1"      # relxill logxi
    m1(13).values = "3. -1"      # relxill Afe
    m1(14).values = "100. 1"     # relxill kTe
    m1(15).values = "1.0 -1"     # relxill refl_frac
    #m1(16).values = "1. 1"       # relxill norm
    m1(16).values = str(norm)      # relxill norm

    m2(1).values = "1. 1"

    # Fit
    Fit.renorm()  # It appears that the renormalization is not working!
    Fit.perform()

    # Calculate Errors
    Fit.error("maximum 3. 2.706 6")

    # Save values
    parval = m1(6).values[0]
    uplim  = m1(6).error[1]
    lolim  = m1(6).error[0]

    #If pegged at the upper boundary
    if (uplim==0.):
       uplim=0.998
       parval=lolim

    #If pegged at the lower boundary
    if (lolim==0.):
       lolim=-0.998
       parval=uplim

    cval.append(parval)
    wval.append(1./((uplim-lolim)/2.)**2.) # There should be a better way!

    # Unload data
    AllData.clear()

    # Unload the model
    AllModels.clear()

    #f.write(str(spin)+' '+str(parval)+' '+str(lolim)+' '+str(uplim)+'\n')

  # Calculate the weighted average and error
  avgVal = np.average(np.array(cval),weights=np.array(wval),returned=True)

  # Output
  #f.write(' ------- \n')
  f.write(str(spin)+' '+str(avgVal[0])+' '+str(np.sqrt(1./avgVal[1]))+'\n')

  # Output
  #f.write(str(spin)+' '+str(parval)+' '+str(lolim)+' '+str(uplim)+'\n')

f.close()
#
#
sys.exit()
# ------------------------------------------------------------------------------
