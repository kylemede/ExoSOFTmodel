#@Author: Kyle Mede, kylemede@astron.s.u-tokyo.ac.jp

params_dict = {
#####################################
# Special Settings for the models ###
#####################################
# Operate in low eccenctricity mode? [bool]
# Then step through in sqrt(e)sin(omega) and sqrt(e)cos(omega) instead of e & omega directly
'low_ecc'   : True,
# Step through parameter space in Time of Center Transit (Inferior Conjunction)?  [bool]
'vary_tc' : False,
# take the time of center transit (inferior conjunction) into account? [bool]
'tc_equal_to' : True,
# force adding a value in degrees to argument of periapsis used in RV orbit fit [double]
# NOTE: Internally, 180deg is already added to arg_peri_rv.
# This is required as the RV values are for the primary and thus omega+pi
# compared to the value used for the secondary in the DI data.
## Another way to think of this is:
## arg_peri_di = arg_peri_companion = arg_peri
## arg_peri_rv = arg_peri_primary   = arg_peri +180
####$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$  MAYBE KILL BOTH omega_offset_rv AND omega_offset_di !?!?!
'omega_offset_rv' : 0.0,
##################################################
## Special settings DI model:
# force adding a value in degrees to argument of periapsis used in DI orbit fit [double]
'omega_offset_di' : 0.0,

###################################################
# Ranges for acceptable random number inputs ######
###################################################
## For Omega, e, T, inc and omega, None indicates to use default ranges.
## For Omega and omega, values can vary outside ranges, but are shifted by 
## +/-360 befire being stored.
# Minimum/Maximum allowed value for the mass of the primary body [double][Msun]
# NOTE: For DI only cases, use mass1 values as total mass and set mass2 values to zero.
#       This will be done during start up otherwise.
'm1_min' : 0.2,
'm1_max' : 2.55,
# Minimum/Maximum allowed value for the mass of the secondary body [double][Msun]
# Note: 1Mj ~ 0.00095 Msun
'm2_min' : 0.00001,
'm2_max' : 0.005,
# Minimum/Maximum allowed value for the Parallax [double][mas]
'para_min' : 1.0,
'para_max' : 100.0,
# Minimum/Maximum allowed value for the Longitude of the Ascending Node [double][deg]
'long_an_min' : 0.0,
'long_an_max' : 360.0,
# Minimum/Maximum allowed value for the Eccentricity [double]
'ecc_min' : 0.0,
'ecc_max' : 0.98,
# Minimum/Maximum value for the Time of Last Periapsis (or Time of Center Transit) [JD] OR None
# Default: [earliestsEpoch-period,earliestEpoch]
't_min' : 2449000,
't_max' : 2453500,
# Minimum/Maximum allowed value for the Period [double][yrs]
'p_min' : 1.0,
'p_max' : 50.0,
# Minimum/Maximum allowed value for the Inclination [double][deg] OR None
# Default: [0,180].  [0,90] is another popular choice.
'inc_min' : 1,
'inc_max' : 90,
# Minimum/Maximum allowed value for the Argument of Perigee [double][deg] OR None
# Default: [0,360]
'arg_peri_min' : -50,
'arg_peri_max' : 90,
# Minimum/Maximum values of Radial Velocity Offsets.  
# Must be one per set of RV data in same order as data comes in RVdata.dat, or the a single value to be used by all [comma separated list of doubles]
'offset_mins' :[-3],
'offset_maxs' :[3],
}

data_dict = {
             # full path to input astrometry data file. [string]
'di_dataFile': '/mnt/HOME/MEGA/Dropbox/EclipseWorkspaceDB/ExoSOFT/examples/DIdata.dat',                
# full path to input radial velocity data file. [string]
'rv_dataFile': '/mnt/HOME/MEGA/Dropbox/EclipseWorkspaceDB/ExoSOFTmodel/examples/RVdata.dat', 
# data mode, choices {'RV','DI','3D'} [string]
'data_mode' : '3D',
# Is the data in the DIdata.dat in PA,SA format? else, it is in E,N (x,y) format [bool]
'pasa'     : False,
}

priors_dict={
############################
#    System Information    #
# ONLY FOR GAUSSIAN PRIORS #
############################
#best estimate of primary's mass, and error [double][Msun]
'm1_est'  : 0.0,
'm1_err'  : 0.0,
#best estimate of secondary's mass, and error [double][Msun]
'm2_est'  : 0.0,
'm2_err'  : 0.0,
#best estimate of parallax, and error [double][mas]
'para_est': 50,
'para_err': 2.5,
##################################
# Choice of which priors to use  #
##################################
# For ALL: False indicates a flat prior. True indicates to use defaults.
'ecc_prior'  : False,
'p_prior'    : True,
# For the inclination prior, use strings or booleans to inducate the specific 
# function to use.  Either sin(i), cos(i) or flat.  True indicates sin(i).
'inc_prior'  : 'cos',
# For m1 and m2: use strings to indicate specific prior function.
# For m1: ['PDMF', 'IMF', True or False], default is 'PDMF'.
'm1_prior'   : True,
# for m2: ['PDMF', 'IMF', 'CMF', True or False], default is 'CMF'
'm2_prior'   : True,
'para_prior' : True,
}

######################
# Merge All dicts#
######################
settings = {}
for key in params_dict:
    settings[key]=params_dict[key]
for key in data_dict:
    settings[key]=data_dict[key]
for key in priors_dict:
    settings[key]=priors_dict[key]
