import numpy as np
import sys, os
from numba import njit

def arg_parse():
    ret = dict()
    for key in sys.argv:
        if key[0] == '--':
            key = key[1:]
            if "=" in key:
                key, value = key.split("=")
                ret[key] = value
            else:
                ret[key] = None
    return ret

def initialize(N, L):
    pass

@njit
def lj_potential(L, r_cut, conf, i, pos): #pos must have shape (1,3)
    #pbc check
    #                                   L/2 - (L/2 - |rj-ri|)%L/2  -> |rj - ri| if |rj-ri|<L/2  else  L/2 - (rj - ri)%L/2
    former = np.sqrt( .5 - (.5 - np.abs( np.power( conf[:i]    - pos , 2 ) ) % .5 ).sum(axis = 1) ) * L   # particle 0 to i   % .5 means PBC
    later  = np.sqrt( .5 - (.5 - np.abs( np.power( conf[i+1:]  - pos , 2 ) ) % .5 ).sum(axis = 1) ) * L   # particle i+1 to N-1
    # cut 
    former = former[former<r_cut]
    later = later[later<r_cut]
    
    r12 = np.power(former,12).sum()+np.power(later,12).sum()
    r6 = np.power(former,6).sum()+np.power(later,6).sum()
    
    pot = r12 - r6
    return pot


#default setting 
settings = {"initial_prefix" : 'cnf_init.npy', "L" : 1, "N" : 10, "overlap" : 0.17, "r_cut" : 5}
output = "cnf_output.npy"
nml = arg_parse()

# Set default values, check keys and typecheck values
defaults = {"nblock":10, "nstep":1000, "temperature":1.0, "r_cut":2.5, "dr_max":0.15}
for key, val in nml.items():
    if key in defaults:
        assert type(val) == type(defaults[key]), key+" has the wrong type"
    else:
        print('Warning', key, 'not in ',list(defaults.keys()))

if 'seed' in nml:
    np.random.seed(nml['seed'])
else:
    np.random.seed()

# Set parameters to input values or defaults
nblock      = nml["nblock"]      if "nblock"      in nml else defaults["nblock"]
nstep       = nml["nstep"]       if "nstep"       in nml else defaults["nstep"]
temperature = nml["temperature"] if "temperature" in nml else defaults["temperature"]
r_cut       = nml["r_cut"]       if "r_cut"       in nml else defaults["r_cut"]
dr_max      = nml["dr_max"]      if "dr_max"      in nml else defaults["dr_max"]
L           = nml["L"]           if "L"           in nml else settings["L"]
N           = nml["N"]           if "N"           in nml else settings["N"]
initial_prefix = nml["initial_prefix"]           if "initial_prefix"           in nml else settings["initial_prefix"]

rho = N/L**3
# read npy configuration
# position nomalize to L
try:
    conf = np.load(settings["initial_prefix"])

except:
    conf = initialize(N, L)

total_E = lj_potential(L, r_cut, conf)

trial_move(conf, )

