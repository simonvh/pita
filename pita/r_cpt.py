#!/usr/bin/env python
from rpy2.robjects.packages import importr
from rpy2.rinterface import RRuntimeError
import rpy2.robjects as robjects

# Try to load changepoint library, install it if it's not installed
try:
    cp = importr('changepoint')
except RRuntimeError:
    utils = importr('utils')
    # Choose first available mirror from the list
    utils.chooseCRANmirror(ind=1)
    utils.install_packages('changepoint')
    cp = importr('changepoint')

def cpt(l):
    """
    Call the mean.cpt function from the R changepoint package and 
    return the changepoint as a float.
    """
    
    vector = robjects.FloatVector([x for x in l])
    result = cp.cpt_mean(vector)

    try: 
        return float(cp.cpts(result)[0])
    except:
        return None
