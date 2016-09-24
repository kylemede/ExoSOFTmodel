#@Author: Kyle Mede, kylemede@astron.s.u-tokyo.ac.jp or kylemede@gmail.com
from __future__ import absolute_import
from __future__ import print_function
import numpy as np
from .cytools import orbit, model_input_pars
import KMlogger
from six.moves import range

log = KMlogger.getLogger('main.model',lvl=100,addFH=False)

class ExoSOFTmodel(object):
    """
    """
    def __init__(self):
        
        ####################
        ## member variables
        ####################    
        #resulting fit values   
        self.chi_squared_3d = 0
        self.chi_squared_di = 0
        self.chi_squared_rv = 0
        self.prior = 0
        ## TRACK BEST CHI SQUAREDS FOUND SO FAR IN HERE?
        ## makes more sense to change this to 'ExoSOFTresults' and name the object 'Results'??!!
        
        ## prior functions??
    
class ExoSOFTparams(object):
    """
    
    
    
    +---+--------------------+---------------+-------------------+-------+
    |   |  Directly Varried  | Model Inputs  | Stored Parameters |       |
    +---+--------------------+---------------+-------------------+-------+
    |   |    direct_pars     | model_in_pars |   stored_pars     |       |
    +---+--------------------+---------------+-------------------+-------+
    | i |     Parameter      |   Parameter   |     Parameter     | units |
    +===+====================+===============+===================+=======+
    | 0 |Mass of Primary (m1)|      m1       |        m1         |  Msun |
    +---+--------------------+---------------+-------------------+-------+
    .
    .
    .$$ FILL THIS OUT!!!!
    
    """
    def __init__(self, omega_offset_di, omega_offset_rv, vary_tc, tc_equal_to, 
                 di_only, low_ecc, range_maxs, range_mins, num_offsets):
        # params that effect calculating the full list of params from the directly varied one
        self.omega_offset_di = omega_offset_di
        self.omega_offset_rv = omega_offset_rv
        self.vary_tc = vary_tc
        self.tc_equal_to = tc_equal_to
        self.di_only = di_only
        self.low_ecc = low_ecc
        ## max/min ranges
        self.maxs = range_maxs 
        self.mins = range_mins
        ## prep versions of all param arrays
        self.num_offsets = num_offsets
        self.direct_pars = np.zeros((9+num_offsets),dtype=np.dtype('d'))
        # model_in_params: [m1,m2,parallax,long_an,e,to,tc,p,inc,arg_peri,arg_peri_di,arg_peri_rv,a_tot_au,K]
        self.model_in_pars = np.zeros((14),dtype=np.dtype('d'))
        # stored_pars: [m1,m2,parallax,long_an,e,to,tc,p,inc,arg_peri,a_tot_au,chi_sqr,K,v1,v2...]
        self.stored_pars = np.zeros((12+num_offsets),dtype=np.dtype('d'))
        self.offsets = np.zeros((num_offsets),dtype=np.dtype('d'))
        #check_pars: [m1, m2, parallax, long_an, e, to/tc, p, inc, arg_peri]
        self.check_pars = np.zeros((9+num_offsets),dtype=np.dtype('d'))
        self.taea = np.zeros((2),dtype=np.dtype('d'))
        
    def make_model_in(self):
        """
        Convert directly varied parameters into a comprehensive list
        of those used ans inputs to during model calculations.
        """        
        model_input_pars(self.direct_pars, self.low_ecc, self.tc_equal_to, 
                   self.vary_tc, self.di_only, self.omega_offset_di, self.omega_offset_rv,
                   self.offsets, self.model_in_pars)
        ## Wrap periodic params into allowed ranges.  ie. long_an and arg_peri
        m_par_ints = [3,9]
        min_max_ints = [3,8]
        for i in [0,1]:
            if self.mins[min_max_ints[i]] > self.model_in_pars[m_par_ints[i]]:
                #print('par was '+str(model_input_pars[m_par_ints[i]]))
                self.model_in_pars[m_par_ints[i]]+=360.0
                #print('now '+str(model_input_pars[m_par_ints[i]]))
            elif self.model_in_pars[m_par_ints[i]] > self.maxs[min_max_ints[i]]:
                #print('par was '+str(model_input_pars[m_par_ints[i]]))
                self.model_in_pars[m_par_ints[i]]-=360.0
                #print('now '+str(model_input_pars[m_par_ints[i]]))
        #print(repr(self.model_in_pars))
        
    def make_stored(self,chi_squared):
        """ 
        Push values in model_in_params, offsets and the resulting 
        chi_squared_3d into an array to be stored on disk during ExoSOFT.  
        Not sure how to make this work with emcee or other tools...
        """
        # model_in_params: [m1,m2,parallax,long_an,e,to,tc,p,inc,arg_peri,arg_peri_di,arg_peri_rv,a_tot_au,K]
        # stored_pars: [m1,m2,parallax,long_an,e,to,tc,p,inc,arg_peri,a_tot_au,chi_sqr,K,v1,v2...]
        self.stored_pars[0:10] = self.model_in_pars[0:10]
        self.stored_pars[10] = self.model_in_pars[12]
        self.stored_pars[11] = chi_squared
        self.stored_pars[11:] = self.offsets[:]
         
    def check_range(self):
        """Determine if all parameters in the full list are within their 
        allowed ranges.
        Range arrays corrispond to parameters in:
        [m1, m2, parallax, long_an, e, to/tc, p, inc, arg_peri, v1,v2,...]
        """        
        self.check_pars[0:5] = self.model_in_pars[0:5]
        if self.vary_tc:
            self.check_pars[5] = self.model_in_pars[6]
        else:
            self.check_pars[5] = self.model_in_pars[5]
        self.check_pars[6:9] = self.model_in_pars[7:10]
        self.check_pars[9:] = self.offsets[:]
        
        if len(self.check_pars)!=len(self.maxs)!=len(self.mins):
            print("LENGTH OF CHECK_PARAMS IS NOT EQUAL TO LENGTH OF MINS OR MAXS!!!")
        in_range = True
        for i in range(9):
            if (self.check_pars[i]>self.maxs[i]) or (self.check_pars[i]<self.mins[i]):
                in_range = False
        return in_range

class ExoSOFTdata(object):
    """
    An object to contain all the necessary data arrays and parameters to 
    calculate matching predicted data with the model.  All member variables 
    will remain constant throughout.
    
    Notes:
    -Except for rv_inst_num array, all other arrays must be ndarrays of double 
     precision floating point numbers (dtype=np.dtype('d')).
    -Arrays, epochs_di, rapa, rapa_err, decsa, and decsa_err must all have same length.
    -Arrays, epochs_rv, rv, rv_err and rv_inst_num must all have same length.
    
    Inputs:
    rv_inst_num = ndarray of positive signed or unsigned integers, of same length
                  as epochs_rv, rv, and rv_err.
    """
    def __init__(self, epochs_di, epochs_rv, rapa, rapa_err, decsa, decsa_err,
                 rv, rv_err, rv_inst_num, data_mode, pasa=False):
        
        self.epochs_di = epochs_di
        self.epochs_rv = epochs_rv
        # x/RA/PA
        self.rapa = rapa
        self.rapa_err = rapa_err
        self.rapa_model = np.zeros((len(epochs_di)),dtype=np.dtype('d'))
        # y/Dec/SA
        self.decsa = decsa
        self.decsa_err = decsa_err
        self.decsa_model = np.zeros((len(epochs_di)),dtype=np.dtype('d'))
        # RV
        self.rv = rv
        self.rv_err = rv_err
        self.rv_model = np.zeros((len(epochs_rv)),dtype=np.dtype('d'))
        # dataset/instrument number
        self.rv_inst_num = rv_inst_num
        self.data_mode = data_mode
        self.pasa = pasa

def ln_posterior(pars, Model, Data, Params, Priors):
    """
    Calculates the likelihood for a given set of inputs.
    Then calculate the natural logarithm of the posterior probability.
    
    -Model is of type ExoSOFTmodel.  Currently just holds resulting fit values.
    -Data is of type ExoSOFTdata, containing all input data and params to 
    produce predicted values of matching units, and arrays for predicted values.
    -Params is of type ExoSOFTparams, an class containing functions for 
    calculating versions of the 'pars' used as model inputs, and a version 
    that would be for storing to disk when ran in ExoSOFT.
    -Priors is of type ExoSOFTpriors, containing funtions for each parameter's
    prior, a function calculate combined prior given list of params, and any 
    variables necessary for those calculations.
    
    """    
    ## convert params from raw values
    Params.direct_pars = pars
    Params.make_model_in()
        
    ## Range check on proposed params, set ln_post=zero if outside ranges.
    ln_post = -np.inf
    in_range = Params.check_range()
    if in_range:         
        ## Call Cython func to calculate orbit. ie. -> predicted x,y,rv values.
        orbit(Params.model_in_pars, Params.offsets, Data.pasa, 
                      Data.data_mode, Data.epochs_di, Data.epochs_rv, 
                      Data.rv_inst_num, Data.rapa_model, Data.decsa_model, 
                      Data.rv_model,Params.taea)
                
        #print('measured rv ',repr(Data.rv))
        #print('model rv ',repr(Data.rv_model))
        ## Calculate chi squareds and then the 3d log likelihood
        #  NOTE: The log(2*pi*sigma**2) term is not included here as the 
        #        likelihood is always used as a likelihood ratio.
        chi_sqr_rv, chi_sqr_rapa, chi_sqr_decsa = 0, 0, 0
        if (len(Data.epochs_rv)>0) and (Data.data_mode!='DI'):
            chi_sqr_rv = np.sum((Data.rv-Data.rv_model)**2 / Data.rv_err**2)
        if (len(Data.epochs_di)>0) and (Data.data_mode!='RV'):
            chi_sqr_rapa = np.sum((Data.rapa-Data.rapa_model)**2 / Data.rapa_err**2)
            chi_sqr_decsa = np.sum((Data.decsa-Data.decsa_model)**2 / Data.decsa_err**2)
        chi_sqr_3d = chi_sqr_rv + chi_sqr_rapa + chi_sqr_decsa
        #print('chi_sqr_rv',chi_sqr_rv)
        #print('chi_sqr_rapa',chi_sqr_rapa)
        #print('chi_sqr_decsa',chi_sqr_decsa)
        #print('chi_sqr_3d',chi_sqr_3d)
        # Remember that chisqr = -2*log(Likelihood).  OR,
        ln_lik = -0.5*chi_sqr_3d
        #print('ln_lik',ln_lik)
        ## Make version of params with chi_sqr_3d for storing during ExoSOFT
        Params.make_stored(chi_sqr_3d)
        ## store the chi sqr values in model object for printing in ExoSOFT.
        Model.chi_squared_3d = chi_sqr_3d
        Model.chi_squared_di = chi_sqr_rapa + chi_sqr_decsa
        Model.chi_squared_rv = chi_sqr_rv
        
        ## Calculate priors
        prior = Priors.priors(Params.model_in_pars)
        Model.prior = prior
        #print('np.log(prior)',np.log(prior))
        #print('prior ',prior)
        ## calculate lnpost
        ln_post = np.log(prior) + ln_lik
        #print('ln_post ',ln_post)
        
    return ln_post

#EOF