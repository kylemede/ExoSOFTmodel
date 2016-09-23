from __future__ import absolute_import
import numpy as np
import multiprocessing
import matplotlib.pyplot as plt
import ExoSOFTmodel
import emcee
from six.moves import range


def main():
    ## load up default settings dictionary
    sd = ExoSOFTmodel.load_settings_dict('./examples/settings.py')
    
    ## instantiate main objects/classes: ExoSOFTpriors, ExoSOFTdata and ExoSOFTparams.  And minor class, ExoSOFTmodel
    Model = ExoSOFTmodel.ExoSOFTmodel()  ###$$$$$$$$$ Kill this by merging into another one? OR merging priors into here?
    
    Params = ExoSOFTmodel.ExoSOFTparams(sd['omega_offset_di'], 
             sd['omega_offset_rv'], sd['vary_tc'], sd['tc_equal_to'], 
             sd['data_mode'], sd['low_ecc'], sd['range_maxs'], sd['range_mins'], 
             sd['num_offsets'])
    
    (epochs_di, rapa, rapa_err, decsa, decsa_err) = ExoSOFTmodel.load_di_data(sd['di_dataFile'])
    (epochs_rv, rv, rv_err, rv_inst_num) = ExoSOFTmodel.load_rv_data(sd['rv_dataFile'])
    
    Data = ExoSOFTmodel.ExoSOFTdata(epochs_di, epochs_rv, rapa, rapa_err, decsa, decsa_err,
                 rv, rv_err, rv_inst_num,sd['data_mode'], sd['pasa'])
    
    Priors = ExoSOFTmodel.ExoSOFTpriors(ecc_prior=sd['ecc_prior'], 
             p_prior=sd['p_prior'], inc_prior=sd['inc_prior'], 
             m1_prior=sd['m1_prior'], m2_prior=sd['m2_prior'], 
             para_prior=sd['para_prior'], inc_min=sd['inc_min'],
             inc_max=sd['inc_max'], p_min=sd['p_min'], p_max=sd['p_max'],
             para_est=sd['para_est'], para_err=sd['para_err'], 
             m1_est=sd['m1_est'], m1_err=sd['m1_err'], m2_est=sd['m2_est'], 
             m2_err=sd['m2_err'])
    
    ln_posterior = ExoSOFTmodel.ln_posterior
    
    
    ## define a set of starting parameters
    # The user can use any reasonable guess here.  For this example, 
    # we will use previously determined best-fit values for simplicity. 
    # pars: [m1,m2,parallax,long_an, e/sqrte_sinomega,to/tc,p,inc,arg_peri/sqrte_cosomega,v1,v2...]
    m2 = ExoSOFTmodel.kg_per_mjup/ExoSOFTmodel.kg_per_msun
    sqrte_sinomega = np.sqrt(0.048)*np.sin((np.pi/180.0)*14.8)
    sqrte_cosomega = np.sqrt(0.048)*np.cos((np.pi/180.0)*14.8)
    start_params = [1.0,m2,50,100.6,sqrte_sinomega,2450639.0,11.9,45.0,sqrte_cosomega,0]
    
    ncpu = multiprocessing.cpu_count()-1
    ndim = len(sd['range_maxs']) # number of parameters in the model 
    nwalkers = 50 # number of MCMC walkers 
    nburn = 10 # "burn-in" to stabilize chains 
    nsteps = 200 # number of MCMC steps to take 
    starting_guesses = []
    for i in range(nwalkers):
        starting_guesses.append(start_params)
    starting_guesses = np.array(starting_guesses,dtype=np.dtype('d'))
    
    sampler = emcee.EnsembleSampler(nwalkers, ndim, ln_posterior, 
                                    args=[Model, Data, Params, Priors], threads=ncpu)
    sampler.run_mcmc(starting_guesses, nsteps)
    
    # chain is of shape (nwalkers, nsteps, ndim)
    # discard burn-in points and reshape
    trace = sampler.chain[:, nburn:, :] 
    trace = trace.reshape(-1, ndim)
    
    ## Show walkers during burn-in (NOTE: by starting at the best-fit, this will look the same as post-burn-in)
    labels = ['m2', 'period', 'inclination']
    fig = plt.figure(figsize=(10,5))
    j=0
    for i, chain in enumerate(sampler.chain[:, :nburn, :].T): 
        if i in [1,6,7]:
            plt.subplot(3, 1, j+1)
            plt.plot(chain, drawstyle='steps', color='k', alpha=0.2)
            plt.ylabel(labels[j])
            j+=1
    plt.show()
    
    ## Show walkers after burn-in
    fig = plt.figure(figsize=(10,5))
    j=0
    for i, chain in enumerate(sampler.chain[:, nburn:, :].T):
        if i in [1,6,7]:
            plt.subplot(3, 1, j+1)
            plt.plot(chain, drawstyle='steps', color='k', alpha=0.2)
            plt.ylabel(labels[j])
            j+=1
    plt.show()
        
    ## inspect posteriors
    fig = plt.figure(figsize=(12,3))
    j=0
    for i in range(ndim):
        if i in [1,6,7]:
            plt.subplot(1,3,j+1)
            plt.hist(trace[:,i], 100, color="k", histtype="step")
            yl = plt.ylim()
            plt.vlines(start_params[i], yl[0], yl[1], color='blue', lw=3, alpha=0.25, label='true')
            plt.title("{}".format(labels[j]))
            plt.legend()
            j+=1
    plt.show()
    
    ## make a corner plot
    #import corner
    #fig = corner.corner(trace, labels=labels, quantiles=[0.16, 0.5, 0.84], truths=start_params)
    #plt.show()

if __name__ == '__main__':
    main()