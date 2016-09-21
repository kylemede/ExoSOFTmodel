from __future__ import absolute_import
import numpy as np
import multiprocessing
import matplotlib.pyplot as plt
import ExoSOFTmodel
import emcee


def main():
    ## load up default settings dictionary
    sd = ExoSOFTmodel.load_settings_dict('.examples/settings.py')
    
    ## instantiate main objects/classes: ExoSOFTpriors, ExoSOFTdata and ExoSOFTparams.  And minor class, ExoSOFTmodel
    Model = ExoSOFTmodel.ExoSOFTmodel()  ###$$$$$$$$$ Kill this by merging into another one? OR merging priors into here?
    
    Params = ExoSOFTmodel.ExoSOFTparams(sd['omega_offset_di'], 
             sd['omega_offset_rv'], sd['vary_tc'], sd['tc_equal_to'], 
             sd['di_only'], sd['low_ecc'], sd['range_maxs'], sd['range_mins'], 
             sd['num_offsets'])
    
    (epochs_di, rapa, rapa_err, decsa, decsa_err) = ExoSOFTmodel.load_di_data(sd['di_dataFile'])
    (epochs_rv, rv, rv_err, rv_inst_num) = ExoSOFTmodel.load_rv_data(sd['rv_dataFile'])
    
    Data = ExoSOFTmodel.ExoSOFTdata(epochs_di, epochs_rv, rapa, rapa_err, decsa, decsa_err,
                 rv, rv_err, rv_inst_num, sd['pasa'])
    
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
    start_params = [1.00580433581,0.000989368473617,49.2248464085,101.526288438,0.0654089809493,2450678.50145,11.9537449022,43.8878838314,0.208918200633,5.23954778104,46.837839368,8.9044445724,-0.00666695395887] 
    #### Need to convert the ecc and arg_peri to low_ecc values!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    ncpu = multiprocessing.cpu_count()
    ndim = len(sd['range_maxs']) # number of parameters in the model 
    nwalkers = 50 # number of MCMC walkers 
    nburn = 1000 # "burn-in" to stabilize chains 
    nsteps = 20000 # number of MCMC steps to take 
    starting_guesses = np.random.rand(nwalkers, ndim)  ### need to get start_params in here!!!!!!!!!!!!
    
    sampler = emcee.EnsembleSampler(nwalkers, ndim, ln_posterior, 
                                    args=[Model, Data, Params, Priors], threads=ncpu)
    sampler.run_mcmc(starting_guesses, nsteps)
    
    # chain is of shape (nwalkers, nsteps, ndim)
    # discard burn-in points and reshape
    trace = sampler.chain[:, nburn:, :] 
    trace = trace.reshape(-1, ndim)
    
    ####### modify below to work for my larger example by only considering a sampling of my parameters, say m2 period and inc. !!!!!!!!!!!!!!!!!!!!!!!!!
    
    ## Show walkers during burn-in (NOTE: by starting at the best-fit, this will look the same as post-burn-in)
    labels = ['intercept', 'slope', 'sigma']
    fig = plt.figure(figsize=(10,5))
    for i, chain in enumerate(sampler.chain[:, :nburn, :].T):
        plt.subplot(3, 1, i+1)
        plt.plot(chain, drawstyle='steps', color='k', alpha=0.2)
        plt.ylabel(labels[i])
        
    ## Show walkers after burn-in
    fig = plt.figure(figsize=(10,5))
    for i, chain in enumerate(sampler.chain[:, nburn:, :].T):
        plt.subplot(3, 1, i+1)
        plt.plot(chain, drawstyle='steps', color='k', alpha=0.2)
        plt.ylabel(labels[i])
        
    ## inspect posteriors
    fig = plt.figure(figsize=(12,3))
    for i in range(ndim):
        plt.subplot(1,3,i+1)
        plt.hist(trace[:,i], 100, color="k", histtype="step")
        yl = plt.ylim()
        plt.vlines(start_params[i], yl[0], yl[1], color='blue', lw=3, alpha=0.25, label='true')
        plt.title("{}".format(labels[i]))
        plt.legend()

    ## make a corner plot
    import corner
    fig = corner.corner(trace, labels=labels, quantiles=[0.16, 0.5, 0.84], truths=start_params)
    

if __name__ == '__main__':
    main()