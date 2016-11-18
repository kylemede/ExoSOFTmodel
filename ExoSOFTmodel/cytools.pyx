# cython: embedsignature=True
import constants as const

#python setup.py build_ext --inplace

#List of available C funcs in math.h:
#https://en.wikipedia.org/wiki/C_mathematical_functions
"""
Docstring at top of cytools.pyx
"""

cdef anomalies(double p, double tc, double to, double ecc, 
               double epoch, double [:] taea):
    """ 
    anomalies(double p, double tc, double to, double ecc, 
               double epoch, double [:] taea)
               
    Calculate the True and Eccentric Anomalies.
    Remember that in RV, there is occasionally a phase shift due to 
    the Tc!=To, that doesn't exist in DI.
    So, for RV pass in both values, for DI just set both to To.
    
    TA and EA (ta,ea) are calculated in radians.
    """
    cdef double ma, multiples, e_prime, ta, ea
    cdef int newton_count
    cdef bint warnings_on
    
    cdef extern from "math.h":
        double sin(double _x)
        double cos(double _x)
        double acos(double _x)
        double fabs(double _x)
        double floor(double _x)
        
    #$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    warnings_on = True  #$$$$ for debugging, kill this after finished testing?
    #$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
        
    ## Calc Mean Anomaly
    ma = (2.0*const.pi*(epoch-2.0*to+tc))/(p*const.days_per_year)
    #print('ma ',ma)
    #Find if over one full orbit and shift into [-2pi,2pi]
    #print("ma/2.0*const.pi ", repr(ma/(2.0*const.pi)))
    multiples = floor(ma/(2.0*const.pi))
    #print('multiples ',multiples)
    ma -= multiples*2.0*const.pi
    #print('ma ',ma)
    # shift into [0,2pi]
    if ma<0:
        ma += 2.0*const.pi
    # set E and TA to M by default for circular orbits
    ea, ta = ma, ma
    
    ## if not circular, calculate their specific values
    if ecc>0:
        #print('ecc ',ecc)
        # double check M is not zero or 2pi
        if ma not in [0.0,2.0*const.pi]:
            #print("starting newton")
            ea_prime = ma + ecc*sin(ma) + ((ecc**2.0)/2.0*ma)*sin(2.0*ma)
            newton_count = 0
            while (fabs(ea-ea_prime)>1e-10) and (newton_count<50):
                ea = ea_prime
                ea_prime = ea - ( (ea-ecc*sin(ea)-ma))/(1.0-ecc*cos(ea))
                newton_count +=1
            # double check E found satisfies original equation
            if (fabs((ea-ecc*sin(ea))-ma)>1e-5) or (newton_count>49):
                if warnings_on:
                    print "PROBLEM!! resulting E from Newton's loop isn't within error limit!!!"
                    if True:
                        print "M = "+str(ma)
                        print "E = "+str(ea)
                        print "T0 = "+str(to)
                        print "Tc = "+str(tc)
                        print "P = "+str(p)
                        print "Eprime = "+str(ea_prime)
                        print "NewtonCount = "+str(newton_count)
            #print('newton_count ',newton_count)
            # calculate TA from E
            ta = acos( (cos(ea)-ecc)/(1.0-ecc*cos(ea)) )
            #print('before sign reversal ta ',ta)
            #print('ea ',ea)
            # both increase properly 0->180, but as E properly continues to 
            # increase from 180->360, TA decreases.  So correcting for that here.
            if ea>const.pi:
                #print('reversing ta sign')
                ta = (2.0*const.pi) - ta
    #print('last in anomalies ta ',ta)
    #print('last in anomalies ea ',ea)
    #print ''
    taea[0],taea[1] = ta, ea
                
cdef orbit_rv(double [:] epochs, double ecc, double to, double tc, double p, 
              double inc, double arg_peri_rv, double k, double [:] offsets,
               int [:] rv_inst_num, double [:] rv_model, double [:] taea):
    """ 
    Calculates the predicted rv for a each epoch in epochs array.
     
    model value = calculated rv + the instrument specific offset.
    """
    cdef double top, ea, ta, rv, arg_peri_rad
    cdef int npts
    cdef extern from "math.h":
        double cos(double _x)
        
    #set ea,ta as zeros for now, they get corrected in anomalies function.
    npts = epochs.shape[0]
    arg_peri_rad = arg_peri_rv*(const.pi/180.0)
    
    #print('taea ',repr(taea))
    
    ## calculate rv for each epoch in the data
    for i in range(npts):
        ## Call anomalies     
        #print('before anomalies')   
        anomalies(p, tc, to, ecc, epochs[i], taea)
        #print('back from anomalies')
        #print('taea[0] ',repr(taea[0]))
        #print('taea[1] ',repr(taea[1]))
        ta, ea = taea[0], taea[1]
        #print('inside orb rv ta ',ta)
        #print('inside orb rv ea ',ea)
        ## Calc predicted RV + the instrument specific offset
        rv = k*( cos(ta+arg_peri_rad) + ecc*cos(arg_peri_rad) )
        #print('inside orb rv ta rv ',rv)
        #print ''
        rv_model[i] = rv + offsets[rv_inst_num[i]]
        
         
cdef orbit_di(double [:] epochs, double long_an, double ecc, double to, 
              double tc, double p, double inc, double arg_peri_di,
              double a_tot_au, double parallax, double [:] rapa_model,
              double [:] decsa_model, bint pasa, double [:] taea):                 
    """
    Calculates the predited x/RA/PA and y/Dec/SA for a single epoch.
    """
    cdef double ea, ta, a_thi, b_thi, f_thi, g_thi, x_thi, y_thi, ra, dec, pa, sa
    cdef double long_an_rad, arg_peri_rad, inc_rad, sep_arcsec
    cdef int npts
    cdef extern from "math.h":
        double sin(double _x)
        double cos(double _x)
        double sqrt(double _x)
        double atan2(double _x,double _y)
        
    #set ea,ta as zeros for now, they get corrected in anomalies function.
    npts = epochs.shape[0]
    long_an_rad = long_an*(const.pi/180.0)
    arg_peri_rad = arg_peri_di*(const.pi/180.0)
    inc_rad = inc*(const.pi/180.0)
    sep_arcsec = a_tot_au*(parallax/1000.0)
    
    ## Calculate predicted x/RA/PA and y/Dec/SA values for each epoch in the data
    for i in range(npts):
        ## Call anomalies        
        anomalies(p, tc, to, ecc, epochs[i], taea)
        ta, ea = taea[0], taea[1]
        #print('inside orb di ta ',ta)
        #print('inside orb di ea ',ea)
        ## Use Thiele-Innes Method to calculate predicted x/RA/PA and y/Dec/SA
        a_thi = sep_arcsec * ( cos(long_an_rad)*cos(arg_peri_rad) - sin(long_an_rad)*sin(arg_peri_rad)*cos(inc_rad))
        b_thi = sep_arcsec * ( sin(long_an_rad)*cos(arg_peri_rad) + cos(long_an_rad)*sin(arg_peri_rad)*cos(inc_rad))
        f_thi = sep_arcsec * (-cos(long_an_rad)*sin(arg_peri_rad) - sin(long_an_rad)*cos(arg_peri_rad)*cos(inc_rad))
        g_thi = sep_arcsec * (-sin(long_an_rad)*sin(arg_peri_rad) + cos(long_an_rad)*cos(arg_peri_rad)*cos(inc_rad))
        # The coordinates of the unit orbital ellipse in the true plane (Binnendijk)
        x_thi = cos(ea)-ecc
        y_thi = sqrt(1.0-ecc*ecc)*sin(ea)
        ## Calculate the predicted x&y in ["], or PA[deg], SA["]
        #  KEY NOTE: x_TH-I = y_plot = North = Dec = A*X +F*Y
        #            y_TH-I = x_plot = East  = RA  = B*X +G*Y
        dec = a_thi*x_thi + f_thi*y_thi
        ra  = b_thi*x_thi + g_thi*y_thi
        ## check which and then store RA/Dec or PA/SA accordingly
        if pasa:
            # convert RA and Dec to PA and SA
            pa = atan2(ra,dec)
            if (pa<0.0):
                pa+=2.0*const.pi
            sa = sqrt(ra*ra+dec*dec)
            rapa_model[i] = pa*(180.0/const.pi)
            decsa_model[i] = sa
        else:
            rapa_model[i] = ra
            decsa_model[i] = dec
        
def model_input_pars(double [:] pars, bint low_ecc, bint tc_equal_to, 
                     bint vary_tc, str data_mode, double omega_offset_di, 
                     double omega_offset_rv, double [:] model_in_pars):
    # pars: [m1,m2,parallax,long_an, e/sqrte_sinomega,to/tc,p,inc,arg_peri/sqrte_cosomega,v1,v2...]
    # model_in_pars: [m1,m2,parallax,long_an,e,to,tc,p,inc,arg_peri,arg_peri_di,arg_peri_rv,a_tot_au,K]
    
    cdef double m1, m2, parallax, long_an, ecc, p, inc, arg_peri, tc
    cdef double sqrte_sinomega, sqrte_cosomega
    cdef double a_tot, top, arg_peri_di, arg_peri_rv
    cdef double ta_temp, half_ea, m_t_tc, delta_t
    
    cdef extern from "math.h":
        double sin(double _x)
        double cos(double _x)
        double atan2(double _x, double _y)
        double sqrt(double _x)
        double pow(double _x, double _y)
    
    ## push pars in array into meaningful variables.
    ########### special conversions for specialty parametrizations ############
    ## Convert sqrt(e)sin(omega)&sqrt(e)cos(omega) => e and omega if required.
    ecc, arg_peri = pars[4], pars[8]
    if low_ecc:
        sqrte_sinomega, sqrte_cosomega = pars[4],pars[8]
        if 0.0 not in [sqrte_sinomega, sqrte_cosomega]:
            ecc = sqrte_sinomega**2  +  sqrte_cosomega**2
            arg_peri = (180.0/const.pi)*atan2(sqrte_sinomega, sqrte_cosomega)        
    ###########################################################################
    [m1, m2, parallax, long_an] = pars[0:4]
    [to, p, inc] = pars[5:8]
    #offsets = pars[9:]
    
    ## update varied omega values to omegaDI and omegaRV here
    # Get the model version of each omega.
    # This is required as the RV values are for the primary and thus omega+pi
    # compared to the value used for the secondary in the DI data.
    ## Another way to think of this is arg_peri_di = arg_peri_companion
    ##                                 arg_peri_rv = arg_peri_primary
    arg_peri_di = arg_peri + omega_offset_di
    arg_peri_rv = arg_peri + omega_offset_rv +180.0
    
    ## Calculate a_tot and a_tot_au
    a_tot = 0.0
    if m1!=0:
        top = p*p*pow(const.sec_per_year,2.0)*const.Grav*const.kg_per_msun*(m1+m2)
        a_tot =pow( top/(4.0*const.pi*const.pi) , (1.0/3.0))
    a_tot_au = a_tot/const.m_per_au
    
    ## set K, Tc/To to defaults, then check if they need to be calculated
    #  properly for use in the RV model.
    k = 0.0
    tc = to
    #if vary_tc:
    #    tc = to
    #else:
    #    to = tc
    if data_mode!="DI":
        ## Calculate K
        #  NOTE: both of the below versions produce identical values of K.
        #        Using 'semi-major version' by default for simplicity/speed.
        # semi-major axis version
        top =  2.0*const.pi * (a_tot/(1.0+(m1/m2))) * sin(inc*(const.pi/180.0)) 
        k =  top / (p*const.sec_per_year*sqrt(1.0-ecc*ecc))
        # masses version
        #cdef double part1, part2, part3
        #part1 = pow( (2.0*const.pi*const.Grav)/(p*const.sec_per_year) , (1.0/3.0) )
        #part2 = pow( (m2*const.kg_per_msun)/((m1+m2)*const.kg_per_msun) , (2.0/3.0) )
        #part3 = sin(inc*(const.pi/180.0)) / sqrt(1.0-ecc*ecc)
        #k = part1 * part2 * part3 
        
        ## Calculate Tc <-> T if needed.  Note, only relevant to RV data.
        if tc_equal_to==False:
            ta_temp = (const.pi/2.0)-arg_peri_rv*(const.pi/180.0)
            half_ea = atan2( sqrt(1.0-ecc)*sin(ta_temp/2.0) , sqrt(1.0+ecc)*cos(ta_temp/2.0) )
            m_t_tc = 2.0*half_ea-ecc*sin(2.0*half_ea);
            delta_t = (m_t_tc*p*const.days_per_year)/(2.0*const.pi);
            if vary_tc:
                # If directly varying Tc, rather than To, then 
                # value in the pars array is Tc, so exchange vars
                to = tc - delta_t
            else:
                tc = to + delta_t
    ## push all calculated values into full list
    #print('m1',type(m1))
    #print('m2',type(m2))
    #print('parallax',type(parallax))
    #print('long_an',type(long_an))
    #print('ecc',type(ecc))
    #print('to',type(to))
    #print('tc',type(tc))
    #print('p',type(p))
    #print('inc',type(inc))
    #print('arg_peri',type(arg_peri))
    #print('arg_peri_di',type(arg_peri_di))
    #print('arg_peri_rv',type(arg_peri_rv))
    #print('a_tot_au',type(a_tot_au))
    #print('k',repr(type(k))+" =  "+str(k))
    
    #for i in range(13):
    #    print(i,' '+repr(model_in_pars[i]))
    #model_in_pars[0:9] = m1, m2, parallax, long_an, ecc, to, tc, p, inc
    #model_in_pars[9:14] = arg_peri, arg_peri_di, arg_peri_rv, a_tot_au, k
    #### Below is a hacky way to do it, but I kept getting the error  ###$$$$$$$$$$$$
    #### due to indexing being wrong.  Fix this !!! #$$$$$
    model_in_pars[0],model_in_pars[1],model_in_pars[2] = m1, m2, parallax 
    model_in_pars[3],model_in_pars[4],model_in_pars[5] = long_an, ecc, to
    model_in_pars[6],model_in_pars[7],model_in_pars[8] = tc, p, inc
    model_in_pars[9],model_in_pars[10] = arg_peri, arg_peri_di
    model_in_pars[11],model_in_pars[12],model_in_pars[13] = arg_peri_rv, a_tot_au, k
    #for i in range(model_in_pars.shape[0]):
    #    print(i,' '+repr(model_in_pars[i]))
    
def orbit(double [:] model_in_pars, double [:] offsets, bint pasa, str data_model,
          double [:] epochs_di, double [:] epochs_rv, int [:] rv_inst_num,
          double [:] rapa_model, double [:] decsa_model, double [:] rv_model, double [:] taea):
    """
    The version of orbit that is called from Python.
    Pass in the parameters and output data array to have loaded up in place.
    This will call fast cdef functions to do the actual calculations.
    """
    cdef double m1, m2, parallax, long_an, ecc, to, tc, p, inc
    cdef double arg_peri, arg_peri_di, arg_peri_rv, a_tot_au, k
    cdef int npts_di, npts_rv
    
    npts_di = epochs_di.shape[0]
    npts_rv = epochs_rv.shape[0]
    
    ## push model input pars in array into meaningful variables.
    [m1, m2, parallax, long_an, ecc, to, tc, p, inc] = model_in_pars[0:9]
    [arg_peri, arg_peri_di, arg_peri_rv, a_tot_au, k] = model_in_pars[9:14]
    
    ## run through each epoch and calc predicted x,y,rv as requested
    if (npts_rv>0) and (data_model!='DI'):       
        ## calc predicted rv if necessary
        orbit_rv(epochs_rv, ecc, to, tc, p, inc, arg_peri_rv, k, offsets, 
                  rv_inst_num,rv_model, taea)
            
    if (npts_di>0) and (data_model!='RV'):
        ## calc predicted DI if necessary
        orbit_di(epochs_di, long_an, ecc, to, tc, p, inc, arg_peri_di, 
                 a_tot_au, parallax, rapa_model, decsa_model, pasa, taea)

#EOF