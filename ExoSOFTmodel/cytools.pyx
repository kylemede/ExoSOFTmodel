import constants as const

#List of available C funcs in math.h:
#https://en.wikipedia.org/wiki/C_mathematical_functions

cdef anomalies(double p, double tc, double to, double ecc, 
               double epoch, double ea, double ta):
    """ 
    Calculate the True and Eccentric Anomalies.
    Remember that in RV, there is occasionally a phase shift due to 
    the Tc!=To, that doesn't exist in DI.
    So, for RV pass in both values, for DI just set both to To.
    
    TA and EA (ta,ea) are calculated in radians.
    """
    cdef double ma, multiples, e_prime
    cdef int newton_count
    cdef bool warnings_on
    
    cdef extern from "math.h":
        double sin(double _x)
        double cos(double _x)
        double acos(double _x)
        double fabs(double _x)
        
    #$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    warnings_on = True  #$$$$ for debugging, kill this after finished testing?
    #$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
        
    ## Calc Mean Anomaly
    ma = (2.0*const.pi*(epoch-2.0*to+tc))/(p*const.days_per_year)
    #Find if over one full orbit and shift into [-2pi,2pi]
    multiples = double(int(ma/2.0*const.pi))
    ma -= multiples*2.0*const.pi
    # shift into [0,2pi]
    if ma<0:
        ma += 2.0*const.pi
    # set E and TA to M by default for circular orbits
    ea, ta = ma, ma
    
    ## if not circular, calculate their specific values
    if ecc>0:
        # double check M is not zero or 2pi
        if ma not in [0.0,2.0*const.pi]:
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
            # calculate TA from E
            ta = acos((cos(ea)-ecc)/(1.0-ecc*cos(ea)))
            # both increase properly 0->180, but as E properly continues to 
            # increase from 180->360, TA decreases.  So correcting for that here.
            if ea>const.pi:
                ta = 2.0*const.pi - ta

cdef orbit_rv(double [:] epochs, double p, double tc, double to, double ecc, 
              double [:] offsets, int [:] rv_inst_num, double k, double arg_peri_rv,
              double [:] rv_model):
    """ 
    Calculates the predicted rv for a each epoch in epochs array.
     
    model value = calculated rv + the instrument specific offset
    """
    cdef double top, ea, ta, ta_temp, half_ea, m_t_tc, delta_t, rv, arg_peri_rad
    cdef int npts
    cdef extern from "math.h":
        double cos(double _x)
        
    npts = epochs.shape[0]
    arg_peri_rad = arg_peri_rv*(const.pi/180.0)
    
    ## calculate rv for each epoch in the data
    for i in range(npts):
        ## Call anomalies        
        anomalies(p, tc, to, ecc, epochs[i], ea, ta)
        ## Calc predicted RV + the instrument specific offset
        rv = k*( cos(ta+arg_peri_rad) + ecc*cos(arg_peri_rad) )
        rv_model[i] = rv + offsets[rv_inst_num[i]]
        
cdef orbit_di(double [:] epochs, double p, double ecc, double to, double tc,
              double arg_peri_di, double long_an,double a_tot_au, double parallax,
              double [:] rapa_model, double [:] decsa_model, bool pasa):                 
    """
    Calculates the predited x/RA/PA and y/Dec/SA for a single epoch.
    """
    cdef double ea, ta, a_thi, b_thi, f_thi, g_thi, x_thi, y_thi, ra, dec, pa, sa
    cdef doubl long_an_rad, arg_peri_rad, inc_rad, sep_arcsec
    cdef int npts
    cdef extern from "math.h":
        double sin(double _x)
        double cos(double _x)
        double sqrt(double _x)
        double atan2(double _x,double _y)
    
    npts = epochs.shape[0]
    long_an_rad = long_an*(const.pi/180.0)
    arg_peri_rad = arg_peri_di*(const.pi/180.0)
    inc_rad = inc*(const.pi/180.0)
    sep_arcsec = a_tot_au*(parallax/1000.0)
    
    ## Calculate predicted x/RA/PA and y/Dec/SA values for each epoch in the data
    for i in range(npts):
        ## Call anomalies        
        anomalies(p, tc, to, ecc, epochs[i], ea, ta)
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
        
def model_input_pars(double [:] pars, bool low_ecc, bool tc_equal_to, 
                     bool di_only, double omega_offset_di, double omega_offset_rv,
                     double [:] offsets, double [:] model_in_pars):
    # pars: [m1,m2,parallax,long_an, e/sqrte_sinomega,to/tc,p,inc,arg_peri/sqrte_cosomega,v1,v2...]
    # model_in_pars: [m1,m2,parallax,long_an,e,to,tc,p,inc,arg_peri,arg_peri_di,arg_peri_rv,a_tot_au,K]
    
    cdef double m1, m2, parallax, long_an, ecc, p, inc, arg_peri
    cdef double sqrte_sinomega, sqrte_cosomega
    cdef double  a_tot, top, arg_peri_di, arg_peri_rv
    
    cdef extern from "math.h":
        double sin(double _x)
        double atan2(double _x)
        double sqrt(double _x)
        double pow(double _x, double _y)
        
    ## push pars in array into meaningful variables.
    ########### special conversions for specialty parametrizations ############
    ## Convert sqrt(e)sin(omega)&sqrt(e)cos(omega) => e and omega if required.
    if low_ecc:
        sqrte_sinomega, sqrte_cosomega = pars[4],pars[8]
        if 0.0 not in [sqrte_sinomega, sqrte_cosomega]:
            ecc = sqrte_sinomega**2  +  sqrte_cosomega**2
            arg_peri = (180.0/const.pi)*atan2(sqrte_sinomega, sqrte_cosomega)
    else:
        ecc, arg_peri = pars[4], pars[8]
    ###########################################################################
    [m1, m2, parallax, long_an] = pars[0:4]
    [to, p, inc] = pars[5:8]
    offsets = pars[9:]
    
    ## update varied omega values to omegaDI and omegaRV here
    # Get the model version of each omega.
    # This is required as the RV values are for the primary and thus omega+pi
    # compared to the value used for the secondary in the DI data.
    arg_peri_di = arg_peri + omega_offset_di
    arg_peri_rv = arg_peri + omega_offset_rv
    
    ## Calculate a_tot and a_tot_au
    if m1!=0:
        top = p*p*pow(const.sec_per_year,2.0)*const.Grav*const.kg_per_msun*(m1+m2)
        a_tot =pow( top/(4.0*const.pi*const.pi) , (1.0/3.0))
    a_tot_au = a_tot/const.m_per_au
    
    ## set K, Tc/To to defaults, then check if they need to be calculated
    #  properly for use in the RV model.
    k = 0.0
    if varytc:
        tc = to
    else:
        to = tc
    if di_only==False:
        ## Calculate K
        #  NOTE: both of the below versions produce identical values of K.
        #        Using 'semi-major version' by default for simplicity/speed.
        # semi-major axis version
        top =  2.0*const.pi * (a_tot/(1.0+(m1/m2))) * sin(inc*(const.pi/180.0) 
        k =  top / (p*const.sec_per_year*sqrt(1.0-ecc*ecc)
        # masses version
        #cdef double part1, part2, part3
        #part1 = pow( (2.0*const.pi*const.Grav)/(p*const.sec_per_year) , (1.0/3.0) )
        #part2 = pow( (m2*const.kg_per_msun)/((m1+m2)*const.kg_per_msun) , (2.0/3.0) )
        #part3 = sin(inc*(const.pi/180.0)) / sqrt(1.0-ecc*ecc)
        #k = part1 * part2 * part3 
        
        ## Calculate Tc <-> T if needed.  Note, only relevant to RV data.
        if tc_equal_to==False:
            cdef double ta_temp, half_ea, m_t_tc, delta_t
            ta_temp = (const.pi/2.0)-arg_peri_rv*(const.pi/180.0)
            half_ea = atan2( sqrt(1.0-ecc)*sin(ta_temp/2.0) , sqrt(1.0+ecc)*cos(ta/2.0) )
            m_t_tc = 2.0*half_ea-ecc*sin(2.0*half_ea);
            delta_t = (m_t_tc*p*const.days_per_year)/(2.0*const.pi);
            if vary_tc:
                # If directly varying Tc, rather than To, then 
                # value in the pars array is Tc, so exchange vars
                tc = to
                to = tc - delta_t
            else:
                tc = to + delta_t
    ## push all calculated values into full list
    model_in_pars[0:9] = [m1, m2, parallax, long_an, ecc, to, tc, p, inc]
    model_in_pars[9:14] = [arg_peri, arg_peri_di, arg_peri_rv, a_tot_au, k]
    

def orbit(double [:] model_in_pars, double [:] offsets, bool pasa,
          double [:] epochs_di, double [:] epochs_rv, int [:] rv_inst_num,
          double [:] rapa_model, double [:] decsa_model, double [:] rv_model):
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
    for epoch in epochs:
        if npts_rv>0:       
            ## calc predicted rv if necessary
            orbit_rv(epochs_rv, p, tc, to, ecc, offsets, rv_inst_num, k,
                      arg_peri_rv, rv_model)
                
        if npts_di>0:
            ## calc predicted DI if necessary
            orbit_di(epochs_di, p, ecc, to, tc, arg_peri_di, long_an, 
                     a_tot_au, parallax, rapa_model, decsa_model, pasa)

#EOF