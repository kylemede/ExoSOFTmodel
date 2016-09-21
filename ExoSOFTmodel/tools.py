#@Author: Kyle Mede, kylemede@astron.s.u-tokyo.ac.jp or kylemede@gmail.com
from __future__ import absolute_import
import numpy as np
import KMlogger

log = KMlogger.getLogger('main.tools',lvl=100,addFH=False)

def load_di_data(filename):
    """
    Load the astrometry data into a numpy array.
    
    file format:
    title 
    column headers
    data
    .
    .
    .
    
    Data must be in the columns:
    obsDate[JD] x["] x_error["] y["] y_error["]
    OR if pasa key in settingsAdvanced ==True, then:
    obsDate[JD] PA[deg] PA_error[deg] SA["] SA_error["]
    """
    epochs_di = []
    rapa = []
    rapa_err = []
    decsa = []
    decsa_err = []
    try:
        if filename[-4:]!='.dat':
            filename = filename+'.dat'
        fl = open(filename, 'r')
        lines = fl.readlines()
        fl.close()
        for line in lines:
            #print "line was:'"+line+"'"
            #log.debug("line was:'"+line+"'")#$$$$$$$$$$$$$$$$$$$$$$$$
            if len(line.split())>2:
                if line.split()[0].replace('.','',1).isdigit() and line.split()[3].replace('.','',1).replace('-','',1).isdigit():
                    epochs_di.append(float(line.split()[0]))
                    rapa.append(float(line.split()[1]))
                    rapa_err.append(float(line.split()[2]))
                    decsa.append(float(line.split()[3]))
                    decsa_err.append(float(line.split()[4]))
                    #print repr([float(line.split()[0]),float(line.split()[1]),float(line.split()[2]),float(line.split()[3]),float(line.split()[4])])
        epochs_di = np.array(epochs_di,dtype=np.dtype('d'))
        rapa = np.array(rapa,dtype=np.dtype('d'))
        rapa_err = np.array(rapa_err,dtype=np.dtype('d'))
        decsa = np.array(decsa,dtype=np.dtype('d'))
        decsa_err = np.array(decsa_err,dtype=np.dtype('d'))
    except:
        log.critical("a problem occured while trying to load DI data. \nPlease check it is formatted correctly.")
    return (epochs_di, rapa, rapa_err, decsa, decsa_err)


def load_rv_data(filename):
    """
    Load the radial velocity data into a numpy array.  Provided jitter values 
    will be added in quadrature with the errors.
    
    file format:
    title 
    column headers
    data
    .
    .
    Empty line between data sets
    data
    .
    .
    
    Data must be in the columns:
    obsDate[JD] RV[m/s] RV_error[m/s] jitter[m/s] datasetNumber[int]
    NOTE: datasetNumber is optional, if not provided they will be automatically set to 0,1,2... following the order of the data in the file.
          If jitter is not provided, it will be assumed zero.
    """
    epochs_rv = []
    rv = []
    rv_err = []
    rv_inst_num = []
    try:
        if filename[-4:]!='.dat':
            filename = filename+'.dat'
        fl = open(filename, 'r')
        lines = fl.readlines()
        fl.close()
        datasetNumLast = 0
        jitterLast = 0
        lastWasDataLine=False
        thisIsDataLine = False
        for line in lines:
            #print "line = "+line  #$$$$$$$$$$$$$$$$$
            lastWasDataLine=thisIsDataLine
            thisIsDataLine=False
            if len(line.split())>2:
                if line.split()[0].replace('.','',1).isdigit() and line.split()[1].replace('.','',1).replace('-','',1).isdigit():
                    thisIsDataLine=True
                    epochs_rv.append(float(line.split()[0]))
                    rv.append(float(line.split()[1]))
                    #if jitter was provided on first line of data set
                    if len(line.split())>3:
                        try:
                            jitterLast = float(line.split()[3])
                        except:
                            log.error("could not convert 4th element of split into jitter.  4th element was: "+str(line.split()[3]))
                    rv_err.append(np.sqrt(float(line.split()[2])**2+jitterLast**2))
                    #if datasetNum was provided on first line of data set
                    if len(line.split())>4:
                        try:
                            datasetNumLast = float(line.split()[4])
                        except:
                            log.error("could not convert 5th element of split into datasetNum.  5th element was: "+str(line.split()[4]))
                    rv_inst_num.append(datasetNumLast)
            if lastWasDataLine and (thisIsDataLine==False):
                jitterLast = 0
                datasetNumLast+=1
                #print 'incrementing datasetNum'  #$$$$$$$$$$$$$$$$$
        epochs_rv = np.array(epochs_rv,dtype=np.dtype('d'))
        rv = np.array(rv,dtype=np.dtype('d'))
        rv_err = np.array(rv_err,dtype=np.dtype('d'))
        rv_inst_num = np.array(rv_inst_num,dtype=np.dtype('i'))
    except:
        log.critical("a problem occured while trying to load RV data.  \nPlease check it is formatted correctly.")
    return (epochs_rv, rv, rv_err, rv_inst_num)

#EOF