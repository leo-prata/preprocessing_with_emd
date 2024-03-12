import os, sys, copy
project_directory = os.getcwd()
sys.path.append(project_directory)

import numpy as np
import emd
from Models.Preprocessing._Generic import generic
from scipy.interpolate import CubicSpline

def _resample_mb(eegdata_, resample):
    duration = len(eegdata_['times'])/eegdata_['sfreq']
    new_times = np.arange(0, duration, 1./resample) if resample is not None else eegdata_['times']

    X = []
    for trial_ in range(eegdata_['X'].shape[0]):
        X.append([])
        for band_ in range(eegdata_['X'].shape[1]):
            X[-1].append([])
            for electrode_ in range(eegdata_['X'].shape[2]):
                cubic_spline = CubicSpline(eegdata_['times'], eegdata_['X'][trial_][band_][electrode_])
                new_signal = cubic_spline(new_times)
                X[-1][-1].append(new_signal)

    eegdata_ = copy.deepcopy(eegdata_)
    eegdata_['sfreq'] = resample
    eegdata_['times'] = new_times
    eegdata_['X'] = np.array(X)
    return eegdata_

class EMD(generic):
        
        def __init__(self, max_imfs=9, resample=128, verbose=False):
            self.verbose = verbose
            self.max_imfs = max_imfs
            self.resample = resample
        
        def transform(self, eegdata):
            imfs = []
            max_imfs = np.inf

            for trial_ in range(eegdata['X'].shape[0]):
                imfs.append([])
                for electrode_ in range(eegdata['X'].shape[1]): #vai passar em cada eletrodo
                    imfs[-1].append(self._emd(eegdata['X'][trial_][electrode_])) #vai mandar pra funcao emd o tempo do eletrodo em cada iteracao
                    if imfs[-1][-1].shape[1] < max_imfs:
                        max_imfs = imfs[-1][-1].shape[1]
                #v = [[imf][imf2][imf3]]

            imfs_ = []
            max_imfs -= 1
            if self.max_imfs is not None:
                if self.max_imfs > max_imfs:
                    self.max_imfs = max_imfs 
            
            for trial_ in range(len(imfs)):
                imfs_.append([])
                for electrode_ in range(len(imfs[trial_])):
                     imfs_[-1].append(imfs[trial_][electrode_][:, :self.max_imfs])      
            
            imfs = np.array(imfs_)
            imfs = np.transpose(imfs, (0, 3, 1, 2))
            
            eegdata = copy.deepcopy(eegdata)
            eegdata['X'] = imfs

            if self.resample is not None:
                eegdata = _resample_mb(eegdata, self.resample)
            return eegdata

        def fit_transform(self, X):
            newX = self.transform(X)
            return newX
        
        def get_params(self, deep):
            return {'verbose': self.verbose, 'inplace': self.inplace}
        
        def _emd(self, X_eletrodo):
            imfs = emd.sift.sift(X_eletrodo, max_imfs=None)
            return imfs

if __name__  == "__main__":
    
    from Datasets.PhysionetMI import Physionet
    from Datasets.ReadFile import write, read

    for sub in range(1, 109+1):
        try:
            print(f'Starting subject number {sub}.', end='\r')
            #eegdata = Physionet().load(subject=sub, events_dict={'left-hand': 0, 'right-hand': 1})
            eegdata = read('__temp__\\eafirst\\phy_%03d.ea'%sub)
            emd_ = EMD(resample=128)
            eegdata = emd_.transform(eegdata)
            #write(eegdata, '__temp__/phy_%03d'%sub)
            write(eegdata, '__temp__\\eafirst\\phy_%03d.emd'%sub)
            print(f'Subject number {sub} saved.    ')
        except:
            a=0


#    from Datasets.CBCIC import CBCIC
#
#    eegdata = CBCIC().load()
#    emd_ = EMD()
#    eegdata = emd_.transform(eegdata)
#    print(eegdata['X'].shape)
#
#        #[trial][eletrodo][tempo]