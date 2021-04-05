# -*- coding: utf-8 -*-
"""
Created on Tue Sep 19 22:56:58 2017

@author: jaehyuk
"""
import numpy as np
import scipy.stats as ss
import scipy.optimize as sopt

def normal_formula(strike, spot, vol, texp, intr=0.0, divr=0.0, cp=1):
    div_fac = np.exp(-texp*divr)
    disc_fac = np.exp(-texp*intr)
    forward = spot / disc_fac * div_fac

    if( texp<0 or vol*np.sqrt(texp)<1e-8 ):
        return disc_fac * np.fmax( cp*(forward-strike), 0 )

    vol_std = np.fmax(vol * np.sqrt(texp), 1.0e-16)
    d = (forward - strike) / vol_std

    price = disc_fac * (cp * (forward - strike) * ss.norm.cdf(cp * d) + vol_std * ss.norm.pdf(d))
    return price

class NormalModel:
    vol, intr, divr = None, None, None
    
    def __init__(self, vol, intr=0, divr=0):
        self.vol = vol
        self.intr = intr
        self.divr = divr
    
    def price(self, strike, spot, texp, cp=1):
        return normal_formula(strike, spot, self.vol, texp, intr=self.intr, divr=self.divr, cp=cp)
    
    def delta(self, strike, spot, vol, texp, intr=0.0, divr=0.0, cp=1):
        ''' 
        <-- PUT your implementation here
        
        '''
        vol_std = np.fmax(vol * np.sqrt(texp), 1.0e-16)
        dn = (spot-strike)/vol_std
        delta_re = cp* ss.norm.cdf(cp * dn)
        return delta_re

    def vega(self, strike, spot, vol, texp, intr=0.0, divr=0.0, cp=1):
        ''' 
        <-- PUT your implementation here
        '''
        vol_std = np.fmax(vol * np.sqrt(texp), 1.0e-16)
        dn = (spot-strike)/vol_std
        vega_re = np.sqrt(texp) * ss.norm.pdf(dn)
        return vega_re

    def gamma(self, strike, spot, vol, texp, intr=0.0, divr=0.0, cp=1):
        ''' 
        <-- PUT your implementation here
        '''
        vol_std = np.fmax(vol * np.sqrt(texp), 1.0e-16)
        dn = (spot-strike)/vol_std
        gamma_re = ss.norm.pdf(dn)/(vol_std *np.sqrt(texp))
        return gamma_re

    def impvol(self, price, strike, spot, texp, cp=1):
        ''' 
        <-- PUT your implementation here
        '''
        func = lambda _vol: normal_formula(strike, spot, _vol, texp, self.intr, self.divr, cp) - price
        imp_vol = sopt.brentq(func, 0, 10)
        return imp_vol