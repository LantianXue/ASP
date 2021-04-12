    # -*- coding: utf-8 -*-
"""
Created on Tue Oct 10

@author: jaehyuk
"""

import numpy as np
import scipy.stats as ss
import scipy.optimize as sopt
from . import normal
from . import bsm
import pyfeng as pf

'''
MC model class for Beta=1
'''
class ModelBsmMC:
    beta = 1.0   # fixed (not used)
    vov, rho = 0.0, 0.0
    sigma, intr, divr = None, None, None
    bsm_model = None
    '''
    You may define more members for MC: time step, etc
    '''
    
    def __init__(self, sigma, vov=0, rho=0.0, beta=1.0, intr=0, divr=0):
        self.sigma = sigma
        self.vov = vov
        self.rho = rho
        self.intr = intr
        self.divr = divr
        self.bsm_model = pf.Bsm(sigma, intr=intr, divr=divr)
        
    def bsm_vol(self, strike, spot, texp=None, sigma=None):
        ''''
        From the price from self.price() compute the implied vol
        this is the opposite of bsm_vol in ModelHagan class
        use bsm_model
        '''
        return 0
    
    def price(self, strike, spot, texp=None, sigma=None, cp=1):
        '''
        Your MC routine goes here
        Generate paths for vol and price first. Then get prices (vector) for all strikes
        You may fix the random number seed
        '''
        n_inter = 100000
        delta_t = texp/n_inter
        np.random.seed(12345)
        z1 = np.random.normal(0,1,n_inter)
        np.random.seed(4567)
        z2 = np.random.normal(0,1,n_inter)
        w1 = self.rho*z1 + np.sqrt(1-np.power(self.rho,2))*z2
        prices = []
        for s in strike:
            p = np.log(spot)
            cp = []
            si = self.sigma
            for i in range(1,n_inter+1):
                p = p + si * w1[i-1]* np.sqrt(delta_t) - 0.5*np.power(si,2)*delta_t
                si = si * np.exp(self.vov*np.sqrt(delta_t)*z1[i-1]-0.5*np.power(self.vov,2)*delta_t)
                c = max(np.exp(p)-s,0)
                cp.append(c)
            price = np.mean(cp)
            prices.append(price)
        prices = np.array(prices)
        return prices

'''
MC model class for Beta=0
'''
class ModelNormalMC:
    beta = 0.0   # fixed (not used)
    vov, rho = 0.0, 0.0
    sigma, intr, divr = None, None, None
    normal_model = None
    
    def __init__(self, sigma, vov=0, rho=0.0, beta=0.0, intr=0, divr=0):
        self.sigma = sigma
        self.vov = vov
        self.rho = rho
        self.intr = intr
        self.divr = divr
        self.normal_model = pf.Norm(sigma, intr=intr, divr=divr)
        
    def norm_vol(self, strike, spot, texp=None, sigma=None):
        ''''
        From the price from self.price() compute the implied vol
        this is the opposite of normal_vol in ModelNormalHagan class
        use normal_model 
        '''
        return 0
        
    def price(self, strike, spot, texp=None, sigma=None, cp=1):
        '''
        Your MC routine goes here
        Generate paths for vol and price first. Then get prices (vector) for all strikes
        You may fix the random number seed
        '''
        n_inter = 10000
        delta_t = texp/n_inter
        np.random.seed(12345)
        z1 = np.random.normal(0,1,n_inter)
        np.random.seed(4567)
        z2 = np.random.normal(0,1,n_inter)
        w1 = self.rho*z1 + np.sqrt(1-np.power(self.rho,2))*z2
        prices = []
        for s in strike:
            p = spot
            cp = []
            si = self.sigma
            for i in range(n_inter+1):
                p = p + si * w1[i-1]* np.sqrt(delta_t)
                si = si * np.exp(self.vov*np.sqrt(delta_t)*z1[i-1]-0.5*np.power(self.vov,2)*delta_t)
                c = max(p-s,0)
                cp.append(c)
            price = np.mean(cp)
            prices.append(price)
        prices = np.array(prices)
        return prices

'''
Conditional MC model class for Beta=1
'''
class ModelBsmCondMC:
    beta = 1.0   # fixed (not used)
    vov, rho = 0.0, 0.0
    sigma, intr, divr = None, None, None
    bsm_model = None
    '''
    You may define more members for MC: time step, etc
    '''
    
    def __init__(self, sigma, vov=0, rho=0.0, beta=1.0, intr=0, divr=0):
        self.sigma = sigma
        self.vov = vov
        self.rho = rho
        self.intr = intr
        self.divr = divr
        self.bsm_model = pf.Bsm(sigma, intr=intr, divr=divr)
        
    def bsm_vol(self, strike, spot, texp=None):
        ''''
        From the price from self.price() compute the implied vol
        this is the opposite of bsm_vol in ModelHagan class
        use bsm_model
        should be same as bsm_vol method in ModelBsmMC (just copy & paste)
        '''
        return 0
    
    def price(self, strike, spot, texp=None, cp=1):
        '''
        Your MC routine goes here
        Generate paths for vol only. Then compute integrated variance and BSM price.
        Then get prices (vector) for all strikes
        You may fix the random number seed
        '''
        n_inter = 100000
        delta_t = texp/n_inter
        np.random.seed(12345)
        z1 = np.random.normal(0,1,n_inter)
        sis = []
        sis.append(self.sigma)
        si = self.sigma
        for i in range(1,n_inter+1):
            si =  si * np.exp(self.vov*np.sqrt(delta_t)*z1[i-1]-0.5*np.power(self.vov,2)*delta_t)
            sis.append(si)
        it = 0
        for i in sis:
            it += 2*np.power(i,2)
        it = it-sis[0]-sis[-1]+2*sis[1]+2*sis[-2]
        it = it/6/n_inter
        spot_now = spot*np.exp(self.rho*(sis[-1]-sis[0])/self.vov-0.5*np.power(self.rho*sis[0],2)*texp*it)
        vol_now = sis[0]*np.sqrt((1-np.power(self.rho,2))*it)
        prices = bsm.price(strike,spot_now,texp,vol_now)
        return prices
'''
Conditional MC model class for Beta=0
'''
class ModelNormalCondMC:
    beta = 0.0   # fixed (not used)
    vov, rho = 0.0, 0.0
    sigma, intr, divr = None, None, None
    normal_model = None
    
    def __init__(self, sigma, vov=0, rho=0.0, beta=0.0, intr=0, divr=0):
        self.sigma = sigma
        self.vov = vov
        self.rho = rho
        self.intr = intr
        self.divr = divr
        self.normal_model = pf.Norm(sigma, intr=intr, divr=divr)
        
    def norm_vol(self, strike, spot, texp=None):
        ''''
        From the price from self.price() compute the implied vol
        this is the opposite of normal_vol in ModelNormalHagan class
        use normal_model
        should be same as norm_vol method in ModelNormalMC (just copy & paste)
        '''
        return 0
        
    def price(self, strike, spot, texp=None, cp=1):
        '''
        Your MC routine goes here
        Generate paths for vol only. Then compute integrated variance and normal price.
        You may fix the random number seed
        '''
        n_inter = 10000
        delta_t = texp/n_inter
        np.random.seed(12345)
        z1 = np.random.normal(0,1,n_inter)
        sis = []
        sis.append(self.sigma)
        si = self.sigma
        for i in range(1,n_inter+1):
            si =  si * np.exp(self.vov*np.sqrt(delta_t)*z1[i-1]-0.5*np.power(self.vov,2)*delta_t)
            sis.append(si)
        it = 0
        for i in sis:
            it += 2*np.power(i,2)
        it = it-sis[0]-sis[-1]
        it = it/2/n_inter
        spot_now = spot+self.rho*(sis[-1]-sis[0])/self.vov
        vol_now = sis[0]*np.sqrt((1-np.power(self.rho,2))*it)
        prices = normal.price(strike,spot_now,texp,vol_now)
        return prices