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
    
    def price(self, strike, spot, texp=None, sigma=None, cp=1,random_seed = 12345):
        '''
        Your MC routine goes here
        Generate paths for vol and price first. Then get prices (vector) for all strikes
        You may fix the random number seed
        '''
        n_interval = 100
        n_iter =10000
        delta_t = texp/n_interval
        prices = []
        np.random.seed(random_seed)
        # get the whole price path and sigma path every iteration
        for i in range(n_iter):
            z1 = np.random.randn(n_interval)
            z2 = np.random.randn(n_interval)
            w1 = self.rho*z1 + np.sqrt(1-np.power(self.rho,2))*z2
            sis = np.exp(self.vov*np.sqrt(delta_t)*z1-0.5*np.power(self.vov,2)*delta_t)
            sis[1:]=sis[:-1]
            sis[0] = self.sigma
            sis = np.cumprod(sis)
            deltap = np.exp(sis*np.sqrt(delta_t)*w1-0.5*np.power(sis,2)*delta_t)
            deltap[0]*=spot
            pts = np.cumprod(deltap)
            prices.append(pts[-1])
        strikes = np.array([strike]*n_iter).T
        callp = -strikes+prices
        callp = np.where(callp>0,callp,0)
        # record the call price among all our MC
        self.cprice_paths = callp
        finalp = callp.mean(axis = 1)
        return finalp

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
        
    def price(self, strike, spot, texp=None, sigma=None, cp=1, random_seed = 12345):
        '''
        Your MC routine goes here
        Generate paths for vol and price first. Then get prices (vector) for all strikes
        You may fix the random number seed
        '''
        n_interval = 100
        n_iter =10000
        delta_t = texp/n_interval
        prices = []
        np.random.seed(random_seed)
        # get the whole price path and sigma path every iteration
        for i in range(n_iter):
            z1 = np.random.randn(n_interval)
            z2 = np.random.randn(n_interval)
            w1 = self.rho*z1 + np.sqrt(1-np.power(self.rho,2))*z2
            sis = np.exp(self.vov*np.sqrt(delta_t)*z1-0.5*np.power(self.vov,2)*delta_t)
            sis[1:]=sis[:-1]
            sis[0] = self.sigma
            sis = np.cumprod(sis)
            deltap = sis*np.sqrt(delta_t)*w1
            deltap[0]+=spot
            pts = np.cumsum(deltap)
            prices.append(pts[-1])
        strikes = np.array([strike]*n_iter).T
        callp = -strikes+prices
        callp = np.where(callp>0,callp,0)
         # record the call price among all our MC
        self.cprice_paths = callp
        finalp = callp.mean(axis = 1)
        return finalp

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
    
    def price(self, strike, spot, texp=None, cp=1,random_seed = 12345):
        '''
        Your MC routine goes here
        Generate paths for vol only. Then compute integrated variance and BSM price.
        Then get prices (vector) for all strikes
        You may fix the random number seed
        '''
        n_interval = 100
        n_iter =10000
        delta_t = texp/n_interval
        prices = []
        np.random.seed(random_seed)
         # get a whole sigma path every iteration
        for i in range(n_iter):
            z1 = np.random.randn(n_interval)
            sis = np.exp(self.vov*np.sqrt(delta_t)*z1-0.5*np.power(self.vov,2)*delta_t)
            sis[1:]=sis[:-1]
            sis[0] = self.sigma
            sis = np.cumprod(sis)
            var = np.power(sis,2)/sis[0]**2
            it = var.mean()
            s0 = spot*np.exp(self.rho/self.vov*(sis[-1]-sis[0])-0.5*np.power(self.rho*sis[0],2)*texp*it)
            sigma_bs = sis[0]*np.sqrt((1-self.rho**2)*it)
            prices.append(bsm.price(strike,s0,texp,sigma_bs))
        prices = np.array(prices)
         # record the call price among our CMC
        self.cprice_paths = prices
        finalp = prices.mean(axis = 0)
        return finalp
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
        
    def price(self, strike, spot, texp=None, cp=1,random_seed =12345):
        '''
        Your MC routine goes here
        Generate paths for vol only. Then compute integrated variance and normal price.
        You may fix the random number seed
        '''
        n_interval = 100
        n_iter =10000
        delta_t = texp/n_interval
        prices = []
        np.random.seed(random_seed)
         # get a whole sigma path every iteration
        for i in range(n_iter):
            z1 = np.random.randn(n_interval)
            sis = np.exp(self.vov*np.sqrt(delta_t)*z1-0.5*np.power(self.vov,2)*delta_t)
            sis[1:]=sis[:-1]
            sis[0] = self.sigma
            sis = np.cumprod(sis)
            var = np.power(sis,2)/sis[0]**2
            it = var.mean()
            s0 = spot+self.rho/self.vov*(sis[-1]-sis[0])
            sigma_nm = sis[0]*np.sqrt((1-self.rho**2)*it)
            prices.append(normal.price(strike,s0,texp,sigma_nm))
        prices = np.array(prices)
        # record all the prices among our CMC
        self.cprice_paths = prices
        finalp = prices.mean(axis = 0)
        return finalp
       