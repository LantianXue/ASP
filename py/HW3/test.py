# -*- coding: utf-8 -*-
"""
Created on Tue Apr 13 13:27:38 2021

@author: Lantian
"""
import numpy as np
import matplotlib.pyplot as plt
import pyfeng as pf
import option_models as opt


# Parameters
strike = np.linspace(75,125,num=25)
forward = 100
sigma = 0.2
texp = 1
vov = 0.5
rho = 0.25
beta = 1

# Create model
sabr_bsm = pf.SabrHagan2002(sigma, vov=vov, rho=rho, beta=beta)
#sabr_bsm.__dict__
# sabr_bsm_mc = opt.sabr.ModelBsmMC(sabr_bsm.sigma, vov=sabr_bsm.vov, rho=sabr_bsm.rho, beta=1)

# price_hagan = sabr_bsm.price(strike, forward, texp)
# price_mc = sabr_bsm_mc.price(strike, forward, texp)

# # make sure the two prices are similar
# print(price_hagan)
# print(price_mc)

# sabr_nm = pf.SabrHagan2002(sigma, vov=vov, rho=rho, beta=0)
# #sabr_bsm.__dict__
# sabr_nm_mc = opt.sabr.ModelNormalMC(sabr_bsm.sigma, vov=sabr_bsm.vov, rho=sabr_bsm.rho, beta=0)

# price_hagan = sabr_nm.price(strike, forward, texp)
# price_mc = sabr_nm_mc.price(strike, forward, texp)

# # make sure the two prices are similar
# print(price_hagan)
# print(price_mc)

# instantiate mc model from the hagan model's parameters
sabr_bsm_cmc = opt.sabr.ModelBsmCondMC(sabr_bsm.sigma, vov=sabr_bsm.vov, rho=sabr_bsm.rho, beta=1)
price_hagan = sabr_bsm.price(strike, forward, texp)
price_mc = sabr_bsm_cmc.price(strike, forward, texp)

# make sure the two prices are similar
print(price_hagan)
print(price_mc)

