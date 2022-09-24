# Modify the correlation of risky return and labor income

import sys
from xml.dom.minicompat import NodeList
from scipy.stats import multivariate_normal as mvn
import numpy as np
from tauchen import *
# from useful import *
# from policy18 import *
# from simulation18 import *
import argparse

import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.rcParams["lines.linewidth"] = 1
mpl.rcParams["savefig.bbox"] = "tight"
mpl.rcParams["figure.figsize"] = (26, 15)
mpl.rcParams["font.size"] = 15
mpl.rcParams["legend.frameon"] = False


def MUc(x):
    return x**(-rhopar)


def inv_MU(u):
    return u**(-1/rhopar)


def U(x):
    return x**(1-rhopar)/(1-rhopar)


# Grid

def inter1d(x1, y1, y2):
    return ((1.0 - x1) * (y1) + (x1) * (y2))


def getwage(rrate):
    return (1.0 - alphapar) * (alphapar / (rrate + deltapar))**(alphapar / (1.0 - alphapar))


def getlevel(x):
    return (scale1 * (np.exp(exponen * (x)) + grmin))


def getomega(x):
    return ((x + 0.0) / size_portfoliochoice)


def getgrid(x):
    return (np.log((x) / scale1 - grmin) / exponen)


# EGM derivatives

def nderiv(val1, val2, val3, x1, x2, x3):
    return ((1.0 - (x3 - x2) / (x3 - x1)) * ((val3 - val2) / (x3 - x2)) + ((x3 - x2) / (x3 - x1)) * ((val2 - val1) / (x2 - x1)))


def POLICY(VF_final,  dVF_final,  save_final,  VF,  dVF,  save,  Portfolio,  K,  Omega,  wagerate):

    # INITIALIZATION #

    VFendo = np.zeros((size_asset, size_risk, size_risk, size_laborincome,
                      size_laborincome, size_shock, size_shock, size_shock, size_shock, size_portfoliochoice+1))
    cohendo = np.zeros((size_asset, size_risk, size_risk, size_laborincome,
                       size_laborincome, size_shock, size_shock, size_shock, size_shock, size_portfoliochoice+1))

    VF_final_old = np.zeros((size_asset, size_risk, size_risk, size_laborincome,
                            size_laborincome, size_shock, size_shock, size_shock, size_shock))

    iter = 0

    critV = 10000.0
    # std::cout << "iter\t"
    #           << "critV\n"

    while critV > epsV and iter < 250:
        # we need copy to make a separate object

        VFendo = np.zeros((size_asset, size_risk, size_risk, size_laborincome,
                          size_laborincome, size_shock, size_shock, size_shock, size_shock, size_portfoliochoice+1))
        cohendo = np.zeros((size_asset, size_risk, size_risk, size_laborincome,
                           size_laborincome, size_shock, size_shock, size_shock, size_shock, size_portfoliochoice+1))
        VF_final_old = VF_final

        # std::cout << std::setprecision(16) << VF[(5, 5)] << "\n"
        # std::cout << std::setprecision(16) << VFnew[(5, 5)] << "\n"

        # main EGM computation
        for asset_index in list(range(size_asset)):
            for risk_index in list(range(size_risk)):
                for risk_pre_index in list(range(size_risk)):
                    for laborincome_index in list(range(size_laborincome)):
                        for laborincome_pre_index in list(range(size_laborincome)):
                            for riskshock_index in list(range(size_shock)):
                                for riskshock_pre_index in list(range(size_shock)):
                                    for laborshock_index in list(range(size_shock)):
                                        for laborshock_pre_index in list(range(size_shock)):
                                            for portfoliochoice_index in list(range(size_portfoliochoice+1)):

                                                tempnext = 0
                                                dtempnext = 0

                                                for risk_indexnext in list(range(size_risk)):
                                                    for laborincome_indexnext in list(range(size_laborincome)):
                                                        for riskshock_indexnext in list(range(size_shock)):
                                                            for laborshock_indexnext in list(range(size_shock)):
                                                                shockstate_next = riskshock_indexnext*size_shock+laborshock_indexnext
                                                                tempnext += risk_trans[risk_index][risk_indexnext] * laborincome_trans[laborincome_index][laborincome_indexnext] * risklaborshock_trans[0, shockstate_next] *\
                                                                    VF[asset_index, risk_indexnext, risk_index,
                                                                        laborincome_indexnext, laborincome_index, riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index, portfoliochoice_index]
                                                                dtempnext += risk_trans[risk_index][risk_indexnext] * laborincome_trans[laborincome_index][laborincome_indexnext] * risklaborshock_trans[0, shockstate_next] *\
                                                                    dVF[asset_index, risk_indexnext, risk_index,
                                                                        laborincome_indexnext, laborincome_index, riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index, portfoliochoice_index]

                                                cohendo[asset_index, risk_index, risk_pre_index, laborincome_index,
                                                        laborincome_pre_index, riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index, portfoliochoice_index] = K[asset_index] + inv_MU(betapar * dtempnext)
                                                VFendo[asset_index, risk_index, risk_pre_index, laborincome_index, laborincome_pre_index, riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index, portfoliochoice_index] = U(
                                                    cohendo[asset_index, risk_index, risk_pre_index, laborincome_index, laborincome_pre_index, riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index, portfoliochoice_index] - K[asset_index]) + betapar * tempnext

        # std::cout << "EGM done\n"

        # rescalings

        for portfoliochoice_index in list(range(size_portfoliochoice+1)):
            for risk_index in list(range(size_risk)):
                for risk_pre_index in list(range(size_risk)):
                    for laborincome_index in list(range(size_laborincome)):
                        for laborincome_pre_index in list(range(size_laborincome)):
                            for riskshock_index in list(range(size_shock)):
                                for riskshock_pre_index in list(range(size_shock)):
                                    for laborshock_index in list(range(size_shock)):
                                        for laborshock_pre_index in list(range(size_shock)):

                                            shockstate_current = riskshock_index * size_shock + laborshock_index

                                            threshold_ii = 0

                                            for asset_index in list(range(size_asset)):

                                                # method 1: cash on hand
                                                cohexo = (1.0 + (r_f + pi + risk_states[risk_index])*riskshock_states[shockstate_current] * Omega[portfoliochoice_index] + r_f * (
                                                    1 - Omega[portfoliochoice_index])) * K[asset_index] + wagerate * laborincome_states[laborincome_index]*laborshock_states[shockstate_current]

                                                if cohexo < cohendo[(0, risk_index, risk_pre_index, laborincome_index,
                                                                     laborincome_pre_index, riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index, portfoliochoice_index)]:
                                                    save[(asset_index, risk_index, risk_pre_index, laborincome_index,
                                                          laborincome_pre_index, riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index, portfoliochoice_index)] = K[0]
                                                    VF[(asset_index, risk_index, risk_pre_index, laborincome_index,
                                                        laborincome_pre_index, riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index, portfoliochoice_index)] = U(cohexo - save[(asset_index, risk_index, risk_pre_index, laborincome_index,
                                                                                                                                                                                                        laborincome_pre_index, riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index, portfoliochoice_index)]) + (
                                                        VFendo[(0, risk_index, risk_pre_index, laborincome_index,
                                                                laborincome_pre_index, riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index, portfoliochoice_index)] - U((cohendo[(0, risk_index, risk_pre_index, laborincome_index,
                                                                                                                                                                                                           laborincome_pre_index, riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index, portfoliochoice_index)] - K[0])))

                                                if cohexo >= cohendo[(0, risk_index, risk_pre_index, laborincome_index,
                                                                     laborincome_pre_index, riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index, portfoliochoice_index)]:

                                                    itest = threshold_ii

                                                    while (itest < size_asset) and cohexo > cohendo[((itest, risk_index, risk_pre_index, laborincome_index,
                                                                                                      laborincome_pre_index, riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index, portfoliochoice_index))]:
                                                        itest += 1

                                                    if (itest == size_asset):
                                                        # extrapolation
                                                        vfweight = (cohexo - cohendo[(size_asset - 2, risk_index, risk_pre_index, laborincome_index,
                                                                                      laborincome_pre_index, riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index, portfoliochoice_index)]) / (cohendo[(
                                                                                          size_asset - 1, risk_index, risk_pre_index, laborincome_index,
                                                                                          laborincome_pre_index, riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index, portfoliochoice_index)] - cohendo[(size_asset - 2, risk_index, risk_pre_index, laborincome_index,
                                                                                                                                                                                                                                  laborincome_pre_index, riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index, portfoliochoice_index)])
                                                        igridL = size_asset - 2
                                                        igridH = size_asset - 1

                                                    else:
                                                        # standard interior
                                                        vfweight = (cohexo - cohendo[(itest - 1, risk_index, risk_pre_index, laborincome_index,
                                                                                      laborincome_pre_index, riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index, portfoliochoice_index)]) / (cohendo[(
                                                                                          itest, risk_index, risk_pre_index, laborincome_index, laborincome_pre_index, riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index, portfoliochoice_index)] - cohendo[(itest - 1, risk_index, risk_pre_index, laborincome_index, laborincome_pre_index, riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index, portfoliochoice_index)])
                                                        igridL = itest - 1
                                                        igridH = itest - 0

                                                    VF[(asset_index, risk_index, risk_pre_index, laborincome_index,
                                                        laborincome_pre_index, riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index, portfoliochoice_index)] = inter1d(vfweight, VFendo[(
                                                            igridL, risk_index, risk_pre_index, laborincome_index,
                                                            laborincome_pre_index, riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index, portfoliochoice_index)], VFendo[(igridH, risk_index, risk_pre_index, laborincome_index,
                                                                                                                                                                                                  laborincome_pre_index, riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index, portfoliochoice_index)])
                                                    save[(asset_index, risk_index, risk_pre_index, laborincome_index,
                                                          laborincome_pre_index, riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index, portfoliochoice_index)] = inter1d(vfweight, K[igridL], K[igridH])

                                                    threshold_ii = min(
                                                        size_asset - 2, itest)

        # std::cout << "rescaling done\n"

        for asset_index in list(range(size_asset)):
            for risk_index in list(range(size_risk)):
                for risk_pre_index in list(range(size_risk)):
                    for laborincome_index in list(range(size_laborincome)):
                        for laborincome_pre_index in list(range(size_laborincome)):
                            for riskshock_index in list(range(size_shock)):
                                for riskshock_pre_index in list(range(size_shock)):
                                    for laborshock_index in list(range(size_shock)):
                                        for laborshock_pre_index in list(range(size_shock)):

                                            # shockstate_pre = riskshock_index * size_shock + laborshock_index

                                            for portfoliochoice_index in list(range(size_portfoliochoice+1)):
                                                tempnext = 0.0
                                                for risk_pre_indexnext in list(range(size_risk)):
                                                    for laborincome_pre_indexnext in list(range(size_laborincome)):
                                                        for riskshock_pre_indexnext in list(range(size_shock)):
                                                            for laborshock_pre_indexnext in list(range(size_shock)):

                                                                shockstate_current = riskshock_pre_indexnext * \
                                                                    size_shock + laborshock_pre_indexnext

                                                                tempnext += risk_trans[risk_pre_index][risk_pre_indexnext]*laborincome_trans[laborincome_pre_index][laborincome_pre_indexnext] * risklaborshock_trans[0, shockstate_current] * VF[asset_index,
                                                                                                                                                                                                                                                  risk_pre_indexnext, risk_pre_index, laborincome_pre_indexnext, laborincome_pre_index, riskshock_pre_indexnext, riskshock_pre_index, laborshock_pre_indexnext, laborshock_pre_index, portfoliochoice_index]

                                                # std::cout << tempnext << "\n"
                                                if (portfoliochoice_index == 0):
                                                    temp = tempnext
                                                    itemp = 0

                                                if (tempnext > temp):

                                                    temp = tempnext
                                                    itemp = portfoliochoice_index

                                            VF_final[asset_index, risk_index, risk_pre_index, laborincome_index, laborincome_pre_index, riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index] = VF[
                                                asset_index, risk_index, risk_pre_index, laborincome_index, laborincome_pre_index, riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index, itemp]
                                            save_final[asset_index, risk_index, risk_pre_index, laborincome_index, laborincome_pre_index, riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index] = save[
                                                asset_index, risk_index, risk_pre_index, laborincome_index, laborincome_pre_index, riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index, itemp]
                                            Portfolio[asset_index, risk_index, risk_pre_index, laborincome_index, laborincome_pre_index,
                                                      riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index] = Omega[itemp]

        # std::cout << "port done\n"
        # std::cout << std::setprecision(16) << VF[(5, 5)] << "\n"

        # computing new derivatives and convergence
        critV = 0.0

        for asset_index in list(range(size_asset)):
            for risk_index in list(range(size_risk)):
                for risk_pre_index in list(range(size_risk)):
                    for laborincome_index in list(range(size_laborincome)):
                        for laborincome_pre_index in list(range(size_laborincome)):
                            for riskshock_index in list(range(size_shock)):
                                for riskshock_pre_index in list(range(size_shock)):
                                    for laborshock_index in list(range(size_shock)):
                                        for laborshock_pre_index in list(range(size_shock)):

                                            if (asset_index >= 2):
                                                dVF_final[(asset_index - 1, risk_index, risk_pre_index, laborincome_index, laborincome_pre_index,
                                                           riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index)] = nderiv(VF_final[(asset_index - 2, risk_index, risk_pre_index, laborincome_index, laborincome_pre_index,
                                                                                                                                                             riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index)], VF_final[(
                                                                                                                                                                 asset_index - 1, risk_index, risk_pre_index, laborincome_index, laborincome_pre_index,
                                                                                                                                                                 riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index)], VF_final[(asset_index, risk_index, risk_pre_index, laborincome_index, laborincome_pre_index,
                                                                                                                                                                                                                                                           riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index)], K[asset_index - 2], K[asset_index - 1], K[asset_index])

                                            critV = max(critV, abs(VF_final[(asset_index, risk_index, risk_pre_index, laborincome_index, laborincome_pre_index,
                                                                             riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index)
                                                                            ] - VF_final_old[(asset_index, risk_index, risk_pre_index, laborincome_index, laborincome_pre_index,
                                                                                              riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index)]))

                                            # left corner
                                            dVF_final[(0, risk_index, risk_pre_index, laborincome_index, laborincome_pre_index,
                                                       riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index)] = (VF_final[(
                                                           1, risk_index, risk_pre_index, laborincome_index, laborincome_pre_index,
                                                           riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index)] - VF_final[(0, risk_index, risk_pre_index, laborincome_index, laborincome_pre_index,
                                                                                                                                                      riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index)]) / (K[1] - K[0])
                                            # right corner
                                            dVF_final[(size_asset - 1, risk_index, risk_pre_index, laborincome_index, laborincome_pre_index,
                                                       riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index)] = (VF_final[(size_asset - 1, risk_index, risk_pre_index, laborincome_index,
                                                                                                                                                   laborincome_pre_index,
                                                                                                                                                   riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index)] - VF_final[(size_asset - 2, risk_index, risk_pre_index, laborincome_index, laborincome_pre_index,
                                                                                                                                                                                                                                              riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index)]) / (K[size_asset - 1] - K[size_asset - 2])

        for asset_index in list(range(size_asset)):
            for risk_index in list(range(size_risk)):
                for risk_pre_index in list(range(size_risk)):
                    for laborincome_index in list(range(size_laborincome)):
                        for laborincome_pre_index in list(range(size_laborincome)):
                            for riskshock_index in list(range(size_shock)):
                                for riskshock_pre_index in list(range(size_shock)):
                                    for laborshock_index in list(range(size_shock)):
                                        for laborshock_pre_index in list(range(size_shock)):
                                            VF[(asset_index, risk_index, risk_pre_index, laborincome_index, laborincome_pre_index,
                                                riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index, portfoliochoice_index)] = relaxVF * VF_final[(asset_index, risk_index, risk_pre_index,
                                                                                                                                                                            laborincome_index, laborincome_pre_index,
                                                                                                                                                                            riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index)] + (1 - relaxVF) * VF[(asset_index, risk_index, risk_pre_index, laborincome_index, laborincome_pre_index,
                                                                                                                                                                                                                                                                                 riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index, portfoliochoice_index)]

        for asset_index in list(range(size_asset)):
            for risk_index in list(range(size_risk)):
                for risk_pre_index in list(range(size_risk)):
                    for laborincome_index in list(range(size_laborincome)):
                        for laborincome_pre_index in list(range(size_laborincome)):
                            for riskshock_index in list(range(size_shock)):
                                for riskshock_pre_index in list(range(size_shock)):
                                    for laborshock_index in list(range(size_shock)):
                                        for laborshock_pre_index in list(range(size_shock)):
                                            for portfoliochoice_index in list(range(size_portfoliochoice+1)):

                                                if (asset_index >= 2):
                                                    dVF[(asset_index - 1, risk_index, risk_pre_index, laborincome_index, laborincome_pre_index,
                                                         riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index, portfoliochoice_index)] = nderiv(VF[(asset_index - 2, risk_index, risk_pre_index, laborincome_index, laborincome_pre_index, riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index, portfoliochoice_index)], VF[(asset_index - 1, risk_index, risk_pre_index, laborincome_index, laborincome_pre_index, riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index, portfoliochoice_index)], VF[(asset_index, risk_index, risk_pre_index, laborincome_index, laborincome_pre_index, riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index, portfoliochoice_index)], K[asset_index - 2], K[asset_index - 1], K[asset_index])

                                                # left corner
                                                dVF[(0, risk_index, risk_pre_index, laborincome_index, laborincome_pre_index,
                                                     riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index, portfoliochoice_index)] = (VF[(1, risk_index, risk_pre_index, laborincome_index,
                                                                                                                                                                  laborincome_pre_index,
                                                                                                                                                                  riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index, portfoliochoice_index)] - VF[(0, risk_index, risk_pre_index, laborincome_index, laborincome_pre_index,
                                                                                                                                                                                                                                                                              riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index, portfoliochoice_index)]) / (K[1] - K[0])
                                                # right corner
                                                dVF[(size_asset - 1, risk_index, risk_pre_index, laborincome_index, laborincome_pre_index,
                                                     riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index, portfoliochoice_index)] = (VF[(size_asset - 1, risk_index, risk_pre_index, laborincome_index,
                                                                                                                                                                  laborincome_pre_index,
                                                                                                                                                                  riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index, portfoliochoice_index)] - VF[(size_asset - 2, risk_index, risk_pre_index, laborincome_index, laborincome_pre_index,
                                                                                                                                                                                                                                                                              riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index, portfoliochoice_index)]) / (K[size_asset - 1] - K[size_asset - 2])

        iter += 1
        print("iteration={:d}, critV={:.8f}".format(iter, critV))
    return VF_final, dVF_final, save_final, VF, dVF, save, Portfolio


p_e_shock1 = 0.0
p_e_shock2 = 0.0
std_e_shock1 = 0.2
std_e_shock2 = 0.2
corr = .8
m_e_shock = 1


p_e_risk = 0.0
p_e_labor = 0.6
std_e_risk = 0.2
std_e_labor = 0.16
m_e_risk = 0
m_e_labor = 1


size_asset = 150  # number of grid points
size_risk = 2*m_e_risk+1   # number of productivity classes
size_shock = 2*m_e_shock+1
size_laborincome = 2*m_e_labor+1
size_portfoliochoice = 100

kmin = 0.0
kmax = 200.0

betapar = 0.9
alphapar = 0.36
deltapar = 0.08
rhopar = 3.0
labor = 1.0219882

epsV = 1.0e-8
epsdist = 1.0e-8
epsK = 1.0e-6
relaxsK = 0.005
relaxVF = 0.000

# grid ants
scale1 = 1.6
grmin = (kmin / scale1) - 1.0
exponen = np.log((kmax / scale1) - grmin) / (size_asset - 1)


pi = 0.00005
corr = 0.999
r_f = 0.03

# Function Definitions:

# Utility


K = np.zeros(size_asset)
Omega = np.zeros(size_portfoliochoice + 1)


riskshock_states, laborshock_states, risklaborshock_trans = tauchenfun2D(
    p_e_shock1, p_e_shock2, m_e_shock, 0.0, 0.0, std_e_shock1, std_e_shock2, corr)

risk_states, risk_trans = tauchenfun1D(
    p_e_risk, m_e_risk, 0.0, std_e_risk)
laborincome_states, laborincome_trans = tauchenfun1D(
    p_e_labor, m_e_labor, 0.0, std_e_labor)

laborincome_states = np.exp(laborincome_states)
laborshock_states = np.exp(laborshock_states)
riskshock_states = np.exp(riskshock_states)
# MARGINAL UTILITY, VALUES FUNCTION AND POLICIES #

# print(risk_trans.shape)
# print(risklaborshock_trans.shape)


# Note for users :: please, always use pointers and save your computer's memory ) == banish all arrays #
VF = np.zeros((size_asset, size_risk, size_risk, size_laborincome,
              size_laborincome, size_shock, size_shock, size_shock, size_shock, size_portfoliochoice+1))   # value function
dVF = np.zeros((size_asset, size_risk, size_risk, size_laborincome,
                size_laborincome, size_shock, size_shock, size_shock, size_shock, size_portfoliochoice+1))  # value function derivative
save = np.zeros((size_asset, size_risk, size_risk, size_laborincome,
                size_laborincome, size_shock, size_shock, size_shock, size_shock, size_portfoliochoice+1))  # value function
# print(VF.shape)

VF_final = np.zeros((size_asset, size_risk, size_risk, size_laborincome,
                    size_laborincome, size_shock, size_shock, size_shock, size_shock))  # value function
dVF_final = np.zeros((size_asset, size_risk, size_risk, size_laborincome,
                      size_laborincome, size_shock, size_shock, size_shock, size_shock))  # value function derivative

save_final = np.zeros((size_asset, size_risk, size_risk, size_laborincome,
                      size_laborincome, size_shock, size_shock, size_shock, size_shock))
Portfolio = np.zeros((size_asset, size_risk, size_risk, size_laborincome,
                      size_laborincome, size_shock, size_shock, size_shock, size_shock))
# cons = ( *)calloc((ARRLLP_dim), sizeof())
distin_final = np.zeros((size_asset, size_risk, size_risk, size_laborincome,
                         size_laborincome, size_shock, size_shock, size_shock, size_shock))
distout_final = np.zeros((size_asset, size_risk, size_risk, size_laborincome,
                          size_laborincome, size_shock, size_shock, size_shock, size_shock))

for asset_index in list(range(size_asset)):
    K[asset_index] = getlevel(asset_index)

for portfoliochoice_index in list(range(size_portfoliochoice+1)):
    Omega[portfoliochoice_index] = getomega(portfoliochoice_index)

# rrate = 0.040237086402090
wagerate = 0.8
distin_final[0] = 1.0
# taxL=0.3

# initializing value function and initial derivatives

for asset_index in list(range(size_asset)):
    for risk_index in list(range(size_risk)):
        for risk_pre_index in list(range(size_risk)):
            for laborincome_index in list(range(size_laborincome)):
                for laborincome_pre_index in list(range(size_laborincome)):
                    for riskshock_index in list(range(size_shock)):
                        for riskshock_pre_index in list(range(size_shock)):
                            for laborshock_index in list(range(size_shock)):
                                for laborshock_pre_index in list(range(size_shock)):
                                    for portfoliochoice_index in list(range(size_portfoliochoice+1)):
                                        shockstate_current = riskshock_index*size_shock+laborshock_index
                                        VF[(asset_index, risk_index, risk_pre_index, laborincome_index, laborincome_pre_index, riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index, portfoliochoice_index)] = U(
                                            wagerate * laborincome_states[laborincome_index]*laborshock_states[shockstate_current] + (1 + r_f + (pi + risk_states[risk_index])*riskshock_states[shockstate_current]) * K[asset_index])  # REQUIERE TO BE INCREASING IN K (the case here)


for asset_index in list(range(size_asset)):
    for risk_index in list(range(size_risk)):
        for risk_pre_index in list(range(size_risk)):
            for laborincome_index in list(range(size_laborincome)):
                for laborincome_pre_index in list(range(size_laborincome)):
                    for riskshock_index in list(range(size_shock)):
                        for riskshock_pre_index in list(range(size_shock)):
                            for laborshock_index in list(range(size_shock)):
                                for laborshock_pre_index in list(range(size_shock)):
                                    for portfoliochoice_index in list(range(size_portfoliochoice+1)):

                                        if (asset_index >= 2):
                                            dVF[(asset_index - 1, risk_index, risk_pre_index, laborincome_index, laborincome_pre_index,  riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index, portfoliochoice_index)] = nderiv(VF[(asset_index - 2, risk_index, risk_pre_index, laborincome_index, laborincome_pre_index,  riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index, portfoliochoice_index)], VF[(
                                                asset_index - 1, risk_index, risk_pre_index, laborincome_index, laborincome_pre_index,  riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index, portfoliochoice_index)], VF[(asset_index, risk_index, risk_pre_index, laborincome_index, laborincome_pre_index, riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index, portfoliochoice_index)], K[asset_index - 2], K[asset_index - 1], K[asset_index])

                                        # left corner
                                        dVF[(0, risk_index, risk_pre_index, laborincome_index, laborincome_pre_index,  riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index, portfoliochoice_index)] = (VF[(1, risk_index, risk_pre_index, laborincome_index,
                                                                                                                                                                                                                                    laborincome_pre_index,  riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index, portfoliochoice_index)] - VF[(0, risk_index, risk_pre_index, laborincome_index, laborincome_pre_index,  riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index, portfoliochoice_index)]) / (K[1] - K[0])
                                        # right corner
                                        dVF[(size_asset - 1, risk_index, risk_pre_index, laborincome_index, laborincome_pre_index,  riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index, portfoliochoice_index)] = (VF[(size_asset - 1, risk_index, risk_pre_index, laborincome_index, laborincome_pre_index, riskshock_index,
                                                                                                                                                                                                                                                 riskshock_pre_index, laborshock_index, laborshock_pre_index, portfoliochoice_index)] - VF[(size_asset - 2, risk_index, risk_pre_index, laborincome_index, laborincome_pre_index,  riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index, portfoliochoice_index)]) / (K[size_asset - 1] - K[size_asset - 2])

# printf("Policy Computation Start")

for asset_index in list(range(size_asset)):
    for risk_index in list(range(size_risk)):
        for risk_pre_index in list(range(size_risk)):
            for laborincome_index in list(range(size_laborincome)):
                for laborincome_pre_index in list(range(size_laborincome)):
                    for riskshock_index in list(range(size_shock)):
                        for riskshock_pre_index in list(range(size_shock)):
                            for laborshock_index in list(range(size_shock)):
                                for laborshock_pre_index in list(range(size_shock)):
                                    shockstate_current = riskshock_index*size_shock+laborshock_index

                                    VF_final[(asset_index, risk_index, risk_pre_index, laborincome_index, laborincome_pre_index, riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index)] = U(
                                        wagerate * laborincome_states[laborincome_index]*laborshock_states[shockstate_current] + (1 + r_f + (pi + risk_states[risk_index])*riskshock_states[shockstate_current]) * K[asset_index])  # REQUIERE TO BE INCREASING IN K (the case here)

# printf("Policy Computation Start")

for asset_index in list(range(size_asset)):
    for risk_index in list(range(size_risk)):
        for risk_pre_index in list(range(size_risk)):
            for laborincome_index in list(range(size_laborincome)):
                for laborincome_pre_index in list(range(size_laborincome)):
                    for riskshock_index in list(range(size_shock)):
                        for riskshock_pre_index in list(range(size_shock)):
                            for laborshock_index in list(range(size_shock)):
                                for laborshock_pre_index in list(range(size_shock)):

                                    if (asset_index >= 2):
                                        dVF_final[(asset_index - 1, risk_index, risk_pre_index, laborincome_index, laborincome_pre_index, riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index)] = nderiv(VF_final[(asset_index - 2, risk_index, risk_pre_index, laborincome_index, laborincome_pre_index, riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index)], VF_final[asset_index - 1,
                                                                                                                                                                                                                                                                                                                                                                                                                            risk_index, risk_pre_index, laborincome_index, laborincome_pre_index, riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index], VF_final[(asset_index, risk_index, risk_pre_index, laborincome_index, laborincome_pre_index, riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index)], K[asset_index - 2], K[asset_index - 1], K[asset_index])

                                    # left corner
                                    dVF_final[(0, risk_index, risk_pre_index, laborincome_index, laborincome_pre_index, riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index)] = (VF_final[(
                                        1, risk_index, risk_pre_index, laborincome_index, laborincome_pre_index, riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index)] - VF_final[(0, risk_index, risk_pre_index, laborincome_index, laborincome_pre_index, riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index)]) / (K[1] - K[0])
                                    # right corner
                                    dVF_final[(size_asset - 1, risk_index, risk_pre_index, laborincome_index, laborincome_pre_index, riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index)] = (VF_final[(size_asset - 1, risk_index, risk_pre_index, laborincome_index, laborincome_pre_index, riskshock_index,
                                                                                                                                                                                                                                 riskshock_pre_index, laborshock_index, laborshock_pre_index)] - VF_final[(size_asset - 2, risk_index, risk_pre_index, laborincome_index, laborincome_pre_index, riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index)]) / (K[size_asset - 1] - K[size_asset - 2])

print("Policy Computation Start")
VF_final, dVF_final, save_final, VF, dVF, save, Portfolio = POLICY(VF_final, dVF_final, save_final, VF,
                                                                   dVF, save, Portfolio, K, Omega, wagerate)
print("Policy Computation Done")


file_name = "test"
for risk_index in list(range(size_risk)):
    for risk_pre_index in list(range(size_risk)):
        for laborincome_index in list(range(size_laborincome)):
            for laborincome_pre_index in list(range(size_laborincome)):
                for riskshock_index in list(range(size_shock)):
                    for riskshock_pre_index in list(range(size_shock)):
                        for laborshock_index in list(range(size_shock)):
                            for laborshock_pre_index in list(range(size_shock)):

                                plt.plot(K, Portfolio[:, risk_index, risk_pre_index, laborincome_index, laborincome_pre_index,
                                         riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index])
                                plt.ylim((0, 1))
                                plt.xlim((0, kmax))
# plt.show()
plt.savefig("./figure/"+file_name+".pdf")
plt.close()

# SIMULATION(save_final, distin_final, & capital1, K)

# for (asset_index = 0 asset_index < size_asset asset_index++)
# {
#     std::cout << asset_index << "," << getlevel(asset_index) << ","

#     #  fprintf(dfilecsv, "%5d\t,%20.15f\t,", asset_index, phi(asset_index))
#     for (risk_index = 0 risk_index < size_risk risk_index++)
#     { # fprintf(dfilecsv, "%20.15f,", VF[(asset_index, risk_index)])
#         std::cout << distout[(asset_index, risk_index)] << ","
#     }
#     std::cout << "\n"
# }

# CreateFolder(".\\csv\\")
# CreateFolder(".\\figure\\")
# std: : string common = "14,pe=e-9,std=0.01,premium=" + std: : to_string(pi) + ",wage=" + std: : to_string(wagerate) + ",rf=" + std: : to_string(r_f) + ",Psize=" + std: : to_string(size_portfoliochoice) + ",rho_c=" + std: : to_string(rhopar) + ",Ksize=" + std: : to_string(size_asset) + ",Kmax=" + std: : to_string(kmax) + ",relaxVF=" + std: : to_string(relaxVF) + ",beta=" + std: : to_string(betapar) + ",corr=" + std: : to_string(corr) + ",Ssize=" + std: : to_string(size_risk) + ".csv "
# std: : string filename_dist = "csv\\dist" + common
# std: : string filename_policy = "csv\\policy" + common
# std: : string filename_VF = "csv\\VF" + common
# std: : string filename_Port = "csv\\Portfolio" + common

# # std::string var = "sometext" + std::to_string(pi)
# # std::cout << var

# std: : ofstream dfilecsv
# dfilecsv.open(filename_dist)
# dfilecsv << "gridnumber,"
#             << "capital,"
# tempcount = 0

# for (risk_index=0 risk_index < size_risk risk_index++)
# {
#     for (risk_pre_index=0 risk_pre_index < size_risk risk_pre_index++)
#     {
#         for (laborincome_index=0 laborincome_index < size_laborincome laborincome_index++)
#         {
#             for (laborincome_pre_index=0 laborincome_pre_index < size_laborincome laborincome_pre_index++)
#             {

#                 dfilecsv << "dist[" << risk_index << " " << risk_pre_index << " " << laborincome_index << laborincome_pre_index << "]"
#                 if (tempcount < size_laborincome * size_laborincome * size_risk * size_risk - 1)
#                 {
#                     dfilecsv << ","
#                 }
#                 else
#                 {
#                     dfilecsv << "\n"
#                 }
#                 tempcount++
#             }
#         }
#     }
# }

# for (asset_index=0 asset_index < size_asset asset_index++)
# {
#     dfilecsv << asset_index << "," << getlevel(asset_index) << ","

#     #  fprintf(dfilecsv, "%5d\t,%20.15f\t,", asset_index, phi(asset_index))
#     for (risk_index=0 risk_index < size_risk risk_index++)
#     {  # fprintf(dfilecsv, "%20.15f,", VF[(asset_index, risk_index)])
#         for (risk_pre_index=0 risk_pre_index < size_risk risk_pre_index++)
#         {
#             for (laborincome_index=0 laborincome_index < size_laborincome laborincome_index++)
#             {
#                 for (laborincome_pre_index=0 laborincome_pre_index < size_laborincome laborincome_pre_index++)
#                 {

#                     if (risk_index * size_risk * size_laborincome * size_laborincome + risk_pre_index * size_laborincome * size_laborincome + laborincome_index * size_laborincome + laborincome_pre_index < size_risk * size_risk * size_laborincome * size_laborincome - 1)
#                     {
#                         dfilecsv << distin_final[(
#                             asset_index, risk_index, risk_pre_index, laborincome_index, laborincome_pre_index)] << ","
#                     }
#                     else
#                     {
#                         dfilecsv << distin_final[(
#                             asset_index, risk_index, risk_pre_index, laborincome_index, laborincome_pre_index)]
#                     }
#                 }
#             }
#         }
#     }

#     dfilecsv << "\n"
# }

# dfilecsv.close()

# std: : ofstream policyfilecsv
# policyfilecsv.open(filename_policy)
# policyfilecsv << "gridnumber,"
# << "capital,"
# tempcount = 0

# for (risk_index=0 risk_index < size_risk risk_index++)
# {
#     for (risk_pre_index=0 risk_pre_index < size_risk risk_pre_index++)
#     {
#         for (laborincome_index=0 laborincome_index < size_laborincome laborincome_index++)
#         {
#             for (laborincome_pre_index=0 laborincome_pre_index < size_laborincome laborincome_pre_index++)
#             {
#                 policyfilecsv << "policy[" << risk_index << " " << risk_pre_index << " " << laborincome_index << laborincome_pre_index << "]"
#                 if (tempcount < size_laborincome * size_laborincome * size_risk * size_risk - 1)
#                 {
#                     policyfilecsv << ","
#                 }
#                 else
#                 {
#                     policyfilecsv << "\n"
#                 }
#                 tempcount++
#             }
#         }
#     }
# }

# for (asset_index=0 asset_index < size_asset asset_index++)
# {
#     policyfilecsv << asset_index << "," << getlevel(asset_index) << ","

#     #  fprintf(dfilecsv, "%5d\t,%20.15f\t,", asset_index, phi(asset_index))
#     for (risk_index=0 risk_index < size_risk risk_index++)
#     {  # fprintf(dfilecsv, "%20.15f,", VF[(asset_index, risk_index)])
#         for (risk_pre_index=0 risk_pre_index < size_risk risk_pre_index++)
#         {
#             for (laborincome_index=0 laborincome_index < size_laborincome laborincome_index++)
#             {
#                 for (laborincome_pre_index=0 laborincome_pre_index < size_laborincome laborincome_pre_index++)
#                 {
#                     if (risk_index * size_risk * size_laborincome * size_laborincome + risk_pre_index * size_laborincome * size_laborincome + laborincome_index * size_laborincome + laborincome_pre_index < size_risk * size_risk * size_laborincome * size_laborincome - 1)
#                     {
#                         policyfilecsv << save_final[(
#                             asset_index, risk_index, risk_pre_index, laborincome_index, laborincome_pre_index)] << ","
#                     }
#                     else
#                     {
#                         policyfilecsv << save_final[(
#                             asset_index, risk_index, risk_pre_index, laborincome_index, laborincome_pre_index)]
#                     }
#                 }
#             }
#         }
#     }

#     policyfilecsv << "\n"
# }

# policyfilecsv.close()

# std: : ofstream VFfilecsv
# VFfilecsv.open(filename_VF)
# VFfilecsv << "gridnumber,"
#             << "capital,"
# tempcount = 0

# for (risk_index=0 risk_index < size_risk risk_index++)
# {
#     for (risk_pre_index=0 risk_pre_index < size_risk risk_pre_index++)
#     {
#         for (laborincome_index=0 laborincome_index < size_laborincome laborincome_index++)
#         {
#             for (laborincome_pre_index=0 laborincome_pre_index < size_laborincome laborincome_pre_index++)
#             {
#                 VFfilecsv << "VF[" << risk_index << " " << risk_pre_index << " " << laborincome_index << laborincome_pre_index << "]"
#                 if (tempcount < size_laborincome * size_laborincome * size_risk * size_risk - 1)
#                 {
#                     VFfilecsv << ","
#                 }
#                 else
#                 {
#                     VFfilecsv << "\n"
#                 }
#                 tempcount++
#             }
#         }
#     }
# }

# for (asset_index=0 asset_index < size_asset asset_index++)
# {
#     VFfilecsv << asset_index << "," << getlevel(asset_index) << ","

#     #  fprintf(dfilecsv, "%5d\t,%20.15f\t,", asset_index, phi(asset_index))
#     for (risk_index=0 risk_index < size_risk risk_index++)
#     {  # fprintf(dfilecsv, "%20.15f,", VF[(asset_index, risk_index)])
#         for (risk_pre_index=0 risk_pre_index < size_risk risk_pre_index++)
#         {
#             for (laborincome_index=0 laborincome_index < size_laborincome laborincome_index++)
#             {
#                 for (laborincome_pre_index=0 laborincome_pre_index < size_laborincome laborincome_pre_index++)
#                 {
#                     if (risk_index * size_risk * size_laborincome * size_laborincome + risk_pre_index * size_laborincome * size_laborincome + laborincome_index * size_laborincome + laborincome_pre_index < size_risk * size_risk * size_laborincome * size_laborincome - 1)
#                     {
#                         VFfilecsv << VF_final[(
#                             asset_index, risk_index, risk_pre_index, laborincome_index, laborincome_pre_index)] << ","
#                     }
#                     else
#                     {
#                         VFfilecsv << VF_final[(
#                             asset_index, risk_index, risk_pre_index, laborincome_index, laborincome_pre_index)]
#                     }
#                 }
#             }
#         }
#     }

#     VFfilecsv << "\n"
# }

# VFfilecsv.close()

# std: : ofstream Portfilecsv
# Portfilecsv.open(filename_Port)
# Portfilecsv << "gridnumber,"
#             << "capital,"
# tempcount = 0

# for (risk_index=0 risk_index < size_risk risk_index++)
# {
#     for (risk_pre_index=0 risk_pre_index < size_risk risk_pre_index++)
#     {
#         for (laborincome_index=0 laborincome_index < size_laborincome laborincome_index++)
#         {
#             for (laborincome_pre_index=0 laborincome_pre_index < size_laborincome laborincome_pre_index++)
#             {
#                 Portfilecsv << "Portfolio[" << risk_index << " " << risk_pre_index << " " << laborincome_index << laborincome_pre_index << "]"

#                 if (tempcount < size_laborincome * size_laborincome * size_risk * size_risk - 1)
#                 {
#                     Portfilecsv << ","
#                 }
#                 else
#                 {
#                     Portfilecsv << "\n"
#                 }
#                 tempcount++
#             }
#         }
#     }
# }

# for (asset_index=0 asset_index < size_asset asset_index++)
# {
#     Portfilecsv << asset_index << "," << getlevel(asset_index) << ","

#     #  fprintf(dfilecsv, "%5d\t,%20.15f\t,", asset_index, phi(asset_index))
#     for (risk_index=0 risk_index < size_risk risk_index++)
#     {  # fprintf(dfilecsv, "%20.15f,", VF[(asset_index, risk_index)])
#         for (risk_pre_index=0 risk_pre_index < size_risk risk_pre_index++)
#         {
#             for (laborincome_index=0 laborincome_index < size_laborincome laborincome_index++)
#             {
#                 for (laborincome_pre_index=0 laborincome_pre_index < size_laborincome laborincome_pre_index++)
#                 {
#                     if (risk_index * size_risk * size_laborincome * size_laborincome + risk_pre_index * size_laborincome * size_laborincome + laborincome_index * size_laborincome + laborincome_pre_index < size_risk * size_risk * size_laborincome * size_laborincome - 1)
#                     {
#                         Portfilecsv << Portfolio[(
#                             asset_index, risk_index, risk_pre_index, laborincome_index, laborincome_pre_index)] << ","
#                     }
#                     else
#                     {
#                         Portfilecsv << Portfolio[(
#                             asset_index, risk_index, risk_pre_index, laborincome_index, laborincome_pre_index)]
#                     }
#                 }
#             }
#         }
#     }

#     Portfilecsv << "\n"
# }

# Portfilecsv.close()
