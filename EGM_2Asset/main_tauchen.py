import sys
from xml.dom.minicompat import NodeList
from scipy.stats import multivariate_normal as mvn
import numpy as np

# mean = np.array([0, 0])
# corr = 1
# covariance = np.array([[1, corr], [corr, 1]])
# dist = mvn(mean=mean, cov=covariance, allow_singular=True)
# print("CDF:", dist.cdf(np.array([0, 0])))

p_e_shock1 = 0.0
p_e_shock2 = 0.0
std_e_shock1 = 0.2
std_e_shock2 = 0.2
corr = 1
m_e_shock = 2


p_e_risk = 0.0
p_e_labor = 0.6
std_e_risk = 0.01
# std_e_labor = 0.16
std_e_labor = 0.25
m_e_risk = 1
m_e_labor = 2


def CDFSTDNormal_1D(x1):
    mean = np.array(0)
    # corr = 1
    covariance = np.array(1)
    dist = mvn(mean=mean, cov=covariance, allow_singular=True)
    value = dist.cdf(np.array([x1]))

    return value


def CDFSTDNormal_2D_r0y0(x1, x2, corr):
    mean = np.array([0, 0])
    # corr = 1
    covariance = np.array([[1, corr], [corr, 1]])
    dist = mvn(mean=mean, cov=covariance, allow_singular=True)
    value = dist.cdf(np.array([x1, x2]))

    return value


def CDFSTDNormal_2D_r0y1(x1, x2, corr):

    value = CDFSTDNormal_1D(x1)-CDFSTDNormal_2D_r0y0(x1, x2, corr)

    return value


def CDFSTDNormal_2D_r1y0(x1, x2, corr):

    value = CDFSTDNormal_1D(x2)-CDFSTDNormal_2D_r0y0(x1, x2, corr)

    return value


def CDFSTDNormal_2D_r1y1(x1, x2, corr):

    value = 1-CDFSTDNormal_2D_r0y1(x1, x2, corr) - CDFSTDNormal_2D_r0y0(
        x1, x2, corr) - CDFSTDNormal_2D_r1y0(x1, x2, corr)

    return value


# print(CDFSTDNormal_1D(2))
# print(CDFSTDNormal_2D_r0y0(0, 0, 1))
# print(CDFSTDNormal_2D_r1y0(2, 0, 1) == 0)
# print(CDFSTDNormal_2D_r0y1(0, 2, 1) == 0)
# print(CDFSTDNormal_2D_r1y1(0, 0, 1))


def tauchenfun2D(rho1,  rho2,  m,  amu1,  amu2,  sigma1,  sigma2,  corr):

    Nodes = 2*m+1

    AR = np.matrix([[rho1, 0], [0, rho2]])
    Sigma = np.matrix([[sigma1**2, corr*sigma1*sigma2],
                       [corr*sigma1*sigma2, sigma2**2]])

    ARAR = np.kron(AR, AR)
    I = np.eye(4)
    Sigma_Vec = Sigma.ravel(order='A')
    Var_Vec = np.linalg.inv(I-ARAR)*np.transpose(Sigma_Vec)


# 	First lets compute unconditiontal variance of yt

# 	Compute stddev of yt

    # std1t = sqrt(var1t)

    std1t = np.sqrt(Var_Vec[0])

# 	 std2t = sqrt(var2t)
    std2t = np.sqrt(Var_Vec[3])

# 	std: : cout << "value of std_e_shock1=" << std1t << "\n"
# 	std: : cout << "value of std_e_shock2=" << std2t << "\n"
# 	Define maximum and minimum grid point

    y1nodes = np.zeros(Nodes)
    y2nodes = np.zeros(Nodes)

    y1nodes[Nodes - 1] = m * std1t
    y1nodes[0] = -y1nodes[Nodes - 1]

# 	Define interior nodes

    y1nodesinterval = (y1nodes[Nodes - 1] - y1nodes[0]) / ((Nodes - 1) * 1.0)

    for i in list(range(1, Nodes)):
        y1nodes[i] = y1nodes[i - 1] + y1nodesinterval

    for i in list(range(Nodes)):
        y1nodes[i] = y1nodes[i] + amu1

    y2nodes[Nodes - 1] = m * std2t
    y2nodes[0] = -y2nodes[Nodes - 1]

# 	Define interior nodes

    y2nodesinterval = (y2nodes[Nodes - 1] - y2nodes[0]) / ((Nodes - 1) * 1.0)

    for i in list(range(1, Nodes)):
        y2nodes[i] = y2nodes[i - 1] + y2nodesinterval

    for i in list(range(Nodes)):
        y2nodes[i] = y2nodes[i] + amu2

    grid1 = np.zeros(Nodes**2)
    grid2 = np.zeros(Nodes**2)
    for i in list(range(Nodes)):
        for j in list(range(Nodes)):
            grid1[i * Nodes + j] = y1nodes[i]
            grid2[i * Nodes + j] = y2nodes[j]

    # Computing transition probability matrix

    transitionMat = np.zeros((Nodes**2, Nodes**2))

    for ir in list(range(Nodes)):
        for iy in list(range(Nodes)):
            tempr = (1 - rho1) * amu1 + rho1 * y1nodes[ir]
            tempy = (1 - rho2) * amu2 + rho2 * y2nodes[iy]
            index_ry = ir * Nodes + iy

# 								        jr = 0,jy =0				(1-rho1)*mu_r + rho1*r_i+epsilon <  r_j+delta_r/2 &&					(1-rho2)*mu_y + rho2*y_i+epsilon <  y_j+delta_y/2
# 								        jr = 0,0 <jy<n2-1			(1-rho1)*mu_r + rho1*r_i+epsilon <  r_j+delta_r/2 &&	y_j-delta_y/2 < (1-rho2)*mu_y + rho2*y_i+epsilon <  y_j+delta_y/2
# 										jr = 0,jy =n2-1			(1-rho1)*mu_r + rho1*r_i+epsilon <  r_j+delta_r/2 &&	y_j-delta_y/2 < (1-rho2)*mu_y + rho2*y_i+epsilon

            jr = 0

            for jy in list(range(1, Nodes)):
                temp = 0
                index_ry_next = jr * Nodes + jy
                # print(jy)
                # print((y1nodes[jr] + y1nodesinterval / 2 - tempr) /
                #       std1t[0, 0])
                # print((y2nodes[jy] + y2nodesinterval / 2 - tempy) / std2t)
                temp += CDFSTDNormal_2D_r0y0((y1nodes[jr] + y1nodesinterval / 2 - tempr) /
                                             std1t[0, 0], (y2nodes[jy] + y2nodesinterval / 2 - tempy) / std2t[0, 0], corr)
                temp -= CDFSTDNormal_2D_r0y0((y1nodes[jr] + y1nodesinterval / 2 - tempr) /
                                             std1t[0, 0], (y2nodes[jy] - y2nodesinterval / 2 - tempy) / std2t[0, 0], corr)
                transitionMat[index_ry][index_ry_next] = temp
                temp = 0

            jy = 0

            temp = 0
            index_ry_next = jr * Nodes + jy
            temp += CDFSTDNormal_2D_r0y0((y1nodes[jr] + y1nodesinterval / 2 - tempr) /
                                         std1t[0, 0], (y2nodes[jy] + y2nodesinterval / 2 - tempy) / std2t[0, 0], corr)
            transitionMat[index_ry][index_ry_next] = temp
            temp = 0

            jy = (Nodes - 1)

            temp = 0
            index_ry_next = jr * Nodes + jy
            temp += CDFSTDNormal_2D_r0y1((y1nodes[jr] + y1nodesinterval / 2 - tempr) /
                                         std1t[0, 0], (y2nodes[jy] - y2nodesinterval / 2 - tempy) / std2t[0, 0], corr)
            transitionMat[index_ry][index_ry_next] = temp
            temp = 0


# 			P_{(ir,iy),(jr,jy)} = {    0 <jr<n1-1,jy=0          r_j-delta_r/2 < (1-rho2)*mu_r + rho2*r_i+epsilon <  r_j+delta_r/2 &&	(1-rho2)*mu_y + rho2*y_i+epsilon <  y_j+delta_y/2						 }
# 			//							0 < jr <n1-1,0<jy<n2-1		r_j-delta_r/2 < (1-rho2)*mu_r + rho2*r_i+epsilon <  r_j+delta_r/2 &&	y_j-delta_y/2 < (1-rho2)*mu_y + rho2*y_i+epsilon <  y_j+delta_y/2
# 			//							0 < jr <n1-1,jy=n2-1		r_j-delta_r/2 < (1-rho2)*mu_r + rho2*r_i+epsilon <  r_j+delta_r/2 &&	y_j-delta_y/2 < (1-rho2)*mu_y + rho2*y_i+epsilon

            for jr in list(range(1, Nodes)):

                for jy in list(range(1, Nodes)):
                    temp = 0
                    index_ry_next = jr * Nodes + jy
                    temp += CDFSTDNormal_2D_r0y0((y1nodes[jr] + y1nodesinterval / 2 - tempr) /
                                                 std1t[0, 0], (y2nodes[jy] + y2nodesinterval / 2 - tempy) / std2t[0, 0], corr)
                    temp -= CDFSTDNormal_2D_r0y0((y1nodes[jr] - y1nodesinterval / 2 - tempr) /
                                                 std1t[0, 0], (y2nodes[jy] + y2nodesinterval / 2 - tempy) / std2t[0, 0], corr)
                    temp -= CDFSTDNormal_2D_r0y0((y1nodes[jr] + y1nodesinterval / 2 - tempr) /
                                                 std1t[0, 0], (y2nodes[jy] - y2nodesinterval / 2 - tempy) / std2t[0, 0], corr)
                    temp += CDFSTDNormal_2D_r0y0((y1nodes[jr] - y1nodesinterval / 2 - tempr) /
                                                 std1t[0, 0], (y2nodes[jy] - y2nodesinterval / 2 - tempy) / std2t[0, 0], corr)
                    transitionMat[index_ry][index_ry_next] = temp
                    temp = 0

                jy = 0

                temp = 0
                index_ry_next = jr * Nodes + jy
                temp += CDFSTDNormal_2D_r0y0((y1nodes[jr] + y1nodesinterval / 2 - tempr) /
                                             std1t[0, 0], (y2nodes[jy] + y2nodesinterval / 2 - tempy) / std2t[0, 0], corr)
                temp -= CDFSTDNormal_2D_r0y0((y1nodes[jr] - y1nodesinterval / 2 - tempr) /
                                             std1t[0, 0], (y2nodes[jy] + y2nodesinterval / 2 - tempy) / std2t[0, 0], corr)
                transitionMat[index_ry][index_ry_next] = temp
                temp = 0

                jy = (Nodes - 1)

                temp = 0
                index_ry_next = jr * Nodes + jy
                temp += CDFSTDNormal_2D_r0y1((y1nodes[jr] + y1nodesinterval / 2 - tempr) /
                                             std1t[0, 0], (y2nodes[jy] - y2nodesinterval / 2 - tempy) / std2t[0, 0], corr)
                temp -= CDFSTDNormal_2D_r0y1((y1nodes[jr] - y1nodesinterval / 2 - tempr) /
                                             std1t[0, 0], (y2nodes[jy] - y2nodesinterval / 2 - tempy) / std2t[0, 0], corr)
                transitionMat[index_ry][index_ry_next] = temp
                temp = 0


# 			//							jr = n1-1,jy =0			r_j-delta_r/2 <	(1-rho1)*mu_r + rho1*r_i+epsilon	&&	(1-rho2)*mu_y + rho2*y_i+epsilon <  y_j+delta_y/2
# 			//							jr = n1-1,0 <jy<n2-1		r_j-delta_r/2 <	(1-rho1)*mu_r + rho1*r_i+epsilon	&&	y_j-delta_y/2 < (1-rho2)*mu_y + rho2*y_i+epsilon <  y_j+delta_y/2
# 										jr = n1-1,jy =n2-1			r_j-delta_r/2 <	(1-rho1)*mu_r + rho1*r_i+epsilon	&&	y_j-delta_y/2 < (1-rho2)*mu_y + rho2*y_i+epsilon

            jr = Nodes - 1

            for jy in list(range(1, Nodes)):
                temp = 0
                index_ry_next = jr * Nodes + jy
                temp += CDFSTDNormal_2D_r1y0((y1nodes[jr] - y1nodesinterval / 2 - tempr) /
                                             std1t[0, 0], (y2nodes[jy] + y2nodesinterval / 2 - tempy) / std2t[0, 0], corr)
                temp -= CDFSTDNormal_2D_r1y0((y1nodes[jr] - y1nodesinterval / 2 - tempr) /
                                             std1t[0, 0], (y2nodes[jy] - y2nodesinterval / 2 - tempy) / std2t[0, 0], corr)
                transitionMat[index_ry][index_ry_next] = temp
                temp = 0

            jy = 0

            temp = 0
            index_ry_next = jr * Nodes + jy
            temp += CDFSTDNormal_2D_r1y0((y1nodes[jr] - y1nodesinterval / 2 - tempr) /
                                         std1t[0, 0], (y2nodes[jy] + y2nodesinterval / 2 - tempy) / std2t[0, 0], corr)
            transitionMat[index_ry][index_ry_next] = temp
            temp = 0

            jy = (Nodes - 1)

            temp = 0
            index_ry_next = jr * Nodes + jy
            temp += CDFSTDNormal_2D_r1y1((y1nodes[jr] - y1nodesinterval / 2 - tempr) /
                                         std1t[0, 0], (y2nodes[jy] - y2nodesinterval / 2 - tempy) / std2t[0, 0], corr)
            transitionMat[index_ry][index_ry_next] = temp
            temp = 0

    if rho1 == 0 and rho2 == 0:
        transitionMat_vec = np.reshape(
            transitionMat[0, :], (1, len(transitionMat)))
        return grid1, grid2, transitionMat_vec
    else:
        return grid1, grid2, transitionMat


def tauchenfun1D(rho, m,  amu, sigma):

    Nodes = 2*m+1
    varyt = (sigma**2)/(1-rho*rho)

    stdyt = np.sqrt(varyt)

    ynodes = np.zeros(Nodes)

    ynodes[Nodes-1] = m*stdyt
    ynodes[0] = -ynodes[Nodes-1]
    ynodesinterval = (ynodes[Nodes-1]-ynodes[0])/((Nodes-1)*1.0)

    for i in list(range(1, Nodes)):
        ynodes[i] = ynodes[i-1]+ynodesinterval

    for i in list(range(Nodes)):
        ynodes[i] = ynodes[i]+amu

    grid = np.zeros(Nodes)
    grid = ynodes

    transitionMat = np.zeros((Nodes, Nodes))

    for j in list(range(Nodes)):
        for k in list(range(1, Nodes)):

            transitionMat[j][k] = CDFSTDNormal_1D((ynodes[k]-(1-rho)*amu-rho*ynodes[j]+ynodesinterval/2.0)/sigma) - CDFSTDNormal_1D(
                (ynodes[k]-(1-rho)*amu-rho*ynodes[j]-ynodesinterval/2.0)/sigma)

        transitionMat[j][0] = CDFSTDNormal_1D(
            (ynodes[0]-(1-rho)*amu-rho*ynodes[j]+ynodesinterval/2.0)/sigma)
        transitionMat[j][Nodes-1] = 1.0-CDFSTDNormal_1D(
            (ynodes[Nodes-1]-(1-rho)*amu-rho*ynodes[j]-ynodesinterval/2.0)/sigma)

    if rho == 0:
        transitionMat_vec = np.reshape(transitionMat[0, :], (1, Nodes))
        return grid, transitionMat
    else:
        return grid, transitionMat


grid1, grid2, ytrans = tauchenfun2D(
    p_e_shock1, p_e_shock2, m_e_shock, 0.0, 0.0, std_e_shock1, std_e_shock2, corr)

grid_risk, ytrans_risk = tauchenfun1D(p_e_risk, m_e_risk, 0.0, std_e_risk)
grid_labor, ytrans_labor = tauchenfun1D(p_e_labor, m_e_labor, 0.0, std_e_labor)

print(ytrans_risk.shape)

print('tauchenfun2D:Shock.hpp')

original_stdout = sys.stdout  # Save a reference to the original standard output

with open('./txt/Shock.hpp', 'w') as f:
    sys.stdout = f  # Change the standard output to the file we created.
    # print("#include <cmath>\n")

    print("// For Copy: PRODUCTIVITY SHOCKS\n")
    print("// For Copy: TRANSITION MATRIX FOR PRODUCTIVITY SHOCKS\n")

    print("const  double  risk_shock_states[{:d}] ={{".format(
        len(grid1)), end="")

    for y in range(len(grid1)):

        # print("{:.15f}\t".format(grid1[y]))

        if (y < len(grid1) - 1):
            print("exp({:.15f}), ".format(grid1[y]), end="")
        else:
            print("exp({:.15f})".format(grid1[y]), end="")

    print("};\n")

    print("const  double  laborincome_shock_states[{:d}] ={{".format(
        len(grid2)), end="")
    for y in range(len(grid2)):

        # print("{:.15f}\t".format(grid2[y]))
        # grid2[y] = np.exp(grid2[y])

        if (y < len(grid2) - 1):
            print("exp({:.15f}), ".format(grid2[y]), end="")
        else:
            print("exp({:.15f})".format(grid2[y]), end="")

    print("};\n")

    print("const  double  risk_labor_shock_trans[{:d}][{:d}] = {{\n".format(
        ytrans.shape[0], ytrans.shape[1]), end="")

    for y in list(range(ytrans.shape[0])):
        print("{", end="")
        for k in list(range(ytrans.shape[1])):

            if (k < ytrans.shape[1] - 1):
                print("{:.15f}, ".format(ytrans[y, k]), end="")
            else:
                print("{:.15f}".format(ytrans[y, k]), end="")

        if (y < ytrans.shape[0] - 1):
            print("},\n", end="")
        else:
            print("}};", end="")
    print("\n")

    print("const  double  risk_labor_shock_pre_trans[{:d}][{:d}] = {{\n".format(
        ytrans.shape[0], ytrans.shape[1]), end="")

    for y in list(range(ytrans.shape[0])):
        print("{", end="")
        for k in list(range(ytrans.shape[1])):

            if (k < ytrans.shape[1] - 1):
                print("{:.15f}, ".format(ytrans[y, k]), end="")
            else:
                print("{:.15f}".format(ytrans[y, k]), end="")

        if (y < ytrans.shape[0] - 1):
            print("},\n", end="")
        else:
            print("}};", end="")
    print("\n")

    # print("};\n", end="")

    print("const  int  size_shock={:d}; \n".format(
        2*m_e_shock+1), end="")

    sys.stdout = original_stdout  # Reset the standard output to its original value


print('tauchenfun1D:Risk_Labor.hpp')

original_stdout = sys.stdout  # Save a reference to the original standard output

with open('./txt/Risk_Labor.hpp', 'w') as f:
    sys.stdout = f  # Change the standard output to the file we created.
    # print("#include <cmath>\n")

    print("// For Copy: PRODUCTIVITY SHOCKS\n")
    print("// For Copy: TRANSITION MATRIX FOR PRODUCTIVITY SHOCKS\n")

    print("const double  risk_states[{:d}] ={{".format(len(grid_risk)), end="")

    # for y in range(len(grid_risk)):

    #     # print("{:.15f}\t".format(grid1[y]))

    #     if (y < len(grid_risk) - 1):
    #         print("{:.15f}, ".format(grid_risk[y]), end="")
    #     else:
    #         print("{:.15f}".format(grid_risk[y]), end="")

    # print("};\n")

    for y in range(len(grid_risk)):

        # print("{:.15f}\t".format(grid1[y]))

        if (y < len(grid_risk) - 1):
            print("exp({:.15f}), ".format(grid_risk[y]), end="")
        else:
            print("exp({:.15f})".format(grid_risk[y]), end="")

    print("};\n")

    print("const  double laborincome_states[{:d}] ={{".format(
        len(grid_labor)), end="")
    for y in range(len(grid_labor)):

        # print("{:.15f}\t".format(grid2[y]))
        # grid2[y] = np.exp(grid2[y])

        if (y < len(grid_labor) - 1):
            print("exp({:.15f}), ".format(grid_labor[y]), end="")
        else:
            print("exp({:.15f})".format(grid_labor[y]), end="")

    print("};\n")

    print("const  double risk_trans[{:d}][{:d}] = {{\n".format(
        ytrans_risk.shape[0], ytrans_risk.shape[1]), end="")

    for y in list(range(ytrans_risk.shape[0])):
        print("{", end="")
        for k in list(range(ytrans_risk.shape[1])):

            if (k < ytrans_risk.shape[1] - 1):
                print("{:.15f}, ".format(ytrans_risk[y, k]), end="")
            else:
                print("{:.15f}".format(ytrans_risk[y, k]), end="")

        if (y < ytrans_risk.shape[0] - 1):
            print("},\n", end="")
        else:
            print("}};", end="")
    print("\n")

    print("const  double laborincome_trans[{:d}][{:d}] = {{\n".format(
        len(grid_labor), len(grid_labor)), end="")

    for y in list(range(len(grid_labor))):
        print("{", end="")
        for k in list(range(len(grid_labor))):

            if (k < len(grid_labor) - 1):
                print("{:.15f}, ".format(ytrans_labor[y, k]), end="")
            else:
                print("{:.15f}".format(ytrans_labor[y, k]), end="")

        if (y < len(grid_labor) - 1):
            print("},\n", end="")
        else:
            print("}};", end="")
    print("\n")

    # print("};\n", end="")

    print("const double  risk_pre_states[{:d}] ={{".format(
        len(grid_risk)), end="")

    for y in range(len(grid_risk)):

        # print("{:.15f}\t".format(grid1[y]))

        if (y < len(grid_risk) - 1):
            print("{:.15f}, ".format(grid_risk[y]), end="")
        else:
            print("{:.15f}".format(grid_risk[y]), end="")

    print("};\n")

    print("const double  laborincome_pre_states[{:d}] ={{".format(
        len(grid_labor)), end="")

    for y in range(len(grid_labor)):

        # print("{:.15f}\t".format(grid2[y]))
        # grid2[y] = np.exp(grid2[y])

        if (y < len(grid_labor) - 1):
            print("exp({:.15f}), ".format(grid_labor[y]), end="")
        else:
            print("exp({:.15f})".format(grid_labor[y]), end="")

    print("};\n")

    print("const  double risk_pre_trans[{:d}][{:d}] = {{\n".format(
        ytrans_risk.shape[0], ytrans_risk.shape[1]), end="")

    for y in list(range(ytrans_risk.shape[0])):
        print("{", end="")
        for k in list(range(ytrans_risk.shape[1])):

            if (k < ytrans_risk.shape[1] - 1):
                print("{:.15f}, ".format(ytrans_risk[y, k]), end="")
            else:
                print("{:.15f}".format(ytrans_risk[y, k]), end="")

        if (y < ytrans_risk.shape[0] - 1):
            print("},\n", end="")
        else:
            print("}};", end="")
    print("\n")

    print("const double  laborincome_pre_trans[{:d}][{:d}] = {{\n".format(
        len(grid_labor), len(grid_labor)), end="")

    for y in list(range(len(grid_labor))):
        print("{", end="")
        for k in list(range(len(grid_labor))):

            if (k < len(grid_labor) - 1):
                print("{:.15f}, ".format(ytrans_labor[y, k]), end="")
            else:
                print("{:.15f}".format(ytrans_labor[y, k]), end="")

        if (y < len(grid_labor) - 1):
            print("},\n", end="")
        else:
            print("}};", end="")
    print("\n")

    print("const int  size_laborincome = {:d}; \n".format(
        len(grid_labor)), end="")

    print("const int  size_risk = {:d}; \n".format(
        len(grid_risk)), end="")

    print("const double  std_labor = {:f}; \n".format(
        std_e_labor), end="")

    sys.stdout = original_stdout  # Reset the standard output to its original value
