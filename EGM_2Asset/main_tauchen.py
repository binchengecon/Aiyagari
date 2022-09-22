import sys
from scipy.stats import multivariate_normal as mvn
import numpy as np

# mean = np.array([0, 0])
# corr = 1
# covariance = np.array([[1, corr], [corr, 1]])
# dist = mvn(mean=mean, cov=covariance, allow_singular=True)
# print("CDF:", dist.cdf(np.array([0, 0])))

p_e1 = 0.000000001
p_e2 = 0.6
std_e1 = 0.2
std_e2 = 0.16
corr = .9
m_e = 3


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


print(CDFSTDNormal_1D(2))
print(CDFSTDNormal_2D_r0y0(0, 0, 1))
print(CDFSTDNormal_2D_r1y0(2, 0, 1) == 0)
print(CDFSTDNormal_2D_r0y1(0, 2, 1) == 0)
print(CDFSTDNormal_2D_r1y1(0, 0, 1))


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

# 	std: : cout << "value of std_e1=" << std1t << "\n"
# 	std: : cout << "value of std_e2=" << std2t << "\n"
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

    return grid1, grid2, transitionMat


grid1, grid2, ytrans = tauchenfun2D(
    p_e1, p_e2, m_e, 0.0, 0.0, std_e1, std_e2, corr)

# print(grid1[2])
# print(np.sum(ytrans >= 0))


print('This message will be displayed on the screen.')

original_stdout = sys.stdout  # Save a reference to the original standard output

with open('./txt/AR1.hpp', 'w') as f:
    sys.stdout = f  # Change the standard output to the file we created.
    # print("#include <cmath>\n")

    print("// For Copy: PRODUCTIVITY SHOCKS\n")
    print("// For Copy: TRANSITION MATRIX FOR PRODUCTIVITY SHOCKS\n")

    print("const double risk_states[{:d}] ={{".format(len(grid1)), end="")

    for y in range(len(grid1)):

        # print("{:.15f}\t".format(grid1[y]))

        if (y < len(grid1) - 1):
            print("{:.15f}, ".format(grid1[y]), end="")
        else:
            print("{:.15f}".format(grid1[y]), end="")

    print("};\n")

    print("const double laborincome_states[{:d}] ={{".format(
        len(grid2)), end="")
    for y in range(len(grid2)):

        # print("{:.15f}\t".format(grid2[y]))
        # grid2[y] = np.exp(grid2[y])

        if (y < len(grid2) - 1):
            print("exp({:.15f}), ".format(grid2[y]), end="")
        else:
            print("exp({:.15f})".format(grid2[y]), end="")

    print("};\n")

    print("const double risk_labor_trans[{:d}][{:d}] = {{\n".format(
        len(grid1), len(grid1)), end="")

    for y in list(range(len(grid1))):
        print("{", end="")
        for k in list(range(len(grid1))):

            if (k < len(grid1) - 1):
                print("{:.15f}, ".format(ytrans[y, k]), end="")
            else:
                print("{:.15f}".format(ytrans[y, k]), end="")

        if (y < len(grid1) - 1):
            print("},\n", end="")
        else:
            print("}};", end="")
    print("\n")

    # print("};\n", end="")

    print("const double risk_pre_states[{:d}] ={{".format(len(grid1)), end="")

    for y in range(len(grid1)):

        # print("{:.15f}\t".format(grid1[y]))

        if (y < len(grid1) - 1):
            print("{:.15f}, ".format(grid1[y]), end="")
        else:
            print("{:.15f}".format(grid1[y]), end="")

    print("};\n")

    print("const double laborincome_pre_states[{:d}] ={{".format(
        len(grid2)), end="")
    for y in range(len(grid2)):

        # print("{:.15f}\t".format(grid2[y]))
        # grid2[y] = np.exp(grid2[y])

        if (y < len(grid2) - 1):
            print("exp({:.15f}), ".format(grid2[y]), end="")
        else:
            print("exp({:.15f})".format(grid2[y]), end="")

    print("};\n")

    print("const double risk_labor_pre_trans[{:d}][{:d}] = {{\n".format(
        len(grid1), len(grid1)), end="")

    for y in list(range(len(grid1))):
        print("{", end="")
        for k in list(range(len(grid1))):

            if (k < len(grid1) - 1):
                print("{:.15f}, ".format(ytrans[y, k]), end="")
            else:
                print("{:.15f}".format(ytrans[y, k]), end="")

        if (y < len(grid1) - 1):
            print("},\n", end="")
        else:
            print("}};", end="")
    print("\n")

    sys.stdout = original_stdout  # Reset the standard output to its original value
