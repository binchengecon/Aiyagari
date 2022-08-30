import matplotlib.pyplot as plt
import numpy as np
import pickle

betapar = 0.96
alphapar = 0.36
deltapar = 0.08
rhopar = 3.0
labor = 1.0219882

size_k = 500
size_m = 7

kmin = 0.0
kmax = 500.0


scale1 = 1.6
grmin = (kmin / scale1) - 1.0
exponen = np.log((kmax / scale1) - grmin) / (size_k - 1)

sstates = np.array([np.exp(-0.600000000000000), np.exp(-0.400000000000000), np.exp(-0.200000000000000),
                   np.exp(0.000000000000000), np.exp(0.200000000000000), np.exp(0.400000000000000), np.exp(0.600000000000000)])


strans = np.array([
    [0.046746218637144, 0.217937777267117, 0.397822606398702, 0.266386738072197,
        0.065169922261456, 0.005754191945237, 0.000182545418147],
    [0.023199661746751, 0.149524091076020, 0.369020347246402, 0.333823905199677,
        0.110578117872631, 0.013276146769082, 0.000577730089437],
    [0.010548958644399, 0.093657511915497, 0.312761268311836, 0.382193227897354,
        0.171253064028981, 0.027919224002876, 0.001666745199056],
    [0.004387354018187, 0.053538402796357, 0.242163972572887, 0.399820541225137,
        0.242163972572887, 0.053538402796357, 0.004387354018187],
    [0.001666745199056, 0.027919224002876, 0.171253064028981, 0.382193227897354,
        0.312761268311837, 0.093657511915497, 0.010548958644399],
    [0.000577730089436, 0.013276146769082, 0.110578117872631, 0.333823905199677,
        0.369020347246403, 0.149524091076020, 0.023199661746751],
    [0.000182545418147, 0.005754191945237, 0.065169922261456, 0.266386738072197, 0.397822606398702, 0.217937777267117, 0.046746218637144]])


def U(x):
    return x**(1-rhopar)/(1-rhopar)


def phi(x):
    return (scale1 * (np.exp(exponen * (x)) + grmin))


def getwage(r):
    return (1.0 - alphapar) * (alphapar / (rrate + deltapar))**(alphapar / (1.0 - alphapar))


def HJB_update(res, V):
    crit = 1
    eps = 1e-8

    while crit > eps:

        Vold = V.copy()
        # print(Vold[0, 0])
        for i in range(size_k):
            for y in range(size_m):
                temp = np.zeros((size_k, 1))

                for ifuture in range(size_k):
                    afuture = phi(ifuture)
                    # print(afuture)
                    cfuture = res[i, y]-afuture
                    if cfuture > 0:
                        tempfuture = 0
                        for yfuture in range(size_m):
                            tempfuture += strans[y, yfuture] * \
                                Vold[ifuture, yfuture]

                        temp[ifuture] = U(cfuture)+betapar*tempfuture
                    else:
                        temp[ifuture] = -np.inf

                V[i, y] = np.max(temp)

        crit = np.max(abs(V-Vold))
        print(crit)

    return V


VF = np.zeros((500, 7))
dVF = np.zeros((500, 7))
resource = np.zeros((500, 7))
grid = np.zeros((500, 1))


rrate = 0.040237086402090
wagerate = getwage(rrate)

for i in range(500):
    for y in range(7):
        # tempnext = 0
        # dtempnext = 0
        # for ynext in range(7):
        #     tempnext = tempnext + strans[y, ynext] * VF[i, ynext]
        #     dtempnext = dtempnext + strans[y, ynext] * dVF[i, ynext]
        grid[i] = phi(i)
        resource[i, y] = (1+rrate)*phi(i)+wagerate*sstates[y]
        VF[i, y] = U(sstates[y] * wagerate + (1 + rrate) * phi(i))


VF = HJB_update(resource, VF)


with open("./EGM_value_function_Bin/pickle_data/VF1e-8", "wb") as f:
    pickle.dump(VF, f)

plt.plot(grid[:], VF[:, :])
plt.savefig("./EGM_value_function_Bin/figure/VF_VFI.pdf")
plt.close()
