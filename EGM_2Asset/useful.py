
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
