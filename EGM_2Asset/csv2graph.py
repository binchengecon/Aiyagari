import csv
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.rcParams["lines.linewidth"] = 1
mpl.rcParams["savefig.bbox"] = "tight"
mpl.rcParams["figure.figsize"] = (26, 15)
mpl.rcParams["font.size"] = 15
mpl.rcParams["legend.frameon"] = False


def graph(file_name):
    path_name = "./csv/"
    file = open(path_name+file_name+'.csv', 'r')
    reader = csv.reader(file, delimiter=',')
    file_header = next(reader)
    file_varnum = len(file_header)
    data = np.array(list(reader)).astype(float)
    file_length = len(data[:, 1])

    plt.plot(data[:, 1], data[:, 2:])
    plt.ylim((0, 1))
    plt.xlim((0,50))
    # plt.show()
    plt.savefig("./figure/"+file_name+".pdf")
    plt.close()


file_name_string = ["Portfolio,pi=0.000000,Psize=100", "Portfolio,pi=0.000000,Psize=300", "Portfolio,pi=0.000000,Psize=500", "Portfolio,pi=0.000000,Psize=2000", "Portfolio,pe=e-4,premium=0.000000,wage=0.340000,rf=0.030000,Psize=100",
                    "Portfolio,pe=e-5,premium=0.000000,wage=0.340000,rf=0.030000,Psize=100", "Portfolio,premium=0.000000,wage=0.340000,rf=0.030000", "Portfolio,pe=e-6,premium=0.000000,wage=0.340000,rf=0.030000,Psize=100", "Portfolio,pe=e-6,premium=0.050000,wage=0.340000,rf=0.030000,Psize=100"]
file_name_string.append(
    "Portfolio,pe=e-6,premium=0.010000,wage=0.340000,rf=0.030000,Psize=100")
file_name_string.append(
    "Portfolio,pe=e-6,premium=0.000050,wage=0.340000,rf=0.030000,Psize=100")
file_name_string.append(
    "Portfolio,pe=e-6,std=0.015,premium=0.000050,wage=0.640000,rf=0.060000,Psize=100")
file_name_string.append(
    "Portfolio,pe=e-6,std=0.2,premium=0.010000,wage=0.640000,rf=0.060000,Psize=100")
file_name_string.append(
    "Portfolio,pe=e-6,std=0.01,premium=0.000050,wage=0.340000,rf=0.030000,Psize=100")
file_name_string.append(
    "Portfolio,pe=e-6,std=0.01,premium=0.000050,wage=0.340000,rf=0.030000,Psize=2")
file_name_string.append(
    "Portfolio,pe=e-6,std=0.2,premium=0.010000,wage=0.340000,rf=0.030000,Psize=2")
file_name_string.append(
    "Portfolio,pe=e-6,std=0.2,premium=0.010000,wage=0.640000,rf=0.030000,Psize=2")
file_name_string.append(
    "Portfolio,pe=e-6,std=0.01,premium=0.000050,wage=0.640000,rf=0.030000,Psize=2")
file_name_string.append(
    "Portfolio7adjust,pe=e-6,std=0.01,premium=0.000050,wage=0.340000,rf=0.030000,Psize=100")
file_name_string.append(
    "Portfolio10,pe=e-6,std=0.01,premium=0.000050,wage=0.340000,rf=0.030000,Psize=100")
file_name_string.append(
    "Portfolio10,pe=e-6,std=0.01,premium=0.000050,wage=0.340000,rf=0.040237,Psize=100,rho_c=2.000000,rho_w=3.000000")
file_name_string.append(
    "Portfolio10,pe=e-6,std=0.01,premium=0.000050,wage=0.340000,rf=0.040237,Psize=100,rho_c=1.500000,rho_w=4.000000")
file_name_string.append(
    "Portfolio12,pe=e-9,std=0.01,premium=0.000050,wage=0.800000,rf=0.030000,Psize=100,rho_c=3.000000,Ksize=300,relaxVF=0.000000")
file_name_string.append(
    "Portfolio12,pe=e-9,std=0.01,premium=0.000000,wage=0.800000,rf=0.030000,Psize=20,rho_c=3.000000")
file_name_string.append(
    "Portfolio12,pe=e-9,std=0.01,premium=0.000000,wage=0.800000,rf=0.030000,Psize=100,rho_c=3.000000")

file_name_string.append(
    "Portfolio11,pe=e-9,std=0.01,premium=0.000050,wage=0.600000,rf=0.030000,Psize=100,rho_c=3.000000")
file_name_string.append(
    "Portfolio12ad,pe=e-9,std=0.01,premium=0.000050,wage=0.800000,rf=0.030000,Psize=100,rho_c=3.000000,Ksize=300,relaxVF=0.000000")
file_name_string.append(
    "Portfolio12ad,pe=e-9,std=0.2,premium=0.000050,wage=0.800000,rf=0.030000,Psize=100,rho_c=3.000000,Ksize=300,relaxVF=0.000000,beta=0.800000")
file_name_string.append(
    "Portfolio12ad,pe=e-9,std=0.2,premium=0.000500,wage=0.800000,rf=0.030000,Psize=100,rho_c=3.000000,Ksize=300,relaxVF=0.000000,beta=0.800000")
file_name_string.append(
    "Portfolio12ad,pe=e-9,std=0.2,premium=0.005000,wage=0.800000,rf=0.030000,Psize=100,rho_c=3.000000,Ksize=300,relaxVF=0.000000,beta=0.800000")
file_name_string.append(
    "Portfolio12ad,pe=e-9,std=0.2,premium=0.050000,wage=0.800000,rf=0.030000,Psize=100,rho_c=3.000000,Ksize=300,relaxVF=0.000000,beta=0.800000")
file_name_string.append(
    "Portfolio12ad,pe=e-9,std=0.2,premium=0.010000,wage=0.800000,rf=0.030000,Psize=100,rho_c=3.000000,Ksize=300,relaxVF=0.000000,beta=0.800000")
file_name_string.append(
    "Portfolio12ad,pe=e-9,std=0.2,premium=0.010000,wage=0.800000,rf=0.030000,Psize=100,rho_c=3.000000,Ksize=300,relaxVF=0.000000,beta=0.900000")
file_name_string.append(
    "Portfolio13,pe=e-9,std=0.2,premium=0.010000,wage=0.800000,rf=0.030000,Psize=100,rho_c=3.000000,rho_w=3.000000,Ksize=300,relaxVF=0.500000,beta=0.900000")
file_name_string.append("Portfolio13,pe=e-9,std=0.2,premium=0.010000,wage=0.800000,rf=0.030000,Psize=100,rho_c=2.000000,rho_w=5.000000,Ksize=300,relaxVF=0.500000,beta=0.900000")
file_name_string.append("Portfolio13,pe=e-9,std=0.2,premium=0.010000,wage=0.800000,rf=0.030000,Psize=100,rho_c=2.000000,rho_w=10.000000,Ksize=300,relaxVF=0.500000,beta=0.900000")
file_name_string.append("Portfolio13,pe=e-9,std=0.2,premium=0.010000,wage=0.800000,rf=0.030000,Psize=100,rho_c=2.000000,rho_w=20.000000,Ksize=300,relaxVF=0.500000,beta=0.900000")



for i in range(len(file_name_string)):
    graph(file_name_string[i])


# path_name = "./csv/"
# file_name = "VF7"
# file = open(path_name+file_name+'.csv', 'r')
# reader = csv.reader(file, delimiter=',')
# file_header = next(reader)
# file_varnum = len(file_header)

# data = np.array(list(reader)).astype(float)
# file_length = len(data[:, 1])

# plt.plot(data[:, 1], data[:, 2:])
# # plt.show()
# plt.savefig("./figure/"+file_name+".pdf")
# plt.close()


# path_name = "./csv/"
# file_name = "policy7"
# file = open(path_name+file_name+'.csv', 'r')
# reader = csv.reader(file, delimiter=',')
# file_header = next(reader)
# file_varnum = len(file_header)
# data = np.array(list(reader)).astype(float)
# file_length = len(data[:, 1])

# plt.plot(data[:, 1], data[:, 2:])
# # plt.show()
# plt.savefig("./figure/"+file_name+".pdf")
# plt.close()

# path_name = "./csv/"
# file_name = "dist7"
# file = open(path_name+file_name+'.csv', 'r')
# reader = csv.reader(file, delimiter=',')
# file_header = next(reader)
# file_varnum = len(file_header)
# data = np.array(list(reader)).astype(float)
# file_length = len(data[:, 1])

# plt.plot(data[:, 1], data[:, 2:])
# # plt.show()
# plt.savefig("./figure/"+file_name+".pdf")
# plt.close()


# path_name = "./csv/"
# file_name = "Portfolio,premium=0.000000,wage=0.340000,rf=0.030000"
# file = open(path_name+file_name+'.csv', 'r')
# reader = csv.reader(file, delimiter=',')
# file_header = next(reader)
# file_varnum = len(file_header)
# data = np.array(list(reader)).astype(float)
# file_length = len(data[:, 1])

# plt.plot(data[:, 1], data[:, 2:])
# # plt.show()
# plt.savefig("./figure/"+file_name+".pdf")
# plt.close()

# path_name = "./csv/"
# file_name = "Portfolio,pi=0.000000"
# file = open(path_name+file_name+'.csv', 'r')
# reader = csv.reader(file, delimiter=',')
# file_header = next(reader)
# file_varnum = len(file_header)
# data = np.array(list(reader)).astype(float)
# file_length = len(data[:, 1])

# plt.plot(data[:, 1], data[:, 2:])
# # plt.show()
# plt.savefig("./figure/"+file_name+".pdf")
# plt.close()
