import csv
from sys import api_version
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
    # plt.ylim((0, 1))
    plt.xlim((0.0, 50))
    # plt.show()
    plt.savefig("./figure/"+file_name+".pdf")
    plt.close()


file_name = []
file_name.append("dist21old,pi=0.030000,wage=1.800000,rf=0.040000,Psize=50,rho_c=1.500000,rho_w=3.000000,Ksize=200,Kmax=500.000000,relaxVF=0.000000,beta=0.750000,corr=1.000000,Ssize=5")
file_name.append("dist21old,pi=0.030000,wage=1.800000,rf=0.040000,Psize=100,rho_c=1.500000,rho_w=3.000000,Ksize=500,Kmax=500.000000,relaxVF=0.000000,beta=0.750000,corr=1.000000,Ssize=5")
file_name.append("policy21old,pi=0.030000,wage=1.800000,rf=0.040000,Psize=100,rho_c=1.500000,rho_w=3.000000,Ksize=500,Kmax=500.000000,relaxVF=0.000000,beta=0.750000,corr=1.000000,Ssize=5")
file_name.append("dist21noextra,pi=0.030000,wage=1.800000,rf=0.040000,Psize=50,rho_c=1.500000,rho_w=3.000000,Ksize=200,Kmax=500.000000,relaxVF=0.000000,beta=0.750000,corr=1.000000,Ssize=5")
file_name.append("dist21noextra_save,pi=0.030000,wage=1.800000,rf=0.040000,Psize=50,rho_c=1.500000,rho_w=3.000000,Ksize=200,Kmax=500.000000,relaxVF=0.000000,beta=0.750000,corr=1.000000,Ssize=5")
file_name.append("dist21noextra_save,pi=0.030000,wage=1.000000,std_l=0.250000,rf=0.040000,Psize=50,rho_c=1.500000,rho_w=3.000000,Ksize=200,Kmax=500.000000,relaxVF=0.000000,beta=0.920000,corr=1.000000,Ssize=5")
file_name.append("dist21noextra_save,pi=0.000000,wage=1.000000,std_l=0.250000,rf=0.040000,Psize=50,rho_c=1.500000,rho_w=3.000000,Ksize=200,Kmax=500.000000,relaxVF=0.000000,beta=0.920000,corr=1.000000,Ssize=5")
file_name.append("dist21nouti,pi=0.000000,wage=1.000000,std_l=0.250000,rf=0.040000,Psize=50,rho_c=1.500000,rho_w=3.000000,Ksize=200,Kmax=500.000000,relaxVF=0.000000,beta=0.920000,corr=1.000000,Ssize=1")
file_name.append("dist21nouti,pi=0.030000,wage=1.000000,std_l=0.250000,rf=0.040000,Psize=50,rho_c=1.500000,rho_w=3.000000,Ksize=200,Kmax=500.000000,relaxVF=0.000000,beta=0.920000,corr=1.000000,Ssize=1")
file_name.append("dist21nouti,pi=0.030000,wage=1.000000,std_l=0.250000,rf=0.040000,Psize=50,rho_c=1.500000,rho_w=3.000000,Ksize=200,Kmax=500.000000,relaxVF=0.000000,beta=0.920000,corr=1.000000,Ssize=3")

graph(file_name[-1])


# for i in range(len(file_name_string)):
#     graph(file_name_string[i])


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
