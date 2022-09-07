import csv
import numpy as np
import matplotlib.pyplot as plt


path_name = "./csv/"
file_name = "VF5"
file = open(path_name+file_name+'.csv', 'r')
reader = csv.reader(file, delimiter=',')
file_header = next(reader)
file_varnum = len(file_header)

data = np.array(list(reader)).astype(float)
file_length = len(data[:, 1])

plt.plot(data[:, 1], data[:, 2:8])
# plt.show()
plt.savefig("./figure/"+file_name+".pdf")
plt.close()


path_name = "./csv/"
file_name = "policy5"
file = open(path_name+file_name+'.csv', 'r')
reader = csv.reader(file, delimiter=',')
file_header = next(reader)
file_varnum = len(file_header)
data = np.array(list(reader)).astype(float)
file_length = len(data[:, 1])

plt.plot(data[:, 1], data[:, 2:8])
# plt.show()
plt.savefig("./figure/"+file_name+".pdf")
plt.close()

# path_name = "./csv/"
# file_name = "dist4"
# file = open(path_name+file_name+'.csv', 'r')
# reader = csv.reader(file, delimiter=',')
# file_header = next(reader)
# file_varnum = len(file_header)
# data = np.array(list(reader)).astype(float)
# file_length = len(data[:, 1])

# plt.plot(data[:, 1], data[:, 2:8])
# # plt.show()
# plt.savefig("./figure/"+file_name+".pdf")
# plt.close()
