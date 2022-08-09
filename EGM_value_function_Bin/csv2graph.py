import csv
import numpy as np
import matplotlib.pyplot as plt


path_name = "./EGM_value_function_Bin/csv/"
file_name = "dist_WF"
file = open(path_name+file_name+'.csv', 'r')
reader = csv.reader(file, delimiter=',')
file_header = next(reader)
file_varnum = len(file_header)
data = np.array(list(reader)).astype(float)
file_length = len(data[:, 1])

plt.plot(data[:, 1], data[:, 2:8])
# plt.show()
plt.savefig("./EGM_value_function_Bin/figure/dist_WF.pdf")
plt.close()


path_name = "./EGM_value_function_Bin/csv/"
file_name = "policy_WF"
file = open(path_name+file_name+'.csv', 'r')
reader = csv.reader(file, delimiter=',')
file_header = next(reader)
file_varnum = len(file_header)
data = np.array(list(reader)).astype(float)
file_length = len(data[:, 1])

plt.plot(data[:, 1], data[:, 2:8])
# plt.show()
plt.savefig("./EGM_value_function_Bin/figure/policy_WF.pdf")
plt.close()
