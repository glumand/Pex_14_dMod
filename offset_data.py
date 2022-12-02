import numpy as np

data = np.loadtxt('geometric_mean.csv', skiprows=1, delimiter=',')
data[:, 0] = data[:, 0] - data[0, 0]
np.savetxt('geometric_mean_offset.csv', data, delimiter=',',
           header='seconds,B_au,A_au,BB_au,AB_au,AA_au,BBB_au,ABB_au,AAB_au,AAA_au')
