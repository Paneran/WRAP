import csv
import matplotlib.pyplot as plt
import numpy as np

with open('data/stm_costa.csv', newline='') as csvfile:
    data = list(csv.reader(csvfile))

dim = np.shape(data)
data = np.reshape(data, dim[0]*dim[1])
data = data[500:1000]
data = np.asarray(data, dtype=float)

plt.plot(data)
plt.show()
