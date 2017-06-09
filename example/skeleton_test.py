from pypcalg import skeleton
import numpy as np

data = np.genfromtxt('input_data.csv', delimiter=',')

adj, sepsets = skeleton(data.tolist(), alpha=0.01, return_sepset=True)

print np.array(adj)

sepsets_dict = {}
for i, r in enumerate(sepsets):
	for j, c in enumerate(r):
		sepsets_dict[(j, i+1)] = c

print sepsets_dict