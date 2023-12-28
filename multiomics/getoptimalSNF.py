# %%
import snf
import sklearn
from snf import datasets 
from snf import designeddatasets   
import numpy as np
from sklearn.cluster import spectral_clustering
from sklearn.metrics import v_measure_score

multi = designeddatasets.load_digits()
affinity_networks = snf.make_affinity(multi.data, metric='euclidean', K=20, mu=0.5)
fused_network = snf.snf(affinity_networks, K=20)
best, second = snf.get_n_clusters(fused_network)
labels = spectral_clustering(fused_network, n_clusters=5)
optima =  v_measure_score(labels, multi.labels)



# %%
