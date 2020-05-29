#!/usr/bin/env python

import math
import sys
import os
from sklearn.cluster import KMeans as kmeans
import numpy as np
import matplotlib.pyplot as plt

contig_file = sys.argv[1] # HaploFlow fasta file
#duplication_ratio_file = sys.argv[2] # quast file pointing to the duplication ratio (will be num_clusters)
min_len = 500 # minimal length of contigs to be clustered
seed = 1234

#with open(duplication_ratio_file, 'r') as dr:
#    for line in dr:
#        if (line.startswith("Assemblies")):
#            continue
#        name, duplication_ratio = line.strip().split()

nr_clusters = 2#math.ceil(float(duplication_ratio))

flow_map = {} # map contigs to their flow
with open(contig_file, 'r') as c:
    for line in c:
        if (not line.startswith('>')):
            length = len(line)
            if (length < min_len):
                del flow_map[nr]
            continue
        contig, nr, f, flow, cc, ccnr = line.strip().split("_")
        flow_map[nr] = flow
       
data = np.fromiter(flow_map.values(), dtype=float)
reshaped = data.reshape(-1,1)

clusters = kmeans(n_clusters = int(nr_clusters), random_state = seed).fit(reshaped)
labels = clusters.labels_
centers = clusters.cluster_centers_

to_write = {}
with open(contig_file, 'r') as c:
    for line in c:
        if (line.startswith('>')):
            contig, nr, f, flow, cc, ccnr = line.strip().split("_")
            label = clusters.predict(np.array([flow]).reshape(-1,1))
            outfile = contig_file.rstrip(".fa") + ("/contig%s.fa") % int(centers[label][0])
        if (outfile not in to_write):
            to_write[outfile] = [line.strip()]
        else:
            to_write[outfile].append(line.strip())

os.makedirs(outfile.rsplit("/",1)[0])
for out in to_write:
    write = "\n".join(to_write[out])
    with open(out, 'a+') as o:
        o.write(write)
        

#points = [[] for i in range(int(nr_clusters))]
#i = 0
#for l in labels:
#    points[l].append(data[i])
#    i += 1
#
#plt.hist(points)
#plt.savefig("test.png")

