#Written by: Jack Zhan
#Exploration of Kmeans on Python

#Create random clustered data

d = 2  # feature space dimension
k = 3  # number of clusters
npercluster = 15
n = k*npercluster
r = .1
# generate random points in unit square
from numpy import *
centers = random.rand(k,d)
print(centers)

F = empty((n,d))
for i in range(k):
    # create a cluster
    start =     i*npercluster
    stop  = (i+1)*npercluster
    F[start:stop,:] = centers[i] + r*(2*random.rand(npercluster,d)-1)

import matplotlib.pyplot as pl

pl.subplot(111,aspect=1)
colors = 'rgb'
for i in range(k):
    start =     i*npercluster
    stop  = (i+1)*npercluster
    pl.plot(F[start:stop,0],F[start:stop,1],'o',color=colors[i]);

pl.savefig('kmeans1.jpg')
pl.close()

# initial choice of means
#means = random.rand(k,d)
means = F[random.choice( range(n), k, replace=False )]
#print(means)

oldassignments = k*ones(n,dtype=int)
count = 0
maxcount = 10
while(True):
    count += 1
    if count>maxcount: break
    # compute the cluster assignments
    displacements = empty((n,d,k))
    for i in range(k):
        displacements[:,:,i] = F - means[i]
    sqdistances = (displacements**2).sum(axis=1)
    assignments = argmin( sqdistances, axis=1 )
    print(assignments)
    if all( assignments == oldassignments ): break
    oldassignments[:] = assignments    
        
    # update the means as the centroids of the clusters

    for i in range(k):
        means[i] = F[assignments==i].mean(axis=0)

pl.subplot(111,aspect=1)
colors = 'rgb'
for i in range(k):
    cluster = assignments==i
    pl.plot(F[cluster,0],F[cluster,1],'o',color=colors[i],alpha=0.75);
    pl.plot(means[i][0],means[i][1],'*',color=colors[i],markersize=10,alpha=0.75)
pl.savefig('kmeans2.jpg')
pl.close()