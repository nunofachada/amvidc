Introduction
============

The AMVIDC algorithm is presented in detail in the following 
publication:

"Spectrometric differentiation of yeast strains using Minimum Volume 
Increase and Minimum Direction Change clustering criteria", N. Fachada, 
M.T. Figueiredo, V.V. Lopes, R.C. Martins and A.C. Rosa. Pattern 
Recognition Letters, 2014 (IN PRESS)

Data format
-----------

Typically, data is presented as a set of samples (or points), each with
a constant number of dimensions. As such, for the rest of this guide,
data matrices are considered to be in the following format:

-   *m* x *n*, with *m* samples (points) and *n* dimensions (variables)

Many times the number of dimensions is too high, making clustering
inefficient. When this occurs, one can reduce data dimensionality using
a number of techniques. In this work, PCA and SVD (which are very
similar) are used via the Matlab native `princomp` and `svd`/`svds`
functions.

Generating data
---------------

This code was inspired on the differentiation of spectrometric data. 
However, to further validate the clustering algorithms, synthetic
data sets can be generated with the generateData function. This function 
generates data in the *m* x *n* format, with *m* samples (points) and 
*n* dimensions (variables) according to a set of parameters, which are 
explained in the source code.

Running the algorithm
=====================

This algorithm is based on AHC, using Minimum Volume Increase (MVI) and 
Minimum Direction Change (MDC) clustering criteria. This algorithm can 
be tested using the clusterdata_amvidc.m function:

    idx = clusterdata_amvidc(X, k, idx_init);

where **X**, **k** and **idx\_init** are the typical data matrix,
maximum number of clusters and initial clustering, respectively. Initial
clustering is required so that all possible new clusters have volume, a
requirement for MVI. `clusterdata_amvidc` function has many optional
parameters, with reasonable defaults, as specified in the following
table:

  Parameter    | Default                |  Options/Description
  ------------ | ---------------------- | ------------------------------------------------------------------------------------------------------
  *volume*     | ‘convhull’             |  Volume type: ‘ellipsoid’ or ‘convhull’
  *tol*        | 0.01                   |  Tolerance for minimum volume ellipse calculation (‘ellipsoid’ volume only)
  *dirweight*  | 0                      |  Direction weight in last iteration (0 means MDC linkage is ignored)
  *dirpower*   | *dirweight* \> 0       |  Convergence power to dirweight (higher values make convergence steeper and occurring more to the end)
  *dirtype*    | ‘svd’                  |  Direction type: ‘pca’, ‘svd’
  *nvi*        | true                   |  Allow negative volume increase?
  *loglevel*   | 3 (show warnings only) |  Log level: 0 (show all messages) to 4 (only show critical errors), default is 3 (show warnings)

For example, to perform clustering using ellipsoid volume taking into
account direction change, where cluster direction is determined using
PCA, one would do:

    idx = clusterdata_mvidc(X, k, idx_init, 'volume', 'ellipsoid', 'dirweight',0.5, 'dirpower', 4, 'dirtype', 'pca');

As specified, the `clusterdata_amvidc` function requires initial clusters
which, if joined, produce new clusters with volume. There are two
clustering functions appropriate for this (but others can be used):

-   initClust.m - Performs very simple initial clustering based
    on AHC with single linkage (nearest neighbor) and user defined
    distance. Each sample is associated with the same cluster of its
    nearest point. Allows to define a minimum size for each cluster,
    distance type (as supported by Matlab `pdist`) and the number of
    clusters which are allowed to have less than the minimum size.
-   pddp.m - Perform PDDP (principal direction divisive
    clustering) on input data. This implementation always selects the
    largest cluster for division, with the algorithm proceeding while
    the division of a cluster yields sub-clusters which can have a
    volume.

Analysis of results
===================

F-score
-------

In this work, the [F-score](http://en.wikipedia.org/wiki/F1_score)
measure was used to evaluate clustering results. The source:fscore.m
function was developed for this purpose. To run this function, do:

    eval = fscore(idx, numclasses, numclassmembers);

where:

-   **idx** - *m* x *1* vector containing the cluster indices of each
    point (as returned by the clustering functions)
-   **numclasses** - Correct number of clusters
-   **numclassmembers** - Vector with the correct size of each cluster
    (or a scalar if all clusters are of the same size)

The `fscore` function returns:

-   **eval** - Value between 0 (worst case) and 1 (perfect clustering)

Plotting clusters
-----------------

Sometimes visualizing how an algorithm grouped clusters can provide
important insight on its effectiveness. Also, it may be important to
visually compare an algorithm’s clustering result with the correct
result. These are the goals of the plotClusters.m function, which
can show two clustering results in the same image (e.g. the correct one
and one returned by an algorithm). You can run `plotClusters` in the
following way:

    h_out = plotClusters(X, dims, idx_marker, idx_encircle, encircle_method, h_in);

where:

-   **X** - Data matrix, *m* x *n*, with m samples (points) and n
    dimensions (variables)
-   **dims** - Number of dimensions (2 or 3)
-   **idx_marker** - Clustering result ^1^ to be shown directly in
    points using markers
-   **idx_encircle** - Clustering result ^1^ to be shown using
    encirclement/grouping of points
-   **encircle_method** - How to encircle the **idx*encircle**
    result: ‘convhull’ (default), ‘ellipsoid’ or ‘none’
-   **h_in** - (Optional) Existing figure handle where to create
    plot

^1^ *m* x *1* vector containing the cluster indices of each point

The `plotClusters` function returns:

-   **h_out** - Figure handle of plot


