# 2016-ICML-gromov-wasserstein

Matlab code to reproduce some of the results of the paper

> Gabriel Peyré, Marco Cuturi, Justin Solomon
> Gromov-Wasserstein Averaging of Kernel and Distance Matrices
> ICML 2016

Only implements the shape interpolation application.

- compute_gw_barycenters.m: implements the computation of the barycenters.
- perform_gw_sinkhorn.m: implement the computation of a coupling solving the GW problem.
- test_distmat_interpol.m: main script to launch the computation. Save figures as .eps files.
- batch_barycenter_distances.m: launch computation in batch mode.
- data/: contains binary image of shapes.
- mesh2d/: the amazing code of Darren Engwirda for meshing a 2-D domain.
- toolbox/: various helper functions.
