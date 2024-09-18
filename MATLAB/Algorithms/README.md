Algorithm Classification
===

There are many algorithms in TIGRE and it can be hard to chose were to start. 
So let me introduce a short classification of the algorithms in TIGRE and a small suggestion on where to start. 
In any case, always check the demos, as all features and algorithms of TIGRE should be showcased there. 


There are two main ways we can classify the algorithms. First is by the way they solve the data minimization problem, 
i.e. the way they make sure the image reocnstructed matches the measured data. 
The second is by their regularization, i.e. by the additional constraints that we add to the problem. A common one, and the most
used one in TIGRE is Total Variation (TV), that tries to make the images "piecewise flat". 


# Gradient descend

Tradititonal iterative algorithms use gradient descend-like algorithm (kaczman method, ART), to solve the image. 

Non regularized:
   - SIRT
   - OS-SART
   - SART
Not based on the maths of gradient descend, but similar (solves the proximal of the problem, using gradient descend) 
   - FISTA
   
Regularized with TV (check demo Algorithms04)
   - ASD-POCS
   - B-ASD-POCS-beta
   - OS-ASD-POCS
   - AwASDS-POCS
   - PCSD
   - AwPCSD
   - SART-TV (this is practically FISTA-TV)
   
**Were do I start?**: For any problem, the first iterative algorithm to try is possibly OS-SART. Good convergence rate with relatively good speed per iteartion. 
For TV, start with ASD-POCS, but caution, the hyperparameters have massive influence on the quality of the reconstruction. 

# Maximum likelihood

This assumes the noise and the data follows Poisson distribution, so mostly useful for very low dose scans. 

   - MLEM
   
# Krylov Subspace algorithms

These are fast coverging iterative algorithms. The have some issues with regards of semiconvergence and loss of ortogonality, so in some cases they may not produce best results,
but the main advantage is that few iterations should produce a good image. 

Non regularized:
     - CGLS
	 - LSQR
	 - LSMR
	 - BA-GMRES
	 - AB-GMRES
Tikhonov regularization
	 - hybrid LSQR 
	 - LSMR with non-zero lambda
TV regularization
	 - IRN-CGLS-TV
	 - hybrid-fLSQR-TV
	 
**Were do I start?**: LSQR. If LSQR doesn't provide stable solutions, BA-GMRES and AB-GMRES are supposed to fix that, but they require many copies (one per iteration) of the image, so they are very memory consuming. For regularized solutions, IRN-CGLS-TV, as the hybrid-fLSQR-TV only works for very small images,
it requires large computational and memory resources. 
   
 