# REO: A Relative-Entropy Optimization Solver
This repository is authored and maintained by Parikshit Shah and Venkat Chandrasekaran.

The relative entropy function in a jointly convex function that appears in a diverse set of applications. Many optimization problems of interest involve the relative entropy cone (or equivalently, the exponential cone). In this code base we implement an interior-point method for solving relative entropy optimization problems.

An important special case of interest is solving signomial optimization problems, for which lower bounds can be obtained via a relative-entropy based convex relaxation. We implement a wrapper for the same in this codebase.

For further details, we refer the reader to the following papers:

V. Chandrasekaran and P. Shah, Relative Entropy Optimization and its Applications, Mathematical Programming, accepted. 
[link: http://users.cms.caltech.edu/~venkatc/cs_rep_preprint.pdf]

V. Chandrasekaran and P. Shah, Relative Entropy Relaxations for Signomial Optimization, SIAM Journal on Optimization, accepted. 
[link: http://users.cms.caltech.edu/~venkatc/cs_signomial_preprint.pdf]
