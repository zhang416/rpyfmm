# RPYFMM 

RPYFMM computes the Rotne-Prager-Yamakawa tensor using the adaptive Fast Multipole Method (FMM). The implementation leverages the Dynamic Adaptive System for Hierarchical Multipole Methods (DASHMM) infrastructure and operates on both shared and distributed memory architectures. 

RPYFMM was created at the Center for Research in Extreme Scale Technologies (CREST) and Indiana University, and was supported by the National Science Foundation. 

RPYFMM package contains the following files: 

1. src/rpy.cc: Source file implementing routines for computing the gradient and Hessian of the spherical harmonics based multipole and local expansions for Laplace potential. 
2. include/rpy.h: Header file implementing the Rotne-Prager-Yamakawa kernel template with DASHMM library. 
3. demo/demo.cc: A program demonstrating the use of RPYFMM. 
4. cmake/Modules: CMake module file that automatically downloads the prerequsite DASHMM library during build. 
5. INSTALL: Installation instruction. 
6. LICENSE: BSD-3 license file. 
7. AUTHORS: Author list. 


