# ProjectWIT
This repository contains our solution to the assignment of the course [Project Mathematical Engineering](https://onderwijsaanbod.kuleuven.be/syllabi/e/H0T46AE.htm) in the master of [Mathematical Engineering at the KU Leuven](https://onderwijsaanbod.kuleuven.be/opleidingen/e/CQ_52354207.htm). The complete assignement is included in the repository.

Makes use of the [Triangle mesh generator/Show me](http://www.cs.cmu.edu/~quake/triangle.html) (files included in this project).
Makes use of the [Eigen framework](http://eigen.tuxfamily.org/) for matrices, matrix operations and system solving in C++.


## Assignment
The aim of the course [Project Mathematical Engineering](https://onderwijsaanbod.kuleuven.be/opleidingen/e/CQ_52354207.htm) (Spring 2017) is to solve a real-life mathematical problem in groups of two. There were 3 different topics to choose from: "Finite element model of respiration of pears", "Bayesian inversion of insulin model" and "Design optimization of cooling device". We worked out the first option.

### Finite element model of respiration of pears
The goal was to simulate the respiration metabolism of pears. To safely store pears and prevent degradation, the oxygen and carbon dioxide concentrations inside the pear must be controlled. This respiration metabolism can be modeled using two coupled non-linear differential equations defined on a three-dimensional bounded domain. To solve these differential equations, we had to use a finite element method (simple finite differential methods would have problems with irregularily shaped objects anyway, so FEM would alwyas be a solid choice). The FEM method had to be implemented in Fortran or C++.

The complete assignement can be found in [pear-wit-project-2016.pdf](https://github.com/PieterAppeltans/ProjectWIT/blob/master/pear-wit-project-2016.pdf), realistic values for the involved constants are given in [pear-wit-project-2016-additions.pdf](https://github.com/PieterAppeltans/ProjectWIT/blob/master/pear-wit-project-2016-additions.pdf).

### How we tackled the problem
To gain insight in the problem, we started with a literature study. We found and used folowing papers, books, and websites:
* The Finite element method for Engineers (K.H. Huebner et al.)
* Remarks around 50 lines of Matlab: short finite element implementation (J. Alberty, C.Carstensen and S.A. Funken)
* http://www.cs.cmu.edu/~quake/triangle.html: Triangle, A Two-Dimensional Quality Mesh Generator and Delaunay Triangulator. J. R. Shewchuk

#### FEM derivation, mesh generation (first 2-3 weeks)
After the literature study, we started with the derivation of the FEM method for this problem. First we simplified the problem by assumming axial (cylindrical) symmetry. To exploit this property, we transformed the differential equations to cylindrical coordinates: x,y,z to r,z,theta. We assumed the variation around the angle theta zero, and thus scrapped all terms differentiating in respect to theta. Next we applied partial integration to get a weak formulation of the problem. For each differential equation, the weak formulation was reformulated to the traditional Ax=b system form, where A and b consisted of several integrals. Then, the integrals were further worked out (giving rise to several coefficient matrices and). Finally, the Jacobian of the system was needed for implementing Newton-Raphson. (2 - 3 weeks)

Simultaneously we worked out mesh generation, by splitting the workload efficiently. We were allowed to use an existing package, and chose for Triangle. It is an easy to use 2D meshing package, with the possibility to set quality constraints (maximum triangle area and minimum angle). Triangle meshing is one of the simplest forms of 2D meshing, but solid in this case because e.g. squares would have more trouble fitting the irregular outside boundary of the pear.

#### Implementation (1-2 weeks)

Next we started implementing the elements in Matlab and C++ (dense version) simultaneously. While implementing the FEM method in Matlab was not required, it allowed for checking the results between Matlab and C++ and Matlab is a language we were very familiar with because of past courses. Not to mention writing and plotting the FEM method in Matlab is way shorter. 

Initially we used the ublas wrapper of [boost](http://www.boost.org/) (this code can be found in [cpp/deprecated](https://github.com/PieterAppeltans/ProjectWIT/tree/master/cpp/deprecated). But because the [Eigen](http://eigen.tuxfamily.org/) ublas interface has linear system solvers built-in, and Eigen allows for neater code (it has operators like * defined for matrices), we switched to Eigen.  

#### Checking and debugging (1-2 weeks)

Ofcourse the implementation wasn't immediately correct. The Matlab version was implemented first, and showed numerical instability. So to check the coefficient matrices (associated with the integrals) we worked out all the matrices on paper for a mesh consisting of a single triangle, and compared it to the coefficient matrices in Matlab and C++. Furthermore, we calculated some of the values for the respiration kinetic functions and their derivatives to r and z on paper, and again compared them to the values of Matlab and C++.

Once those were all correct, there still were problems with numerical stability. However, up until then we'd always worked with less than 1000 triangles in the mesh. By driving it up to about 1500 triangles and choosing the best conditioned case (low temperature good oxygen and carbon dioxide concentrations), our implementation seemed to be working as intended.

But because this is numerical software, seeming is not always being. Solutions that seem about right can be very wrong. A way to check whether the implementation is correct, is comparing to an analytical solution. We using a sphere as the problem shape, which allows further simplification. Using Mupad, the analytical solution for a fixed theta and z was worked out and compared to the numerical solution. Also, the solution for half a circle (half a slice of the sphere) was numerically determined. This solution had to show symmetry around the center point. Our FEM implementation passed both checks!

#### Further improvements and GUI (1 week)

Now, those coefficient matrices are very sparse. Especially for very fine meshes, which give rise to large coefficient matrices, a sparse implementation can be of big help. That's why we added a sparse implementation, which speeds the code up with a factor 5-10.

Additionally, instead of Newton-Raphson, a quasi-Newton-Rapshon method can be used. By calculating the Jacobian only once and never recalculating it, the calculation time is reduced by another factor of 5. This comes with a big but: if the mesh isn't fine enough, it will have trouble converging. A smarter implementation would combine both and refrain from calculating the Jacobian until it detects trouble converging. However, this falls out of the scope of the assignment.

For easier demonstration, a graphical user interface was made using tkinter and matplotlib. It allows setting the temperature and concentration conditions, mesh generation requirements and the figure to mesh, and which FEM method to use. It shows a contourplot. 

## The repository

### C++

The cpp folder contains a completed implementation of the assignment written in c++ (a dense and a sparse version).
The sparse version can use the normal Newton-Raphson method, or a quasi-Newton-Raphson method (only calculating
the Jacobian once).

### GUI

The gui folder contains all Python code and bash scripts needed to run the GUI (python tk_interface.py).

### Matlab

The matlab folder contains all Matlab code to solve the assignment in Matlab.
The correctness of the code is analytically verified for a linearized version of the problem,
using a circle/sphere as figure and fem.m as solver. 

### Mesh generation (triangle)

The mesh folder contains the triangle package, bash scripts to generate meshes, Python scripts to generate figures
to be meshed, and calculated meshes (\*.1.\*).

### Theory (FEM)

Derivation of the FEM method for this specific problem, on which all implemented methods are based. Pictures of a few relevant pages in The Finite element method for Engineers (K.H. Huebner et al.).

## Run it on your own
* Install prerequisites: Eigen, Python (matplotlib, numpy, tkinter) and Matlab
* Clone this repository
* $ python gui/tk_interface.py 


Pieter Appeltans & Lennart Bulteel, MIT License
