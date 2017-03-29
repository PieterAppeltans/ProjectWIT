# ProjectWIT
This repository contains our solution to the assignment of the course [Project Mathematical Engineering](https://onderwijsaanbod.kuleuven.be/syllabi/e/H0T46AE.htm) in the master of Mathematical Engineering at the KU Leuven [https://onderwijsaanbod.kuleuven.be/opleidingen/e/CQ_52354207.htm#activetab=diploma_omschrijving]. The complete assignement is included in the repository.

Makes use of the [Triangle mesh generator/Show me](http://www.cs.cmu.edu/~quake/triangle.html) (files included in this project).
Makes use of the [Eigen framework](http://eigen.tuxfamily.org/) for matrices, matrix operations and system solving in C++.


## Assignement
The aim of the course [Project Mathematical Engineering] (Spring 2017) is to solve a real-life mathematical problem in groups of two. There were 3 different topics to chose from: "Finite element model of respiration of pears", "Bayesian inversion of insulin model" and "Design optimization of cooling device". We worked out the first option.
### Finite element model of respiration of pears
The goal of this assignement was to simulate the respiration metabolism of pears. This is useful because to store the pears in such a way that quality is maintained, the oxygen and carbon dioxide concentrations inside the pear must be controlled. This respiration metabolism can be modeled using two coupled non-linear differential equations defined on a three-dimensional bounded domain. To solve this differential equations we had to use a finite element method. Code must be in Fortran or C++.
The complete assignement can be found in pear-wit-project-2016.pdf, additional information over realistic values of the constants can be foun  in pear-wit-project-2016-additions.pdf.

### How did we tackled the problem
To get insight in the problem, we started with a literature study. We found and used folowing papers, books, and websites:
* The Finite element method for Engineers (K.H. Huebner et al.)
* Remarks around 50 lines of Matlab: short finite element implementation (J. Alberty, C.Carstensen and S.A. Funken)
* http://www.cs.cmu.edu/~quake/triangle.html: Triangle, A Two-Dimensional Quality Mesh Generator and Delaunay Triangulator. J. R. Shewchuk

In a following stap we simplified the problem by assumming that the pear is cyclindrical symmetric. To uitbuiten this property we transformed the differential equations to cyclindrical coordinates(r,z,theta) and making the concentrations independent of theta. Next we applied partial integration to get a weak formulation of the problem. TODO ... (2 - 3 weeks)

For the triangulation, we were allowed to use an existing package. We choose for Triangle because ... .

Next we started implementing the elements in Matlab ... (1-2 weeks).

Implementation in C++. (1-2 weeks)

Checking + Debugging (1-2 weeks)

Improvements (1 week)

GUI (1 week)


## C++

The cpp folder contains a completed implementation of the assignment written in c++ (a dense and a sparse version).
The sparse version can use the normal Newton-Raphson method, or a quasi-Newton-Raphson method (only calculating
the Jacobian once).

## GUI

The gui folder contains all Python code (tkinter) and bash scripts needed to run the GUI (python tk_interface.py).

## Matlab

The matlab folder contains all Matlab code to solve the assignment in Matlab.
The correctness of the code is analytically verified for a linearized version of the problem,
using a circle/sphere as figure and fem.m as solver. 

## Mesh generation (triangle)

The mesh folder contains the triangle package, bash scripts to generate meshes, Python scripts to generate figures
to be meshed, and calculated meshes (\*.1.\*).

## Theory (FEM)

Derivation of the FEM method for this specific problem, on which all implemented methods are based.

## Run it your own
* Download Triangle
* Download Eigen
* Edit code


Pieter Appeltans & Lennart Bulteel
