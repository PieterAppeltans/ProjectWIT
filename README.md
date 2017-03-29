# ProjectWIT

Makes use of the [Triangle mesh generator/Show me](http://www.cs.cmu.edu/~quake/triangle.html) (files included in this project).
Makes use of the [Eigen framework](http://eigen.tuxfamily.org/) for matrices, matrix operations and system solving in C++.

Solution for the task of the course
[Project Mathematical Engineering](https://onderwijsaanbod.kuleuven.be/syllabi/e/H0T46AE.htm)
in the master of Mathematical Engineering at the KU Leuven. Assignment included.

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
to be meshed, and calculated meshes (*.1.*).

## Theory (FEM)

Derivation of the FEM method for this specific problem, on which all implemented methods are based.
