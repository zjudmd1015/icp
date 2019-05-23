# ICP-mini-project
Small collaborative mini-project focused on implementing the Iterative Closest Point algorithm.

## Introduction
C++ implementation of 3-dimensional ICP (Iterative Closest Point) algorithm. A Python implementation of ICP by [Clay Flannigan](https://github.com/ClayFlannigan/icp) was referred and rewritten into a C++ version in this project.

- ICP finds a best fit rigid body transformation between two point sets. Correspondence between the points is not assumed. Included is an SVD-based least-squared best-fit algorithm for corresponding point sets.
- In this version, exhaustive search method is used to find the nearest neighbor for each point.
- Eigen library is used for matrices operations.

> [ICP Wiki](https://en.wikipedia.org/wiki/Iterative_closest_point) | [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page)

Use command `$ g++ icp.cpp test.cpp` to make, and use command `$ ./a.out` to run it.
