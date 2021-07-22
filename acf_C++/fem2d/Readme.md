## fem2d

**C++** implementation of **acf/fem2d**, using **Eigen** as linear algebra and linear system solver. ([Remarks around 50 lines of Matlab: short finite element implementation](https://www.math.hu-berlin.de/~cc/cc_homepage/download/1999-AJ_CC_FS-50_Lines_of_Matlab.pdf))

This program employs the finite element method to calculate a numerical solution $U$ which approximates the solution $u$ to the 2D Laplace problem ($P$) with mixed boundary conditions.

![image-20210721200809937](C:\Users\陈威\AppData\Roaming\Typora\typora-user-images\image-20210721200809937.png)

### Build

***

```shell
mkdir build
cd build
cmake ..
make
```

#### Dependencies

* [Eigen](https://eigen.tuxfamily.org/index.php?title=Main_Page): we use Eigen's dense, sparse linear algebra, linear solver (SimplicialLDLT).  **To use custom linear solvers, you can implement a new interface (subclass) to out** `LinSysSolver` **class**.

### How to use

***

you can modify implementation of function: `stima3`, `stima4`, `f`, `g` and `u_d` to define the problem to be solved. These five functions are all defined in `src/Mesh.cpp`

* `stima3`: calculate local stiffness matrix of triangle elements
* `stima4`: calculate local stiffness matrix of parallelogram elements
* `f`: this function is called with the coordinates of points in $\Omega$ and it returns the volume forces at these locations
* `g`: this function is called with the coordinates of points in $\Gamma_{\mathrm{N}}$ and returns the corresponding stresses
* `u_d`: this function is called with the coordinates of points in $\Gamma_{\mathrm{D}}$ and returns the values at the corresponding locations

#### Input

fem2d load data from `input` directory.

* `coordinates.dat`: vertices for mesh
* `elements3.dat`: triangle elements for mesh
* `elements4.dat`: parallelogram elements for mesh
* `dirichlet.dat`: edges in dirichlet boundary condition
* `nuemann.dat`: edges in nuemann boundary condition

#### Output

fem2d compute $u$ and write it to `output/u.txt`