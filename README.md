# NumericalMethods
An implementation in C of the numerical methods from S.D. Conte's "Elementary Numerical Analysis" (1980 Brazilian version published by "Editora Globo"). All source code is located in the **/source** folder.

**THIS PROJECT IS UNDER DEVELOPMENT AND HAS NOT BEEN THOROUGHLY TESTED YET**

*Everything seems to be working so far though*

## Nonlinear Equation Solvers
The following is a list of the iterative root finding methods implemented for the purpose of solving an equation of the form ![equation](https://latex.codecogs.com/png.latex?f%28x%29%20%3D%200):
- Bisection Method
- Linear Iteration
- Aitken's Δ squared process
- Newton's Method
- Secant Method
- False Position Method (*Regula Falsi*)
- Müller's Process

## Nonlinear Equation System Solvers
The following is a list of the iterative methods implemented for the purpose of solving a system of the form ![equation](https://latex.codecogs.com/png.latex?%5Cbegin%7Balign*%7D%20f_1%28x_1%2C%20%5Cdots%2C%20x_n%29%20%26%3D%200%20%5C%5C%20f_2%28x_1%2C%20%5Cdots%2C%20x_n%29%20%26%3D%200%20%5C%5C%20%5Cvdots%20%5C%5C%20f_n%28x_1%2C%20%5Cdots%2C%20x_n%29%20%26%3D%200%20%5C%5C%20%5Cend%7Balign*%7D)

- Linear Iteration
- Aitken's Δ squared process

The following methods are planned to be added in the future one matrix methods are implemented:
- Newton's Method
- Accelerated Pseudo-Newton's Method (Algorithm 2.17)

## Interpolation Methods
The following is a list of the interpolation methods implemented:
- Lagrange Polynomials
- Aitken Interpolation
- Ascending Newton's Method of Finite Differences
- Descending Newton's Method of Finite Differences
- Stirling's Formula
- Everett's Formula

