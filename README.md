# NQL - Numerical Quadrature Library

This library implements in MATLAB an automatic numerical integration routine using __Newton-Cotes__ and __Gauss-Legendre__ formulae.

## Main structure
The library consists of
- [x] `simulation.m` which is the interaction script from which the user sets the desired input parameters s.a.:
	* interval of integration _I = [a,b]_ specified through its endpoints _a,b_;
	* integrand, or target function f(x), specified as `kernel`;
	* number of discretisation nodes _N_ in _I_ for the numerical quadrature;
- [x] `main.m` computes the necessary parameters according to the choice of formulae provided by the user (nodes, affine map coefficients, jacobian etc...);
- [x] `legendre` if Gauss-Legendre formulae is selected it computes the nodes of the _N-th_ degree Legendre's polynomial while displaying the recursive polynomials used for its construction;
- [x] `quadrature.m` defines the workflow for the computation of the weigths associated to the Lagrangian basis' functions;
- [x] `lagrangian_weights.m` computes the aforementioned weights;
- [x] `lagrange.m` computes iteratively the _NxN_ coefficient matrix associated to the _N_ Lagrange's polynomials of degree _N-1_;
- [x] `extract_basis.m` for each node computes the coefficients of all monomials associated recursively to Lagrange's polynomial associated to that node;

## Accessories
To support data management and visualisation the following functions have been provided:
- [x] `cmap.m` defines a colormap based on the interpolation of 8 colors and according to the number of Legendre's and Lagrange's functions are defined;
- [x] `doubledot.m` computes the inner product between two vectors but in ascending order of values so to avoid numerical cancellation;

## Suggested benchmarks
With educational purpose in mind one might show that Netwon-Cotes formula has a _order_-of-precision d=N-1_ in accordance with the analysis.  To do that the user must specify
```matlab
integrand = @(x) x.^k;
```
where `k` is a number less or equal to _N-1_. Gauss-Legendre formulae on the other hand shows a numerical convergence also for polynomials of degree _d=2N-1_ which disagrees with the theoretical result. The reason for this is the low accuracy of the calculated roors of Legendre's polynomial computed using linear interpolation. 