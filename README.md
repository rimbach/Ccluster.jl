# Ccluster.jl

Ccluster.jl is a julia wrapper for Ccluster (https://github.com/rimbach/Ccluster.git)
Ccluster implements a clustering algorithm for univariate polynomials whose
coefficients are algebraic complex numbers.

## Brief description

The main function provided by Ccluster.jl is **ccluster**.
It takes as input
a polynomial *P*, 
a square complex box *B*, 
a bound *eps*, 
a strategy and
a verbosity flag.

It outputs a set of *natural clusters* of roots together with the sum of multiplicities
of the roots in each cluster.
A cluster is a complex disc *D* containing at least one root, 
and it is natural when *3D* contains the same roots
than *D*.
Each root of *P* in *B* is in exactly one cluster of the output, and clusters may contain
roots of *P* in *2B*

The implemented algorithm is described here:
https://dl.acm.org/citation.cfm?id=2930939

Please cite:
https://link.springer.com/chapter/10.1007/978-3-319-96418-8_28
if you use it in your research.

## Installation

You can install it within Julia with:

```
Pkg.clone("https://github.com/rimbach/Ccluster.jl.git")
Pkg.build("Ccluster")
```
Ccluster depends on Nemo that will be automatically installed.

For graphical outputs, install the package CclusterPlot with
```
Pkg.clone("https://github.com/rimbach/CclusterPlot.jl.git")
```
CclusterPlot depends on PyCall and PyPlot, and requires that matplotlib is installed
on your system.
It is heavy both to install and to load.

## Usage

### Simple example: clustering the roots of a Mignotte-like polynomial
see the file examples/mignotte.jl
```
using Nemo

R, x = PolynomialRing(QQ, "x")

d=64
a=14
P = x^d - 2*((2^a)*x-1)^2 #mignotte polynomial

using Ccluster

bInit = [fmpq(0,1),fmpq(0,1),fmpq(4,1)] #box centered in 0 + sqrt(-1)*0 with width 4
eps = Nemo.fmpq(1,100)                  #eps = 1/100
verbosity = 0                           #nothing printed

Res = ccluster(P, bInit, eps, verbosity)
```
Res in an array of couples (sum of multiplicity, disc):
```
63-element Array{Any,1}:
 Any[1, Nemo.fmpq[975//1024, 1025//1024, 15//2048]]      
 ⋮                                                      
 Any[1, Nemo.fmpq[-2995//4096, 4805//4096, 15//8192]] 
 Any[2, Nemo.fmpq[0, 0, 15//16384]]                     # the cluster with sum of multiplicity 2
 Any[1, Nemo.fmpq[6935//8192, -8955//8192, 15//16384]]
 Any[1, Nemo.fmpq[6935//8192, 8955//8192, 15//16384]]
```
If you have installed CclusterPlot.jl, you can plot the clusters with:
```
using CclusterPlot

plotCcluster(Res, bInit, false)
```
when replacing false with true, the graphical output is focus on the clusters.

### Defining a polynomial
**ccluster** takes as input a function prototyped as:
```
function getApproximation( dest::Ptr{Nemo.acb_poly}, p::Int )
```

Here is an example for a polynomial with rational coefficients:
```
using Nemo
R, x = PolynomialRing(Nemo.QQ, "x")

function getApproximation( dest::Ptr{acb_poly}, p::Int )
    P = R(fmpq(1,3)*x^3 + fmpq(1,2)*x^2 + fmpq(1,1)*x + 1)
    ccall((:acb_poly_set_fmpq_poly, :libarb), Void,
                (Ptr{acb_poly}, Ptr{fmpq_poly}, Int), dest, &P, prec)
end
```

### Defining an initial box  
The initial box is an array of three Nemo.fmpq defining respectively
the real part of the center,
the imaginary part of the center and
the width of the box.

The following code:
```
bInit = [fmpq(0,1),fmpq(0,1),fmpq(150,1)]
```
defines a box centered in 0+*i*0 with width 150.

### The bound *eps* and the strategy
The strategy is an integer that can take its value in:
* 7: in which case the clusters are natural clusters with radius at most *eps*.
* 15: in which case the clusters are either 
  natural clusters with radius at most *eps* 
  or natural clusters with exactly one root of multiplicity 1. 

The *eps* is a rational number:
```
eps = fmpq(1,100)
```
or
```
eps = fmpq(1, fmpz(2)^(-53))
```

Unless you know what you are doing, setting *eps* to 0 is a very bad idea.

### The verbosity flag
Use 0 unless you want statistics on the solving process.

### Example
Below is a minimal example of use of **ccluster** for computing the 
roots of Bernoulli polynomials.

```
using Nemo
using Ccluster

R, x = PolynomialRing(Nemo.QQ, "x")

n = 8 #degree
P = zero(R)
Nemo.bernoulli_cache(n)
for k = 0:n
    coefficient = (Nemo.binom(n,k))*(Nemo.bernoulli(n-k))
    P = P + coefficient*x^k
end #P is now the Bernoulli polynomial of degree 8

function getAppBern( dest::Ptr{acb_poly}, prec::Int )
    ccall((:acb_poly_set_fmpq_poly, :libarb), Void,
                (Ptr{acb_poly}, Ptr{fmpq_poly}, Int), dest, &P, prec)
end

bInit = [fmpq(0,1),fmpq(0,1),fmpq(150,1)] #the initial box: an array [cRe, cIm, w]
eps = fmpq(1,100) #an escape bound
#compute the roots    
Res = ccluster(getAppBern, bInit, eps, 15, 0)
``` 

The output of this code is the array of clusters of roots of P:
```
8-element Array{Any,1}:
 Any[1, Nemo.fmpq[1125//4096, 0, 1125//16384]]           
 Any[1, Nemo.fmpq[6375//8192, 0, 1125//16384]]           
 Any[1, Nemo.fmpq[375//256, -1125//4096, 1125//16384]]   
 Any[1, Nemo.fmpq[375//256, 1875//8192, 1125//16384]]    
 Any[1, Nemo.fmpq[-1875//4096, -1125//4096, 1125//16384]]
 Any[1, Nemo.fmpq[-1875//8192, 0, 1125//16384]]          
 Any[1, Nemo.fmpq[20625//16384, 0, 1125//32768]]         
 Any[1, Nemo.fmpq[-1875//4096, 1875//8192, 1125//32768]]
```
each element of the latter array being an array which
* second element is a complex disk (defined by the real and
imaginary parts of its center and its radius)
* first element is the sum of multiplicities of the roots in the disk.

If you care about geometry, so do we. You can plot the clusters with:
```
plotCcluster(Res, bInit, false)
```
The last argument is a flag telling the function wether to focus 
on clusters (when *true*) or not (when *false*).
