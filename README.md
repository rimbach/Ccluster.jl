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

## Installation
It depends on Nemo, PyCall and PyPlot.
You can install it within Julia with:

```
Pkg.clone("https://github.com/rimbach/Ccluster.jl.git")
Pkg.build("Ccluster")
```

## Usage

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
