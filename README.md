# Ccluster.jl

Ccluster.jl is a Julia wrapper for Ccluster (https://github.com/rimbach/Ccluster.git)
Ccluster implements a clustering algorithm for univariate polynomials whose
coefficients are complex numbers.

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
See the file examples/mignotte.jl
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
each element of Res being an array which
* second element is a complex disc (defined by the real and
imaginary parts of its center and its radius)
* first element is the sum of multiplicities of the roots in the disk.

If you care about geometry, so do we.
If you have installed CclusterPlot.jl, you can plot the clusters with:
```
using CclusterPlot

plotCcluster(Res, bInit, false)
```
The last argument is a flag telling the function wether to focus 
on clusters (when *true*) or not (when *false*).

### Other example: clustering the roots of a polynomial whose coefficients are roots of polynomials
See the file examples/coeffsBernoulli.jl
#### Find the 64 roots of the Bernoulli polynomial of degree 64
```
using Nemo

R, x = PolynomialRing(Nemo.QQ, "x")

n = 64 #degree
P = zero(R)
bernoulli_cache(n)
for k = 0:n
    coefficient = (binom(n,k))*(bernoulli(n-k))
    P = P + coefficient*x^k
end #P is now the Bernoulli polynomial of degree 8

using Ccluster

bInit = [fmpq(0,1),fmpq(0,1),fmpq(100,1)] #box centered in 0 + sqrt(-1)*0 with width 100
eps = fmpq(1, fmpz(2)^10)               #eps = 2^-10
verbosity = 0                           #nothing printed
Coeffs = ccluster(P, bInit, eps, verbosity)
```
#### Define an approximation function for the polynomial whose coefficients are the found roots
```
function getApproximation( dest::Ptr{acb_poly}, prec::Int )
    eps = fmpq(1, fmpz(2)^prec)
    Qre = zero(R)
    Qim = zero(R)
    for i=1:n
        btemp = [ Coeffs[i][2][1], Coeffs[i][2][2], 2*Coeffs[i][2][3] ]
        temp = ccluster(P, btemp, eps, 0)
        Qre = Qre + temp[1][2][1]*x^(i-1)
        Qim = Qim + temp[1][2][2]*x^(i-1)
    end
    Ccluster.ptr_set_2fmpq_poly( dest, Qre, Qim, prec )
end
```
#### Cluster the roots
```
bInit = [fmpq(0,1),fmpq(0,1),fmpq(100,1)] #box centered in 0 + sqrt(-1)*0 with width 100
eps = fmpq(1, 100)                      #eps = 1/100
verbosity = 0                           #nothing printed
Roots = ccluster(getApproximation, bInit, eps, 1)
```
Output (total time in s on a Intel(R) Core(TM) i7-7600U CPU @ 2.80GHz):
```
 -------------------Ccluster: ----------------------------------------
 -------------------Input:    ----------------------------------------
|box: cRe: 0                cIm: 0                wid: 100            |
|eps: 1/100                                                           |
|strat: newton tstarOpt predPrec anticip                              |
 -------------------Output:   ----------------------------------------
|number of clusters:                                 63               |
|number of solutions:                                63               |
 -------------------Stats:    ----------------------------------------
|tree depth:                                         20               |
|tree size:                                        2096               |
|total time:                                   3.540887               |
 ---------------------------------------------------------------------
63-element Array{Any,1}:
 Any[1, Nemo.fmpq[-3125//32768, 5125//8192, 375//65536]]    
 ⋮                                                              
 Any[1, Nemo.fmpq[211625//262144, -105125//262144, 375//524288]]
 ```

### Defining an approximation function
**ccluster** takes as input a function prototyped as:
```
function getApproximation( dest::Ptr{Nemo.acb_poly}, p::Int )
```
Here is an example for a polynomial with complex coefficients (see also the file examples/spiral.jl)
```
degr = 64

function getApproximation( dest::Ptr{acb_poly}, prec::Int )
    
    CC = ComplexField(prec)
    R2, y = PolynomialRing(CC, "y")
    res = R2(1)
    for k=1:degr
        modu = fmpq(k,degr)
        argu = fmpq(4*k,degr)
        root = modu*Nemo.exppii(CC(argu))
        res = res * (y-root)
    end
    Ccluster.ptr_set_acb_poly(dest, res)
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

### The bound *eps*

*eps* is a rational number:
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
