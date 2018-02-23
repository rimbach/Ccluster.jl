# Ccluster.jl

Ccluster.jl is a julia wrapper for Ccluster (https://github.com/rimbach/Ccluster.git)

You can install it with:

julia> Pkg.clone("https://github.com/rimbach/Ccluster.jl.git")

Here is a minimal exemple:

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

#in the current version, ccluster takes a function in input
#this function fill the acb_poly dest with a prec-bit approximation of P
function getAppBern( dest::Ptr{acb_poly}, prec::Int )
    ccall((:acb_poly_set_fmpq_poly, :libarb), Void,
                (Ptr{acb_poly}, Ptr{fmpq_poly}, Int), dest, &P, prec)
end

#the initial box: an array [cRe, cIm, w]
#defining the square complex box centered in (cRe,cIm)
#with width w

bInit = [fmpq(0,1),fmpq(0,1),fmpq(150,1)]

#an escape bound
eps = fmpq(1,100)
    
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

```
plotCcluster(Res, bInit, false)
```
will draw the clusters of roots.