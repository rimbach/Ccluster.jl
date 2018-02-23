#Ccluster.jl

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
end

function getAppBern( dest::Ptr{acb_poly}, prec::Int )
    ccall((:acb_poly_set_fmpq_poly, :libarb), Void,
                (Ptr{acb_poly}, Ptr{fmpq_poly}, Int), dest, &P, prec)
end

bInit = [fmpq(0,1),fmpq(0,1),fmpq(150,1)]
eps = fmpq(1,100)
    
Res = ccluster(getAppBern, bInit, eps, 15, 0);
plotCcluster(Res, bInit, false)
``` 
