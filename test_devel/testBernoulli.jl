using Nemo

RR, x = PolynomialRing(Nemo.QQ, "x")

n = 64 #degree
P = zero(RR)
bernoulli_cache(n)
for k = 0:n
    global P
    coefficient = (Nemo.binomial(n,k))*(bernoulli(n-k))
    P = P + coefficient*x^k
end #P is now the Bernoulli polynomial of degree 64

using Ccluster

#compute the roots in [-1,1] + i[-1,1]
bInit = [fmpq(0,1),fmpq(0,1),fmpq(2,1)] #box centered in 0 + sqrt(-1)*0 with width 100
precision = 53                          #get clusters of size 2^-53
Res = ccluster(P, bInit, precision);
#compute all the roots
Res = ccluster(P, precision);

using CclusterPlot #only if you have installed CclusterPlot.jl

plotCcluster(Res)

