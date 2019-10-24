using Nemo
using Ccluster

R, x = PolynomialRing(Nemo.QQ, "x")

n = 63 #degree
k=Int(floor(log2(n+1)))
P = one(R)
for i = 1:k
    global P
    P = x*P*P +1
end

bInit = [fmpq(-1,1),fmpq(0,1),fmpq(10,1)] #box centered in 0 + sqrt(-1)*0 with width 4
precision = 53                          #get clusters of size 2^-53
    
Res = ccluster(P, bInit, precision, verbosity="silent");

using CclusterPlot #only if you have installed CclusterPlot.jl

plotCcluster(Res, bInit, focus=true) #use true instead of false to focus on clusters
