using Nemo
using Ccluster

R, x = PolynomialRing(Nemo.QQ, "x")

# n = 63 #degree
# k=Int(floor(log2(n+1)))
k=5
r=2
Pm2 = one(R)
Pm1 = x
P = one(R)
for i = 2:k
    global P, Pm2, Pm1
    P = Pm1^r + x*Pm2^(r^2)
    Pm2 = Pm1
    Pm1 = P
end

print("P: $P\n")
print("degree: $(degree(P))\n\n");

bInit = [fmpq(-1,1),fmpq(0,1),fmpq(10,1)] #box centered in 0 + sqrt(-1)*0 with width 4
precision = 53                          #get clusters of size 2^-53
    
Res = ccluster(P, bInit, precision, verbosity="brief");

using CclusterPlot #only if you have installed CclusterPlot.jl

plotCcluster(Res, bInit, focus=true) #use true instead of false to focus on clusters
