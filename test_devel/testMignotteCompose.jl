using Nemo

R, x = PolynomialRing(QQ, "x")

n = 64 #degree
k = 10
a = 140
# P = x^n - (fmpq(a,1)x -fmpq(1,1))^k
P = x^n - ((fmpq(a,1)x -fmpq(1,1))^k)*((fmpq(a,1)x +fmpq(1,1))^k)

# print("P: $P\n");
precision = 53

using Ccluster

Res = ccluster(P, precision, verbosity="silent");

using CclusterPlot #only if you have installed CclusterPlot.jl

plotCcluster(Res) #use true instead of false to focus on clusters

