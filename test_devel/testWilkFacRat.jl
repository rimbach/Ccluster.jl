using Nemo

R, x = PolynomialRing(Nemo.QQ, "x")

n = 16 #degree
P = one(R)
for k = 1:n
    global P
    Ptemp = R(x - fmpq(Nemo.factorial(big(k)),Nemo.factorial(big(n))+1))
    P = P*Ptemp
end

# print("P: $P\n\n")
precision = 53

using Ccluster

Res = ccluster(P, precision, verbosity="silent");

using CclusterPlot #only if you have installed CclusterPlot.jl

plotCcluster(Res) #use true instead of false to focus on clusters
