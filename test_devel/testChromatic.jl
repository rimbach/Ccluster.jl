using Nemo
using Ccluster

R, x = PolynomialRing(Nemo.QQ, "x")

k=3
A = one(R)
C = x-fmpq(2,1)
for i = 1:k
    global A, C
    D = (x-fmpq(2,1))*C*C + fmpq(2,1)*(x-fmpq(1,1))*A*A
    A = (x-fmpq(1,1))^3*A*A*A*A + C*C*C*C
    C = C*C*D
end
P = A*x*(x-1)
print("degree(P): $(degree(P))\n")
# print("P: $P\n")

#compute all the roots
Res = ccluster(P, 53, verbosity="stats"); 
using CclusterPlot #only if you have installed CclusterPlot.jl
plotCcluster(Res)
