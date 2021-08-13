using Nemo
using Ccluster

R, x = PolynomialRing(Nemo.QQ, "x")

n = 16 #degree
P0 = one(R)
P1 = (fmpq(1,1)-x)
P=zero(R)
for i = 1:n-1
    global P0, P1, P
    P = (2*i+1-x)*P1 - (i^2)*P0
    P0=P1
    P1=P
end
print("degree(P): $(degree(P))\n")
# print("P: $P\n")
    
Res = ccluster(P, 53);
