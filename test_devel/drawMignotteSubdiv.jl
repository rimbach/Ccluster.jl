using Nemo

R, x = PolynomialRing(QQ, "x")

d=64
a=14
P = x^d - 2*((2^a)*x-1)^2 #mignotte polynomial

function getApproximation( dest::Ptr{acb_poly}, precision::Int )
    
    global P
    
    CC = ComplexField(precision)
    R2, y = PolynomialRing(CC, "y")
    poly = R2(P)
    
    Ccluster.ptr_set_acb_poly(dest, poly)

end

using Ccluster

bInit = [fmpq(0,1),fmpq(0,1),fmpq(4*4,5)] #box centered in 0 + sqrt(-1)*0 with width 100
precision = 53                          #get clusters of size 2^-53
Roots = ccluster(getApproximation, bInit, precision, verbosity="silent")
# 
using CclusterPlot #only if you have installed CclusterPlot.jl
# 
# plotCcluster(Roots, bInit, focus=true) #use true instead of false to focus on clusters

# Res = ccluster(getApproximation, bInit, precision, strat=23, verbosity="brief");

res, dis = Ccluster.ccluster_draw( getApproximation, bInit, precision, strategy="V4", verbosity="silent")
# 
CclusterPlot.plotCcluster_subdiv( res, dis, bInit)

res, dis = Ccluster.ccluster_draw( getApproximation, bInit, precision, verbosity="silent")

CclusterPlot.plotCcluster_subdiv( res, dis, bInit)
