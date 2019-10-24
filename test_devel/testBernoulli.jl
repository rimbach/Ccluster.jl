#
#  Copyright (C) 2018 Remi Imbach
#
#  This file is part of Ccluster.
#
#  Ccluster is free software: you can redistribute it and/or modify it under
#  the terms of the GNU Lesser General Public License (LGPL) as published
#  by the Free Software Foundation; either version 2.1 of the License, or
#  (at your option) any later version.  See <http://www.gnu.org/licenses/>.
#

using Nemo
using Ccluster

R, x = PolynomialRing(Nemo.QQ, "x")

n = 64 #degree
P = zero(R)
Nemo.bernoulli_cache(n)
for k = 0:n
    global P
    coefficient = (Nemo.binom(n,k))*(Nemo.bernoulli(n-k))
    P = P + coefficient*x^k
end

P2 = zero(R)
Nemo.bernoulli_cache(n)
for k = 0:15
    global P2
    coefficient = (Nemo.binom(n,k))*(Nemo.bernoulli(n-k))
    P2 = P2 + coefficient*x^k
end

function getAppBern( dest::Ref{acb_poly}, prec::Int )
    ccall((:acb_poly_set_fmpq_poly, :libarb), Cvoid,
                (Ref{acb_poly}, Ref{fmpq_poly}, Int), dest, P, prec)
end

print("P: $P\n");

bInit = [fmpq(0,1),fmpq(0,1),fmpq(1,1)] #box centered in 0 + sqrt(-1)*0 with width 4
precision = 53                          #get clusters of size 2^-53
    
# local tests
Res = ccluster(getAppBern, bInit, precision, strategy="default", verbosity="brief");
# Res = ccluster(P,          bInit, precision, strategy="default", verbosity="brief");
# Res = ccluster(P, P2,      bInit, precision, strategy="default", verbosity="brief");

# global tests
# Res = ccluster(getAppBern, precision, strategy="default", verbosity="brief");
# Res = ccluster(P,          precision, strategy="default", verbosity="brief");
# Res = ccluster(P, P2,      precision, strategy="default", verbosity="brief");

# using CclusterPlot #only if you have installed CclusterPlot.jl

# plotCcluster(Res, bInit, focus=true) #use true instead of false to focus on clusters

# plotCcluster(Res) #use true instead of false to focus on clusters

