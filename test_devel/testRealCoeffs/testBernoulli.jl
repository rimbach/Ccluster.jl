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
using CclusterPlot

R, x = PolynomialRing(Nemo.QQ, "x")

n = 128 #degree
P = zero(R)
Nemo.bernoulli_cache(n)
for k = 0:n
    coefficient = (Nemo.binom(n,k))*(Nemo.bernoulli(n-k))
    P = P + coefficient*x^k
end
# N = round(Int,4+ceil(log2(1+log2(degree(P)))))
# print("N: $N, 2^$N: $(2^N), N-log2(degree(P)): $(N-log2(degree(P)))\n")


function getAppBern( dest::Ptr{acb_poly}, prec::Int )
    ccall((:acb_poly_set_fmpq_poly, :libarb), Void,
                (Ptr{acb_poly}, Ptr{fmpq_poly}, Int), dest, &P, prec)
end

bInit = [fmpq(0,1),fmpq(0,1),fmpq(100,1)]
eps = fmpq(1,10)

strategy = 23
# strategy = 31 # stop when compact
strategy = 55
# strategy = 63 # stop when compact
verbosity = 2
Res = ccluster(getAppBern, bInit, eps, strategy, verbosity);
# plotCcluster(Res, bInit, true)

