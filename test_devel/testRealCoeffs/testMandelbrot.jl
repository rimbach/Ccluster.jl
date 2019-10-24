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

n = 127 #degree
k=Int(floor(log2(n+1)))
P = one(R)
for i = 1:k
    P = x*P*P +1
end
# degree(P)
# N = round(Int,4+ceil(log2(1+log2(degree(P)))))
# print("N: $N, 2^$N: $(2^N), N-log2(degree(P)): $(N-log2(degree(P)))\n")
# print("P: $P\n")

function getAppMan( dest::Ptr{acb_poly}, prec::Int )
    ccall((:acb_poly_set_fmpq_poly, :libarb), Void,
                (Ptr{acb_poly}, Ptr{fmpq_poly}, Int), dest, &P, prec)
end

# bInit = [fmpq(-1,1),fmpq(0,1),fmpq(4,1)]
# bInit = [fmpq(-16,10),fmpq(0,1),fmpq(1,10)]
bInit = [fmpq(0,1),fmpq(0,1),fmpq(10,1)]
# eps = fmpq(1,10)
eps = fmpq(1,2^(53))
    
# strategy = 23
strategy = 31 # stop when compact
# strategy = 55
# strategy = 63 # stop when compact
verbosity = 2

Res = ccluster(getAppMan, bInit, eps, strategy, verbosity);
plotCcluster(Res, bInit, true) 
