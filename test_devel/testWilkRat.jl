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

R, x = PolynomialRing(Nemo.QQ, "x")

n = 128 #degree
P = one(R)
for k = 1:n
    Ptemp = R(x - fmpq(k,n+1))
    P = P*Ptemp
end
N = round(Int,5+ceil(log2(1+log2(degree(P)))))
print("N: $N, 2^$N: $(2^N), N-log2(degree(P)): $(N-log2(degree(P)))\n")

function getAppWilkRat( dest::Ptr{acb_poly}, prec::Int )
    ccall((:acb_poly_set_fmpq_poly, :libarb), 
     Void, (Ptr{acb_poly}, Ptr{fmpq_poly}, Int), 
            dest,          &P,             prec)
end

include("../julia/ccluster.jl")
# bInit = box(Nemo.fmpq(0,1),Nemo.fmpq(0,1),Nemo.fmpq(2,1))
bInit = [Nemo.fmpq(1,6),Nemo.fmpq(1,4),Nemo.fmpq(2,1)]
# bInit = box(Nemo.fmpq(1,3),Nemo.fmpq(1,3),Nemo.fmpq(2,1))
eps = Nemo.fmpq(1,100)
    
Res = Ccluster(getAppWilkRat, bInit, eps, 15, 1);
plotCcluster(Res, bInit, false)