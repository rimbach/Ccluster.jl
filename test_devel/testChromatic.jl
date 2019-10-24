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

n = 63 #degree
k=3
A = one(R)
C = x-fmpq(2,1)
for i = 1:k
    D = (x-fmpq(2,1))*C*C + fmpq(2,1)*(x-fmpq(1,1))*A*A
    A = (x-fmpq(1,1))^3*A*A*A*A + C*C*C*C
    C = C*C*D
end
P = A*x*(x-1)
print("degree(P): $(degree(P))\n")
# N = round(Int,4+ceil(log2(1+log2(degree(P)))))
# print("N: $N, 2^$N: $(2^N), N-log2(degree(P)): $(N-log2(degree(P)))\n")
# print("P: $P\n")

function getAppChrom( dest::Ptr{acb_poly}, prec::Int )
    ccall((:acb_poly_set_fmpq_poly, :libarb), Void,
                (Ptr{acb_poly}, Ptr{fmpq_poly}, Int), dest, &P, prec)
end

# bInit = [fmpq(-1,1),fmpq(0,1),fmpq(4,1)]
bInit = [fmpq(3,2),fmpq(0,1),fmpq(10,1)]
eps = fmpq(1,10)
# eps = fmpq(1,2^(53))
    
Res = ccluster(getAppChrom, bInit, eps, 23, 2);
plotCcluster(Res, bInit, false) 
