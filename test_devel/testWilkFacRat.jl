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

n = 64 #degree
P = one(R)
for k = 1:n
    Ptemp = R(x - fmpq(Nemo.fac(k),Nemo.fac(n)+1))
    P = P*Ptemp
end

# print("P: $P\n\n")

function getAppWilkFacRat( dest::Ptr{acb_poly}, prec::Int )
    ccall((:acb_poly_set_fmpq_poly, :libarb), 
      Void, (Ptr{acb_poly}, Ptr{fmpq_poly}, Int), 
             dest,          &P,             prec)
end

include("../julia/ccluster.jl")
# bInit = box(Nemo.fmpq(0,1),Nemo.fmpq(0,1),Nemo.fmpq(2,1))
bInit = [Nemo.fmpq(1,6),Nemo.fmpq(1,4),Nemo.fmpq(2,1)]
eps = Nemo.fmpq(0,1)
    
Res = Ccluster(getAppWilkFacRat, bInit, eps, 15, 1);
plotCcluster(Res, bInit, false)