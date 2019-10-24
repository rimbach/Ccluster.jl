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

# include("functionsTstar.jl")

R, x = PolynomialRing(QQ, "x")

n = 64 #degree
k = 10
a = 140
# P = x^n - (fmpq(a,1)x -fmpq(1,1))^k
P = x^n - ((fmpq(a,1)x -fmpq(1,1))^k)*((fmpq(a,1)x +fmpq(1,1))^k)

 print("P: $P\n");
function getAppMign( dest::Ptr{acb_poly}, prec::Int )
    ccall((:acb_poly_set_fmpq_poly, :libarb), 
      Void, (Ptr{acb_poly}, Ptr{fmpq_poly}, Int), 
             dest,          &P,             prec)
end
# print("P: $P \n")
include("../julia/ccluster.jl")

bInit = [Nemo.fmpq(0,1),Nemo.fmpq(0,1),Nemo.fmpq(100,1)]
eps = Nemo.fmpq(1,10)


Res = Ccluster(getAppMign, bInit, eps, 15, 0);
plotCcluster(Res, bInit, false)
