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

degr = 64 #degree


function getApp( dest::Ptr{acb_poly}, prec::Int )
    
    CC = ComplexField(prec)
    R2, y = PolynomialRing(CC, "y")
    res = R2(1)
    for k=1:degr
        modu = fmpq(k,degr)
        argu = fmpq(4*k,degr)
        root = modu*Nemo.exppii(CC(argu))
        res = res * (y-root)
    end
    ccall((:acb_poly_set, :libarb), Void,
                (Ptr{acb_poly}, Ptr{acb_poly}, Int), 
                 dest,         &res,          prec)

end

# CCt = ComplexField(53)
# R2t, y = PolynomialRing(CCt, "y")
# p = R2t(1)
# getApp(&p, 53)
# print("res: $(getApproximation(53))\n\n")
# print("degree(res): $(degree( getApproximation(53) )) \n\n")

include("../julia/ccluster.jl")
# bInit = box(Nemo.fmpq(0,1),Nemo.fmpq(0,1),Nemo.fmpq(10,1))
bInit = [Nemo.fmpq(1,6),Nemo.fmpq(1,4),Nemo.fmpq(3,1)]

# withNewton, stopWhenCompact, useTstarOptimised, predictPrec, anticipate, countSols
quo = fmpz(2)
quo = quo^53
print("quo: $quo\n")
eps = Nemo.fmpq(1,quo)
eps = fmpq(1,100)
    
Res = Ccluster(getApp, bInit, eps, 15, 1);
plotCcluster(Res, bInit, false)