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

# t = [(x-0)]
nbiter = 4

function iterate( table, n, CC )
    
    index = 1
    bound = length(table)
    while index<=bound
        p = splice!(table,1)
#         modu = fmpq(1,3^(2*n-1))
        modu = fmpq(1,4^(2*(n-1)))
        argu = fmpq(1,3)
        root1 = modu*Nemo.exppii(CC(2*argu))
        root2 = modu*Nemo.exppii(CC(4*argu))
        root3 = modu*Nemo.exppii(CC(6*argu))
        append!(table, [p-root1, p-root2, p-root3])
        index+=1
    end
end

function getAppClusters( dest::Ptr{acb_poly}, prec::Int )

    CC = ComplexField(prec)
    R2, y = PolynomialRing(CC, "y")
    
    table = [(y-0)]
    for i=1:nbiter
        iterate( table, i, CC )
    end
    
    res = R2(1)
    for i=1:length(table)
        res = res*table[i]
    end
    
    ccall((:acb_poly_set, :libarb), Void,
                (Ptr{acb_poly}, Ptr{acb_poly}, Int), 
                 dest,         &res,          prec)

end

include("../julia/ccluster.jl")
bInit = [fmpq(0,1),fmpq(0,1),fmpq(3,1)]

eps = fmpq(1,2000)

Res = Ccluster(getAppClusters, bInit, eps, 15, 2);
plotCcluster(Res, bInit, false)
