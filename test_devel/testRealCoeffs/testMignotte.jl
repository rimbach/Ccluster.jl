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

R, x = PolynomialRing(QQ, "x")

n = 181 #degree
bitsize = 14
bitsize = Int64(floor(bitsize/2))
print("bitsize: $bitsize\n");
P = x^n - (fmpq(fmpz(2)^(bitsize-1),1)*x -1)^2
# P = x^n - (fmpq(fmpz(2)^(bitsize)-1,1)*x -1)^2
 print("P: $P\n");

function getAppMign( dest::Ptr{acb_poly}, prec::Int )
    ccall((:acb_poly_set_fmpq_poly, :libarb), 
      Void, (Ptr{acb_poly}, Ptr{fmpq_poly}, Int), 
             dest,          &P,             prec)
end
# print("P: $P \n")

# bInit = box(Nemo.fmpq(0,1),Nemo.fmpq(0,1),Nemo.fmpq(1,1))
bInit = [fmpq(0,1),fmpq(0,1),fmpq(100,1)]
# bInit = box(Nemo.fmpq(0,1),Nemo.fmpq(0,1),Nemo.fmpq(1,2))

eps = Nemo.fmpq(1,100)
# eps = Nemo.fmpq(1,1000)

# strategy = 23
strategy = 31 # stop when compact
# strategy = 55
# strategy = 63 # stop when compact
verbosity = 2

Res = ccluster(getAppMign, bInit, eps, strategy, verbosity)
plotCcluster(Res, bInit, true)
