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

R, x = PolynomialRing(QQ, "x")

d = 64 #degree
bitsize = 14
bitsize = Int64(floor(bitsize/2))
P=x^d - 2*((2^bitsize)*x-1)^2

bInit = [fmpq(-1,1),fmpq(0,1),fmpq(100,1)] #box centered in 0 + sqrt(-1)*0 with width 4
precision = 53                          #get clusters of size 2^-53
    
Res = ccluster(P, bInit, precision, verbosity="silent");

using CclusterPlot #only if you have installed CclusterPlot.jl

plotCcluster(Res, bInit, focus=true) #use true instead of false to focus on clusters
