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
using Printf

# import Ccluster: printClusters

# function printClusters(out, clusters, prec::Int)
#     R::ArbField = RealField(prec)
#     C::AcbField = ComplexField(prec)
#     write(out, "-------------------------------------\n")
#     for index in 1:length(clusters)
#         mult = clusters[index][1]
#         B    = clusters[index][2]
#         if isa( clusters[index][2], Array{fmpq, 1})
#             if length(clusters[index][2])==3
#                 B = Ccluster.box( clusters[index][2][1], clusters[index][2][2], clusters[index][2][3] )
#             end
#             if length(clusters[index][2])==2
#                 B = Ccluster.box( clusters[index][2][1], fmpq(0,1), clusters[index][2][2] )
#             end
#         end
#         
#         bRe::arb = ball(R(Ccluster.getCenterRe(B)), R(fmpq(1,2)*Ccluster.getWidth(B)))
#         bIm::arb = ball(R(Ccluster.getCenterIm(B)), R(fmpq(1,2)*Ccluster.getWidth(B)))
#         if isa( clusters[index][2], Array{fmpq, 1})
#             if length(clusters[index][2])==2
#                 bIm = ball(R(0), R(0))
#             end
#         end
#         b::acb   = C(bRe, bIm)
#         s = @sprintf("*** cluster with sum of multiplicity %4d *** \n", mult); write(out,s);
#         write(out, "$(b)\n");
#     end
#     write(out, "-------------------------------------\n")
#     return
# end

R, x = PolynomialRing(QQ, "x")

d = 64 #degree
bitsize = 14
bitsize = Int64(floor(bitsize/2))
P=x^d - 2*((2^bitsize)*x-1)^2

bInit = [fmpq(-1,1),fmpq(0,1),fmpq(100,1)] #box centered in 0 + sqrt(-1)*0 with width 4
precision = 53                          #get clusters of size 2^-53
    
# Res = ccluster(P, bInit, precision, verbosity="silent");
Res = ccluster(P, bInit, precision, verbosity="results");

# using CclusterPlot #only if you have installed CclusterPlot.jl

# plotCcluster(Res, bInit, focus=true) #use true instead of false to focus on clusters

# printClusters(stdout, Res, 53)

# Res = risolate(P, precision, verbosity="silent");
Res = risolate(P, 2, verbosity="results");

# printClusters(stdout, Res, 53)
