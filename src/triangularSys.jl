#
#  Copyright (C) 2019 Remi Imbach
#
#  This file is part of Ccluster.
#
#  Ccluster is free software: you can redistribute it and/or modify it under
#  the terms of the GNU Lesser General Public License (LGPL) as published
#  by the Free Software Foundation; either version 2.1 of the License, or
#  (at your option) any later version.  See <http://www.gnu.org/licenses/>.
#

# utility functions for triangular systems
#

function getDeg(pol, ind::Int)::Array{Int,1} #get the vector of degrees of pol
    if ind==1
        return [Nemo.degree(pol)]
    else
        d::Int = Nemo.degree(pol)
        degrees = Int[]
        if d==-1
            for i=1:ind-1
                push!(degrees,0)
            end
        else
            global degrees
            degrees = getDeg( pol.coeffs[1], ind-1)
            for i = 2:d+1
                degtemp::Array{Int,1} = getDeg( pol.coeffs[i], ind-1)
                for j=1:ind-1
                    if degtemp[j]>degrees[j]
                        degrees[j] = degtemp[j]
                    end
                end
            end
        end
        push!(degrees, d)
        return degrees
    end
end

function evalUniFMPQPol(P::Nemo.fmpq_poly, b::Nemo.acb, prec::Int)::Nemo.acb
    CC::Nemo.AcbField = Nemo.ComplexField(prec)
    R::Nemo.AcbPolyRing, dummy = Nemo.PolynomialRing(CC, "dummy")
    res::Nemo.acb_poly = R(0)
    ccall((:acb_poly_set_fmpq_poly, :libarb), 
      Cvoid, (Ref{acb_poly}, Ref{fmpq_poly}, Int), 
              res,           P,             prec)
    return evaluate(res, CC(b))
end
    
function evalPolAt(P, b::Array{Nemo.acb, 1}, prec::Int)::Nemo.acb
    if (length(b)==1)
        res = evalUniFMPQPol(P, b[1], prec)
    else
        btemp::Nemo.acb=pop!(b)
        CC::Nemo.AcbField = Nemo.ComplexField(prec)
        R::Nemo.AcbPolyRing, dummy = Nemo.PolynomialRing(CC, "dummy")
        pol::Nemo.acb_poly = R(0);
        for index=0:degree(P)
            pol = pol + evalPolAt(coeff(P,index), b, prec)*dummy^index
        end
        res::Nemo.acb = evaluate(pol, CC(btemp))
        push!(b,btemp)
    end
    return res
end

function getPolAt(P, b::Array{Nemo.acb, 1}, prec::Int)::Nemo.acb_poly
#     print("--------------------------\n")
    CC::Nemo.AcbField = Nemo.ComplexField(prec)
    R::Nemo.AcbPolyRing, dummy = Nemo.PolynomialRing(CC, "dummy")
    res::Nemo.acb_poly = R(0)
    for index=0:Nemo.degree(P)
        res = res + evalPolAt(Nemo.coeff(P,index), b, prec)*dummy^index
    end
#     print("--------------------------\n")
#     print("res $(length(b)): $res\n")
#     print("--------------------------\n")
    return res
end
