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

mutable struct algClus                      # represents an algebraic cluster \alpha
    _nbSols::Int                            # the sum of multiplicity of roots in the cluster
    _prec::Int                              # the precision of the cluster; i.e. w() < 2*2^{-prec}
    _isolatingBox::Ccluster.box             # a box isolating \alpha st width(_isolatingBox) < 2*2^(-prec)
    _approx::Nemo.acb                       # satisfies |\alpha - a |<2^{-prec}
    
    _CC::Ptr{Ccluster.connComp}             # a pointer on the connected component of boxes
    _initBox::Ccluster.box                  # the initial box i.e. root in the subdivision tree
    
    function algClus( objCC::Ccluster.connComp, ptrCC::Ptr{Ccluster.connComp}, bInit::Ccluster.box, prec::Int )
        z = new()
        
        z._nbSols = Ccluster.getNbSols(objCC)
        z._prec = prec
        z._initBox = Ccluster.box( Ccluster.getCenterRe(bInit), Ccluster.getCenterIm(bInit), Ccluster.getWidth(bInit) )
        z._CC = ptrCC
        
        z._isolatingBox = Ccluster.getComponentBox(objCC,z._initBox)
        R = RealField(z._prec)
        C = ComplexField(z._prec)
        bRe::Nemo.arb = Nemo.ball(R(Ccluster.getCenterRe(z._isolatingBox)), R(fmpq(1,2)*Ccluster.getWidth(z._isolatingBox)))
        bIm::Nemo.arb = Nemo.ball(R(Ccluster.getCenterIm(z._isolatingBox)), R(fmpq(1,2)*Ccluster.getWidth(z._isolatingBox)))
        z._approx = C(bRe, bIm)
        return z
    end
    
    function algClus( a::algClus )
        z = new()
        
        z._nbSols = a._nbSols
        z._prec = a._prec
        z._initBox = Ccluster.box( Ccluster.getCenterRe(a._initBox), Ccluster.getCenterIm(a._initBox), Ccluster.getWidth(a._initBox) )
        z._CC = Ccluster.copy_Ptr(a._CC)
        
        z._isolatingBox = Ccluster.box( Ccluster.getCenterRe(a._isolatingBox), Ccluster.getCenterIm(a._isolatingBox), Ccluster.getWidth(a._isolatingBox) )
        R = RealField(z._prec)
        C = ComplexField(z._prec)
        bRe::Nemo.arb = Nemo.ball(R(Ccluster.getCenterRe(z._isolatingBox)), R(fmpq(1,2)*Ccluster.getWidth(z._isolatingBox)))
        bIm::Nemo.arb = Nemo.ball(R(Ccluster.getCenterIm(z._isolatingBox)), R(fmpq(1,2)*Ccluster.getWidth(z._isolatingBox)))
        z._approx = C(bRe, bIm)
        return z
    end
end

function toStr(a::algClus)
    res = ""
    res = res * "algebraic cluster: prec: $(a._prec), nbsols: $(a._nbSols)\n"
    res = res * Ccluster.toStr(a._isolatingBox)
    res = res * "\n approx: $(a._approx)"
    res = res * "\n"
    res 
    return res
end

#copying
function copyIn( dest::algClus, src::algClus )
    dest._nbSols = src._nbSols
    dest._prec = src._prec
    dest._isolatingBox = Ccluster.box( Ccluster.getCenterRe(src._isolatingBox), 
                                       Ccluster.getCenterIm(src._isolatingBox), 
                                       Ccluster.getWidth(src._isolatingBox) )
    dest._initBox = Ccluster.box( Ccluster.getCenterRe(src._initBox), Ccluster.getCenterIm(src._initBox), Ccluster.getWidth(src._initBox) )
    dest._CC = Ccluster.copy_Ptr(src._CC)
    R = RealField(src._prec)
    C = ComplexField(src._prec)
    bRe::Nemo.arb = Nemo.ball(R(Ccluster.getCenterRe(dest._isolatingBox)), R(fmpq(1,2)*Ccluster.getWidth(dest._isolatingBox)))
    bIm::Nemo.arb = Nemo.ball(R(Ccluster.getCenterIm(dest._isolatingBox)), R(fmpq(1,2)*Ccluster.getWidth(dest._isolatingBox)))
    dest._approx = C(bRe, bIm);
end

function copyIn( dest::Array{algClus,1}, src::Array{algClus,1} )
    for index in 1:length(src)
        copyIn( dest[index], src[index] )
    end
end

function clusCopy(a::algClus)
    return algClus(a)
end

function clusCopy(a::Array{algClus,1})
    res=[]
    for index in 1:length(a)
        push!(res, algClus(a[index]))
    end
    return res
end
