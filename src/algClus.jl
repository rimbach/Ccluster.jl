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
    _approx::acb                            # satisfies |\alpha - a |<2^{-prec}
    
    _CC::Ptr{Ccluster.connComp}             # a pointer on the connected component of boxes
    _initBox::Ccluster.box                  # the initial box i.e. root in the subdivision tree
    
    function algClus( objCC::Ccluster.connComp, ptrCC::Ptr{Ccluster.connComp}, bInit::Ccluster.box, prec )
        z = new()
        
        z._nbSols = Ccluster.getNbSols(objCC)
        z._prec = prec
        z._initBox = Ccluster.box( Ccluster.getCenterRe(bInit), Ccluster.getCenterIm(bInit), Ccluster.getWidth(bInit) )
        z._CC = ptrCC
        
        z._isolatingBox = Ccluster.getComponentBox(objCC,z._initBox)
        R = RealField(z._prec)
        C = ComplexField(z._prec)
        bRe = ball(R(Ccluster.getCenterRe(z._isolatingBox)), R(fmpq(1,2)*Ccluster.getWidth(z._isolatingBox)))
        bIm = ball(R(Ccluster.getCenterIm(z._isolatingBox)), R(fmpq(1,2)*Ccluster.getWidth(z._isolatingBox)))
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
        bRe = ball(R(Ccluster.getCenterRe(z._isolatingBox)), R(fmpq(1,2)*Ccluster.getWidth(z._isolatingBox)))
        bIm = ball(R(Ccluster.getCenterIm(z._isolatingBox)), R(fmpq(1,2)*Ccluster.getWidth(z._isolatingBox)))
        z._approx = C(bRe, bIm)
        return z
    end
end
