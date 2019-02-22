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

type disk
    _center_real_den::Int
    _center_real_num::Int
    _center_imag_den::Int
    _center_imag_num::Int
    _radius_den::Int
    _radius_num::Int
    
    function disk()
        z = new()
        ccall( (:compDsk_init, :libccluster), 
             Void, (Ptr{disk},), 
                    &z)
        finalizer(z, _disk_clear_fn)
        return z
    end
    
    function disk(re::fmpq, im::fmpq, rad::fmpq)
        z = new()
        
        ccall( (:compDsk_init, :libccluster), 
             Void, (Ptr{disk},), 
                    &z)
        ccall( (:compDsk_set_3realRat, :libccluster), 
             Void, (Ptr{disk}, Ptr{fmpq}, Ptr{fmpq}, Ptr{fmpq}), 
                    &z,        &re,       &im,       &rad)
        finalizer(z, _disk_clear_fn)
        return z
    end
end

function _disk_clear_fn(d::disk)
    ccall( (:compDsk_clear, :libccluster), 
         Void, (Ptr{disk},), 
                &d)
end

function getCenterRe(d::disk)
    res = fmpq(0,1)
    ccall( (:compDsk_get_centerRe, :libccluster), 
             Void, (Ptr{fmpq}, Ptr{disk}), 
                    &res,      &d)
    return res
end

function getCenterIm(d::disk)
    res = fmpq(0,1)
    ccall( (:compDsk_get_centerIm, :libccluster), 
             Void, (Ptr{fmpq}, Ptr{disk}), 
                    &res,      &d)
    return res
end

function getRadius(d::disk)
    res = fmpq(0,1)
    ccall( (:compDsk_get_radius, :libccluster), 
             Void, (Ptr{fmpq}, Ptr{disk}), 
                    &res,      &d)
    return res
end

function inflateDisk(d::disk, ratio::fmpq)
    res = disk()
    ccall( (:compDsk_inflate_realRat, :libccluster), 
             Void, (Ptr{disk}, Ptr{disk}, Ptr{fmpq}), 
                    &res,      &d,        &ratio)
    return res
end

function isSeparated(d::disk, qMainLoop::listConnComp, qResults::listConnComp, qAllResults::listConnComp, discardedCcs::listConnComp )
    res = ccall( (:ccluster_compDsk_is_separated_DAC, :libccluster), 
                   Cint, (Ptr{disk}, Ptr{listConnComp}, Ptr{listConnComp}, Ptr{listConnComp}, Ptr{listConnComp}), 
                          &d,        &qMainLoop,       &qResults,           &qAllResults,       &discardedCcs)
    return Bool(res)
end
    
function toStr(d::disk)
    res = ""
    res = res * "Disk: center: $(getCenterRe(d))" 
    res = res * " + i*$(getCenterIm(d))"
    res = res * ", radius: $(getRadius(d))"
    return res
end

# function isDiskInDisk(d1::disk, d2::disk)
#     RR = RealField(64)
#     #compute the square of the distance between the center of d1 and the center of d2
#     dist::fmpq = (getCenterRe(d1)-getCenterRe(d2))^2 + (getCenterIm(d1)-getCenterIm(d2))^2
#     if dist>= (getRadius(d2))^2
#         return false
#     end
#     sqrtdist::arb = sqrt(RR(dist))+RR(getRadius(d1))
#     if sqrtdist <= RR(getRadius(d2))
#         return true
#     else 
#         return false
#     end
#     
#     
# end
    
