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

# saveBoxForDebug = []

type connComp
#   list of boxes 
    _boxes_begin::Ptr{Void}
    _boxes_end::Ptr{Void}
    _boxes_size::Cint
    _boxes_clear::Ptr{Void}
#     _list_of_boxes::listBox
#   width: fmpq 
    _width_den::Int
    _width_num::Int
#   infRe: fmpq
    _infRe_den::Int
    _infRe_num::Int
#   supRe: fmpq
    _supRe_den::Int
    _supRe_num::Int
#   infIm: fmpq
    _infIm_den::Int
    _infIm_num::Int
#   supIm: fmpq
    _supIm_den::Int
    _supIm_num::Int
#   nSols: int
    _nSols::Cint
#   nwSpd: fmpz
    _nwSpd::Int
#   appPr: int
    _appPr::Int
#   newSu: int
    _newSu::Cint
    
    function connComp()
        z = new()
        ccall( (:connCmp_init, :libccluster), 
             Void, (Ptr{connComp},), 
                    &z)
        finalizer(z, _connComp_clear_fn)
        return z
    end
    
#     function connComp(b::box)
#         push!(saveBoxForDebug, b)
#         z = new()
#         ccall( (:connCmp_init_compBox, :libccluster), 
#              Void, (Ptr{connComp}, Ptr{box}), 
#                     &z,            &b)
#         finalizer(z, _connComp_clear_fn)
#         return z
#     end
    
end

function _connComp_clear_fn(cc::connComp)
    ccall( (:connCmp_clear, :libccluster), 
         Void, (Ptr{connComp},), 
                &cc)
end

function getNbSols(cc::connComp)
    return Int(cc._nSols)
end

function getComponentBox(cc::connComp, initialBox::box)
    
    res = box()
    ccall( (:connCmp_componentBox, :libccluster), 
             Void, (Ptr{box}, Ptr{connComp}, Ptr{box}), 
                    &res,      &cc,           &initialBox);
    return res
    
end

function pop( cc::connComp )
    res = ccall( (:connCmp_pop, :libccluster), 
                  Ptr{box}, (Ptr{connComp},), 
                                 &cc)                        
    resobj::box = unsafe_load(res)
    return resobj
end

function isEmpty(cc::connComp)
    res = ccall( (:connCmp_is_empty, :libccluster), 
                  Cint, (Ptr{connComp},), 
                        &cc )
    return Bool(res)
end