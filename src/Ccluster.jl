VERSION >= v"0.4.0-dev+6521" && __precompile__()

module Ccluster

import Nemo: fmpq, acb_poly, fmpq_poly, QQ

using PyCall
using PyPlot

###############################################################################
#
#   Set up environment / load libraries
#
###############################################################################

const pkgdir = realpath(joinpath(dirname(@__FILE__), ".."))
const libdir = joinpath(pkgdir, "local", "lib")

const libccluster = joinpath(pkgdir, "local", "lib", "libccluster")

const matplotlib_patches = PyNULL()

function __init__()
    push!(Libdl.DL_LOAD_PATH, libdir)
    if is_linux() 
            Libdl.dlopen(libccluster)
    end
    println("")
    println("Welcome to Ccluster version 0.0.1")
    println("")
#     println("Nemo comes with absolutely no warranty whatsoever")
#     println("")

    copy!(matplotlib_patches, pyimport_conda("matplotlib.patches", "matplotlib"))
end

# using Nemo
# using PyCall
# using PyPlot

include("box.jl")
include("listBox.jl")
include("connComp.jl")
include("listConnComp.jl")
include("plotCcluster.jl")

# global PGLOBALCCLUSTERFMPQ = fmpq_poly(0);

__init__()
   
PGLOBALCCLUSTERFMPQ = fmpq_poly(0);

function getApp_FMPQ( dest::Ptr{acb_poly}, prec::Int )
    ccall((:acb_poly_set_fmpq_poly, :libarb), Void,
                (Ptr{acb_poly}, Ptr{fmpq_poly}, Int), dest, &PGLOBALCCLUSTERFMPQ, prec)
end

function ccluster( P_FMPQ::fmpq_poly, initialBox::Array{fmpq,1}, eps::fmpq, verbose::Int)
    
    return ccluster( getApp_FMPQ, initialBox, eps, 23, verbose)
    
end

function ccluster( P_FMPQ::fmpq_poly, initialBox::Array{fmpq,1}, eps::fmpq, strat::Int, verbose::Int)
    
    ccall((:fmpq_poly_set, :libflint), Void,
                (Ptr{fmpq_poly}, Ptr{fmpq_poly}), &PGLOBALCCLUSTERFMPQ, &P_FMPQ)

    return ccluster( getApp_FMPQ, initialBox, eps, strat, verbose)
    
end

function ccluster( getApprox::Function, initialBox::Array{fmpq,1}, eps::fmpq, verbose::Int)
    
    return ccluster( getApprox, initialBox, eps, 23, verbose)
end

function ccluster( getApprox::Function, initialBox::Array{fmpq,1}, eps::fmpq, strat::Int, verbose::Int)
    
    initBox::box = box(initialBox[1],initialBox[2],initialBox[3])

    queueResults = ccluster( getApprox, initBox, eps, strat, verbose )
    
    for sol in queueResults
        tempBO = sol[2]
        sol[2] = [getCenterRe(tempBO),getCenterIm(tempBO),fmpq(3,4)*getWidth(tempBO)]
    end
    
    return queueResults
end

function ccluster( getApprox::Function, initBox::box, eps::fmpq, strat::Int, verbose::Int)
    
    const getApp_c = cfunction(getApprox, Void, (Ptr{acb_poly}, Int))
    
    lccRes = listConnComp()
    ccall( (:ccluster_interface_forJulia, :libccluster), 
             Void, (Ptr{listConnComp}, Ptr{Void},    Ptr{box}, Ptr{fmpq}, Int,   Int), 
                    &lccRes,           getApp_c,   &initBox, &eps,      strat, verbose )
     
    queueResults = []
    while !isEmpty(lccRes)
        tempCC = pop(lccRes)
        tempBO = getComponentBox(tempCC,initBox)
        push!(queueResults, [getNbSols(tempCC),tempBO])
    end
    
    return queueResults
    
end

function ccluster_refine(qRes::listConnComp, 
                         getApprox::Function, 
                         CC::listConnComp, 
                         initBox::box,
                         eps::fmpq, 
                         strat::Int, 
                         verbose::Int = 0 )
    
#     initBox::box = box(initialBox[1],initialBox[2],initialBox[3])
    const getApp_c = cfunction(getApprox, Void, (Ptr{acb_poly}, Int))
    
#     lccRes = listConnComp()
    ccall( (:ccluster_refine_forJulia, :libccluster), 
             Void, (Ptr{listConnComp}, Ptr{listConnComp}, Ptr{Void}, Ptr{box}, Ptr{fmpq}, Int,   Int), 
                    &qRes,             &CC,               getApp_c,  &initBox, &eps,      strat, verbose )
     
#     queueResults = []
#     while !isEmpty(lccRes)
#         tempCC = pop(lccRes)
#         tempBO = getComponentBox(tempCC,initBox)
#         push!(queueResults, [getNbSols(tempCC),tempBO])
#     end
#     
#     return queueResults
    
end

function ccluster_DAC_first(qRes::listConnComp, 
                            qMainLoop::listConnComp, 
                            discardedCcs::listConnComp,
                            getApprox::Function, 
                            nbSols::Int,
                            initBox::box, eps::fmpq, strat::Int, verbose::Int = 0 )
    
    const getApp_c = cfunction(getApprox, Void, (Ptr{acb_poly}, Int))
    
    ccall( (:ccluster_DAC_first_interface_forJulia, :libccluster), 
             Void, (Ptr{listConnComp}, Ptr{listConnComp}, Ptr{listConnComp}, 
                    Ptr{Void},Int, Ptr{box}, Ptr{fmpq}, Int,   Int), 
                    &qRes,             &qMainLoop,        &discardedCcs,
                    getApp_c,  nbSols, &initBox, &eps,      strat, verbose )
                    
    return
    
end

function ccluster_DAC_next(qRes::listConnComp, 
                            qMainLoop::listConnComp, 
                            discardedCcs::listConnComp,
                            getApprox::Function,
                            nbSols::Int,
                            initBox::box, eps::fmpq, strat::Int, verbose::Int = 0 )
    
    const getApp_c = cfunction(getApprox, Void, (Ptr{acb_poly}, Int))
    
    ccall( (:ccluster_DAC_next_interface_forJulia, :libccluster), 
             Void, (Ptr{listConnComp}, Ptr{listConnComp}, Ptr{listConnComp}, 
                    Ptr{Void}, Int, Ptr{box}, Ptr{fmpq}, Int,   Int), 
                    &qRes,             &qMainLoop,        &discardedCcs,
                    getApp_c,  nbSols, &initBox, &eps,      strat, verbose )
    
    return 
    
end

function ccluster_draw( getApprox::Function, initialBox::Array{fmpq,1}, eps::fmpq, strat::Int, verbose::Int = 0 )
    
    initBox::box = box(initialBox[1],initialBox[2],initialBox[3])
    const getApp_c = cfunction(getApprox, Void, (Ptr{acb_poly}, Int))
    
    lccRes = listConnComp()
    lcbDis = listBox()
    
    ccall( (:ccluster_interface_forJulia_draw, :libccluster), 
             Void, (Ptr{listConnComp},Ptr{listBox}, Ptr{Void},    Ptr{box}, Ptr{fmpq}, Int,   Int), 
                    &lccRes, &lcbDis,          getApp_c,   &initBox, &eps,      strat, verbose )
     
    queueResults = []
    while !isEmpty(lccRes)
        tempCC = pop(lccRes)
#         tempBO = getComponentBox(tempCC,initBox)
#         push!(queueResults, [getNbSols(tempCC),[getCenterRe(tempBO),getCenterIm(tempBO),fmpq(3,4)*getWidth(tempBO)]])
        push!(queueResults, tempCC)
    end
    
    queueDiscarded = []
    while !isEmpty(lcbDis)
        tempB = pop(lcbDis)
        push!(queueDiscarded, tempB)
#         push!(queueDiscarded, [getCenterRe(tempB), getCenterIm(tempB), getWidth(tempB)])
#         push!(queueDiscarded, [getCenterIm(tempB)])
    end
    
    return queueResults, queueDiscarded
    
end

export ccluster, plotCcluster

# POLY_GLOBAL = fmpq_poly()
# 
# function GETAPP_GLOBAL( dest::Ptr{acb_poly}, prec::Int )
#     ccall((:acb_poly_set_fmpq_poly, :libarb), Void,
#                 (Ptr{acb_poly}, Ptr{fmpq_poly}, Int), dest, &POLY_GLOBAL, prec)
# end
# 
# function Ccluster( poly::fmpq_poly, initialBox::Array{fmpq,1}, eps::fmpq, strat::Int, verbose::Int = 0 )
#     
#     POLY_GLOBAL = poly
#     
#     initBox::box = box(initialBox[1],initialBox[2],initialBox[3])
#     const getApp_c = cfunction(GETAPP_GLOBAL, Void, (Ptr{acb_poly}, Int))
#     
#     lccRes = listConnComp()
#     
#     ccall( (:ccluster_interface_forJulia, :libccluster), 
#              Void, (Ptr{listConnComp}, Ptr{Void},    Ptr{box}, Ptr{fmpq}, Int,   Int), 
#                     &lccRes,           getApp_c,   &initBox, &eps,      strat, verbose )
#      
#      
#     queueResults::Array{Array{fmpq,1},1} = []
#     while !isEmpty(lccRes)
#         temp = getComponentBox(pop(lccRes),initBox)
#         push!(queueResults, [getCenterRe(temp),getCenterIm(temp),fmpq(3,4)*getWidth(temp)])
#     end
#     
#     return queueResults
#     
# end  

end # module
