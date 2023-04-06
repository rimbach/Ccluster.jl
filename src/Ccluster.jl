__precompile__()

module Ccluster

using Libdl

import Nemo: fmpq, fmpz, acb_poly, fmpq_poly, QQ, prec, parent,
             degree, coeff, accuracy_bits

using Ccluster_jll

###############################################################################
#
#   Set up environment / load libraries
#
###############################################################################

#DEPRECATED
#const pkgdir = realpath(joinpath(dirname(@__FILE__), ".."))
#const libdir = joinpath(pkgdir, "local", "lib")

#const libccluster = joinpath(pkgdir, "local", "lib", "libccluster")
#libccluster is now defined by Ccluster_jll

function __init__()

#    if "HOSTNAME" in keys(ENV) && ENV["HOSTNAME"] == "juliabox"
#       push!(Libdl.DL_LOAD_PATH, "/usr/local/lib")
#    elseif Sys.islinux()
#        push!(Libdl.DL_LOAD_PATH, libdir)
#    else
#        push!(Libdl.DL_LOAD_PATH, libdir)
#   end
#     println("")
#     println("Welcome to Ccluster version 0.0.1")
#     println("")
#       print("libccluster: $libccluster \n")
#       Libdl.dlopen(libccluster)
#     push!(Libdl.DL_LOAD_PATH, "/home/imbach/.julia/artifacts/7524922475b968e0b4750f25dfc030c1fc43d65e/lib/")
end

include("box.jl")
include("listBox.jl")
include("connComp.jl")
include("listConnComp.jl")
include("disk.jl")
include("algClus.jl")
include("triangularSys.jl")
include("Tcluster.jl")
include("rand_generator.jl")

__init__()
   
function ptr_set_acb_poly( dest::Ref{acb_poly}, src::acb_poly )
#     ccall((:acb_poly_set, libarb), Nothing,
    ccall((:compApp_poly_set, libccluster), Nothing,
                (Ref{acb_poly}, Ref{acb_poly}, Int), 
                 dest,         src,          precision(parent(src)))
end

function ptr_set_2fmpq_poly( dest::Ref{acb_poly}, re::fmpq_poly, im::fmpq_poly, prec::Int )
#     ccall((:acb_poly_set2_fmpq_poly, libarb), Nothing,
    ccall((:compApp_poly_set2_fmpq_poly, libccluster), Nothing,
                (Ref{acb_poly}, Ref{fmpq_poly}, Ref{fmpq_poly}, Int), 
                 dest,         re,             im,        prec)
end

# check the RELATIVE ACCURACY in bits
function checkAccuracy( pol::acb_poly, prec::Int )
    
    for i=0:degree(pol)
        if accuracy_bits( coeff(pol, i) ) < prec
            return 0
        end
    end
    return 1
    
end

# returns the MINIMUM RELATIVE ACCURACY in bits
function getAccuracy( pol::acb_poly )::Int
    res::Int = accuracy_bits( coeff(pol, 0) )
    for i=1:degree(pol)
        if accuracy_bits( coeff(pol, i) ) < res
            res = accuracy_bits( coeff(pol, i) )
        end
    end
    return res
    
end

function parseVerbosity( verbosity::String )::Int 
    if verbosity=="silent"
        return 0
    elseif verbosity=="brief"
        return 1
    elseif verbosity=="stats"
        return 2
    elseif verbosity=="results"
        return 1
    else
        return 0
    end
end

#print clusters
function printClusters(out, clusters, prec::Int)
    R::ArbField = RealField(prec)
    C::AcbField = ComplexField(prec)
    write(out, "-------------------------------------\n")
    for index in 1:length(clusters)
        mult = clusters[index][1]
        B    = clusters[index][2]
        if isa( clusters[index][2], Array{fmpq, 1})
            if length(clusters[index][2])==3
                B = Ccluster.disk( clusters[index][2][1], clusters[index][2][2], clusters[index][2][3] )
            end
            if length(clusters[index][2])==2
                B = Ccluster.disk( clusters[index][2][1], fmpq(0,1), clusters[index][2][2] )
            end
        end
        if isa( clusters[index][2], Ccluster.box)
            B = Ccluster.disk( Ccluster.getCenterRe(B), Ccluster.getCenterIm(B), fmpq(3,4)*Ccluster.getWidth(B) )
        end
        
        bRe::arb = ball(R(Ccluster.getCenterRe(B)), R(Ccluster.getRadius(B)))
        bIm::arb = ball(R(Ccluster.getCenterIm(B)), R(Ccluster.getRadius(B)))
        if isa( clusters[index][2], Array{fmpq, 1})
            if length(clusters[index][2])==2
                bIm = ball(R(0), R(0))
            end
        end
        b::acb   = C(bRe, bIm)
        s = @sprintf("*** cluster with sum of multiplicity %4d *** \n", mult); write(out,s);
        write(out, "$(b)\n");
    end
    write(out, "-------------------------------------\n")
    return
end

function ccluster( getApprox::Function, 
                   initialBox::Array{fmpq,1}, 
                   precision::Int;
                   strategy="default", #a strategy: Int
                   verbosity="silent" )#a verbosity flag; by defaults, nothing is printed
                                       #options are "silent", "brief" and "results"
      
      initBox::box = box(initialBox[1],initialBox[2],initialBox[3]); 
      
      queueResults = ccluster( getApprox, initBox, precision, strategy=strategy, verbosity=verbosity )
    
      for sol in queueResults
          tempBO = sol[2]
          sol[2] = [getCenterRe(tempBO),getCenterIm(tempBO),fmpq(3,4)*getWidth(tempBO)]
      end
      
      return queueResults
      
end

function ccluster( getApprox::Function, 
                   initBox::box, 
                   precision::Int;
                   strategy="default", #a strategy: Int
                   verbosity="silent" )#a verbosity flag; by defaults, nothing is printed
                                       #options are "silent", "brief" and "results"
    getApp_c = @cfunction( $getApprox, Cvoid, (Ptr{acb_poly}, Int))
    lccRes = listConnComp()
    verbose::Int = parseVerbosity(verbosity)
    eps = fmpq(1, fmpz(2)^precision)
                     
    ccall( (:ccluster_forJulia_func, libccluster), 
             Nothing, (Ref{listConnComp}, Ptr{Cvoid}, Ref{box}, Ref{fmpq}, Cstring,  Int, Int), 
                       lccRes,            getApp_c,   initBox,  eps,       strategy, 1,   verbose )
                     
    queueResults = []
    
    while !isEmpty(lccRes)
        tempCC = pop(lccRes)
        tempBO = getComponentBox(tempCC,initBox)
        push!(queueResults, [getNbSols(tempCC),tempBO])
    end
    
    if verbosity=="results"
        printClusters(stdout, queueResults, precision)
    end
    
    return queueResults                                 
end

function ccluster( getApprox::Function, 
                   precision::Int;
                   strategy="default", #a strategy: Str
                   verbosity="silent" )#a verbosity flag; by defaults, nothing is printed
                                       #options are "silent", "brief" and "results"
    getApp_c = @cfunction( $getApprox, Cvoid, (Ptr{acb_poly}, Int))
    lccRes = listConnComp()
    verbose::Int = parseVerbosity(verbosity)
    eps = fmpq(1, fmpz(2)^precision)
    
    initBox::box = box(fmpq(0,1),fmpq(0,1),fmpq(0,1));
    
    ccall( (:ccluster_global_forJulia_func, libccluster), 
             Nothing, (Ref{listConnComp}, Ptr{Cvoid}, Ref{box}, Ref{fmpq}, Cstring,  Int, Int), 
                       lccRes,            getApp_c,   initBox,  eps,       strategy, 1,   verbose )
                     
    queueResults = []
    
    while !isEmpty(lccRes)
        tempCC = pop(lccRes)
        tempBO = getComponentBox(tempCC,initBox)
        tempBO = [getCenterRe(tempBO),getCenterIm(tempBO),fmpq(3,4)*getWidth(tempBO)]
        push!(queueResults, [getNbSols(tempCC),tempBO])
    end
    
    if verbosity=="results"
        printClusters(stdout, queueResults, precision)
    end
    
    return queueResults                                 
end

function ccluster( P_FMPQ::fmpq_poly, 
                   initialBox::Array{fmpq,1}, 
                   precision::Int;
                   strategy="default", #a strategy: Int
                   verbosity="silent" )#a verbosity flag; by defaults, nothing is printed
                                       #options are "silent", "brief" and "results"
                                    
      initBox::box = box(initialBox[1],initialBox[2],initialBox[3]); 
      
      queueResults = ccluster( P_FMPQ, initBox, precision, strategy=strategy, verbosity=verbosity )
    
      for sol in queueResults
          tempBO = sol[2]
          sol[2] = [getCenterRe(tempBO),getCenterIm(tempBO),fmpq(3,4)*getWidth(tempBO)]
      end
      
      return queueResults
end

function ccluster( P_FMPQ::fmpq_poly, 
                   initBox::box, 
                   precision::Int;
                   strategy="default", #a strategy: Int
                   verbosity="silent" )#a verbosity flag; by defaults, nothing is printed
                                       #options are "silent", "brief" and "results"
     
    lccRes = listConnComp()
    verbose::Int = parseVerbosity(verbosity)
    eps = fmpq(1, fmpz(2)^precision)
                       
    ccall( (:ccluster_forJulia_realRat_poly, libccluster), 
             Nothing, (Ref{listConnComp}, Ref{fmpq_poly}, Ref{box}, Ref{fmpq}, Cstring, Int, Int), 
                       lccRes,            P_FMPQ,         initBox,  eps,       strategy, 1,  verbose )
                     
    queueResults = []
    
    while !isEmpty(lccRes)
        tempCC = pop(lccRes)
        tempBO = getComponentBox(tempCC,initBox)
        push!(queueResults, [getNbSols(tempCC),tempBO])
    end
    
    if verbosity=="results"
        printClusters(stdout, queueResults, precision)
    end
    
    return queueResults                                   
end

function ccluster( P_FMPQ::fmpq_poly, 
                   precision::Int;
                   strategy="default", #a strategy: Int
                   verbosity="silent" )#a verbosity flag; by defaults, nothing is printed
                                       #options are "silent", "brief" and "results"
     
    lccRes = listConnComp()
    verbose::Int = parseVerbosity(verbosity)
    eps = fmpq(1, fmpz(2)^precision)
    
    initBox::box = box(fmpq(0,1),fmpq(0,1),fmpq(0,1));
    
    ccall( (:ccluster_global_forJulia_realRat_poly, libccluster), 
             Nothing, (Ref{listConnComp}, Ref{fmpq_poly}, Ref{box}, Ref{fmpq}, Cstring, Int, Int), 
                       lccRes,            P_FMPQ,         initBox,  eps,       strategy, 1,  verbose )
                     
    queueResults = []
    
    while !isEmpty(lccRes)
        tempCC = pop(lccRes)
        tempBO = getComponentBox(tempCC,initBox)
        tempBO = [getCenterRe(tempBO),getCenterIm(tempBO),fmpq(3,4)*getWidth(tempBO)]
        push!(queueResults, [getNbSols(tempCC),tempBO])
    end
    
    if verbosity=="results"
        printClusters(stdout, queueResults, precision)
    end
    
    return queueResults                                   
end

function risolate( P_FMPQ::fmpq_poly, 
                   initBox::Array{fmpq,1}, 
                   precision::Int;
                   strategy="default", #a strategy: Int
                   verbosity="silent" )#a verbosity flag; by defaults, nothing is printed
                                       #options are "silent", "brief" and "results"
    
     initialBox::box = box(initBox[1],fmpq(0,1),initBox[2]); 
    
     queueResults = risolate( P_FMPQ, initialBox, precision, strategy=strategy, verbosity=verbosity )
    
      for sol in queueResults
          tempBO = sol[2]
          sol[2] = [getCenterRe(tempBO),getWidth(tempBO)]
      end
      
      return queueResults                                   
end

function risolate( P_FMPQ::fmpq_poly, 
                   initBox::box,
                   precision::Int;
                   strategy="default", #a strategy: Int
                   verbosity="silent" )#a verbosity flag; by defaults, nothing is printed
                                       #options are "silent", "brief" and "results"
     
    lccRes = listConnComp()
    verbose::Int = parseVerbosity(verbosity)
    eps = fmpq(1, fmpz(2)^precision)
    
#     ccall( (:risolate_global_forJulia_realRat_poly, libccluster), 
#              Nothing, (Ref{listConnComp}, Ref{fmpq_poly}, Ref{box}, Ref{fmpq}, Cstring, Int, Int), 
#                        lccRes,            P_FMPQ,         initBox,  eps,       strategy, 1,  verbose )
                    
    ccall( (:risolate_forJulia_realRat_poly, libccluster), 
             Nothing, (Ref{listConnComp}, Ref{fmpq_poly}, Ref{box}, Ref{fmpq}, Cstring, Int, Int), 
                       lccRes,            P_FMPQ,         initBox,  eps,       strategy, 1,  verbose )
                       
    queueResults = []
    
    while !isEmpty(lccRes)
        tempCC = pop(lccRes)
        tempBO = getComponentBox(tempCC,initBox)
        push!(queueResults, [getNbSols(tempCC),tempBO])
    end
    
    if verbosity=="results"
        printClusters(stdout, queueResults, precision)
    end
    
    return queueResults                                   
end

function risolate( P_FMPQ::fmpq_poly, 
                   precision::Int;
                   strategy="default", #a strategy: Int
                   verbosity="silent" )#a verbosity flag; by defaults, nothing is printed
                                       #options are "silent", "brief" and "results"
     
    lccRes = listConnComp()
    verbose::Int = parseVerbosity(verbosity)
    eps = fmpq(1, fmpz(2)^precision)
    
    initBox::box = box(fmpq(0,1),fmpq(0,1),fmpq(0,1));
    
    ccall( (:risolate_global_forJulia_realRat_poly, libccluster), 
             Nothing, (Ref{listConnComp}, Ref{fmpq_poly}, Ref{box}, Ref{fmpq}, Cstring, Int, Int), 
                       lccRes,            P_FMPQ,         initBox,  eps,       strategy, 1,  verbose )
                     
    queueResults = []
    
    while !isEmpty(lccRes)
        tempCC = pop(lccRes)
        tempBO = getComponentBox(tempCC,initBox)
        tempBO = [getCenterRe(tempBO),getWidth(tempBO)]
        push!(queueResults, [getNbSols(tempCC),tempBO])
    end
    
    if verbosity=="results"
        printClusters(stdout, queueResults, precision)
    end
    
    return queueResults                                   
end

function ccluster( Preal_FMPQ::fmpq_poly, Pimag_FMPQ::fmpq_poly, 
                   initialBox::Array{fmpq,1}, 
                   precision::Int;
                   strategy="default", #a strategy: Int
                   verbosity="silent" )#a verbosity flag; by defaults, nothing is printed
                                       #options are "silent", "brief" and "results"
      initBox::box = box(initialBox[1],initialBox[2],initialBox[3]); 
      
      queueResults = ccluster( Preal_FMPQ, Pimag_FMPQ, initBox, precision, strategy=strategy, verbosity=verbosity )
    
      for sol in queueResults
          tempBO = sol[2]
          sol[2] = [getCenterRe(tempBO),getCenterIm(tempBO),fmpq(3,4)*getWidth(tempBO)]
      end
      
      return queueResults                                      
end

function ccluster( Preal_FMPQ::fmpq_poly, Pimag_FMPQ::fmpq_poly, 
                   initBox::box, 
                   precision::Int;
                   strategy="default", #a strategy: Int
                   verbosity="silent" )#a verbosity flag; by defaults, nothing is printed
                                       #options are "silent", "brief" and "results"
    lccRes = listConnComp()
    verbose::Int = parseVerbosity(verbosity)
    eps = fmpq(1, fmpz(2)^precision)
    
#     ccall( (:ccluster_interface_forJulia_realRat_poly_real_imag, libccluster), 
#              Nothing, (Ref{listConnComp}, Ref{fmpq_poly}, Ref{fmpq_poly}, Ref{box}, Ref{fmpq}, Int,   Int), 
#                        lccRes,            Preal_FMPQ,     Pimag_FMPQ,     initBox,  eps,       55, verbose )
                       
    ccall( (:ccluster_forJulia_realRat_poly_real_imag, libccluster), 
             Nothing, (Ref{listConnComp}, Ref{fmpq_poly}, Ref{fmpq_poly}, Ref{box}, Ref{fmpq}, Cstring, Int, Int), 
                       lccRes,            Preal_FMPQ,     Pimag_FMPQ,     initBox,  eps,       strategy,   1,   verbose )
                     
    queueResults = []
    
    while !isEmpty(lccRes)
        tempCC = pop(lccRes)
        tempBO = getComponentBox(tempCC,initBox)
        push!(queueResults, [getNbSols(tempCC),tempBO])
    end
    
    if verbosity=="results"
        printClusters(stdout, queueResults, precision)
    end
    
    return queueResults                                   
end

function ccluster( Preal_FMPQ::fmpq_poly, Pimag_FMPQ::fmpq_poly, 
                   precision::Int;
                   strategy="default", #a strategy: Int
                   verbosity="silent" )#a verbosity flag; by defaults, nothing is printed
                                       #options are "silent", "brief" and "results"
    lccRes = listConnComp()
    verbose::Int = parseVerbosity(verbosity)
    eps = fmpq(1, fmpz(2)^precision)
    
    initBox::box = box(fmpq(0,1),fmpq(0,1),fmpq(0,1));
    
    ccall( (:ccluster_global_forJulia_realRat_poly_real_imag, libccluster), 
             Nothing, (Ref{listConnComp}, Ref{fmpq_poly}, Ref{fmpq_poly}, Ref{box}, Ref{fmpq}, Cstring, Int, Int), 
                       lccRes,            Preal_FMPQ,     Pimag_FMPQ,     initBox,  eps,       strategy,   1,   verbose )
                     
    queueResults = []
    
    while !isEmpty(lccRes)
        tempCC = pop(lccRes)
        tempBO = getComponentBox(tempCC,initBox)
        tempBO = [getCenterRe(tempBO),getCenterIm(tempBO),fmpq(3,4)*getWidth(tempBO)]
        push!(queueResults, [getNbSols(tempCC),tempBO])
    end
    
    if verbosity=="results"
        printClusters(stdout, queueResults, precision)
    end
    
    return queueResults                                   
end

function ccluster_solve(getApprox::Function, 
                        initBox::box, 
                        eps::fmpq, 
                        strategy::String, 
                        verbose::Int)::listConnComp

    getApp_c = @cfunction( $getApprox, Cvoid, (Ptr{acb_poly}, Int))
    
    lccRes = listConnComp()
                     
        ccall( (:ccluster_forJulia_func, libccluster), 
             Nothing, (Ref{listConnComp}, Ptr{Cvoid}, Ref{box}, Ref{fmpq}, Cstring,  Int, Int), 
                       lccRes,            getApp_c,   initBox,  eps,       strategy, 1,   verbose )
    
    return lccRes
    
end

function ccluster_solve(getApprox::Function, 
                        eps::fmpq, 
                        strategy::String, 
                        verbose::Int)::Tuple{listConnComp, box}

    getApp_c = @cfunction( $getApprox, Cvoid, (Ptr{acb_poly}, Int))
    
    initBox::box = box(fmpq(0,1),fmpq(0,1),fmpq(0,1));
    
    lccRes = listConnComp()
                     
        ccall( (:ccluster_global_forJulia_func, libccluster), 
             Nothing, (Ref{listConnComp}, Ptr{Cvoid}, Ref{box}, Ref{fmpq}, Cstring,  Int, Int), 
                       lccRes,            getApp_c,   initBox,  eps,       strategy, 1,   verbose )
    
    return lccRes, initBox
    
end

function ccluster_solve(P_FMPQ::fmpq_poly, 
                              eps::fmpq, 
                              strategy::String, 
                              verbose::Int)::Tuple{listConnComp, box}
    
    initBox::box = box(fmpq(0,1),fmpq(0,1),fmpq(0,1));
    
    lccRes = listConnComp()
                       
        ccall( (:ccluster_global_forJulia_realRat_poly, libccluster), 
             Nothing, (Ref{listConnComp}, Ref{fmpq_poly}, Ref{box}, Ref{fmpq}, Cstring, Int, Int), 
                       lccRes,            P_FMPQ,         initBox,  eps,       strategy, 1,  verbose )
                       
    return lccRes, initBox
    
end

function ccluster_solve_forTcluster(getApprox::Function, 
                        eps::fmpq, 
                        strategy::String, 
                        verbose::Int)::Tuple{listConnComp, box, Int}

    getApp_c = @cfunction( $getApprox, Cvoid, (Ptr{acb_poly}, Int))
    
    initBox::box = box(fmpq(0,1),fmpq(0,1),fmpq(0,1));
    
    lccRes = listConnComp()
                     
    resCall::Int =  ccall( (:ccluster_global_forJulia_forTcluster_func, libccluster), 
             Int, (Ref{listConnComp}, Ptr{Cvoid}, Ref{box}, Ref{fmpq}, Cstring,  Int, Int), 
                       lccRes,            getApp_c,   initBox,  eps,       strategy, 1,   verbose )
#     print("resCall: $resCall\n")
    return lccRes, initBox, resCall
    
end


function ccluster_refine(qRes::listConnComp, 
                         getApprox::Function, 
                         CC::listConnComp, 
                         initBox::box,
                         eps::fmpq, 
                         strategy::String, 
                         verbose::Int = 0 )
    
    getApp_c = @cfunction( $getApprox, Cvoid, (Ptr{acb_poly}, Int))
                     
    ccall( (:ccluster_forJulia_refine, libccluster), 
             Nothing, (Ref{listConnComp}, Ref{listConnComp}, Ptr{Cvoid}, Ref{box}, Ref{fmpq}, Cstring,  Int,  Int), 
                      qRes,               CC,                getApp_c,   initBox,  eps,       strategy, 1,    verbose )
                    
end

function ccluster_DAC_first(qRes::listConnComp,
                            qAllRes::listConnComp,
                            qMainLoop::listConnComp, 
                            discardedCcs::listConnComp,
                            getApprox::Function, 
                            nbSols::Int,
#                             initialBox::Array{fmpq,1},
                            initBox::box,
                            eps::fmpq, strat::Int, verbose::Int = 0 )
    
#     initBox::box = box(initialBox[1],initialBox[2],initialBox[3])
    
#     const getApp_c = cfunction(getApprox, Nothing, (Ref{acb_poly}, Int))
    getApp_c = @cfunction( $getApprox, Cvoid, (Ptr{acb_poly}, Int))
    
    ccall( (:ccluster_DAC_first_interface_forJulia, libccluster), 
             Nothing, (Ref{listConnComp}, Ref{listConnComp}, Ref{listConnComp}, Ref{listConnComp}, 
                    Ref{Nothing},Int, Ref{box}, Ref{fmpq}, Int,   Int), 
                     qRes,              qAllRes,           qMainLoop,         discardedCcs,
                    getApp_c,  nbSols,  initBox,  eps,      strat, verbose )
                    
    return
    
end

function ccluster_DAC_next(qRes::listConnComp,
                            qAllRes::listConnComp,
                            qMainLoop::listConnComp, 
                            discardedCcs::listConnComp,
                            getApprox::Function,
                            nbSols::Int,
#                             initialBox::Array{fmpq,1},
                            initBox::box, 
                            eps::fmpq, strat::Int, verbose::Int = 0 )
    
#     initBox::box = box(initialBox[1],initialBox[2],initialBox[3])
    
#     const getApp_c = cfunction(getApprox, Nothing, (Ref{acb_poly}, Int))
    getApp_c = @cfunction( $getApprox, Cvoid, (Ptr{acb_poly}, Int))
    
    ccall( (:ccluster_DAC_next_interface_forJulia, libccluster), 
             Nothing, (Ref{listConnComp}, Ref{listConnComp}, Ref{listConnComp}, Ref{listConnComp}, 
                    Ref{Nothing}, Int, Ref{box}, Ref{fmpq}, Int,   Int), 
                     qRes,              qAllRes,           qMainLoop,         discardedCcs,
                    getApp_c,  nbSols,  initBox,  eps,      strat, verbose )
    
    return 
    
end

function ccluster_draw( getApprox::Function, 
                        initialBox::Array{fmpq,1}, 
                        precision::Int;
                        strategy="default", #a strategy: Int
                        verbosity="silent" )#a verbosity flag; by defaults, nothing is printed
                                            #options are "silent", "brief" and "results"
    
    initBox::box = box(initialBox[1],initialBox[2],initialBox[3])
    getApp_c = @cfunction( $getApprox, Cvoid, (Ptr{acb_poly}, Int))
    
    verbose::Int = parseVerbosity(verbosity)
    eps = fmpq(1, fmpz(2)^precision)
    
    lccRes = listConnComp()
    lcbDis = listBox()
#     
#     ccall( (:ccluster_interface_forJulia_draw, libccluster), 
#              Nothing, (Ref{listConnComp},Ref{listBox}, Ref{Nothing},    Ref{box}, Ref{fmpq}, Int,   Int), 
#                      lccRes,  lcbDis,          getApp_c,    initBox,  eps,      strat, verbose )
                     
    ccall( (:ccluster_forJulia_draw, libccluster), 
             Nothing, (Ref{listConnComp}, Ref{listBox}, Ref{Nothing},    Ref{box}, Ref{fmpq}, Cstring, Int,   Int), 
                       lccRes,            lcbDis,       getApp_c,        initBox,  eps,       strategy, 1,    verbose )
                     
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

#DEPRECATED interfaces

PGLOBALCCLUSTERFMPQ = fmpq_poly(0);

function getApp_FMPQ( dest::Ref{acb_poly}, prec::Int )
# function getApp_FMPQ( dest::acb_poly, prec::Int )
#     ccall((:acb_poly_set_fmpq_poly, libarb), Nothing,
#                 (Ref{acb_poly}, Ref{fmpq_poly}, Int), dest, PGLOBALCCLUSTERFMPQ, prec)
    ccall((:acb_poly_set_fmpq_poly, libarb), Cvoid,
                (Ref{acb_poly}, Ref{fmpq_poly}, Int), dest, PGLOBALCCLUSTERFMPQ, prec)
end

function ccluster( P_FMPQ::fmpq_poly, initialBox::Array{fmpq,1}, eps::fmpq, verbose::Int)
    
    print(" WARNING: syntax ccluster( ::fmpq_poly, ::box, ::fmpq, verb::Int) is deprecated\n")
    print("             use ccluster( ::fmpq_poly, ::box, ::Int,  verbosity=verb) instead\n")
    
    return ccluster( P_FMPQ, initialBox, eps, 23, verbose)
    
end

function ccluster( P_FMPQ::fmpq_poly, initialBox::Array{fmpq,1}, eps::fmpq, strat::Int, verbose::Int)
    
    print(" WARNING: syntax ccluster( ::fmpq_poly, ::box, ::fmpq, strat::Int, verb::Int) is deprecated\n")
    print("             use ccluster( ::fmpq_poly, ::box, ::Int,  strategy=strat, verbosity=verb) instead\n")
    
    ccall((:fmpq_poly_set, :libflint), Nothing,
                (Ref{fmpq_poly}, Ref{fmpq_poly}), PGLOBALCCLUSTERFMPQ, P_FMPQ)

    return ccluster( getApp_FMPQ, initialBox, eps, strat, verbose)
    
end

function ccluster( getApprox::Function, initialBox::Array{fmpq,1}, eps::fmpq, verbose::Int)
    print(" WARNING: syntax ccluster( ::Function, ::Array{fmpq,1}, ::fmpq, verb::Int) is deprecated\n")
    print("             use ccluster( ::Function, ::Array{fmpq,1}, ::Int,  verbosity=verb) instead\n")
    return ccluster( getApprox, initialBox, eps, 23, verbose)
end

function ccluster( getApprox::Function, initialBox::Array{fmpq,1}, eps::fmpq, strat::Int, verbose::Int)
    print(" WARNING: syntax ccluster( ::Function, ::Array{fmpq,1}, ::fmpq, strat::Int, verb::Int) is deprecated\n")
    print("             use ccluster( ::Function, ::Array{fmpq,1}, ::Int,  strategy=strat, verbosity=verb) instead\n")
    initBox::box = box(initialBox[1],initialBox[2],initialBox[3])

    queueResults = ccluster( getApprox, initBox, eps, strat, verbose )
    
    for sol in queueResults
        tempBO = sol[2]
        sol[2] = [getCenterRe(tempBO),getCenterIm(tempBO),fmpq(3,4)*getWidth(tempBO)]
    end
    
    return queueResults
end

function ccluster( getApprox::Function, initBox::box, eps::fmpq, strat::Int, verbose::Int)
    print(" WARNING: syntax ccluster( ::Function, ::box, ::fmpq, strat::Int, verb::Int) is deprecated\n")
    print("             use ccluster( ::Function, ::box, ::Int,  strategy=strat, verbosity=verb) instead\n")
#     const getApp_c = cfunction(getApprox, Nothing, (Ref{acb_poly}, Int))
    getApp_c = @cfunction( $getApprox, Cvoid, (Ptr{acb_poly}, Int))
    lccRes = listConnComp()
    ccall( (:ccluster_interface_forJulia, libccluster), 
             Nothing, (Ref{listConnComp}, Ptr{Cvoid},    Ref{box}, Ref{fmpq}, Int,   Int), 
                     lccRes,           getApp_c,    initBox,  eps,      strat, verbose )
    queueResults = []
    while !isEmpty(lccRes)
        tempCC = pop(lccRes)
        tempBO = getComponentBox(tempCC,initBox)
        push!(queueResults, [getNbSols(tempCC),tempBO])
    end
    return queueResults
    
end

export ccluster, risolate, printClusters

end # module
