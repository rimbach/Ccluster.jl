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

export TCLUSTER_POLS, TCLUSTER_CFEV, TCLUSTER_CLUS, TCLUSTER_DEGS, TCLUSTER_STRA, TCLUSTER_VERB

# global variables
TCLUSTER_POLS = [[]]
TCLUSTER_CFEV = []
TCLUSTER_CLUS = [[]]
TCLUSTER_DEGS = [[]]
TCLUSTER_STRA = [23] 
TCLUSTER_VERB = [0] 

### Main function
function clusterTriSys(b::Array{Ccluster.box,1}, prec::Int)::Array{Array{Ccluster.algClus,1},1}

    global TCLUSTER_STRA, TCLUSTER_VERB
    global TCLUSTER_POLS, TCLUSTER_CFEV, TCLUSTER_CLUS, TCLUSTER_DEGS
    
    actualPol::Int = length(b)
    
    if actualPol==1 #terminal case
        clusters::Array{Array{Ccluster.algClus,1},1} = Ccluster.clusterPol(b[1], prec)
    else            #other cases
        btemp::Ccluster.box = pop!(b)
        clusterstemp::Array{Array{Ccluster.algClus,1},1} = clusterTriSys( b, prec )
        clusters=[]
        while length(clusterstemp)>0
            clus::Array{Ccluster.algClus,1}=pop!(clusterstemp)
            clusterstemp2::Array{Array{Ccluster.algClus,1},1} = Ccluster.clusterPolInFiber(clus, btemp, prec)
            while length(clusterstemp2)>0
                push!(clusters, pop!(clusterstemp2))
            end
        end
        
    end 
    
    return clusters
end

### Compute clusters for first equation: univariate polynomial
# approximation function
function getAppFirst( dest::Ptr{acb_poly}, prec::Int )::Cvoid
    global TCLUSTER_POLS
    ccall((:acb_poly_set_fmpq_poly, :libarb), 
            Cvoid, (Ptr{acb_poly}, Ref{fmpq_poly}, Int), 
                   dest,           TCLUSTER_POLS[1][1],            prec)
end

# version for box
function clusterPol(b::Ccluster.box, prec::Int)::Array{Array{Ccluster.algClus,1},1}

    global TCLUSTER_STRA, TCLUSTER_VERB
    global TCLUSTER_POLS, TCLUSTER_CFEV, TCLUSTER_CLUS, TCLUSTER_DEGS
    
    TCLUSTER_CFEVsave::Array{Ccluster.algClus,1} = TCLUSTER_CFEV[1]
    
    TCLUSTER_CFEV[1]=[]
    clusters=Array{Ccluster.algClus,1}[]
    eps = Nemo.fmpq(1, Nemo.fmpz(2)^(prec-1))
    
#     qRes::Ccluster.listConnComp = Ccluster.ccluster_solve(getAppSys, b, eps, TCLUSTER_STRA[1], TCLUSTER_VERB[1]);
    qRes::Ccluster.listConnComp = Ccluster.ccluster_solve(getAppFirst, b, eps, TCLUSTER_STRA[1], TCLUSTER_VERB[1]);
    while !Ccluster.isEmpty(qRes)
        objCC, ptrCC = Ccluster.pop_obj_and_ptr(qRes)
        push!(clusters, [Ccluster.algClus(objCC, ptrCC, b, prec)] )
    end
    
    TCLUSTER_CFEV[1] = TCLUSTER_CFEVsave
    return clusters
end

# version for Ccluster.algClus
function clusterPol(b::Ccluster.algClus, prec::Int)::Array{Array{Ccluster.algClus,1},1}
    
    global TCLUSTER_STRA, TCLUSTER_VERB
    global TCLUSTER_POLS, TCLUSTER_CFEV, TCLUSTER_CLUS, TCLUSTER_DEGS
    
    TCLUSTER_CFEVsave::Array{Ccluster.algClus,1} = TCLUSTER_CFEV[1]
    TCLUSTER_CFEV[1]=[]
#     clusters::Array{Array{Ccluster.algClus,1},1} = Ccluster.refine_algClus( b, getAppSys, prec, TCLUSTER_STRA[1], TCLUSTER_VERB[1])
    clusters::Array{Array{Ccluster.algClus,1},1} = Ccluster.refine_algClus( b, getAppFirst, prec, TCLUSTER_STRA[1], TCLUSTER_VERB[1])
    TCLUSTER_CFEV[1] = TCLUSTER_CFEVsave
    return clusters
    
end

### find the next floor of a TAC that has to be refined
function nextFloorToRefine( actualPol::Int, clus::Array{Ccluster.algClus,1}, prec::Int)::Int

    global TCLUSTER_CFEV, TCLUSTER_DEGS
    
    degrees::Array{Int,1} = TCLUSTER_DEGS[1][actualPol]
    curPrecs::Array{Int,1} = Ccluster.getPrecs(TCLUSTER_CFEV[1])
    res::Int = actualPol -1
    
    while (res>=1)&&((degrees[res]==0)||(curPrecs[res]>=prec))
        res = res - 1
    end
    
    return res
    
end

### Compute clusters for other equations: recursive
function clusterPolInFiber(a::Array{Ccluster.algClus,1}, b::Ccluster.box, prec::Int)::Array{Array{Ccluster.algClus,1},1}
    
    global TCLUSTER_STRA, TCLUSTER_VERB
    global TCLUSTER_POLS, TCLUSTER_CFEV, TCLUSTER_CLUS, TCLUSTER_DEGS
    
    TCLUSTER_CFEVsave::Array{Ccluster.algClus,1} = TCLUSTER_CFEV[1]
    #floor d of the TAC
    actualPol::Int = length(a) + 1
    #push a in TCLUSTER_CLUS[1][actualPol-1] that should be empty
    push!(TCLUSTER_CLUS[1][actualPol-1], a)
    
    #initialize clusters
    clusters=Array{Ccluster.algClus,1}[]
    while length( TCLUSTER_CLUS[1][actualPol-1] )>0
        
        c::Array{Ccluster.algClus,1} = pop!(TCLUSTER_CLUS[1][actualPol-1])
        TCLUSTER_CFEV[1] = c
        eps::Nemo.fmpq = Nemo.fmpq(1, Nemo.fmpz(2)^(prec-1))
        
        qRes::Ccluster.listConnComp = Ccluster.ccluster_solve(getAppSys, b, eps, TCLUSTER_STRA[1], TCLUSTER_VERB[1]);
        c = TCLUSTER_CFEV[1]
        
        while !Ccluster.isEmpty(qRes)
            cc::Array{Ccluster.algClus,1} = Ccluster.clusCopy(c)
            objCC, ptrCC = Ccluster.pop_obj_and_ptr(qRes)
            push!(cc, Ccluster.algClus( objCC, ptrCC, b, prec ) )
            push!(clusters, cc )
        end
            
    end
    
    TCLUSTER_CFEV[1] = TCLUSTER_CFEVsave
    
    return clusters
end

#version for Ccluster.algClus
function clusterPolInFiber(a::Array{Ccluster.algClus,1}, b::Ccluster.algClus, prec::Int)::Array{Array{Ccluster.algClus,1},1}
    
    global TCLUSTER_STRA, TCLUSTER_VERB
    global TCLUSTER_POLS, TCLUSTER_CFEV, TCLUSTER_CLUS, TCLUSTER_DEGS
    
    TCLUSTER_CFEVsave::Array{Ccluster.algClus,1} = TCLUSTER_CFEV[1]
    #floor d of the TAC
    actualPol::Int = length(a) + 1
    #push a in TCLUSTER_CLUS[1][actualPol-1] that should be empty
    push!(TCLUSTER_CLUS[1][actualPol-1], a)
    
    #initialize clusters
    clusters=Array{Ccluster.algClus,1}[]
    while length( TCLUSTER_CLUS[1][actualPol-1] )>0
        
        c::Array{Ccluster.algClus,1} = pop!(TCLUSTER_CLUS[1][actualPol-1])
        TCLUSTER_CFEV[1] = c
        qRes::Array{Array{Ccluster.algClus,1},1} = Ccluster.refine_algClus(b, getAppSys, prec, TCLUSTER_STRA[1], TCLUSTER_VERB[1])
        c = TCLUSTER_CFEV[1]
        
        for index = 1:length(qRes)
                cc::Array{Ccluster.algClus,1} = Ccluster.clusCopy(c)
                push!(cc, qRes[index][1])
                push!(clusters, cc) 
        end
            
    end
    
    TCLUSTER_CFEV[1] = TCLUSTER_CFEVsave
    
    return clusters
end

function getAppSys( dest::Ptr{acb_poly}, prec::Int )::Cvoid

    global TCLUSTER_STRA, TCLUSTER_VERB
    global TCLUSTER_POLS, TCLUSTER_CFEV, TCLUSTER_CLUS, TCLUSTER_DEGS
    
    actualPol::Int = length(TCLUSTER_CFEV[1]) + 1
    Ptemp = TCLUSTER_POLS[1][actualPol]
    
    if actualPol==1 #should never enter here
        ccall((:acb_poly_set_fmpq_poly, :libarb), 
            Cvoid, (Ptr{acb_poly}, Ref{fmpq_poly}, Int), 
                   dest,           Ptemp,            prec)
    else
        #find the higher floor to be refined
        ind::Int = Ccluster.nextFloorToRefine( actualPol, TCLUSTER_CFEV[1], prec)
        
        while ind>0
        
            #deconstruct the tower
            upperpart::Array{Ccluster.algClus,1} = TCLUSTER_CFEV[1][ind + 1 : length(TCLUSTER_CFEV[1])]
            TCLUSTER_CFEV[1] = TCLUSTER_CFEV[1][1:ind]
            #refine the tower
            if length( TCLUSTER_CFEV[1] ) ==1 # TCLUSTER_CFEV[1] contains just an algebraic cluster
                local clus, clust, clusters
                clus::Ccluster.algClus = TCLUSTER_CFEV[1][1]
                clusters::Array{Array{Ccluster.algClus,1},1} = Ccluster.clusterPol(clus, prec)
                clust::Array{Ccluster.algClus,1} = pop!(clusters)
                Ccluster.copyIn( TCLUSTER_CFEV[1][1], clust[1] )
            else # TCLUSTER_CFEV[1] is a TAC
                local clus
                clus = pop!(TCLUSTER_CFEV[1]) #the last floor of the TAC
                clusters = Ccluster.clusterPolInFiber(TCLUSTER_CFEV[1], clus, prec)
                clust = pop!(clusters)
                push!(TCLUSTER_CFEV[1],clust[length(clust)]) #just to extend the size of TCLUSTER_CFEV[1]
                Ccluster.copyIn( TCLUSTER_CFEV[1], clust )
            end
            #reconstruct the towers
            for j = 1:length(upperpart)
                for i = 1:length(clusters)
                    push!( clusters[i], Ccluster.clusCopy( upperpart[j] ) )
                end
                push!(TCLUSTER_CFEV[1], upperpart[j])
            end
            # push additionnal clusters in queue
            while length(clusters)>0
                push!(TCLUSTER_CLUS[1][actualPol-1], pop!(clusters))
            end
            # check if there is still something to refine
            ind = Ccluster.nextFloorToRefine( actualPol, TCLUSTER_CFEV[1], prec)
            
        end
        
        approx::Array{Nemo.acb,1} = Ccluster.getApproximation(TCLUSTER_CFEV[1],prec)
        
        Ptemp2::Nemo.acb_poly = Ccluster.getPolAt(Ptemp,approx,prec)
        
#         print("actualPol: $actualPol, Ptemp: $Ptemp \n\n")
        
        ccall((:acb_poly_set, :libarb), 
            Cvoid, (Ptr{acb_poly}, Ref{acb_poly}, Int), 
                   dest,          Ptemp2,         prec)
    end   
end
