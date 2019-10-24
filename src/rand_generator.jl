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

using Printf

function rand_gen_uni( var, ring, deg, bitsize )
    res = ring(0)
    inf = -2^(bitsize-1)
    sup =  2^(bitsize-1)
    for index = 0:deg
        monomial = fmpq(rand(inf:sup))*var^index
        res = res + monomial
    end
    return res
end

function rand_gen( vars, rings, degr, bitsize, index = [1] )
#     cur_ring = pop!(rings)
#     cur_var  = pop!(vars)
#     print("var: $(vars[index[1]]), index[1]: $(index[1])\n")
    res = rings[index[1]](0);
    if index[1]==length(vars)
#         print("ici\n")
        res = rand_gen_uni( vars[index[1]], rings[index[1]], degr, bitsize )
    else
        for deg = 0:degr
            cur_var = vars[index[1]]
            index[1] +=1
            res = res + rand_gen( vars, rings, (degr-deg), bitsize, index )*cur_var^deg
            index[1] -=1
        end
    end
    
    #check that res has degree at least 1 in vars[index[1]]
    if (index[1]==1)
        while degree(res)==0
            res = rand_gen( vars, rings, degr, bitsize, [1] )
        end
    end
    return res
end

function rand_gen_multiple( vars, rings, degree, bitsize, index = [1] )
    a = rand_gen( vars, rings, Int(floor(degree/2)), bitsize)
    a = a*a
    b = rand_gen( vars, rings, degree, bitsize, [2])
    c = rand_gen( vars, rings, degree, bitsize, [2])
    res = (b*vars[1] + c)^(Int(floor((degree+1)/2)) - Int(floor(degree/2)))  
    res = a*res
    return res
end

function rand_gen_system( listofdegrees, bitsize )
    #create the variables names
    varnames=[]
    for index in 1:length(listofdegrees)
        push!(varnames, @sprintf("x%d",index))
    end
#     print("$varnames\n");
    #create the first ring
    R, v = PolynomialRing(QQ, varnames[1])
    #create the rings and the vars
    rings=[]
    vars=[]
    push!(rings, R)
    push!(vars, v)
    for index in 2:length(listofdegrees)
        R, v = PolynomialRing(rings[index-1], varnames[index])
        push!(rings, R)
        push!(vars, v)
    end
    reverse!(rings)
    reverse!(vars)
    res = []
    nbpols = length(listofdegrees)
#     print("nbpols: $nbpols\n")
    for index in 1:nbpols
        push!(res, rand_gen(vars[index:nbpols],
                            rings[index:nbpols],
                            listofdegrees[nbpols+1-index],
                            bitsize))
    end
    reverse!(res)
    return res, rings, vars
end

function rand_gen_system_multiple( listofdegrees, bitsize )
    #create the variables names
    varnames=[]
    for index in 1:length(listofdegrees)
        push!(varnames, @sprintf("x%d",index))
    end
#     print("$varnames\n");
    #create the first ring
    R, v = PolynomialRing(QQ, varnames[1])
    #create the rings and the vars
    rings=[]
    vars=[]
    push!(rings, R)
    push!(vars, v)
    for index in 2:length(listofdegrees)
        R, v = PolynomialRing(rings[index-1], varnames[index])
        push!(rings, R)
        push!(vars, v)
    end
    reverse!(rings)
    reverse!(vars)
    res = []
    nbpols = length(listofdegrees)
#     print("nbpols: $nbpols\n")
    for index in 1:nbpols-1
        push!(res, rand_gen_multiple(vars[index:nbpols],
                            rings[index:nbpols],
                            listofdegrees[nbpols+1-index],
                            bitsize))
    end
    push!(res, rand_gen([vars[nbpols]],
                            [rings[nbpols]],
                            listofdegrees[1],
                            bitsize))
    reverse!(res)
    return res, rings, vars
end

function write_HOM4PS_sys( filename, listofpolys )
    open(filename,"w") do f
        write(f, "{\n")
        for index in 1:length(listofpolys)
            write(f," $(listofpolys[index]);\n")
        end
        write(f, "}\n")
    end
end

function write_Bertini_sys( filename, listofpolys, vars )
    open(filename,"w") do f
        write(f, "\nCONFIG\n\nMPTYPE:2;\n\nEND;\n\nINPUT\n\n")
        write(f, "variable_group ")
        reverse!(vars)
        for index in 1:length(vars)-1
            write(f, "$(vars[index]),")
        end
        write(f, "$(vars[length(vars)]);\n")
        reverse!(vars)
        write(f, "function ")
        for index in 1:length(vars)-1
            write(f, "p$(index),")
        end
        write(f, "p$(length(vars));\n\n")
        
        for index in 1:length(vars)
            write(f, "p$(index) = $(listofpolys[index]);\n")
        end
        write(f, "\nEND;\n")
    end
end

function write_julia_sys( filename, listofpolys)
    open(filename,"w") do f
        Rstring = string("R1, x1 = PolynomialRing(QQ, \"x1\")\n")
        write(f, Rstring)
        for index in 2:length(listofpolys)
            Rstring = string("R$(index), x$(index) = PolynomialRing(R$(index-1), \"x$(index)\")\n")
            write(f, Rstring)
        end
        write(f, "polys = [\n")
        for index in 1:length(listofpolys)
            Rstring = string("R$(index)(")
            write(f, Rstring)
            write(f, "$(listofpolys[index]))")
            if index<length(listofpolys)
                write(f, ",\n")
            else 
                write(f, "]\n")
            end
        end
    end
end

# Rx, x = PolynomialRing(QQ, "x")
# a = rand_gen_uni(x, Rx, 5, 10)
# Syx, y = PolynomialRing(Rx, "y")
# r = rand_gen([y,x],[Syx,Rx],5,5)
# Tzyx, z = PolynomialRing(Syx, "z")
# t = rand_gen([z,y,x],[Tzyx,Syx,Rx],5,5)
