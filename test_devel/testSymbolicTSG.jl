using Nemo

R, x = PolynomialRing(Nemo.QQ, "x")

n = 16 #degree
P = zero(R)
Nemo.bernoulli_cache(n)
for k = 0:n
    coefficient = (Nemo.binom(n,k))*(Nemo.bernoulli(n-k))
    P = P + coefficient*x^k
end

# function getAppBern( dest::Ptr{acb_poly}, prec::Int )
#     ccall((:acb_poly_set_fmpq_poly, :libarb), Void,
#                 (Ptr{acb_poly}, Ptr{fmpq_poly}, Int), dest, &P, prec)
# end

function getAppBern( prec::Int )
    CC = ComplexField(prec)
    R2, y = PolynomialRing(CC, "y")
    res = R2(0)
    for k=0:degree(P)
        res = res + CC(coeff(P,k))*y^k
    end
    return res
end

# print("P: $P\n")
# print("Approx of P: $(getAppBern(53))\n\n\n")
Papprox = getAppBern(53)

function symbTS( PA::acb_poly, prec::Int )
    d = degree(PA)
#     ring = Pappox.parent
    CC = ComplexField(prec)
    R3, m = PolynomialRing(CC, "m")
    res = [R3(0)]
    
    for k=0:d
        for i=k:d
            res[k+1] = res[k+1] + Nemo.binom(i,k)*coeff(PA,i)*m^(i-k)
        end
        if k<d
            push!(res, R3(0))
        end
    end
    return res
end

function symbAdd( pols1, pols2, prec::Int)
    d1 = length(pols1) -1
    d2 = length(pols2) -1
    res = []
#     d = min(d1,d2)
    for k=0:min(d1,d2)
        push!(res, pols1[k+1] + pols2[k+1])
#         pols1[k+1] = pols1[k+1] + pols2[k+1]
    end
    if d1<d2
        for k=(d1+1):d2
            push!(res,  pols2[k+1])
#             push!(pols1,  pols2[k+1])
        end
    else
        for k=(d2+1):d1
            push!(res,  pols1[k+1])
        end
    end
    return res
end

function symbSub( pols1, pols2, prec::Int)
    d1 = length(pols1) -1
    d2 = length(pols2) -1
    res = []
    d = min(d1,d2)
    for k=0:min(d1,d2)
        push!(res, pols1[k+1] - pols2[k+1])
#         pols1[k+1] = pols1[k+1] - pols2[k+1]
    end
    if d1<d2
        for k=(d1+1):d2
            push!(res,  -pols2[k+1])
#             push!(pols1,  -pols2[k+1])
        end
    else
        for k=(d2+1):d1
            push!(res,  pols1[k+1])
        end
    end
    return res
end

function symbMonMul( pols, mon::acb_poly, deg::Int, prec::Int)
    d = length(pols) -1
    CC = ComplexField(prec)
    R3, m = PolynomialRing(CC, "m")
    res=[]
#     for k=0:d
#         pols[k+1] = pols[k+1]*mon
#     end
#     
#     pols = reverse(pols)
#     for k=1:deg
#         push!( pols, R3(0) )
#     end
#     pols = reverse(pols)

    for k=1:deg
        push!( res, R3(0) )
    end
    for k=0:d
        push!( res, pols[k+1]*mon )
    end
    return res
end

function symbMul( pols1, pols2, prec::Int )
    res=[]
    d1 = length(pols1) -1
    for k=0:d1
        res = symbAdd( res, symbMonMul( pols2, pols1[k+1], k, prec), prec )
    end
    return res
end

function symbOneGraeffe( pols, prec::Int )
    CC = ComplexField(prec)
    R3, m = PolynomialRing(CC, "m")
    d = length(pols) -1
    polse=[]
    polso=[]
    for k=0:d
        if mod(k,2)==0
            push!(polse, pols[k+1])
        else
            push!(polso, pols[k+1])
        end
    end
    polse = symbMul(polse,polse,prec)
    polso = symbMul(polso,polso,prec)
    insert!(polso, 1, R3(0))
    res = symbSub( polse, polso, prec)
    if mod(d,2)==1
        for k=0:d
            res[k+1] = -res[k+1]
        end
    end
    return res
end

function symbNGraeffe( pols, prec::Int)
    d = length(pols) -1
    N = round(Int,4+ceil(log2(1+log2(d))))
    res = symbOneGraeffe( pols, prec )
    for i=2:N
        res = symbOneGraeffe( res, prec )
    end
    return res
end

function evaluateSymbTS( pols, centerRe::fmpq, centerIm::fmpq, prec::Int )
    d = length(pols) -1
    CC = ComplexField(prec)
    R2, y = PolynomialRing(CC, "y")
    res = R2(0)
    for k=0:d
        res = res + evaluate( pols[k+1], CC(centerRe, centerIm) )*y^k
    end
    return res
end

function evaluateSymbTS_first( pols, centerRe::fmpq, centerIm::fmpq, prec::Int )
    d = length(pols) -1
    CC = ComplexField(prec)
#     R2, y = PolynomialRing(CC, "y")
#     res = R2(0)
#     for k=0:d
#         res = res + evaluate( pols[k+1], CC(centerRe, centerIm) )*y^k
#     end
#     return res
    return evaluate( pols[d], CC(centerRe, centerIm) )
end

function numeTS(Papprox, centerRe::fmpq, centerIm::fmpq, prec )
    CC = ComplexField(prec)
    R2, y = PolynomialRing(CC, "y")
    comp = y+CC(centerRe, centerIm)
    return Nemo.compose(Papprox, comp)
end

function numeOneGraeffe(Papprox, prec)
    CC = ComplexField(prec)
    R2, y = PolynomialRing(CC, "y")
    d=degree(Papprox)
    Pe = R2(0)
    Po = R2(0)
    for k=0:d
        if mod(k,2)==0
            Pe = Pe + coeff(Papprox,k)*y^(Int(floor(k/2)))
        else
            Po = Po + coeff(Papprox,k)*y^(Int(floor(k/2)))
        end
    end
    return ((-1)^d)*(Pe^2 - y*Po^2)
end

function numeNGraeffe( Papprox, prec::Int)
    d = degree(Papprox)
    N = round(Int,4+ceil(log2(1+log2(d))))
    res = numeOneGraeffe( Papprox, prec )
    for i=2:N
        res = numeOneGraeffe( res, prec )
    end
    return res
end

precnum = 53
Papprox = getAppBern(precnum)
Papprox = numeTS(Papprox, fmpq(1,1), fmpq(1,1), precnum)
# Papprox = numeOneGraeffe(Papprox, precnum)
Papprox = numeNGraeffe(Papprox, precnum)

precsymb = 113
P1 = symbTS( getAppBern(precsymb), precsymb )
# P1 = symbOneGraeffe( P1, precsymb )
P1 = symbNGraeffe( P1, precsymb )
# P2 = evaluateSymbTS( P1, fmpq(1,1), fmpq(1,1), precsymb)
temp = evaluateSymbTS_first( P1, fmpq(1,1), fmpq(1,1), precsymb )
print("coeff $(n-1) of Papprox: $(coeff(Papprox,n-1))\n\n")
# print("P2: $P2\n\n")
print("temp: $temp\n\n")

tic()
# P1 = symbTS(Papprox,53)
P1 = symbTS( getAppBern(precsymb), precsymb )
P1 = symbNGraeffe( P1, precsymb )
ellapsedTime = toq()
# 
tic()
for i=1:100
#     P2 = evaluateSymbTS( P1, fmpq(1,1), fmpq(1,1), 53)
    temp = evaluateSymbTS_first( P1, fmpq(1,1), fmpq(1,1), precsymb )
end
ellapsedTime1 = toq()
# 
tic()
for i=1:100
    P3 = numeTS(Papprox, fmpq(1,1), fmpq(1,1), precnum)
#     P3 = numeNGraeffe(P3, precnum)
end
ellapsedTime2 = toq()
# 
# # print("P2: $P2\n\n")
# # print("P3: $P3\n\n")
print("symbolic: $ellapsedTime\n")
print("with evaluations: $ellapsedTime1\n")
print("with numeric: $ellapsedTime2\n")
