println("Testing Ccluster");
using Test
using AbstractAlgebra
using Nemo
using Ccluster

@testset "box" begin
    B = Ccluster.box( fmpq(0,1), fmpq(0,1), fmpq(2,1) )
    @test isa(B, Ccluster.box)
    C = Ccluster.box(B)
    @test isa( Ccluster.getCenterRe(B), fmpq )
    @test isa( Ccluster.getCenterIm(B), fmpq )
    @test isa( Ccluster.getWidth(B), fmpq )
    @test isa( Ccluster.getCenterRe(C), fmpq )
    @test isa( Ccluster.getCenterIm(C), fmpq )
    @test isa( Ccluster.getWidth(C), fmpq )
    @test Ccluster.getCenterRe(B) == Ccluster.getCenterRe(C)
    @test Ccluster.getCenterIm(B) == Ccluster.getCenterRe(C)
    @test Ccluster.getWidth(B) == Ccluster.getWidth(C)
end

@testset "disk" begin
    D = Ccluster.disk( fmpq(0,1), fmpq(0,1), fmpq(1,1) )
    @test isa( D, Ccluster.disk)
    @test isa( Ccluster.getCenterRe(D), fmpq )
    @test isa( Ccluster.getCenterIm(D), fmpq )
    @test isa( Ccluster.getRadius(D), fmpq )
end

@testset "risolate for fmpq poly" begin
    #integer polynomials
    #Bernoulli polynomial
    R, x = PolynomialRing(Nemo.QQ, "x")
    n = 64 #degree
    P = zero(R)
    bernoulli_cache(n)
    for k = 0:n
        coefficient = (Nemo.binomial(n,k))*(bernoulli(n-k))
        P = P + coefficient*x^k
    end #P is now the Bernoulli polynomial of degree 64
    #local solving
    B = Ccluster.box( fmpq(0,1), fmpq(0,1), fmpq(2,1) )
    Res = risolate( P, B, 53, strategy="default", verbosity="silent" )
    @test length(Res)==4
    @test isa(Res, Array{Any, 1})
    @test isa(Res[1], Array{Any, 1})
    @test isa(Res[1][1], Int64)
    @test isa(Res[1][2], Ccluster.box)
    #same but array of fmpq input
    B = [ fmpq(0,1), fmpq(2,1) ]
    Res = risolate( P, B, 53, strategy="default", verbosity="silent" )
    @test length(Res)==4
    @test isa(Res, Array{Any, 1})
    @test isa(Res[1], Array{Any, 1})
    @test isa(Res[1][1], Int64)
    @test isa(Res[1][2], Array{fmpq, 1})
    #global solving
    Res = risolate( P, 53, strategy="default", verbosity="silent" )
    @test length(Res)==16
    @test isa(Res, Array{Any, 1})
    @test isa(Res[1], Array{Any, 1})
    @test isa(Res[1][1], Int64)
    @test isa(Res[1][2], Array{fmpq, 1})
end

@testset "ccluster for fmpq poly" begin
    #integer polynomials
    #Mignotte polynomial
    R, x = PolynomialRing(Nemo.QQ, "x")
    d = 64 #degree
    bitsize = 14
    bitsize = Int64(floor(bitsize/2))
    P=x^d - 2*((2^bitsize)*x-1)^2
    #local solving
    B = Ccluster.box( fmpq(0,1), fmpq(0,1), fmpq(2,1) )
    Res = ccluster( P, B, 53, strategy="default", verbosity="silent" )
    @test length(Res)==17
    @test isa(Res, Array{Any, 1})
    @test isa(Res[1], Array{Any, 1})
    @test isa(Res[1][1], Int64)
    @test isa(Res[1][2], Ccluster.box)
    #same but array of fmpq input
    B = [ fmpq(0,1), fmpq(0,1), fmpq(2,1) ]
    Res = ccluster( P, B, 53, strategy="default", verbosity="silent" )
    @test length(Res)==17
    @test isa(Res, Array{Any, 1})
    @test isa(Res[1], Array{Any, 1})
    @test isa(Res[1][1], Int64)
    @test isa(Res[1][2], Array{fmpq, 1})
    #global solving
    Res = ccluster( P, 53, strategy="default", verbosity="silent" )
    @test length(Res)==63
    @test isa(Res, Array{Any, 1})
    @test isa(Res[1], Array{Any, 1})
    @test isa(Res[1][1], Int64)
    @test isa(Res[1][2], Array{fmpq, 1})
end

#spiral polynomials
degr=64
function getApproximationSpiral( dest::Ptr{acb_poly}, precision::Int )
    
    function getAppSpiral( degree::Int, prec::Int )::Nemo.acb_poly
        CC = ComplexField(prec)
        R2, y = PolynomialRing(CC, "y")
        res = R2(1)
        for k=1:degree
            modu = fmpq(k,degree)
            argu = fmpq(4*k,degree)
            root = modu*Nemo.exppii(CC(argu))
            res = res * (y-root)
        end
        return res
    end
    
    precTemp::Int = 2*precision
    poly = getAppSpiral( degr, precTemp)
    
    Ccluster.ptr_set_acb_poly(dest, poly)

end

@testset "ccluster for black box poly" begin
    
    #local solving
    B = Ccluster.box( fmpq(0,1), fmpq(0,1), fmpq(1,4) )
    Res = ccluster( getApproximationSpiral, B, 53, strategy="default", verbosity="silent" )
    @test length(Res)==8
    @test isa(Res, Array{Any, 1})
    @test isa(Res[1], Array{Any, 1})
    @test isa(Res[1][1], Int64)
    @test isa(Res[1][2], Ccluster.box)
    #same but array of fmpq input
    B = [ fmpq(0,1), fmpq(0,1), fmpq(1,4) ]
    Res = ccluster( getApproximationSpiral, B, 53, strategy="default", verbosity="silent" )
    @test length(Res)==8
    @test isa(Res, Array{Any, 1})
    @test isa(Res[1], Array{Any, 1})
    @test isa(Res[1][1], Int64)
    @test isa(Res[1][2], Array{fmpq, 1})
    #global solving
    Res = ccluster( getApproximationSpiral, 53, strategy="default", verbosity="silent" )
    @test length(Res)==64
    @test isa(Res, Array{Any, 1})
    @test isa(Res[1], Array{Any, 1})
    @test isa(Res[1][1], Int64)
    @test isa(Res[1][2], Array{fmpq, 1})
end

@testset "random generation of triangular systems" begin
    typeOfSyst = [3,3,3]
    heigthOfSyst = 3
    polys, rings, vars = Ccluster.rand_gen_system( typeOfSyst, heigthOfSyst )
    @test length(polys)==length(rings)==length(vars)==3
    @test isa(polys, Array{Any, 1})
    @test isa(polys[1], fmpq_poly)
    @test isa(polys[2], AbstractAlgebra.Generic.Poly{Nemo.fmpq_poly})
    @test isa(polys[3], AbstractAlgebra.Generic.Poly{AbstractAlgebra.Generic.Poly{Nemo.fmpq_poly}})
    @test isa(rings, Array{Any, 1})
    @test isa(vars, Array{Any, 1})
end

@testset "solving triangular systems" begin
    Rx, x = PolynomialRing(Nemo.QQ, "x") #Ring of polynomials in x with rational coefficients
    Syx, y = PolynomialRing(Rx, "y") #Ring of polynomials in y with coefficients that are in Rx
    d1=30
    c = 10
    delta=128
    d2=10
    twotodelta = fmpq(2,1)^(delta)
    f  = Rx( x^d1 - (twotodelta*x-1)^(c) )
    g = Syx( y^d2 - x^d2 )
    precision = 53
    bInitx = [fmpq(0,1),fmpq(0,1),fmpq(10,1)]
    #local solving
    nbSols, clusters, ellapsedTime = tcluster( [f,g], [bInitx], precision, verbosity = "silent" );
    @test isa(nbSols, Int64)
    @test nbSols==100
    @test isa(clusters, Array{Any, 1})
    @test length(clusters)==1
    @test isa(clusters[1], Array{Any, 1})
    @test isa(clusters[1][1], Int64)
    @test clusters[1][1]==100
    @test isa(clusters[1][2], Array{acb, 1})
    #global solving, precision = 100
    precision = 100
    nbSols, clusters, ellapsedTime = tcluster( [f,g], precision, verbosity = "silent" );
    @test nbSols==300
    @test length(clusters)==201
    #global solving, precision = 200
    precision = 200
    nbSols, clusters, ellapsedTime = tcluster( [f,g], precision, verbosity = "silent" );
    @test nbSols==300
    @test length(clusters)==210
    #global solving, precision = 300
    precision = 300
    nbSols, clusters, ellapsedTime = tcluster( [f,g], precision, verbosity = "silent" );
    @test nbSols==300
    @test length(clusters)==300
end

println("End of testing Ccluster");
