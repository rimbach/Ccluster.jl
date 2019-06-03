using Nemo

R, x = PolynomialRing(Nemo.QQ, "x")

n = 64 #degree
P = zero(R)
bernoulli_cache(n)
for k = 0:n
    global P
    coefficient = (binom(n,k))*(bernoulli(n-k))
    P = P + coefficient*x^k
end #P is now the Bernoulli polynomial of degree 8

using Ccluster

bInit = [fmpq(0,1),fmpq(0,1),fmpq(100,1)] #box centered in 0 + sqrt(-1)*0 with width 100
eps = fmpq(1, fmpz(2)^10)               #eps = 2^-10
verbosity = 0                           #nothing printed
Coeffs = ccluster(P, bInit, eps, verbosity)

function getApproximation( dest::Ptr{acb_poly}, prec::Int )
    eps = fmpq(1, fmpz(2)^prec)
    Qre = zero(R)
    Qim = zero(R)
    for i=1:n
        btemp = [ Coeffs[i][2][1], Coeffs[i][2][2], 2*Coeffs[i][2][3] ]
        temp = ccluster(P, btemp, eps, 0)
        Qre = Qre + temp[1][2][1]*x^(i-1)
        Qim = Qim + temp[1][2][2]*x^(i-1)
    end
    Ccluster.ptr_set_2fmpq_poly( dest, Qre, Qim, prec )
end

bInit = [fmpq(0,1),fmpq(0,1),fmpq(100,1)] #box centered in 0 + sqrt(-1)*0 with width 100
eps = fmpq(1, 100)                      #eps = 1/100
verbosity = 0                           #nothing printed
Roots = ccluster(getApproximation, bInit, eps, 1)

# using CclusterPlot #only if you have installed CclusterPlot.jl

# plotCcluster(Roots, bInit, false) #use true instead of false to focus on clusters
