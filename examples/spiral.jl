using Nemo
using Ccluster

bInit = [fmpq(0,1),fmpq(0,1),fmpq(4,1)] #box centered in 0 + sqrt(-1)*0 with width 4
eps = fmpq(1,100)                       #eps = 1/100
verbosity = 0                           #nothing printed

degr=64
function getApproximation( dest::Ptr{acb_poly}, precision::Int )
    
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
    
    while Ccluster.checkAccuracy( poly, precision ) == 0
            precTemp = 2*precTemp
            poly = getAppTemp(precTemp)
    end
    
    Ccluster.ptr_set_acb_poly(dest, poly)

end

Res = ccluster(getApproximation, bInit, eps, verbosity)

using CclusterPlot #only if you have installed CclusterPlot.jl

plotCcluster(Res, bInit, false) #use true instead of false to focus on clusters
