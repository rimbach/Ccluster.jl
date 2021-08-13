using Nemo

R, x = PolynomialRing(Nemo.QQ, "x")

nbIter = 3

function iterate( table, n, CC )
    
    index = 1
    bound = length(table)
    while index<=bound
        p = splice!(table,1)
#         modu = fmpq(1,3^(2*n-1))
        modu = fmpq(1,4^(2*(n-1)))
        argu = fmpq(1,3)
        root1 = modu*Nemo.exppii(CC(2*argu))
        root2 = modu*Nemo.exppii(CC(4*argu))
        root3 = modu*Nemo.exppii(CC(6*argu))
        append!(table, [p-root1, p-root2, p-root3])
        index+=1
    end
end

function getApproximation( dest::Ptr{acb_poly}, precision::Int )

    function getAppClusters( it::Int, prec::Int  )::Nemo.acb_poly
#         print("it: $(it)\n")
        CC = ComplexField(prec)
        R2, y = PolynomialRing(CC, "y")
        table = [(y-0)]
        for i=1:it
            iterate( table, i, CC )
        end
        res = R2(1)
        for i=1:length(table)
            res = res*table[i]
        end
        return res
    end
    precTemp::Int = 2*precision
    poly = getAppClusters( nbIter, precTemp)
#     print("poly: $(poly)\n")
#     print("precision: $precision \n");
#     for i=0:degree(poly)
#         print("$i-th coeff, accuracy_bits: $(accuracy_bits( coeff(poly, i) )) \n");
#     end
#     
#     while Ccluster.checkAccuracy( poly, precision ) == 0
#             precTemp = 2*precTemp
#             poly = getAppClusters( nbIter, precTemp)
# #             print("poly: $(poly)\n")
#             print("poly: $(poly)\n")
#             print("precision: $precision \n");
#             for i=0:degree(poly)
#                 print("$i-th coeff, accuracy_bits: $(accuracy_bits( coeff(poly, i) )) \n");
#             end
#     end
    
    Ccluster.ptr_set_acb_poly(dest, poly)

end

using Ccluster

Res = ccluster(getApproximation, 53, verbosity="silent")

using CclusterPlot #only if you have installed CclusterPlot.jl

plotCcluster(Res) #use true instead of false to focus on clusters
