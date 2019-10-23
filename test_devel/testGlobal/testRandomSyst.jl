using Nemo
using Ccluster

typeOfSyst = [9,9,9]
heigthOfSyst = 10

nbOfSols = 1
for i=1:length(typeOfSyst)
  global nbOfSols
  nbOfSols = nbOfSols*typeOfSyst[i]
end

polys, rings, vars = Ccluster.rand_gen_system( typeOfSyst, heigthOfSyst )

precision = 53

nbSols, clusters, ellapsedTime = tcluster( polys, precision, verbosity = "silent" );
print("time to solve the system: $ellapsedTime \n")
print("number of clusters: $(length(clusters))\n")
print("number of solutions: $(nbSols)\n")

#comparison
bInitx = [fmpq(0,1),fmpq(0,1),fmpq(10,1)^40]
nbSols, clusters, ellapsedTime = tcluster( polys, [bInitx], precision, verbosity = "silent" );
print("time to solve the system: $ellapsedTime \n")
print("number of clusters: $(length(clusters))\n")
print("number of solutions: $(nbSols)\n")

print("          should be: $(nbOfSols)\n")
