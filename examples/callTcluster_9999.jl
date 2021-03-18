using Nemo
using Ccluster

#maximum size of clusters in output
prec = 53

#load the system to be solved
include( "9999_random.jl" );

#solve the system locally
bInit = [fmpq(0,1),fmpq(0,1),fmpq(2,1)]
nbSols, clusters, ellapsedTime = tcluster( polys, [bInit], prec, verbosity = "brief" );

#solve the system globally
nbSols, clusters, ellapsedTime = tcluster( polys, prec, verbosity = "brief" );
