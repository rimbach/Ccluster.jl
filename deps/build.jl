oldwdir = pwd()

@show CCLUSTER_VERSION = "master"

pkgdir = dirname(dirname(@__FILE__))
wdir = joinpath(pkgdir, "deps")
vdir = joinpath(pkgdir, "local")
NemoLibsDir = joinpath(pkgdir, "..", "Nemo", "local")

if !ispath(vdir)

    mkdir(vdir)

    if !ispath(joinpath(vdir, "lib"))
        mkdir(joinpath(vdir, "lib"))
    end
else
    println("Deleting old $vdir")
    rm(vdir, force=true, recursive=true)
    mkdir(vdir)
    mkdir(joinpath(vdir, "lib"))
end

# INSTALL CCLUSTER 

if !is_windows()
#   println("Cloning Ccluster ... ")
#   try
#     run(`git clone https://github.com/rimbach/Ccluster.git`)
#     cd(joinpath("$wdir", "Ccluster"))
#     run(`git checkout CCLUSTER_VERSION`)
#     cd(wdir)
#   catch
# #     if ispath(joinpath("$wdir", "Ccluster"))
# #       open(`patch -R --forward -d arb -r -`, "r", open("../deps-PIE-ftbfs.patch"))
# #       cd(joinpath("$wdir", "arb"))
# #       run(`git fetch`)
# #       run(`git checkout $ARB_VERSION`)
# #       cd(wdir)
# #     end
#   end
# #   open(`patch --forward -d arb -r -`, "r", open("../deps-PIE-ftbfs.patch"))
#   println("DONE")
end

if is_windows()
    if Int == Int32
        println("No binaries for 32 bits windows yet ... ")
    else
    end
else
#     println("Building ccluster ... ")
#     cd(joinpath("$wdir", "arb"))
#    withenv("LD_LIBRARY_PATH"=>"$vdir/lib", "LDFLAGS"=>LDFLAGS) do
#       run(`./configure --prefix=$vdir --disable-static --enable-shared --with-mpir=$vdir --with-mpfr=$vdir --with-flint=$vdir`)
#       run(`make -j4`)
#       run(`make install`)
#    end
#    println("DONE")
end
#TODO
