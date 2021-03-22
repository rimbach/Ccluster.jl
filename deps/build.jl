using Libdl

oldwdir = pwd()

@show CCLUSTER_VERSION = "v1.1.5"

pkgdir = dirname(dirname(@__FILE__))
wdir = joinpath(pkgdir, "deps")
vdir = joinpath(pkgdir, "local")
# NemoLibsDir = joinpath(pkgdir, "../../packages/Nemo/9nu4c/deps/usr")
NemoLibsDir = Base.find_package("Nemo")
NemoLibsDir = Base.Filesystem.dirname(NemoLibsDir)
NemoLibsDir = Base.Filesystem.dirname(NemoLibsDir)
NemoLibsDir = joinpath( NemoLibsDir, "deps/usr" )

if Sys.isapple() && !("CC" in keys(ENV))
   ENV["CC"] = "clang"
   ENV["CXX"] = "clang++"
end

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

LDFLAGS = "-Wl,-rpath,$vdir/lib -Wl,-rpath,\$\$ORIGIN/../share/julia/site/v$(VERSION.major).$(VERSION.minor)/Nemo/local/lib"
DLCFLAGS = "-fPIC -fno-common"

# INSTALL CCLUSTER #temp

cd(wdir)

function download_dll(url_string, location_string)
   try
      run(`curl -o $(location_string) -L $(url_string)`)
   catch
      download(url_string, location_string)
   end
end

cd(wdir)

if !Sys.iswindows()
  println("Cloning Ccluster ... ")
  try
    run(`git clone https://github.com/rimbach/Ccluster.git`)
    cd(joinpath("$wdir", "Ccluster"))
    run(`git checkout $CCLUSTER_VERSION`)
    cd(wdir)
  catch
    if ispath(joinpath("$wdir", "Ccluster"))
      #open(`patch -R --forward -d arb -r -`, "r", open("../deps-PIE-ftbfs.patch"))
      cd(joinpath("$wdir", "Ccluster"))
      run(`git fetch`)
      run(`git checkout $CCLUSTER_VERSION`)
      cd(wdir)
    end
  end
# #   open(`patch --forward -d arb -r -`, "r", open("../deps-PIE-ftbfs.patch"))
  println("DONE")
end

cd(wdir)

if Sys.iswindows()
    if Int == Int32
        println("No binaries for 32 bits windows yet ... ")
    else
        println("downloading binaries ... ")
        download_dll("https://github.com/rimbach/Ccluster/releases/download/v1.1.5/libccluster.dll", joinpath(vdir, "lib", "libccluster.dll"))
        try
            run(`ln -sf $NemoLibsDir\\bin\\libflint.dll $vdir\\lib\\libflint-15.dll`)
        catch
            cp(joinpath(NemoLibsDir, "bin", "libflint.dll"), joinpath(vdir, "lib", "libflint-15.dll"), remove_destination=true)
        end
        try
            run(`ln -sf $NemoLibsDir\\bin\\libarb.dll $vdir\\lib\\libarb-2.dll`)
        catch
            cp(joinpath(NemoLibsDir, "bin", "libarb.dll"), joinpath(vdir, "lib", "libarb-2.dll"), remove_destination=true)
        end
    end
    println("DONE")
else
    println("Building Ccluster ... ")
    cd(joinpath("$wdir", "Ccluster"))
    withenv("LD_LIBRARY_PATH"=>"$vdir/lib", "LDFLAGS"=>LDFLAGS) do
      run(`./configure --prefix=$vdir --disable-static --enable-shared --disable-pthread --with-flint=$NemoLibsDir --with-arb=$NemoLibsDir`)
      run(`make library -j4`)
      run(`make install`)
    end
    println("DONE")
end

push!(Libdl.DL_LOAD_PATH, joinpath(vdir, "lib"))

cd(oldwdir)

#TODO
