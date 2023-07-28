#{ sources ? import ./nix/sources.nix { } }:
let
  #pkgs = import <nixpkgs> { system = "x86_64-darwin"; };
  pkgs = import <nixpkgs> {};
  sources = import ./nix/sources.nix;
  # Use clang
  #stdenv = pkgs.gccStdenv;
  # Allow broken packages
  #allowBroken = true;
  # Use darwin frameworks if on OSX
  frameworks = pkgs.darwin.apple_sdk.frameworks;
  # Fix hdf5 output for python
  hdf5 = pkgs.symlinkJoin {
    name = "hdf5";
    paths = [
      pkgs.hdf5
      pkgs.hdf5.dev
    ];
  };
in pkgs.mkShell {
  name = "decoNascentEnv";
  buildInputs = [
    # Generic dependencies
    pkgs.git
    # C/C++ Dependencies
    pkgs.gcc
    #pkgs.gfortran
    pkgs.ccacheWrapper
    pkgs.clang
    pkgs.bedtools
    pkgs.blas
    pkgs.gettext
    pkgs.libiconv
    pkgs.libsvm
    pkgs.liblinear
    pkgs.pv
    pkgs.mawk
    pkgs.pigz
    # Library Dependencies
    pkgs.netcdf
    pkgs.graphviz
    pkgs.harfbuzz
    pkgs.fribidi
    pkgs.zlib
    pkgs.gsl
    pkgs.curl
    pkgs.bzip2
    pkgs.libpng
    # Python Dependencies
    pkgs.python39
    pkgs.poetry
    pkgs.black
    # R dependencies
    pkgs.R
    pkgs.rPackages.tidyverse
    pkgs.rPackages.plyr
    pkgs.rPackages.ggthemes
    pkgs.rPackages.R_utils
    pkgs.rPackages.e1071
    pkgs.rPackages.LiblineaR	
    pkgs.rPackages.progress
    pkgs.rPackages.optparse
    pkgs.rPackages.MASS
    pkgs.rPackages.corrplot
    pkgs.rPackages.Matrix
    pkgs.rPackages.distances
    pkgs.rPackages.rstan
    pkgs.rPackages.codetools
    pkgs.rPackages.languageserver
    # System dependencies
    frameworks.Accelerate
    frameworks.Foundation
  ];
  shellHook = ''
    # Update library paths
    export PKG_CONFIG_PATH=$PKG_CONFIG_PATH:${pkgs.graphviz}/lib/pkgconfig
    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${pkgs.stdenv.cc.cc.lib}/lib
    export HDF5_DIR=${hdf5}
    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${pkgs.hdf5}/lib
    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${pkgs.netcdf}/lib
    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${pkgs.graphviz}/lib
    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${pkgs.gettext}/lib
    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${pkgs.libiconv}/lib
    mkdir -p "$(pwd)/_libs"
    export R_LIBS_USER="$(pwd)/_libs"
    # Rscript -e "install.packages(c('nnls'), repos='https://cloud.r-project.org/', dependencies=TRUE, type='source')"
    # source $(poetry env info --path)/bin/activate
  '';

}
