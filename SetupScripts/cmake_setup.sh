#### CMAKE environment setup ####
#### Instead of /home/jongho/Software/Cmake, please change the right path depend on your machine ####

if test "x$LD_LIBRARY_PATH" = "x" ; then
    export LD_LIBRARY_PATH="/home/jongho/Software/Cmake/cmake-3.9.5-Linux-x86_64/lib"
else
    export LD_LIBRARY_PATH="/home/jongho/Software/Cmake/cmake-3.9.5-Linux-x86_64/lib":${LD_LIBRARY_PATH}
fi
if test "x$C_INCLUDE_PATH" = "x" ; then
    export C_INCLUDE_PATH="/home/jongho/Software/Cmake/cmake-3.9.5-Linux-x86_64/include"
else
    export C_INCLUDE_PATH="/home/jongho/Software/Cmake/cmake-3.9.5-Linux-x86_64/include":${C_INCLUDE_PATH}
fi
if test "x$CXX_INCLUDE_PATH" = "x" ; then
    export CXX_INCLUDE_PATH="/home/jongho/Software/Cmake/cmake-3.9.5-Linux-x86_64/include"
else
    export CXX_INCLUDE_PATH="/home/jongho/Software/Cmake/cmake-3.9.5-Linux-x86_64/include":${CXX_INCLUDE_PATH}
fi
if test "x$PATH" = "x" ; then
    export PATH="/home/jongho/Software/Cmake/cmake-3.9.5-Linux-x86_64/bin"
else
    export PATH="/home/jongho/Software/Cmake/cmake-3.9.5-Linux-x86_64/bin":${PATH}
fi
if test "x$MANPATH" = "x" ; then
    export MANPATH="/home/jongho/Software/Cmake/cmake-3.9.5-Linux-x86_64/man"
else
    export MANPATH="/home/jongho/Software/Cmake/cmake-3.9.5-Linux-x86_64/man":${MANPATH}
fi
