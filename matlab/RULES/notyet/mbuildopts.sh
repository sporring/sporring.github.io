#
# mbuildopts.sh Shell script for configuring mbuild, the C/C++ stand alone 
#   applications creation script.
#
# usage:        Do not call this file directly; it is sourced by the
#               mbuild shell script.  Modify only if you don't like the
#               defaults after running mbuild.  No spaces are allowed
#               around the '=' in the variable assignment.
#
# SELECTION_TAGs occur in template option files and are used by MATLAB
# tools, such as mex and mbuild, to determine the purpose of the contents
# of an option file. These tags are only interpreted when preceded by '#'
# and followed by ':'.
#
#SELECTION_TAG_ML_OPT: Build and link with MATLAB C/C++ Math Library
#
# Copyright (c) 1984-1998 by The MathWorks, Inc.
# All Rights Reserved.
# $Revision: 1.1 $  $Date: 2003-11-24 17:04:09 $
#----------------------------------------------------------------------------
    case "$Arch" in
        Undetermined)
#----------------------------------------------------------------------------
# Change this line if you need to specify the location of the TMW_ROOT
# root directory.  The mbuild script needs to know where to find utility
# routines so that it can determine the architecture; therefore, this
# assignment needs to be done while the architecture is still
# undetermined.
#----------------------------------------------------------------------------
            TMW_ROOT="$TMW_ROOT"
            ;;
        alpha)
#----------------------------------------------------------------------------
# Alpha C compiler options
            CC='cc'
            CFLAGS="-I$TMW_ROOT/extern/include -ieee"
            CLIBS="-L$TMW_ROOT/extern/lib/$Arch $OTHER_LIBS -lmmfile -lmcc -lmatlb -lmat -lmx -lm"
            COPTIMFLAGS='-O -DNDEBUG'
            CDEBUGFLAGS='-g'
            
            SHLCFLAGS="$CFLAGS"

# Alpha C++ compiler options
	    CPP=cxx
            CPPFLAGS="-I$TMW_ROOT/extern/include/cpp -I$TMW_ROOT/extern/include -x cxx -std1 -DALPHA  -DUNIX -Olimit 5000 -ieee"
            CPPLIBS="-L$TMW_ROOT/extern/lib/$Arch $OTHER_LIBS -lmatpp -lmmfile -lmcc -lmatlb -lmat -lmx -lots -lm"
            CPPOPTIMFLAGS='-O -DNDEBUG'
            CPPDEBUGFLAGS='-g'
#
            FC='f77'
            FFLAGS="-I$TMW_ROOT/extern/include"
            FLIBS="-L$TMW_ROOT/extern/lib/$Arch $OTHER_LIBS -lmmfile -lmcc -lmatlb -lmat -lmx -lm"
            FOPTIMFLAGS='-O'
            FDEBUGFLAGS='-g'
#
            LD="$COMPILER"
            LDFLAGS=''
            LDOPTIMFLAGS=''
            LDDEBUGFLAGS=''

            SHLLD='ld'
            SHLLDFLAGS="-shared -expect_unresolved \* -hidden \`cat $EXPFILE\`"
            SHLMAKEDEF="awk '{printf \"-exported_symbol %s\\n\", \$0;}'"
#----------------------------------------------------------------------------
            ;;
        hp700)
#----------------------------------------------------------------------------
# HP C compiler flags
            CC='cc'
            CFLAGS="-Wp,-H65535 -Aa -I$TMW_ROOT/extern/include +DA1.1"
            CLIBS="-L$TMW_ROOT/extern/lib/$Arch $OTHER_LIBS -lmmfile -lmcc -lmatlb -lmat -lmx -lm"
            COPTIMFLAGS='-O -DNDEBUG'
            CDEBUGFLAGS='-g'

            SHLCFLAGS="+z -D_HPUX_SOURCE -Wp,-H65535 -Aa +DAportable -I$TMW_ROOT/extern/include"

# HP C++ compiler flags
            CPP='aCC'
            CPPFLAGS="+d -I$TMW_ROOT/extern/include/cpp -I$TMW_ROOT/extern/include -D_HPUX_SOURCE -DHP700 -DUNIX +DA1.1"
            CPPLIBS="-L$TMW_ROOT/extern/lib/$Arch $OTHER_LIBS -lmatpp -lmmfile -lmcc -lmatlb -lmat -lmx -lm"
            CPPOPTIMFLAGS='-O -DNDEBUG'
            CPPDEBUGFLAGS='-g'
            CPPLDFLAGS="-I$TMW_ROOT/extern/include/cpp -I$TMW_ROOT/extern/include -D_HPUX_SOURCE -DHP700 -DUNIX +DA1.1"
#
            FC='f77'
            FFLAGS="-I$TMW_ROOT/extern/include +DA1.1"
            FLIBS="-Wl,-L$TMW_ROOT/extern/lib/$Arch $OTHER_LIBS -lmmfile -lmcc -lmatlb -lmat -lmx -lm"
            FOPTIMFLAGS='-O'
            FDEBUGFLAGS='-g'
#
            LD="$COMPILER"
            LDFLAGS=''
            LDOPTIMFLAGS=''
            LDDEBUGFLAGS=''

            SHLLD='ld'
            SHLLDFLAGS="-b -c $EXPFILE"
            SHLMAKEDEF="awk '{printf \"+e %s\\n\", \$0;}'"
#----------------------------------------------------------------------------
            ;;
        ibm_rs)
#----------------------------------------------------------------------------
# IBM RS/6000 C compiler flags
            CC='cc'
            CFLAGS="-I$TMW_ROOT/extern/include"
            CLIBS="-L$TMW_ROOT/extern/lib/$Arch $OTHER_LIBS -lmmfile -lmcc -lmatlb -lmat -lmx -lm"
            COPTIMFLAGS='-O -DNDEBUG'
            CDEBUGFLAGS='-g'

            SHLCFLAGS="$CFLAGS"

# IBM RS/6000 C++ compiler flags
	    CPP='xlC'
            CPPFLAGS="-I$TMW_ROOT/extern/include/cpp -I$TMW_ROOT/extern/include -+ -Q -qfloat=nomaf -qlanglvl=ansi -qspill=2048 -D_ALL_SOURCE -DIBM_RS -DUNIX"
            CPPLIBS="-L$TMW_ROOT/extern/lib/$Arch $OTHER_LIBS -lmatpp -lmmfile -lmcc -lmatlb -lmat -lmx -lm"
            CPPOPTIMFLAGS='-O -DNDEBUG'
            CPPDEBUGFLAGS='-g'
#
            FC='f77'
            FFLAGS="-I$TMW_ROOT/extern/include"
            FLIBS="-L$TMW_ROOT/extern/lib/$Arch $OTHER_LIBS -lmmfile -lmcc -lmatlb -lmat -lmx -lm"
            FOPTIMFLAGS='-O'
            FDEBUGFLAGS='-g'
#
            LD="$COMPILER"
            LDFLAGS=''
            LDOPTIMFLAGS=''
            LDDEBUGFLAGS=''
            
            SHLLD="$LD"
            SHLLDFLAGS="-bE:$EXPFILE -bM:SRE -bnoentry"
            SHLMAKEDEF="cat"
#----------------------------------------------------------------------------
            ;;
        lnx86)
#----------------------------------------------------------------------------
#
# Default to libc5 based development (ie. RedHat4.2)
#
            CC='gcc'
            CPP='g++'
            if [ -f /etc/redhat-release ]; then
		OS=`cat /etc/redhat-release`
		version=`expr "$OS" : '.*\([0-9][0-9]*\)\.'`
#
# Use these compilers for RedHat5.* systems
#
		if [ "$version" = "5" ]; then
		    CC='i486-linuxlibc5-gcc'
		    CPP='i486-linuxlibc5-g++'
		fi
	    elif [ -f /etc/debian_version ]; then
	        OS=`cat /etc/debian_version`
		version=`expr "$OS" : '.*\([0-9][0-9]*\)\.'`
#
# Use these compilers for Debian 2.* systems
#
		if [ "$version" = "2" ]; then
		    CC='i486-linuxlibc1-gcc'
		    CPP='i486-linuxlibc1-g++'
		fi
	    fi
            CFLAGS="-ansi -I$TMW_ROOT/extern/include"
            CLIBS="-L$TMW_ROOT/extern/lib/$Arch $OTHER_LIBS -lmmfile -lmcc -lmatlb -lmat -lmx -lm"
            COPTIMFLAGS='-O -DNDEBUG'
            CDEBUGFLAGS='-g'

            SHLCFLAGS="$CFLAGS"

# Linux C++ compiler flags
            CPPFLAGS="-I$TMW_ROOT/extern/include/cpp -I$TMW_ROOT/extern/include -DLNX86 -DGCC -DUNIX -x c++"
            CPPLIBS="-L$TMW_ROOT/extern/lib/$Arch $OTHER_LIBS -lmatpp -lmmfile -lmcc -lmatlb -lmat -lmx -lm"
            CPPOPTIMFLAGS='-O -DNDEBUG'
            CPPDEBUGFLAGS='-g'
#
# These flags use f2c and gcc for building FORTRAN applications.
# The fort77 script invokes the f2c command transparently,
# so it can be used like a real FORTRAN compiler.
#
            FC='fort77'
            FFLAGS="-I$TMW_ROOT/extern/include"
#
# Use these flags for the Absoft F77 FORTRAN Compiler
#
        #   FC='f77'
        #   FFLAGS="-I$TMW_ROOT/extern/include -f -N9 -N70"

            FLIBS="-L$TMW_ROOT/extern/lib/$Arch $OTHER_LIBS -lmmfile -lmcc -lmatlb -lmat -lmx -lm"
            FOPTIMFLAGS='-O'
            FDEBUGFLAGS='-g'
#
            LD="$COMPILER"
            LDFLAGS=''
            LDOPTIMFLAGS=''
            LDDEBUGFLAGS=''

            SHLLD="$LD"
            SHLLDFLAGS='-shared'
            SHLMAKEDEF="cat"
#----------------------------------------------------------------------------
            ;;
        sgi)
#----------------------------------------------------------------------------
# SGI C compiler flags
            CC='cc'
            CFLAGS="-I$TMW_ROOT/extern/include -o32"
            CLIBS="-L$TMW_ROOT/extern/lib/$Arch $OTHER_LIBS -lmmfile -lmcc -lmatlb -lmat -lmx -lm"
            COPTIMFLAGS='-O -DNDEBUG'
            CDEBUGFLAGS='-g'

            SHLCFLAGS="$CFLAGS"

# SGI C++ compiler flags
	    CPP=CC
            CPPFLAGS="+d -I$TMW_ROOT/extern/include/cpp -I$TMW_ROOT/extern/include -DSGI -DUNIX -mips2 -32 -woff 3262 -exceptions"
            CPPLIBS="-L$TMW_ROOT/extern/lib/$Arch $OTHER_LIBS -lmatpp -lmmfile -lmcc -lmatlb -lmat -lmx -lm"
            CPPOPTIMFLAGS='-O1 Olimit 7000 -DNDEBUG'
            CPPDEBUGFLAGS='-g'
            CPPLDFLAGS='-Wl,-woff,85 -mips2 -exceptions'
#
            FC='f77'
            FFLAGS="-I$TMW_ROOT/extern/include -o32"
            FLIBS="-L$TMW_ROOT/extern/lib/$Arch $OTHER_LIBS -lmmfile -lmcc -lmatlb -lmat -lmx -lm"
            FOPTIMFLAGS='-O'
            FDEBUGFLAGS='-g'
#
            LD="$COMPILER"
            LDFLAGS='-o32'
            LDOPTIMFLAGS=''
            LDDEBUGFLAGS=''

            SHLLD='ld'
            SHLLDFLAGS="-o32 -shared -exports_file $EXPFILE"
            SHLMAKEDEF="cat"
#----------------------------------------------------------------------------
            ;;
        sgi64)
# R8000 only: The default action of mbuild is to generate full MIPS IV
#             (R8000) instruction set.
#----------------------------------------------------------------------------
# SGI64 C compiler flags
            CC='cc'
            CFLAGS="-I$TMW_ROOT/extern/include -ansi -mips4 -64"
            CLIBS="-L$TMW_ROOT/extern/lib/$Arch $OTHER_LIBS -lmmfile -lmcc -lmatlb -lmat -lmx -lm"
            COPTIMFLAGS='-O -DNDEBUG'
            CDEBUGFLAGS='-g'

            SHLCFLAGS="$CFLAGS"
# SGI64 C++ compiler flags
            CPP='CC'
            CPPFLAGS="-I$TMW_ROOT/extern/include/cpp -I$TMW_ROOT/extern/include -DSGI64 -mips4 -64 -exceptions"
            CPPLIBS="-L$TMW_ROOT/extern/lib/$Arch $OTHER_LIBS -lmatpp -lmmfile -lmcc -lmatlb -lmat -lmx -lm"
            CPPOPTIMFLAGS='-O -DNDEBUG'
            CPPDEBUGFLAGS='-g'
            CPPLDFLAGS='-Wl,-woff,85 -mips4 -64'
#
            FC='f77'
            FFLAGS="-I$TMW_ROOT/extern/include"
            FLIBS="-L$TMW_ROOT/extern/lib/$Arch $OTHER_LIBS -lmmfile -lmcc -lmatlb -lmat -lmx -lm"
            FOPTIMFLAGS='-O'
            FDEBUGFLAGS='-g'
#
            LD="$COMPILER"
            LDFLAGS='-ansi -mips4 -64'
            LDOPTIMFLAGS=''
            LDDEBUGFLAGS=''

            SHLLD='ld'
            SHLLDFLAGS="-64 -mips4 -shared -exports_file $EXPFILE"
            SHLMAKEDEF="cat"
#----------------------------------------------------------------------------
            ;;
        sol2)
#----------------------------------------------------------------------------
# Solaris C compiler flags
            CC='gcc'
            CFLAGS="-ansi -I/opt/matlab5/extern/include"
            CLIBS="-L/opt/matlab5/extern/lib/sol2 $OTHER_LIBS -lmmfile -lmcc -lmatlb -lmat -lmx -lm"
            COPTIMFLAGS='-O -DNDEBUG'
            CDEBUGFLAGS='-g'
            
            SHLCFLAGS="-dalign -KPIC -I/opt/matlab5/extern/include"
            
# Solaris C++ compiler flags
            CPP="CC"
	    CCV=`CC -V 2>&1`
	    version=`expr "$CCV" : '.*\([0-9][0-9]*\)\.'`
	    if [ "$version" = "5" ]; then
		    CPP='CC -compat'
	    fi
            CPPFLAGS="-I/opt/matlab5/extern/include/cpp -I/opt/matlab5/extern/include -DSOL2 -DUNIX -DX11 +d"
            CPPLIBS="-L/opt/matlab5/extern/lib/sol2 $OTHER_LIBS -lmatpp -lmmfile -lmcc -lmatlb -lmat -lmx -lnsl -lm"
            CPPOPTIMFLAGS='-DNDEBUG -O3'
            CPPDEBUGFLAGS='-g'
#
            FC='f77'
            FFLAGS="-I/opt/matlab5/extern/include"
            FLIBS="-L/opt/matlab5/extern/lib/$Arch $OTHER_LIBS -lmmfile -lmcc -lmatlb -lmat -lmx -lm"
            FOPTIMFLAGS='-O'
            FDEBUGFLAGS='-g'
#
            LD="$COMPILER"
            LDFLAGS=''
            LDOPTIMFLAGS=''
            LDDEBUGFLAGS=''

            SHLLD=$LD
            SHLLDFLAGS="-G -M $EXPFILE"
            SHLMAKEDEF="awk 'BEGIN {printf \"{\\n\\tglobal:\\n\";} 
                                   {printf \"\\t\\t%s;\\n\", \$0;} 
                             END   {printf \"\\tlocal:\\n\\t\\t*;\\n};\\n\";}'"
#----------------------------------------------------------------------------
            ;;
    esac
#############################################################################
#
# Architecture independent lines:
#
#     Set and uncomment any lines which will apply to all architectures.
#
#----------------------------------------------------------------------------
#           CC="$CC"
#           CFLAGS="$CFLAGS"
#           COPTIMFLAGS="$COPTIMFLAGS"
#           CDEBUGFLAGS="$CDEBUGFLAGS"
#           CLIBS="$CLIBS"
#
#           FC="$FC"
#           FFLAGS="$FFLAGS"
#           FOPTIMFLAGS="$FOPTIMFLAGS"
#           FDEBUGFLAGS="$FDEBUGFLAGS"
#           FLIBS="$FLIBS"
#
#           LD="$LD"
#           LDFLAGS="$LDFLAGS"
#           LDOPTIMFLAGS="$LDOPTIMFLAGS"
#           LDDEBUGFLAGS="$LDDEBUGFLAGS"
#----------------------------------------------------------------------------
#############################################################################
