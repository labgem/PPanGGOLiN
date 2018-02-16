#
#	--- Main  groups of targets to be updated ---
#

# V1.04-a  04-10-1997  MD  add lint
# V1.04-b  18-11-1997  MD  replace cc with gcc, add -lemu option
# V1.04-c  11-01-1998  MD  archive name as variable
# V1.05-a  25-01-1998  MD  add exemain.c -> nem_exe
# V1.05-b  25-01-1998  MD  add mexmain.c and mex target -> nem_exe.mexsol
# V1.05-c  26-01-1998  MD  add exememo.c -> nem_exe and mexmemo.c -> .mexsol
# V1.06-a  17-06-1998  MD  change nem_noi -> nem_mod
# V1.06-b  23-06-1998  MD  flag -Wall
# V1.06-c  23-06-1998  MD  add target t_nem_arg
# V1.06-d  30-06-1998  MD  add target t_nem_mod
# V1.07-a  29-02-1999  MD  archive name : new version

CFLAGS = -Wall

# Flags for MS-DOS compilation using DJGPP's gcc
#CFLAGS = -lemu -Wall

# Select C compiler here
CC = gcc
#CC = cc

# Warning messages


# Name of uncompressed archive file to create
ARC = nem107.tar

#	Optimized application executables (default make)
exe :		txt2hlp nem_exe geo2nei randord err2 tcpu

#	Matlab executable
mex :		nem_exe.mexsol

#	Debuggable executables
dbg :		g_nem_exe g_geo2nei g_randord

#	Archive
arc :		$(ARC).Z

#	Everything
all :		remobj mex exe dbg arc

#	Test version of nem_exe
test :		nem_tmp

#	Tight syntax check
check :
	lint nem_exe.c nem_arg.c nem_alg.c nem_nei.c \
	nem_mod.c nem_rnd.c nem_ver.c nem_hlp.c \
	lib_io.c exemain.c exememo.c
	lint err2.c
	lint tcpu.c

#	Remove all object files
remobj :
	\rm *.o

#
#	--- Single targets ---
#

#	Archive commands to update source file archiving
$(ARC).Z :	nem_exe.c nem_arg.c nem_alg.c nem_nei.c \
		nem_mod.c nem_rnd.c nem_ver.c nem_hlp.c \
		nem_typ.h nem_arg.h nem_alg.h nem_nei.h \
		nem_mod.h nem_rnd.h nem_ver.h nem_hlp.h \
		lib_io.c  lib_io.h genmemo.h exemain.c exememo.c \
		mexmain.c mexmemo.c mainfunc.h \
		geo2nei.c randord.c \
		err2.c tcpu.c txt2hlp.c \
		simimg.str simimg.dat simimg.nei simimg.cr \
		simimg.cf simimg.mf \
		makefile
	tar cvf $(ARC) nem_*.[ch] lib_*.[ch] exe*.c mex*.c gen*h \
	err2.c geo2nei.c randord.c txt2hlp.c tcpu.c mainfunc.h \
	simimg.str simimg.dat simimg.nei simimg.cr simimg.cf simimg.mf \
	nem_user.txt makefile
	gzip -vf $(ARC)
	\cp -p $(ARC).gz /home_f/mdang/public_html/Progs/$(ARC).gz


#	Linkage commands to update the executables
nem_exe.mexsol : nem_exe.o nem_arg.o nem_alg.o nem_nei.o nem_mod.o nem_rnd.o \
		nem_ver.o nem_hlp.o lib_io.o mexmain.c mexmemo.c genmemo.h
	cmex nem_exe.o nem_arg.o nem_alg.o nem_nei.o nem_mod.o nem_rnd.o \
	nem_ver.o nem_hlp.o lib_io.o  mexmain.c mexmemo.c \
	-o nem_exe -lm -O

nem_tmp : 	nem_exe.o nem_arg.o nem_alg.o nem_nei.o nem_mod.o nem_rnd.o \
		nem_ver.o nem_hlp.o lib_io.o exemain.o exememo.o 
	$(CC) nem_exe.o nem_arg.o nem_alg.o nem_nei.o nem_mod.o nem_rnd.o \
	nem_ver.o nem_hlp.o lib_io.o  exemain.o exememo.o \
	-o nem_tmp -lm -O $(CFLAGS)

nem_exe : 	nem_exe.o nem_arg.o nem_alg.o nem_nei.o nem_mod.o nem_rnd.o \
		nem_ver.o nem_hlp.o lib_io.o exemain.o exememo.o 
	$(CC) nem_exe.o nem_arg.o nem_alg.o nem_nei.o nem_mod.o nem_rnd.o \
	nem_ver.o nem_hlp.o lib_io.o  exemain.o exememo.o \
	-o nem_exe -lm -O $(CFLAGS)

g_nem_exe : 	g_nem_exe.o g_nem_arg.o g_nem_alg.o g_nem_nei.o \
		g_nem_mod.o g_nem_rnd.o g_nem_ver.o g_nem_hlp.o g_lib_io.o \
                g_exemain.o g_exememo.o
	$(CC) g_nem_exe.o g_nem_arg.o g_nem_alg.o g_nem_nei.o \
	g_nem_mod.o g_nem_rnd.o g_nem_ver.o g_nem_hlp.o g_lib_io.o \
        g_exemain.o g_exememo.o \
	-o g_nem_exe -lm -g $(CFLAGS)


geo2nei :	geo2nei.o lib_io.o
	$(CC) geo2nei.o lib_io.o \
	-o geo2nei  -lm  -O $(CFLAGS)

g_geo2nei :	g_geo2nei.o g_lib_io.o
	$(CC) g_geo2nei.o g_lib_io.o \
	-o g_geo2nei  -lm  -O $(CFLAGS)

randord :	randord.o nem_rnd.o
	$(CC) randord.o nem_rnd.o \
	-o randord  -lm  -O $(CFLAGS)

g_randord :	g_randord.o g_nem_rnd.o
	$(CC) g_randord.o g_nem_rnd.o \
	-o g_randord  -lm  -g $(CFLAGS)

tcpu :	tcpu.c
	$(CC) tcpu.c -o tcpu  -lm  -O $(CFLAGS)


#	Linkage commands to update the test executables

#V1.06-c
t_nem_arg :	g_nem_arg.o g_lib_io.o g_nem_hlp.o g_exememo.o g_nem_ver.o
	$(CC) g_nem_arg.o g_lib_io.o g_nem_hlp.o g_exememo.o g_nem_ver.o \
	-o $@  -lm  -g $(CFLAGS)

#V1.06-d
t_nem_mod.out :	nem_mod.c t_nem_mod.c t_nem_mod.sh
	$(CC) g_exememo.o t_nem_mod.c -o t_nem_mod -lm  -g $(CFLAGS)
	@t_nem_mod.sh > tmp.out
	@echo "\ndiff tmp.out t_nem_mod.out :\n"
	if ( diff tmp.out t_nem_mod.out ) ; \
	then \
          echo "test of nem_mod.c OK" ; touch t_nem_mod.out ; \
	else \
          echo "\07\nWarning \07: output of test of nem_mod.c changed" ; \
#	  echo "Change t_nem_mod.out ?" ;  ; \
#	  if "$a" = "y"; then \mv tmp.out t_nem_mod.out; fi ; \
	fi
	@\rm tmp.out
	@\rm t_nem_mod

t_nem_mod :	nem_mod.c t_nem_mod.c g_exememo.o 
	$(CC) g_exememo.o t_nem_mod.c \
	-o $@  -lm  -g $(CFLAGS)
#
#	--- Secondary targets (object files) ---
#

#	Compilation commands to update the object files
nem_exe.o :	nem_exe.c nem_typ.h nem_arg.h nem_alg.h nem_rnd.h lib_io.h
	$(CC) -c nem_exe.c -O $(CFLAGS)

g_nem_exe.o :	nem_exe.c nem_typ.h nem_arg.h nem_alg.h nem_rnd.h lib_io.h
	$(CC) -c nem_exe.c -g $(CFLAGS) -o g_nem_exe.o


nem_arg.o :	nem_arg.c nem_typ.h nem_arg.h
	$(CC) -c nem_arg.c -O $(CFLAGS)

g_nem_arg.o :	nem_arg.c nem_typ.h nem_arg.h
	$(CC) -c nem_arg.c -g $(CFLAGS) -o g_nem_arg.o


nem_alg.o :	nem_alg.c nem_typ.h nem_alg.h nem_nei.h nem_mod.h nem_rnd.h
	$(CC) -c nem_alg.c -O $(CFLAGS)

g_nem_alg.o :	nem_alg.c nem_typ.h nem_alg.h nem_nei.h nem_mod.h nem_rnd.h
	$(CC) -c nem_alg.c -g $(CFLAGS) -o g_nem_alg.o


nem_nei.o :	nem_nei.c nem_typ.h nem_nei.h
	$(CC) -c nem_nei.c -O $(CFLAGS)

g_nem_nei.o :	nem_nei.c nem_typ.h nem_nei.h
	$(CC) -c nem_nei.c -g $(CFLAGS) -o g_nem_nei.o


nem_mod.o :	nem_mod.c nem_typ.h nem_mod.h
	$(CC) -c nem_mod.c -O $(CFLAGS)

g_nem_mod.o :	nem_mod.c nem_typ.h nem_mod.h
	$(CC) -c nem_mod.c -g $(CFLAGS) -o g_nem_mod.o


nem_hlp.o :	nem_hlp.c nem_hlp.h
	$(CC) -c nem_hlp.c -O $(CFLAGS)

g_nem_hlp.o :	nem_hlp.c nem_hlp.h
	$(CC) -c nem_hlp.c -g $(CFLAGS) -o g_nem_hlp.o


nem_ver.o :	nem_ver.c nem_typ.h nem_ver.h
	$(CC) -c nem_ver.c -O $(CFLAGS)

g_nem_ver.o :	nem_ver.c nem_typ.h nem_ver.h
	$(CC) -c nem_ver.c -g $(CFLAGS) -o g_nem_ver.o


nem_rnd.o :	nem_rnd.c nem_rnd.h
	$(CC) -c nem_rnd.c -O $(CFLAGS)

g_nem_rnd.o :	nem_rnd.c nem_rnd.h
	$(CC) -c nem_rnd.c -g $(CFLAGS) -o g_nem_rnd.o


lib_io.o :	lib_io.c lib_io.h
	$(CC) -c lib_io.c -O $(CFLAGS)

g_lib_io.o :	lib_io.c lib_io.h
	$(CC) -c lib_io.c -g $(CFLAGS) -o g_lib_io.o

exemain.o :	exemain.c
	$(CC) -c exemain.c -O $(CFLAGS)

g_exemain.o :	exemain.c
	$(CC) -c exemain.c -g $(CFLAGS) -o g_exemain.o

exememo.o :	exememo.c genmemo.h
	$(CC) -c exememo.c -O $(CFLAGS)

g_exememo.o :	exememo.c genmemo.h
	$(CC) -c exememo.c -g $(CFLAGS) -o g_exememo.o

geo2nei.o :	geo2nei.c lib_io.h
	$(CC) -c geo2nei.c -O $(CFLAGS)

g_geo2nei.o :	geo2nei.c lib_io.h
	$(CC) -c geo2nei.c -g $(CFLAGS) -o g_geo2nei.o


randord.o :	randord.c nem_rnd.h
	$(CC) -c randord.c -O $(CFLAGS)

g_randord.o :	randord.c nem_rnd.h
	$(CC) -c randord.c -g $(CFLAGS) -o g_randord.o

#
#	--- Auxiliary targets (automated source files) ---
#

nem_hlp.c :	nem_user.txt
	./txt2hlp nem_user.txt nem_hlp y \
	"Project NEM : display of long help" > /dev/null
