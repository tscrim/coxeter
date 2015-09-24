ACLOCAL_AMFLAGS = -I m4

coxeter_LDADD = libcoxeter.la

lib_LTLIBRARIES = libcoxeter.la
libcoxeter_la_SOURCES = affine.cpp automata.cpp bits.cpp cells.cpp \
                        commands.cpp constants.cpp coxgroup.cpp coxtypes.cpp \
                        error.cpp fcoxgroup.cpp files.cpp general.cpp \
                        graph.cpp help.cpp interactive.cpp interface.cpp \
                        invkl.cpp io.cpp kl.cpp klsupport.cpp main.cpp \
                        memory.cpp minroots.cpp posets.cpp schubert.cpp \
                        special.cpp transducer.cpp typeA.cpp type.cpp \
                        uneqkl.cpp wgraph.cpp
libcoxeter_la_LDFLAGS = $(LIBCOXETER_LDFLAGS) -version-info 3:1:0

coxeterinclude_HEADERS = *.h


