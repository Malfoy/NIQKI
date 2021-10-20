bin_PROGRAMS = nihm

common_SOURCES = \
  common.h \
  zstr.hpp \
  strict_fstream.hpp \
  zstr_make_unique_polyfill.h \
  genome.h genome.cpp \
  nihm_index.h nihm_index.cpp

nihm_SOURCES = \
  $(common_SOURCES) \
  optionparser.h \
  nihm.cpp

MOSTLYCLEANFILES=*~
CLEANFILES=*~ $(noinst_SCRIPTS)
DISTCLEANFILES=*~
MAINTAINERCLEANFILES=*~

