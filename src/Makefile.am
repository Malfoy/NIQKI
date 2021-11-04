bin_PROGRAMS = niqki

common_SOURCES = \
  common.h \
  zstr.hpp \
  strict_fstream.hpp \
  zstr_make_unique_polyfill.h \
  genome.h genome.cpp \
  niqki_index.h niqki_index.cpp

niqki_SOURCES = \
  $(common_SOURCES) \
  optionparser.h \
  niqki.cpp

MOSTLYCLEANFILES=*~
CLEANFILES=*~ $(noinst_SCRIPTS)
DISTCLEANFILES=*~
MAINTAINERCLEANFILES=*~

