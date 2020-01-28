bin_PROGRAMS = nthits

nthits_CPPFLAGS = -I$(top_srcdir)/Common

nthits_SOURCES = \
        Common/Uncompress.cpp \
        Common/Uncompress.h \
        Common/SignalHandler.cpp \
        Common/SignalHandler.h \
        Common/StringUtil.h \
        Common/Fcontrol.cpp \
        Common/Fcontrol.h \
        Common/nthash.hpp \
        Common/ntHashIterator.hpp \
        Common/ntcard.hpp \
        Common/CBFilter.hpp \
        nthits.cpp

dist_doc_DATA = \
	ChangeLog \
	LICENSE \
	README.md

EXTRA_DIST=autogen.sh


# Check the C++ source code for white-space errors with clang-format.
clang-format:
	for i in nthits.cpp Common/ntcard.hpp; do clang-format -style=file $$i >$$i.fixed; done
	for i in nthits.cpp Common/ntcard.hpp; do diff -su $$i $$i.fixed && rm -f $$i.fixed; done
	if ls *.fixed; then exit 1; fi
