bin_PROGRAMS = nthits

nthits_CPPFLAGS = -I$(top_srcdir)/Common -I$(top_srcdir)/btllib/include

nthits_LDADD = -L$(top_srcdir)/btllib/lib -lbtllib

nthits_SOURCES = \
        Common/StringUtil.h \
        Common/ntcard.hpp \
        Common/CBFilter.hpp \
        btllib/include/btllib/nthash.hpp \
        btllib/include/btllib/status.hpp \
        btllib/include/btllib/bloom_filter.hpp \
        btllib/include/btllib/counting_bloom_filter.hpp \
        btllib/include/btllib/seq_reader.hpp \
        btllib/include/external/cpptoml.hpp \
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
