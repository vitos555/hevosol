bin_PROGRAMS = solver test

solver_SOURCES = main.c
solver_LDADD = libhvs-0.1.la

test_SOURCES = test.c
test_LDADD = libhvs-0.1.la

lib_LTLIBRARIES = libhvs-0.1.la

libhvs_0_1_la_SOURCES = libhvs/hevosol.c libhvs/fileutil.c libhvs/errorutil.c libhvs/hermiteutil.c libhvs/factorialutil.c libhvs/gmres.c libhvs/vectorutil.c
libhvs_0_1_la_LDFLAGS = -version-info 1:0:1