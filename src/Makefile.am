CFLAGS = -g -O2 ${OPENMP_CFLAGS}

bin_PROGRAMS = solver0 solver2 solver3 solver4
check_PROGRAMS = testsolver testmatrixutil testgmres

solver0_SOURCES = main.c
solver0_LDADD = libhvs0-0.1.la
solver0_CFLAGS = -DNMOMENTS=0

solver2_SOURCES = main.c
solver2_LDADD = libhvs2-0.1.la
solver2_CFLAGS = -DNMOMENTS=3

solver3_SOURCES = main.c
solver3_LDADD = libhvs3-0.1.la
solver3_CFLAGS = -DNMOMENTS=3

solver4_SOURCES = main.c
solver4_LDADD = libhvs4-0.1.la
solver4_CFLAGS = -DNMOMENTS=4

testsolver_SOURCES = test/test_solver.c
testsolver_LDADD = libhvs-0.1.la
testmatrixutil_SOURCES = test/test_matrixutil.c
testmatrixutil_LDADD = libhvs-0.1.la
testgmres_SOURCES = test/test_gmres.c
testgmres_LDADD = libhvs-0.1.la

lib_LTLIBRARIES = libhvs0-0.1.la libhvs2-0.1.la libhvs3-0.1.la libhvs4-0.1.la 

libhvs0_0_1_la_SOURCES = libhvs/hevosol.c libhvs/fileutil.c libhvs/errorutil.c libhvs/hermiteutil.c libhvs/factorialutil.c libhvs/gmres.c libhvs/vectorutil.c libhvs/matrixutil.c
libhvs0_0_1_la_LDFLAGS = -version-info 1:0:1 -fopenmp
libhvs0_0_1_la_CFLAGS = -DNMOMENTS=0

libhvs2_0_1_la_SOURCES = libhvs/hevosol.c libhvs/fileutil.c libhvs/errorutil.c libhvs/hermiteutil.c libhvs/factorialutil.c libhvs/gmres.c libhvs/vectorutil.c libhvs/matrixutil.c
libhvs2_0_1_la_LDFLAGS = -version-info 1:0:1 -fopenmp
libhvs2_0_1_la_CFLAGS = -DNMOMENTS=2

libhvs3_0_1_la_SOURCES = libhvs/hevosol.c libhvs/fileutil.c libhvs/errorutil.c libhvs/hermiteutil.c libhvs/factorialutil.c libhvs/gmres.c libhvs/vectorutil.c libhvs/matrixutil.c
libhvs3_0_1_la_LDFLAGS = -version-info 1:0:1 -fopenmp
libhvs3_0_1_la_CFLAGS = -DNMOMENTS=3

libhvs4_0_1_la_SOURCES = libhvs/hevosol.c libhvs/fileutil.c libhvs/errorutil.c libhvs/hermiteutil.c libhvs/factorialutil.c libhvs/gmres.c libhvs/vectorutil.c libhvs/matrixutil.c
libhvs4_0_1_la_LDFLAGS = -version-info 1:0:1 -fopenmp
libhvs4_0_1_la_CFLAGS = -DNMOMENTS=4
