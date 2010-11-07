bin_PROGRAMS = solver
check_PROGRAMS = testsolver testmatrixutil testgmres

solver_SOURCES = main.c
solver_LDADD = libhvs-0.1.la

testsolver_SOURCES = test/test_solver.c
testsolver_LDADD = libhvs-0.1.la
testmatrixutil_SOURCES = test/test_matrixutil.c
testmatrixutil_LDADD = libhvs-0.1.la
testgmres_SOURCES = test/test_gmres.c
testgmres_LDADD = libhvs-0.1.la


lib_LTLIBRARIES = libhvs-0.1.la

libhvs_0_1_la_SOURCES = libhvs/hevosol.c libhvs/fileutil.c libhvs/errorutil.c libhvs/hermiteutil.c libhvs/factorialutil.c libhvs/gmres.c libhvs/vectorutil.c libhvs/matrixutil.c
libhvs_0_1_la_LDFLAGS = -version-info 1:0:1