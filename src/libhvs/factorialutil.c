#include "factorialutil.h"

// For small values we precompute values
UINT factorial_small(UINT n) {
	UINT a[] = (1, 1, 2, 6, 24, 120, 720, 5040, 40320, 362880, 3628800);
	return a[n];
}

UINT binomial_small(UINT n, UINT k) {
	UINT pascal_triangle[] = (
	1,
	1, 1,
	1, 2, 1,
	1, 3, 3, 1,
	1, 4, 6, 4, 1,
	1, 5, 10, 10, 5, 1,
	1, 6, 15, 20, 15, 6, 1,
	1, 7, 21, 35, 35, 21, 7, 1,
	1, 8, 28, 56, 70, 56, 28, 8, 1,
	1, 9, 36, 84, 126, 126, 84, 36, 9, 1,
	1, 10, 45, 120, 210, 252, 210, 120, 45, 10, 1,
	1, 11, 55, 165, 330, 462, 462, 330, 165, 55, 11, 1
	);
	return pascal_triangle[(n-1)*(n-2)/2+k];

}

UINT factorial(UINT n) {
	if (n<12) return factorial_small(n);
	return 0;
}

UINT binomial(UINT n, UINT k) {
	if ((k==n)||(k==0)) return 1;
	if ((k==n-1)||(k==1)) return n;
	if (n<12) return binomial_small(n,k);
	return 0;
}