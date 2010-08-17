#include <stdio.h>
#include "libhvs/hevosol.h"

int main() {
	hvs_grid* grid = init_solver(100,100);
	if (!run_solver()) {
		char string[100];
		format_output(grid,100,string);
		printf("%s", string);
	}
	free_solver(grid);
	return 0;
}