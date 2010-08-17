//
// Copyright (C) 2010, Vitalii Ostrovskyi <vitalii@ostrovskyi.org.ua>
// Author: Vitalii Ostrovskyi <vitalii@ostrovskyi.org.ua>
//

#include <stdlib.h>
#include <stdio.h>
#include "hevosol.h"

hvs_grid* init_solver(UINT gridwidth, UINT gridheight) {
	hvs_grid *grid = (hvs_grid *) malloc(sizeof(hvs_grid));
	if (grid == NULL) 
		return NULL;
	// Init grid
	grid->width = gridwidth;
	grid->height = gridheight;
	grid->array = (hvs_grid_point *) malloc(gridwidth*gridheight*sizeof(hvs_grid_point));
	if (grid->array == NULL) {
		free(grid);
		return NULL;
	}
	return grid;
}

void free_solver(hvs_grid *grid) {
	free(grid->array);
	free(grid);
}

int run_solver() {
	return 0;
}

int format_output(hvs_grid *grid, size_t size, char* output) {
	return snprintf(output, 100, "Width: %u, height: %u.\n", grid->width, grid->height);
}