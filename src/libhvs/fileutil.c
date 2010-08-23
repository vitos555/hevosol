#include "fileutil.h"
#include <stdlib.h>

int initfile(const char *filename, hvs_file **file) {
	hvs_file *ifile = (hvs_file *)malloc(sizeof(hvs_file));
	FILE *fh;
	if ((fh = fopen(filename,"r")) == NULL) {
		return HVS_ERR;
	}
	ifile->fh = fh;
	*file = ifile;
	return HVS_OK;
}

void closefile(hvs_file **file) {
	fclose((*file)->fh);
	free(*file);
	*file = NULL;
}

ssize_t readdata(const hvs_file *file, size_t count, hvs_position *pos, hvs_vorticity *vort) {
	int readcount = 0;
	int currread = 0;
	while (!feof(file->fh) && readcount<count) {
		if ((currread = fscanf(file->fh, "%f\t%f\t%f\n", &(pos[readcount].x), &(pos[readcount].y), &(vort[readcount])))==3) {
			readcount++;
		} else if (currread == EOF) {
			continue;
		} else {
#ifdef HVS_DEBUG
			fprintf(stderr,"Currread: %i\n",currread);
#endif
			return HVS_ERR_WRONG_FILE_FORMAT;
		}
	}
	return readcount;
}

ssize_t writedata(const hvs_state *state, const char *filename) {
	int writecount = 0;
	FILE *fh;
	if ((fh = fopen(filename, "w")) == NULL) {
		return HVS_ERR;
	}
	while(writecount<state->size) {
		fprintf(fh, "%f\t%f\t%f\n", state->grid[writecount].x, state->grid[writecount].y, state->vorticity_field[writecount]);
		writecount++;
	}
	fclose(fh);
	return writecount;
}