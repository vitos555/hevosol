#include <stdlib.h>
#include "fileutil.h"

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

ssize_t read_vorticity(const hvs_file *file, size_t count, hvs_position *pos, hvs_vorticity *vort) {
	int readcount = 0;
	int currread = 0;
	while (!feof(file->fh) && readcount<count) {
#if HVS_FLOAT_TYPE==HVS_LONG_DOUBLE
		if ((currread = fscanf(file->fh, "%Lf\t%Lf\t%Lf\n", &(pos[readcount].x), &(pos[readcount].y), &(vort[readcount])))==3) {
#else
		if ((currread = fscanf(file->fh, "%f\t%f\t%f\n", &(pos[readcount].x), &(pos[readcount].y), &(vort[readcount])))==3) {
#endif
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

ssize_t write_vorticity(const hvs_state *state, const char *filename) {
	int writecount = 0;
	FILE *fh;
	if ((fh = fopen(filename, "w+")) == NULL) {
		return HVS_ERR;
	}
	fprintf(fh,"=Time=\n");
#if HVS_FLOAT_TYPE==HVS_LONG_DOUBLE
	fprintf(fh,"%Lf",state->curtime);
#else
	fprintf(fh,"%f",state->curtime);
#endif
	fprintf(fh,"=Vorticity x	y	vorticity=\n");
	while(writecount<state->size) {
#if HVS_FLOAT_TYPE==HVS_LONG_DOUBLE
		fprintf(fh, "%Lf\t%Lf\t%.18Lf\n", state->grid[writecount].x, state->grid[writecount].y, state->vorticity_field[writecount]);
#else
		fprintf(fh, "%f\t%f\t%.12f\n", state->grid[writecount].x, state->grid[writecount].y, state->vorticity_field[writecount]);
#endif
		writecount++;
	}
	fclose(fh);
	return writecount;
}

ssize_t append_vorticity(const hvs_state *state, const char *filename) {
	int writecount = 0;
	FILE *fh;
	if ((fh = fopen(filename, "a+")) == NULL) {
		return HVS_ERR;
	}
	fprintf(fh,"=Time=\n");
#if HVS_FLOAT_TYPE==HVS_LONG_DOUBLE
	fprintf(fh,"%Lf\n",state->curtime);
#else
	fprintf(fh,"%f\n",state->curtime);
#endif
	fprintf(fh,"=Vorticity x	y	vorticity=\n");
	while(writecount<state->size) {
#if HVS_FLOAT_TYPE==HVS_LONG_DOUBLE
		fprintf(fh, "%Lf\t%Lf\t%.18Lf\n", state->grid[writecount].x, state->grid[writecount].y, state->vorticity_field[writecount]);
#else
		fprintf(fh, "%f\t%f\t%.12f\n", state->grid[writecount].x, state->grid[writecount].y, state->vorticity_field[writecount]);
#endif
		writecount++;
	}
	fclose(fh);
	return writecount;
}

ssize_t write_params(const hvs_params *params, const char *filename) {
	int writecount = 0;
	FILE *fh;
	if ((fh = fopen(filename, "w+")) == NULL) {
		return HVS_ERR;
	}
	fprintf(fh,"=Params=\n");
#if HVS_FLOAT_TYPE==HVS_LONG_DOUBLE
	fprintf(fh,"lambda=%Lf,nu=%Lf\n",params->lambda0,params->nu);
#else
	fprintf(fh,"lambda=%f,nu=%f\n",params->lambda0,params->nu);
#endif
	writecount+=2;
	fclose(fh);
	return writecount;
}
ssize_t append_params(const hvs_params *params, const char *filename) {
	int writecount = 0;
	FILE *fh;
	if ((fh = fopen(filename, "a+")) == NULL) {
		return HVS_ERR;
	}
	fprintf(fh,"=Params=\n");
#if HVS_FLOAT_TYPE==HVS_LONG_DOUBLE
	fprintf(fh,"lambda=%Lf,nu=%Lf\n",params->lambda0,params->nu);
#else
	fprintf(fh,"lambda=%f,nu=%f\n",params->lambda0,params->nu);
#endif
	writecount+=2;
	fclose(fh);
	return writecount;
}

ssize_t write_moments(const hvs_state *state, const char *filename) {
	int writecount = 0;
	int i0;
	FILE *fh;
	if ((fh = fopen(filename, "w+")) == NULL) {
		return HVS_ERR;
	}
	fprintf(fh,"=Moments (0,0)	(2,0)	(1,1)	(0,2) =\n");
	for(i0=0;i0<state->ncenters;i0++) {
#if HVS_FLOAT_TYPE==HVS_LONG_DOUBLE
	fprintf(fh,"%Lf\t%Lf\t%Lf\t%Lf\n",
		state->moments[i0][MOM_INDEX(0,0)],
		state->moments[i0][MOM_INDEX(2,0)],
		state->moments[i0][MOM_INDEX(1,1)],
		state->moments[i0][MOM_INDEX(0,2)]);
#else
	fprintf(fh,"%Lf\t%Lf\t%Lf\t%Lf\n",
		state->moments[i0][MOM_INDEX(0,0)],
		state->moments[i0][MOM_INDEX(2,0)],
		state->moments[i0][MOM_INDEX(1,1)],
		state->moments[i0][MOM_INDEX(0,2)]);
#endif
	writecount+=1;
	}
	fclose(fh);
	return writecount;
}
ssize_t append_moments(const hvs_state *state, const char *filename) {
	int writecount = 0;
	int i0;
	FILE *fh;
	if ((fh = fopen(filename, "a+")) == NULL) {
		return HVS_ERR;
	}
	fprintf(fh,"=Moments (0,0)	(2,0)	(1,1)	(0,2) =\n");
	for(i0=0;i0<state->ncenters;i0++) {
#if HVS_FLOAT_TYPE==HVS_LONG_DOUBLE
	fprintf(fh,"%Lf\t%Lf\t%Lf\t%Lf\n",
		state->moments[i0][MOM_INDEX(0,0)],
		state->moments[i0][MOM_INDEX(2,0)],
		state->moments[i0][MOM_INDEX(1,1)],
		state->moments[i0][MOM_INDEX(0,2)]);
#else
	fprintf(fh,"%f\t%f\t%f\t%f\n",
		state->moments[i0][MOM_INDEX(0,0)],
		state->moments[i0][MOM_INDEX(2,0)],
		state->moments[i0][MOM_INDEX(1,1)],
		state->moments[i0][MOM_INDEX(0,2)]);
#endif
	writecount+=1;
	}
	fclose(fh);
	return writecount;
}

ssize_t write_centers(const hvs_state *state, const char *filename) {
	int writecount = 0;
	int i0;
	FILE *fh;
	if ((fh = fopen(filename, "w+")) == NULL) {
		return HVS_ERR;
	}
	fprintf(fh,"=Centers x	y=\n");
	for(i0=0;i0<state->ncenters;i0++) {
#if HVS_FLOAT_TYPE==HVS_LONG_DOUBLE
	fprintf(fh,"%Lf\t%Lf\n",
		state->centers[i0].x,
		state->centers[i0].y);
#else
	fprintf(fh,"%Lf\t%Lf\n",
		state->centers[i0].x,
		state->centers[i0].y);
#endif
	writecount+=1;
	}
	fclose(fh);
	return writecount;
}
ssize_t append_centers(const hvs_state *state, const char *filename) {
	int writecount = 0;
	int i0;
	FILE *fh;
	if ((fh = fopen(filename, "a+")) == NULL) {
		return HVS_ERR;
	}

	fprintf(fh,"=Centers x	y=\n");
	for(i0=0;i0<state->ncenters;i0++) {
#if HVS_FLOAT_TYPE==HVS_LONG_DOUBLE
	fprintf(fh,"%Lf\t%Lf\n",
		state->centers[i0].x,
		state->centers[i0].y);
#else
	fprintf(fh,"%Lf\t%Lf\n",
		state->centers[i0].x,
		state->centers[i0].y);
#endif
	writecount+=1;
	}
	fclose(fh);
	return writecount;
}
