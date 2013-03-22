#include <stdlib.h>
#include "fileutil.h"

int initfile(const char *filename, hvs_file **file) {
	hvs_file *ifile = (hvs_file *)malloc(sizeof(hvs_file));
	FILE *fh;
	if ((fh = fopen(filename,"r")) == NULL) {
		*file = NULL;
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
	FLOAT_TYPE x,y,tvort;
	while (!feof(file->fh) && readcount<count) {
#if HVS_FLOAT_TYPE==HVS_LONG_DOUBLE
		if ((currread = fscanf(file->fh, "%Lf\t%Lf\t%Lf\n", &x, &y, &tvort))==3) {
#elif HVS_FLOAT_TYPE==HVS_FLOAT
		if ((currread = fscanf(file->fh, "%f\t%f\t%f\n", &x, &y, &tvort))==3) {
#else
		if ((currread = fscanf(file->fh, "%lf\t%lf\t%lf\n", &x, &y, &tvort))==3) {
#endif
			pos[readcount].x=(FLOAT_TYPE)x;
			pos[readcount].y=(FLOAT_TYPE)y;
			vort[readcount]=(FLOAT_TYPE)tvort;
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

ssize_t read_centers(const hvs_file *file, size_t count, hvs_centers pos) {
	int readcount = 0;
	int currread = 0;
	FLOAT_TYPE x,y;
	while (!feof(file->fh) && readcount<count) {
#if HVS_FLOAT_TYPE==HVS_LONG_DOUBLE
		if ((currread = fscanf(file->fh, "%Lf\t%Lf\n", &x, &y))==2) {
#elif HVS_FLOAT_TYPE==HVS_FLOAT
		if ((currread = fscanf(file->fh, "%f\t%f\t\n", &x, &y))==2) {
#else
		if ((currread = fscanf(file->fh, "%lf\t%lf\t\n", &x, &y))==2) {
#endif
			pos[readcount].x=x;
			pos[readcount].y=y;
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

ssize_t read_moments(const hvs_file *file, size_t count, hvs_center *pos, hvs_moment *moment) {
	int readcount = 0;
	int currread = 0;
	FLOAT_TYPE x,y,m;
	int i;
	while (!feof(file->fh) && readcount<count) {
#if HVS_FLOAT_TYPE==HVS_LONG_DOUBLE
		if ((currread = fscanf(file->fh, "%Lf\t%Lf", &x, &y))==2) {
#elif HVS_FLOAT_TYPE==HVS_FLOAT
		if ((currread = fscanf(file->fh, "%f\t%f", &x, &y))==2) {
#else
		if ((currread = fscanf(file->fh, "%lf\t%lf", &x, &y))==2) {
#endif
			pos[readcount].x=x;
			pos[readcount].y=y;
		} else if (currread == EOF) {
			continue;
		} else {
#ifdef HVS_DEBUG
			fprintf(stderr,"Currread: %i\n",currread);
#endif
			return HVS_ERR_WRONG_FILE_FORMAT;
		}
		for(i=0;i<NCOMBS;i++) {
			if ((currread = fscanf(file->fh,"\t"))==1) {
			} else if (currread == EOF) {
				continue;
			} else {
#ifdef HVS_DEBUG
				fprintf(stderr,"Currread: %i\n",currread);
#endif
			}
#if HVS_FLOAT_TYPE==HVS_LONG_DOUBLE
			if ((currread = fscanf(file->fh, "%Lf", &m))==1) {
#elif HVS_FLOAT_TYPE==HVS_FLOAT
			if ((currread = fscanf(file->fh, "%f", &m))==1) {
#else
			if ((currread = fscanf(file->fh, "%lf", &m))==1) {
#endif
				if ((i==1)||(i==2)) m=0.0;
				moment[readcount][i]=m;
			} else if (currread == EOF) {
				continue;
			} else {
#ifdef HVS_DEBUG
				fprintf(stderr,"Currread: %i\n",currread);
#endif
			}
		}
		if ((currread = fscanf(file->fh,"\n"))==1) {
		} else if (currread == EOF) {
			continue;
		} else {
#ifdef HVS_DEBUG
			fprintf(stderr,"Currread: %i\n",currread);
#endif
		}
		readcount++;
	}
	return readcount;
}

ssize_t read_time(const hvs_file *file, FLOAT_TYPE *time) {
	int readcount = 0;
	if (!feof(file->fh)) {
#if HVS_FLOAT_TYPE==HVS_LONG_DOUBLE
		readcount = fscanf(file->fh, "%Lf\n", time);
#elif HVS_FLOAT_TYPE==HVS_FLOAT
		readcount = fscanf(file->fh, "%f\n", time);
#else
		readcount = fscanf(file->fh, "%lf\n", time);
#endif
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
	fprintf(fh,"t0=%Lf,t1=%Lf,timestep=%Lf\n",params->t0,params->t1,params->timestep);
#else
	fprintf(fh,"lambda=%f,nu=%f\n",params->lambda0,params->nu);
	fprintf(fh,"t0=%f,t1=%f,timestep=%f\n",params->t0,params->t1,params->timestep);
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
	fprintf(fh,"t0=%Lf,t1=%Lf,timestep=%Lf\n",params->t0,params->t1,params->timestep);
#else
	fprintf(fh,"lambda=%f,nu=%f\n",params->lambda0,params->nu);
	fprintf(fh,"t0=%f,t1=%f,timestep=%f\n",params->t0,params->t1,params->timestep);
#endif
	writecount+=2;
	fclose(fh);
	return writecount;
}

ssize_t write_moments(const hvs_state *state, const char *filename) {
	int writecount = 0;
	int i0,i;
	FILE *fh;
	if ((fh = fopen(filename, "w+")) == NULL) {
		return HVS_ERR;
	}
	fprintf(fh,"=Moments (0,0)	(2,0)	(1,1)	(0,2) =\n");
	for(i0=0;i0<state->ncenters;i0++) {
		for(i=0;i<NCOMBS;i++)
#if HVS_FLOAT_TYPE==HVS_LONG_DOUBLE
			fprintf(fh,"%Lf\t",state->moments[i0][i]);
#else
			fprintf(fh,"%f\t",state->moments[i0][i]);
#endif
		fprintf(fh,"\n");
		writecount+=1;
	}
	fclose(fh);
	return writecount;
}
ssize_t write_tmp_moments(const hvs_state *state, FLOAT_TYPE time, const char *filename) {
	int writecount = 0;
	int i0,i;
	FILE *fh;
	if ((fh = fopen(filename, "w+")) == NULL) {
		return HVS_ERR;
	}
#if HVS_FLOAT_TYPE==HVS_LONG_DOUBLE
	fprintf(fh,"%Lf\n",time);
#else
	fprintf(fh,"%f\n",time);
#endif
	for(i0=0;i0<state->ncenters;i0++) {
#if HVS_FLOAT_TYPE==HVS_LONG_DOUBLE
		fprintf(fh,"%Lf\t%Lf\t",
			state->centers[i0].x,
			state->centers[i0].y);
#else
		fprintf(fh,"%f\t%f\t",
			state->centers[i0].x,
			state->centers[i0].y);
#endif
		for(i=0;i<NCOMBS;i++)
#if HVS_FLOAT_TYPE==HVS_LONG_DOUBLE
			fprintf(fh,"%Lf\t",state->moments[i0][i]);
#else
			fprintf(fh,"%f\t",state->moments[i0][i]);
#endif
		fprintf(fh,"\n");
		writecount+=1;
	}
	fclose(fh);
	return writecount;
}
ssize_t append_moments(const hvs_state *state, const char *filename) {
	int writecount = 0;
	int i0,i;
	FILE *fh;
	if ((fh = fopen(filename, "a+")) == NULL) {
		return HVS_ERR;
	}
	fprintf(fh,"=Moments ");
	for(i=0;i<NCOMBS;i++)
		fprintf(fh,"(%d,%d)\t",COMBS_IND1(i),COMBS_IND2(i));
	fprintf(fh, "=\n");
	for(i0=0;i0<state->ncenters;i0++) {
		for(i=0;i<NCOMBS;i++)
#if HVS_FLOAT_TYPE==HVS_LONG_DOUBLE
			fprintf(fh,"%Lf\t",state->moments[i0][i]);
#else
			fprintf(fh,"%f\t",state->moments[i0][i]);
#endif
		fprintf(fh,"\n");
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
	fprintf(fh,"%f\t%f\n",
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
	fprintf(fh,"%f\t%f\n",
		state->centers[i0].x,
		state->centers[i0].y);
#endif
	writecount+=1;
	}
	fclose(fh);
	return writecount;
}
