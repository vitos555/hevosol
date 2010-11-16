#include <stdio.h>
#include <errno.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <time.h>
#include "libhvs/hevosol.h"
#include "libhvs/fileutil.h"
#include "libhvs/hermiteutil.h"

void usage(FILE *fh,const char **argv);

int main(int argc, char **argv) {
	hvs_state* state;
	hvs_params params;
	int status,i,c,index,quiet=0,reqargs=0;
	char *outfile=NULL;
	float t;
	time_t starttime,endtime;
	memset(&params,0,sizeof(hvs_params));
	// Default values of the parameters
	params.nu=(FLOAT_TYPE)0.1;
	params.lambda0=(FLOAT_TYPE)1.0;
	params.timestep=(FLOAT_TYPE)0.01;
	params.t0=(FLOAT_TYPE)0.0;
	params.t1=(FLOAT_TYPE)1.0;
	params.initvortfile = NULL;
	params.initcentersfile = NULL;
	params.initmomentsfile = NULL;

	opterr = 0; 
	while ((c = getopt (argc, argv, "hqv::c::t:l:n:b:e:m::o:")) != -1)
		switch (c)
		{
		case 'v':
			params.initvortfile = optarg;
			reqargs = 1;
			break;
		case 'm':
			params.initmomentsfile = optarg;
			break;
		case 'o':
			outfile = optarg;
			break;
		case 't':
			sscanf(optarg,"%f",&t);
			params.timestep=t;
			break;
		case 'l':
			sscanf(optarg,"%f",&t);
			params.lambda0=t;
			break;
		case 'n':
			sscanf(optarg,"%f",&t);
			params.nu=t;
			break;
		case 'b':
			sscanf(optarg,"%f",&t);
			params.t0=t;
			break;
		case 'e':
			sscanf(optarg,"%f",&t);
			params.t1=t;
			break;
		case 'h':
			usage(stdout,argv);
			return 0;
		case 'q':
			quiet = 1;
			break;
		case '?':
		default:
			usage(stderr,argv);
			return 1;
	}
	if (!params.timestep) {
		fprintf(stderr,"Timestep is missing.\n");
		usage(stderr,argv);
		return 1;
	} else if (!params.t1) {
		fprintf(stderr,"t1 is missing.\n");
		usage(stderr,argv);
		return 1;
	} else if (!params.lambda0) {
		fprintf(stderr,"Initial lambda is missing.\n");
		usage(stderr,argv);
		return 1;
	} else if (!params.nu) {
		fprintf(stderr,"Nu is missing.\n");
		usage(stderr,argv);
		return 1;
	} else if (!outfile) {
		fprintf(stderr,"Output file name is missing.\n");
		usage(stderr,argv);
		return 1;
	} else if ((params.initvortfile==NULL)&&(params.initmomentsfile==NULL)) {
		fprintf(stderr,"Vorticity data file is missing.\n");
		usage(stderr,argv);
		return 1;
	}

	// Get the time when the solver starts working
	starttime = time(NULL);
	
	write_params(&params,outfile);
	
	// Initialize solver
	if ((status = init_solver(&params, &state)) != HVS_OK) {
		hvserror(status, "Init error");
		return 0;
	}
	
	// Integrate
	update_vorticity_field(state,&params);
	append_vorticity(state,outfile);
	for (i=0;i<10;i++) {
		params.t0 = 0.1*i;
		params.t1 = 0.1*(i+1);
		if ((status = run_solver(&params, state)) == HVS_OK) {
			// Write vorticity to output file
			append_vorticity(state,outfile);
		} else {
			hvserror(status,"Run error");
			break;
		}
	}
	
	// Deinitialize solver
	free_solver(&state);
	endtime = time(NULL);
	if(!quiet)
		printf("Total computation time: %d sec\n",(int)(endtime-starttime));
	return 0;
}

void usage(FILE *fh,const char **argv) {
	fprintf(fh, "Hermitian vorticity solver.\n");
	fprintf(fh, "Usage: %s [-h|-q] -v file -o file [-t float|-l float|-n float|-b float|-e float]\n", argv[0]);
	fprintf(fh, "-h\t\tShow this message.\n");
	fprintf(fh, "-v file\t\tVorticity data file.\n");
	fprintf(fh, "-o file\t\tOutput file name.\n");
	fprintf(fh, "-t float\tIntegration time step. Default: 0.01.\n");
	fprintf(fh, "-n float\tViscosity of the fluid. Default: 0.1.\n");
	fprintf(fh, "-l float\tInitial lambda. Default: 1.0.\n");
	fprintf(fh, "-b float\tInitial time. Default: 0.0\n");
	fprintf(fh, "-e float\tEnd time. Default: 1.0.\n");
	fprintf(fh, "-q\t\tQuiet mode. Output only to the output file.\n");
}