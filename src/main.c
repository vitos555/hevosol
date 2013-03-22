#include <stdio.h>
#include <unistd.h>
#include <errno.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <time.h>
#include "libhvs/hevosol.h"
#include "libhvs/fileutil.h"
#include "libhvs/hermiteutil.h"

void usage(FILE *fh,char **argv);

int main(int argc, char **argv) {
	hvs_state* state;
	hvs_params params;
	int status,i,c,index,quiet=0,reqargs=0,num=0,recovery_mode=0;
	char *outfile=NULL;
	float t,t1,t2;
	time_t starttime,endtime;
	memset(&params,0,sizeof(hvs_params));
	// Default values of the parameters
	params.nu=(FLOAT_TYPE)0.1;
	params.lambda0=(FLOAT_TYPE)1.0;
	params.timestep=(FLOAT_TYPE)0.01;
	params.t0=(FLOAT_TYPE)0.0;
	params.t1=(FLOAT_TYPE)1.0;
	params.xmin=(FLOAT_TYPE)-5.0;
	params.xmax=(FLOAT_TYPE)5.0;
	params.xstep=(FLOAT_TYPE)0.1;
	params.ymin=(FLOAT_TYPE)-5.0;
	params.ymax=(FLOAT_TYPE)5.0;
	params.ystep=(FLOAT_TYPE)0.1;
	params.initvortfile = NULL;
	params.initcentersfile = NULL;
	params.initmomentsfile = NULL;
	params.tmpmomentsfile = NULL;
	params.tmpmomentstime = -1.0;

	opterr = 0; 
	while ((c = getopt(argc, argv, "hqrv:c:t:l:n:b:e:m:o:x:y:p:")) != -1)
		switch (c)
		{
		case 'v':
			params.initvortfile = optarg;
			reqargs = 1;
			break;
		case 'm':
			params.initmomentsfile = optarg;
			break;
		case 'p':
			sscanf(optarg,"%d",&num);
			break;
		case 'o':
			outfile = optarg;
			break;
		case 't':
			sscanf(optarg,"%f",&t);
			params.timestep=(FLOAT_TYPE)t;
			break;
		case 'l':
			sscanf(optarg,"%f",&t);
			params.lambda0=(FLOAT_TYPE)t;
			break;
		case 'n':
			sscanf(optarg,"%f",&t);
			params.nu=(FLOAT_TYPE)t;
			break;
		case 'b':
			sscanf(optarg,"%f",&t);
			params.t0=(FLOAT_TYPE)t;
			break;
		case 'e':
			sscanf(optarg,"%f",&t);
			params.t1=(FLOAT_TYPE)t;
			break;
		case 'x':
			sscanf(optarg,"%f:%f:%f",&t,&t1,&t2);
			params.xmin=(FLOAT_TYPE)t;
			params.xmax=(FLOAT_TYPE)t2;
			params.xstep=(FLOAT_TYPE)t1;
			break;
		case 'y':
			sscanf(optarg,"%f:%f:%f",&t,&t1,&t2);
			params.xmin=(FLOAT_TYPE)t;
			params.ymax=(FLOAT_TYPE)t2;
			params.ystep=(FLOAT_TYPE)t1;
			break;
		case 'h':
			usage(stdout,argv);
			return 0;
		case 'q':
			quiet = 1;
			break;
		// Recovery mode - store moments file after each step.
		case 'r':
			recovery_mode = 1;
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

	// Store moments after each step into tmpmomentsfile
	if (recovery_mode) {
		params.tmpmomentsfile = "hevosol.tmpmomentsfile";
	}

	// Get the time when the solver starts working
	starttime = time(NULL);
	
	write_params(&params,outfile);
	
	// Initialize solver
	if ((status = init_solver(&params, &state)) != HVS_OK) {
		hvserror(status, "Init error");
		return 1;
	}
	
	// Output initial state
	update_vorticity_field(state);
	append_centers(state,outfile);
	append_moments(state,outfile);
	append_vorticity(state,outfile);
	
	//Integrate
	t1=params.t0;
	t2=params.t1;
	if (t1>t2) {
		hvserror(0,"End time should be greater then beginning.");
		return 1;
	}
	if (num>0) {
		for(t=t1;t<t2;t+=num*params.timestep) {
			params.t0=t;
			params.t1=MIN(t+num*params.timestep,t2);
			if ((status = run_solver(&params, state)) == HVS_OK) {
				// Write info to output file
				append_centers(state,outfile);
				append_moments(state,outfile);
				append_vorticity(state,outfile);
			} else if (status == HVS_SKIP_STEP) {
				// Do nothing
			} else {
				hvserror(status,"Run error");
				return 1;
			}
		}
	} else {
		if ((status = run_solver(&params, state)) == HVS_OK) {
			// Write info to output file
			append_centers(state,outfile);
			append_moments(state,outfile);
			append_vorticity(state,outfile);
		} else {
			hvserror(status,"Run error");
			return 1;
		}
	}
	// Deinitialize solver
	free_solver(&state);
	endtime = time(NULL);
	if(!quiet) {
#ifdef HVS_PROFILE
		printf("Initialization moments time: %d sec\n",timings.init_moments);
		printf("Initialization coeficients time: %d sec\n",timings.init_coefs);
		printf("Initialization time: %d sec\n",timings.init);
		printf("Equation evaluation time: %d sec\n",timings.eval_equation);
		printf("RK step time: %d sec\n",timings.rk_step);
		printf("Update vorticity time: %d sec\n",timings.vorticity_update);
		printf("Integration time: %d sec\n",timings.integration);
#endif
		printf("Total computation time: %d sec\n",(int)(endtime-starttime));
	}
	return 0;
}

void usage(FILE *fh,char **argv) {
	fprintf(fh, "Hermitian vorticity solver.\n");
	fprintf(fh, "Usage: %s [-h|-q] -v file -o file [-t float|-l float|-n float|-b float|-e float]\n", argv[0]);
	fprintf(fh, "-h\t\tShow this message.\n");
	fprintf(fh, "-r\t\tRun in safe mode (save recovery moments file after each step).\n");
	fprintf(fh, "-v file\t\tVorticity data file.\n");
	fprintf(fh, "-m file\t\tMoments data file.\n");
	fprintf(fh, "-o file\t\tOutput file name.\n");
	fprintf(fh, "-t float\tIntegration time step. Default: 0.01.\n");
	fprintf(fh, "-p number\tNumber of integration time steps between vorticity\n");
	fprintf(fh, "\t\tfield outputs. Default: 0 (output vorticity only\n");
	fprintf(fh, "\t\tat the beginning and at the end).\n");
	fprintf(fh, "-n float\tViscosity of the fluid. Default: 0.1.\n");
	fprintf(fh, "-l float\tInitial lambda. Default: 1.0.\n");
	fprintf(fh, "-b float\tInitial time. Default: 0.0\n");
	fprintf(fh, "-e float\tEnd time. Default: 1.0.\n");
	fprintf(fh, "-x min:step:max\tGrid size in x direction. Default: -5.0:0.1:5.0.\n");
	fprintf(fh, "\t\tUsed only with moments data file.\n");
	fprintf(fh, "-y min:step:max\tGrid size in y direction. Default: -5.0:0.1:5.0.\n");
	fprintf(fh, "\t\tUsed only with moments data file.\n");
	fprintf(fh, "-q\t\tQuiet mode. Output only to the output file.\n");
}