#define FLOAT_TYPE float
#define UINT unsigned int

typedef struct {
	char	*initvortfile; /* Initial vorticity field file name */
	FLOAT_TYPE	t0;
	FLOAT_TYPE	t1;
	FLOAT_TYPE	timestep;
} hvs_params;

typedef struct {
	FLOAT_TYPE x,y;
} hvs_position;

typedef struct {
	FLOAT_TYPE	xval, yval;
} hvs_vector;

typedef FLOAT_TYPE hvs_moment;

typedef struct {
	hvs_position	position;
	hvs_moment	velocity;
} hvs_grid_point;

typedef struct {
	UINT	width;
	UINT	height;
	hvs_grid_point *array;
} hvs_grid;

typedef struct {
    hvs_position	position;
    hvs_vector		vector;
} hvs_field;

hvs_grid* init_solver(UINT gridwidth, UINT gridheight);
void free_solver(hvs_grid *grid);
int run_solver();
int format_output(hvs_grid *grid, size_t size, char *output);
