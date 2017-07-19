/*****************************************************************************
*
* MODULE:	r.damflood
*
* AUTHOR:	Roberto Marzocchi - roberto.marzocchi[]supsi.ch (2008)
*		Massimiliano Cannata - massimiliano.cannata[]supsi.ch (2008)
*
* PURPOSE:	Estimate the area potentially inundated in case of dam breaking
*
*
* COPYRIGHT:	(C) 2008 by Istituto Scienze della Terra-SUPSI
*
*		This program is free software under the GNU General Public
*		Licence (>=2). Read the file COPYING that cames with GRASS
*		for details.
*
/***************************************************************************/



/* libraries*/
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <grass/gis.h>
#include <math.h>
#include <grass/glocale.h>
#include <grass/gmath.h>
#include <grass/dbmi.h>
#include <grass/linkm.h>
#include <grass/bitmap.h>
//#include <grass/interpf.h>




/* simple functions*/
#define min(A,B) ((A) < (B) ? (A):(B))
#define max(A,B) ((A) > (B) ? (A):(B))
#define min4(A,B,C,D) min(min(A,B),min(C,D))
#define max4(A,B,C,D) max(max(A,B),max(C,D))

#define ALLOC_DIM 10000
#define PI 3.14159265



//#define hmin 0.01
/* lo calcolo dopo in funzione della velocita' massima 
#define timestep 0.1 */

int **G_alloc_imatrix(int rows, int cols)
{
	int **mmm;
   int i;
   mmm = (int **)G_calloc(rows, sizeof(int *));
	mmm[0] = (int *)G_calloc(rows * cols, sizeof(int));
	for (i = 1; i < rows; i++)
		mmm[i] = mmm[i - 1] + cols;
 	return mmm;
}

void G_free_imatrix(int **mmm)
{
	G_free(mmm[0]);
	G_free(mmm);
	mmm = NULL;
	return;
}


float **G_alloc_fmatrix(int rows, int cols)
{
	float **m;
   int i;
   m = (float **)G_calloc(rows, sizeof(float *));
	m[0] = (float *)G_calloc(rows * cols, sizeof(float));
	for (i = 1; i < rows; i++)
		m[i] = m[i - 1] + cols;
 	return m;
}

void G_free_fmatrix(float **m)
{
	G_free(m[0]);
	G_free(m);
	m = NULL;
	return;
}

double *G_alloc_vector(size_t n)
{
	return (double *)G_calloc(n, sizeof(double));
}


double **G_alloc_dmatrix(int rows, int cols)
{
	double **mm;
	int i;

	mm = (double **)G_calloc(rows, sizeof(double *));
	mm[0] = (double *)G_calloc(rows * cols, sizeof(double));
	for (i = 1; i < rows; i++)
	mm[i] = mm[i - 1] + cols;

	return mm;
}

void G_free_dmatrix(double **mm)
{
	G_free(mm[0]);
	G_free(mm);
	mm = NULL;
 	return;
}

void G_free_vector(double *v)
{
	G_free(v);
	v = NULL;
	return;
}



float velocita_breccia(int i,double h)
{
	//double h;
	//int i;
	float g=9.81;
	float v;

	if(i==1){
		v=0.93*sqrt(h);
	}else if (i==2){
		v=0.4*sqrt(2*g*h);
	}
	return v;
}

//*********************************************************************************************
/* main program */
int main(int argc, char *argv[]){

 /* typedef struct
  {
    double sx, sy, sz, dir;
  } startpt;

  startpt *a_start,*tmp_start;*/

  struct Cell_head cellhd;
  struct Cell_head window;
  struct History history;

  /* input-output raster files */
  int infd_ELEV,infd_LAKE,infd_DAMBREAK, infd_MANNING;
  int outfd_H,outfd_VEL,outfd_VEL_DIR,outfd_HMAX,outfd_T_HMAX,outfd_I_HMAX,outfd_VMAX,outfd_T_VMAX,outfd_I_VMAX,outfd_DIR_VMAX,outfd_IMAX,outfd_T_IMAX,outfd_WAVEFRONT;
  float g=9.81;
  /* mapset name locator */
  char *mapset_ELEV,*mapset_LAKE,*mapset_DAMBREAK,*mapset_MANNING;

  /* buffer for input-output raster */
  FCELL *inrast_ELEV;			/* elevation map [m]*/
  DCELL *inrast_LAKE;			/* water elevation in the map [m]*/
  FCELL *inrast_DAMBREAK;		/* break in the dam*/
  FCELL *inrast_MANNING;		/* manning*/
  DCELL *outrast_VEL;		   /* velocity [m/s] */
  DCELL *outrast_VEL_DIR;
  DCELL *outrast_H;	         /* water elevation output [m]*/
  DCELL *outrast_HMAX;
  DCELL *outrast_I_HMAX;
  DCELL *outrast_T_HMAX;
  DCELL *outrast_VMAX;
  DCELL *outrast_I_VMAX;
  DCELL *outrast_T_VMAX;
  DCELL *outrast_DIR_VMAX;
  DCELL *outrast_IMAX;
  DCELL *outrast_T_IMAX;
  DCELL *outrast_WAVEFRONT;
  /* cell counters */
  int nrows, ncols;
  int row, col;
  int num_cell, num_break;
  int method;
  int warn1=0,warn2=0;
  float Q=0.0, vol_res,fall, volume=0.0;
  float res_ew ,res_ns;
  float R_i;
  /* memory matrix */
  double **m_h1, **m_h2, **m_u1, **m_u2, **m_v1, **m_v2;
  float **m_z, **m_DAMBREAK, **m_m;
  float **m_hmax, **m_t_hmax, **m_i_hmax, **m_vmax, **m_t_vmax, **m_i_vmax, **m_dir_vmax, **m_imax, **m_t_imax, **m_wavefront;
  int ** m_lake;
  /*  other variables  */
  double water_elevation, profondita_soglia, Z_piu, Z_meno;
  double hmin=0.1;
  double h_dx, h_sx, h_up, h_dw, Fup, Fdw, Fdx, Fsx, Gup, Gdw, Gdx, Gsx;
  int temp_i;
  double temp_d, v_tot;
  double dZ_dx_down, dZ_dx_up, dZ_dx, dZ_dy_down, dZ_dy_up, dZ_dy;
  double cr_down, cr_up;
  double F, G, S;
  double u, v, V;
  double u_sx, u_dx, v_dx, v_sx, v_up, v_dw, u_up, u_dw;
  char* strcat(char* s, const char* ct);
  int m=1, M=1;
  int i, i_cont;
  double vel_0=0.0, vel_max=0.0, t;
  float timestep, velocity;
  // Parameters for the optimization of timestep using the CFL stability condition
  float timestep_ct, timestep_ct_temp;
  int DELTAT, TSTOP;
  char name1[20],name2[20],name3[20],name4[20],name5[20],name6[20],name7[20],name8[20],name9[20];
  int tmp_int, test;
  char tmp[15];
  int ntimes, pp;
  double *times;
  double opt_t, time;


  /* variables to handle user inputs and outputs */
  char *ELEV, *LAKE, *DAMBREAK, *MANNING, *OUT_VEL, *OUT_H, *OUT_HMAX, *OUT_VMAX, *OUT_IMAX, *OUT_WAVEFRONT;


  /***********************************************************************************************************************/
  /* GRASS structure */
  struct GModule *module;
  struct Option *input_ELEV, *input_LAKE, *input_DAMBREAK,*input_MANNING, *input_DELTAT, *input_TSTOP, *input_TIMESTEP;
  struct
    {
	struct Option *opt_t;
    }
    parm;
  struct Option *output_VEL,*output_H,*output_HMAX, *output_VMAX, *output_IMAX, *output_WAVEFRONT;
  struct Flag *flag_d;
  struct {
		struct Option *met;
	} opt;
 /* initialize GRASS */
  G_gisinit(argv[0]);

  //pgm_name = argv[0];

  module = G_define_module();
  module->keywords = _("raster, dambreak");
  module->description = _("Estimate the area potentially inundated in case of dam break");

  //OPTIONS

  flag_d = G_define_flag();
  flag_d->key = 'd';
  flag_d->description = _("Flow direction additional output (aspect visualization)");
  //flag_d ->guisection  = _("Options");

  opt.met = G_define_option();
  opt.met->key = "method";
  opt.met->type = TYPE_STRING;
  opt.met->required = NO;
  opt.met->description = _("Computational method for initial velocity estimation");
  opt.met->options = "dambreak-without_hypotesis,uniform drop in of lake,small dam breach";
  opt.met->answer = "dambreak-without_hypotesis";
  //opt.met->guisection  = _("Options");

  /* LEGENDA
  total_dambreak-without_hypotesis = nessuna ipotesi sulla velocita' iniziale
  total_dambreak = Hp altezza critica --> velocita' iniziale (h critica) = 0.93*sqrt(h);
  small_dam_breach = Hp stramazzo --> velocita' iniziale = 0.4*sqrt(2*g*h)
  */

  input_TIMESTEP = G_define_option();
  input_TIMESTEP->key	= "timestep";
  input_TIMESTEP->type = TYPE_DOUBLE;
  input_TIMESTEP->required = NO;
  input_TIMESTEP->multiple   = NO;
  input_TIMESTEP->answer     = "0.01";
  input_TIMESTEP->description = _("Initial computational time step [s] - CFL condition");
  //input_TIMESTEP->guisection  = _("Options");


  /* Define different options */
  input_ELEV = G_define_option();
  input_ELEV->key	= "elev";
  input_ELEV->type = TYPE_STRING;
  input_ELEV->required = YES;
  input_ELEV->gisprompt = "old,cell,raster";
  input_ELEV->description = _("Name of elevation raster map (including lake bathymetry and dam)");
  input_ELEV->guisection  = _("Input options");

  input_LAKE = G_define_option();
  input_LAKE->key	= "lake";
  input_LAKE->type = TYPE_STRING;
  input_LAKE->required = YES;
  input_LAKE->gisprompt = "old,cell,raster";
  input_LAKE->description = _("Name of water depth raster map");
  input_LAKE->guisection  = _("Input options");

  input_DAMBREAK = G_define_option();
  input_DAMBREAK->key = "dambreak";
  input_DAMBREAK->type = TYPE_STRING;
  input_DAMBREAK->required = YES;
  input_DAMBREAK->gisprompt = "old,cell,raster";
  input_DAMBREAK->description = _("Name of dam breach width raster map");
  input_DAMBREAK->guisection  = _("Input options");

  input_MANNING = G_define_option();
  input_MANNING->key = "manning";
  input_MANNING->type = TYPE_STRING;
  input_MANNING->required = YES;
  input_MANNING->gisprompt = "old,cell,raster";
  input_MANNING->description = _("Name of Manning's roughness coefficient raster map");
  input_MANNING->guisection  = _("Input options");

  input_TSTOP = G_define_option();
  input_TSTOP->key	= "tstop";
  input_TSTOP->type = TYPE_INTEGER;
  input_TSTOP->required = YES;
  input_TSTOP->multiple   = NO;
  input_TSTOP->description = _("Simulation time lenght [s]");
  input_TSTOP->guisection  = _("Input options");

  //OUTPUTS options
  input_DELTAT = G_define_option();
  input_DELTAT->key	= "deltat";
  input_DELTAT->type = TYPE_INTEGER;
  input_DELTAT->required = NO;
  input_DELTAT->multiple   = NO;
  input_DELTAT->description = _("Time-lag for output generation [s]");
  input_DELTAT->guisection  = _("Output options");

  parm.opt_t = G_define_option();
  parm.opt_t->key	= "opt_t";
  parm.opt_t->type = TYPE_DOUBLE;
  parm.opt_t->required = NO;
  //parm.opt_t->multiple   = NO;
  parm.opt_t->multiple   = YES;
  parm.opt_t->description = _("Additional instants for output map generation [s]");
  parm.opt_t->guisection  = _("Output options");

  output_H = G_define_option();
  output_H ->key = "h";
  output_H ->type = TYPE_STRING;
  output_H ->required = NO;
  output_H ->gisprompt = "new,cell,raster";
  output_H ->description = _("Prefix for water depth output raster maps");
  output_H ->guisection  = _("Output options");

  output_VEL = G_define_option();
  output_VEL ->key = "vel";
  output_VEL ->type = TYPE_STRING;
  output_VEL ->required = NO;
  output_VEL ->gisprompt = "new,cell,raster";
  output_VEL ->description = _("Prefix for water velocity output raster maps");
  output_VEL ->guisection  = _("Output options");

  output_HMAX = G_define_option();
  output_HMAX ->key = "hmax";
  output_HMAX ->type = TYPE_STRING;
  output_HMAX ->required = NO;
  output_HMAX ->gisprompt = "new,cell,raster";
  output_HMAX ->description = _("Name of output maximum water depth raster map; relative intensity [h*v] and time map are also generated");
  output_HMAX ->guisection  = _("Output options");

  output_VMAX = G_define_option();
  output_VMAX ->key = "vmax";
  output_VMAX ->type = TYPE_STRING;
  output_VMAX ->required = NO;
  output_VMAX ->gisprompt = "new,cell,raster";
  output_VMAX ->description = _("Name of output maximum water velocity raster map; relative intensity [h*v] and time map are also generated");
  output_VMAX ->guisection  = _("Output options");

  output_IMAX = G_define_option();
  output_IMAX ->key = "imax";
  output_IMAX ->type = TYPE_STRING;
  output_IMAX ->required = NO;
  output_IMAX ->gisprompt = "new,cell,raster";
  output_IMAX ->description = _("Name of output maximum intensity [h*v] raster map; relative time map are also generated");
  output_IMAX ->guisection  = _("Output options");

  output_WAVEFRONT = G_define_option();
  output_WAVEFRONT ->key = "wavefront";
  output_WAVEFRONT ->type = TYPE_STRING;
  output_WAVEFRONT ->required = NO;
  output_WAVEFRONT ->gisprompt = "new,cell,raster";
  output_WAVEFRONT ->description = _("Name of output wave front time[s] raster map");
  output_WAVEFRONT ->guisection  = _("Output options");
  

  if (G_parser(argc, argv))
    exit(EXIT_FAILURE);


  /***********************************************************************************************************************/
  /* get entered parameters */
  ELEV=input_ELEV->answer;
  LAKE=input_LAKE->answer;
  DAMBREAK=input_DAMBREAK->answer;
  MANNING=input_MANNING->answer;
  sscanf(input_TIMESTEP->answer, "%f", &timestep);
  //timestep=input_TIMESTEP->answer;
  if (parm.opt_t->answer != NULL) {
	  for (i = 0; parm.opt_t->answers[i]; i++) ;
	  ntimes=i;
	  times = G_alloc_vector(ntimes);
     for (i = 0; i < ntimes; i++) {
		  sscanf(parm.opt_t->answers[i], "%lf", &opt_t);
		  times[i]=opt_t;
	  }
	  pp=0;
  }
  	OUT_H=output_H->answer;
  	OUT_VEL=output_VEL->answer;
  	OUT_HMAX=output_HMAX->answer;
  	OUT_VMAX=output_VMAX->answer;
	OUT_IMAX=output_IMAX->answer;
        OUT_WAVEFRONT=output_WAVEFRONT->answer;
  	if (input_DELTAT->answer!= NULL){
		DELTAT = atoi(input_DELTAT->answer);
	}
	//DELTAT = atoi(input_DELTAT->answer);
  	TSTOP = atoi(input_TSTOP->answer);



  /* find maps in mapset */
  mapset_ELEV = G_find_cell2 (ELEV, "");
  if (mapset_ELEV == NULL)
    G_fatal_error (_("cell file [%s] not found"), ELEV);

  mapset_LAKE = G_find_cell2 (LAKE, "");
  if (mapset_LAKE == NULL)
    G_fatal_error (_("cell file [%s] not found"), LAKE);

  mapset_DAMBREAK = G_find_cell2 (DAMBREAK, "");
  if (mapset_DAMBREAK == NULL)
    G_fatal_error (_("cell file [%s] not found"), DAMBREAK);

  mapset_MANNING = G_find_cell2 (MANNING, "");
  if (mapset_MANNING == NULL)
    G_fatal_error (_("cell file [%s] not found"), MANNING);


	/* open input raster files */
  if ( (infd_ELEV = G_open_cell_old (ELEV, mapset_ELEV)) < 0)
    G_fatal_error (_("Cannot open cell file [%s]"), ELEV);
  if ( (infd_LAKE = G_open_cell_old (LAKE, mapset_LAKE)) < 0)
    G_fatal_error (_("Cannot open cell file [%s]"), LAKE);
  if ( (infd_DAMBREAK = G_open_cell_old (DAMBREAK, mapset_DAMBREAK)) < 0)
    G_fatal_error (_("Cannot open cell file [%s]"), DAMBREAK);
  if ( (infd_MANNING = G_open_cell_old (MANNING, mapset_MANNING)) < 0)
    G_fatal_error (_("Cannot open cell file [%s]"), MANNING);


  /* Check for some output map */
    if ((output_H->answer == NULL)
    && (output_VEL->answer == NULL)
    && (output_HMAX->answer == NULL)
    && (output_VMAX->answer == NULL)
    && (output_IMAX->answer == NULL)
    && (output_WAVEFRONT->answer == NULL)){
  		G_fatal_error(_("Sorry, you must choose an output map."));
    }

    if (flag_d->answer && (output_H->answer == NULL) ) {
    	G_fatal_error("You choose flow direction map without flow velocity prefix name");
    }

	/* check legal output name */
	if (OUT_H) {
	 if (G_legal_filename (OUT_H) < 0)
		G_fatal_error (_("[%s] is an illegal name"), OUT_H);
	}
	if (OUT_VEL) {
	 if (G_legal_filename (OUT_VEL) < 0)
		G_fatal_error (_("[%s] is an illegal name"), OUT_VEL);
	}
	if (OUT_HMAX) {
	 if (G_legal_filename (OUT_HMAX) < 0)
		G_fatal_error (_("[%s] is an illegal name"), OUT_HMAX);
	}
	if (OUT_VMAX) {
	 if (G_legal_filename (OUT_VMAX) < 0)
		G_fatal_error (_("[%s] is an illegal name"), OUT_VMAX);
	}
	if (OUT_IMAX) {
	 if (G_legal_filename (OUT_IMAX) < 0)
		G_fatal_error (_("[%s] is an illegal name"), OUT_IMAX);
	}
        if (OUT_WAVEFRONT) {
	 if (G_legal_filename (OUT_WAVEFRONT) < 0)
		G_fatal_error (_("[%s] is an illegal name"), OUT_WAVEFRONT);
	}
  /* type of dam failure*/
  if (strcmp(opt.met->answer, "uniform drop in of lake") == 0)
  	method=1;
  else if (strcmp(opt.met->answer, "small dam breach") == 0)
	method=2;
  else if (strcmp(opt.met->answer, "dambreak-without_hypotesis") == 0)
  	method=3;
  else
	G_fatal_error(_("Unknown method. Please, select a computational method"));

  /* Allocate input buffer */
  inrast_ELEV = G_allocate_f_raster_buf();
  inrast_LAKE = G_allocate_d_raster_buf();
  inrast_DAMBREAK = G_allocate_f_raster_buf();
  inrast_MANNING = G_allocate_f_raster_buf();


  /* get windows rows & cols */
  nrows	= G_window_rows();
  ncols	= G_window_cols();
  G_get_window(&window);
  res_ew = window.ew_res;
  res_ns = window.ns_res;

  /* allocate memory matrix */
  m_DAMBREAK = G_alloc_fmatrix(nrows,ncols);
  m_m = G_alloc_fmatrix(nrows,ncols);
  m_z = G_alloc_fmatrix(nrows,ncols);
  m_h1 = G_alloc_dmatrix(nrows,ncols);
  m_h2 = G_alloc_dmatrix(nrows,ncols);
  m_u1 = G_alloc_dmatrix(nrows,ncols);
  m_u2 = G_alloc_dmatrix(nrows,ncols);
  m_v1 = G_alloc_dmatrix(nrows,ncols);
  m_v2 = G_alloc_dmatrix(nrows,ncols);
  m_lake = G_alloc_imatrix(nrows,ncols);
  if (OUT_HMAX) {
  		m_hmax = G_alloc_fmatrix(nrows,ncols);
  		m_t_hmax = G_alloc_fmatrix(nrows,ncols);
  		m_i_hmax = G_alloc_fmatrix(nrows,ncols);
  }
  if (OUT_VMAX) {
  		m_vmax = G_alloc_fmatrix(nrows,ncols);
  		m_t_vmax = G_alloc_fmatrix(nrows,ncols);
  		m_i_vmax = G_alloc_fmatrix(nrows,ncols);
  		if (flag_d->answer) {
  			m_dir_vmax = G_alloc_fmatrix(nrows,ncols);
  		}
  }
  if (OUT_IMAX) {
  		m_imax = G_alloc_fmatrix(nrows,ncols);
  		m_t_imax = G_alloc_fmatrix(nrows,ncols);
  }
  if (OUT_WAVEFRONT) {
  		m_wavefront = G_alloc_fmatrix(nrows,ncols);
  }
	G_message("Reading input maps");
  	for (row = 0; row < nrows; row++)
  		{
      G_percent (row, nrows, 2);

      /* read a line input maps into buffers*/
      if (G_get_f_raster_row (infd_ELEV, inrast_ELEV, row) < 0)
	G_fatal_error (_("Could not read from <%s>"),ELEV);
      if (G_get_d_raster_row (infd_LAKE, inrast_LAKE, row) < 0)
	G_fatal_error (_("Could not read from <%s>"),LAKE);
      if (G_get_f_raster_row (infd_DAMBREAK, inrast_DAMBREAK, row) < 0)
	G_fatal_error (_("Could not read from <%s>"),DAMBREAK);
      if (G_get_f_raster_row (infd_MANNING, inrast_MANNING, row) < 0)
	G_fatal_error (_("Could not read from <%s>"),MANNING);

      /* read every cell in the line buffers */
      for (col = 0; col < ncols; col++)
			{

	  			/* store values in memory matrix (attenzione valori nulli!)*/

			  	m_DAMBREAK[row][col] = ((FCELL *) inrast_DAMBREAK)[col];
			  	m_m[row][col] = ((FCELL *) inrast_MANNING)[col];
			  	m_z[row][col] = ((FCELL *) inrast_ELEV)[col];
			  	m_h1[row][col] = ((DCELL *) inrast_LAKE)[col];
			  	m_h2[row][col] = m_h1[row][col];
			  	m_u1[row][col] = 0.0;
			  	m_u2[row][col] = 0.0;
			  	m_v1[row][col] = 0.0;
		  	  	m_v2[row][col] = 0.0;
		  		if (OUT_HMAX) {
		  	  		m_hmax[row][col] = m_h1[row][col];
		  	  		m_t_hmax[row][col] = 0.0;
		  	  		m_i_hmax[row][col]=0.0;
		  	  	}
		  	  	if (OUT_VMAX) {
		  	  		m_vmax[row][col] = 0.0;
		  	  		m_t_vmax[row][col] = 0.0;
		  	  		m_i_vmax[row][col]=0.0;
		  	  		if (flag_d->answer) {
		  	  			m_dir_vmax[row][col]=0.0;
		  	  		}
		  	  	}
		  	  	if (OUT_IMAX) {
		  	  		m_imax[row][col] = 0.0;
		  	  		m_t_imax[row][col] = 0.0;
		  	  	}
				if (OUT_WAVEFRONT) {
		  	  		m_wavefront[row][col] = 0.0;
		  	  	}
			}
    	}


  G_free(inrast_ELEV);
  G_free(inrast_LAKE);
  G_free(inrast_DAMBREAK);
  G_free( inrast_MANNING);
  G_close_cell(infd_ELEV);
  G_close_cell(infd_LAKE);
  G_close_cell (infd_DAMBREAK);
  G_close_cell (infd_MANNING);

	if (method==1 || method==2) {
		num_cell=0;
		/* cerco il lago */
		for (row = 0; row < nrows; row++)
			{
			for (col = 0; col < ncols; col++)
				{
				if (m_h1[row][col] != 0 ){
					water_elevation = m_h1[row][col] + m_z[row][col];
					m_lake[row][col]=1;
					num_cell++;
				}else {
					m_lake[row][col]=0;
				}
		}}
		num_break=0;
		G_message("Searching dam breach");
		for (row = 0; row < nrows; row++)
			{
			for (col = 0; col < ncols; col++)
				{
				/*  rottura diga */
				if (m_DAMBREAK[row][col] > 0){
					num_break++;
					G_message("(%d,%d)Cell Dam Breach n° %d",row,col,num_break);
					//printf("ho trovato la rottura diga (%d,%d)=\n", row,col);
					//while(!getchar()){ }
					profondita_soglia = water_elevation-(m_z[row][col] - m_DAMBREAK[row][col]) ;
					vel_0=velocita_breccia(method,profondita_soglia);
					volume = volume + profondita_soglia * res_ew * res_ns;
					//printf("Q=%.2f\n",Q);
					//vel_0 = 0.93 * sqrt(profondita_soglia);
					if (vel_0 > vel_max)
					vel_max=vel_0;
				} else {
					G_fatal_error(_("Don't find the dambreak - Please select a correct map or adjust the computational region"));
				}

				if (m_DAMBREAK[row][col] > 0) {
					//cambio il DTM inserendo la rottura della diga
					m_z[row][col]=m_z[row][col]-m_DAMBREAK[row][col];
					m_lake[row][col]=0;
					m_h1[row][col]=profondita_soglia;
					//printf("la velocità vel_max vale: %f\n",vel_max);
					//while(!getchar()){ }
					}
		  	}
		}
		G_message("the number of lake cell is': %d\n", num_cell);

		/**************************************/
		/* timestep in funzione di V_0 e res */
		/**************************************/
		//timestep=0.01;
		//timestep= ((res_ns+res_ew)/2.0)/(vel_max*50.0);
		// DEVELOPEMENT
		//*****************************************************************************
	 	G_message("vel on the dam break cells is %.2f, and timestep is %.2f",vel_max, timestep);
	 }

      		//*****************************************************************************
		// Uniform drop in of lake (method=1 or method =2) 
		//*****************************************************************************
		if (method==1 || method==2){
			/* calcolo l'abbassamento sul lago*/
			if (num_cell!=0) {
				fall = (volume) / (num_cell * res_ew * res_ns);
				//printf("volume=%f, fall=%f\n",volume,fall);
				//while(!getchar()){ }
			}
			vol_res=0.0;

			for (row = 1; row < nrows-1; row++)
		   {
				for (col = 1; col < ncols-1; col++)
				{
					if (m_DAMBREAK[row][col]>0){
						// ragiona se ha senso
						m_h2[row][col]=m_h1[row][col]-fall;
						if (m_h2[row][col]<=0) {
							m_h2[row][col]=0.0;
							if (m_h1[row][col]>0) {
							// questo warning va modificato perchè vale per ogni cella ---> bisogna metterne uno generico che valga quando tutte le celle sono con h=0
								num_break--;
								if (num_break==0){
									if (warn2==0){
										G_warning("At the time %.0fs no water go out from lake",t);
									}
								}
							}
						}
					}
					if (m_lake[row][col]==1){
						m_h2[row][col]=m_h1[row][col]-fall;
						if (m_h2[row][col]<=0) {
							vol_res = vol_res-m_h2[row][col]*res_ew * res_ns;
							m_lake[row][col]=-1;
							m_h2[row][col]=0.0;
							num_cell--;
						}
					}
			}}//end two for cicles
		} //end if
	// there isn't interest to find where is the lake --> everywhere m_lake[row][col]=0 
	if (method==3){
	 	for (row = 0; row < nrows; row++)
			{
			for (col = 0; col < ncols; col++)
				{
				/*  rottura diga */
				if (m_DAMBREAK[row][col] > 0){
					m_z[row][col]=m_z[row][col]-m_DAMBREAK[row][col];
					m_DAMBREAK[row][col]=-1.0;
					//timestep=0.01;
					m_lake[row][col]=0;
				} else {
					m_lake[row][col]=0;
	 }}}}

  	G_percent(nrows, nrows, 1);	/* finish it */

  	G_message("Model running");

	/* calculate time step loop */
	for(t=0; t<=TSTOP; t+=timestep){
	//printf("************************************************\n");
	//while(!getchar()){ }
	//G_percent(t, TSTOP, timestep);
	// ciclo sui tempi
        //G_message("timestep =%f,t=%f",timestep,t);
	   if (t>M*100){
	   	if (M*100!=(m-1)*DELTAT)
	   		G_percent(ceil(t), TSTOP, 2);
	   		G_message("t:%d",M*100);
	   	M++;
	   }
	  	//printf("t:%lf\n",t);
	
	// DESCRIPTION OF METHOD (italian --> TRASLATE)
   	// primo ciclo: calcolo nuove altezze dell'acqua al tempo t+1
   	// 					- a valle della diga applico l'equazione di continuita' delle shallow water
   	// 					  in pratica la nuova altezza e' valutata attraverso un bilancio dei
   	// 					  flussi in ingresso e in uscita nelle due direzioni principali
   	// 					- a monte delle diga:
	//                                      	- nel metodo 1 e 2 :l'equazione di continuita' e' applicata al volume del lago
   	//					  	fisicamente questo porta a una minore realisticita' ma evita le oscillazioni che
   	// 					  	sono causa di instabilita' numerica
	//						- nel caso piu' generale si applicano le equazioni a tutto il lago


		for (row = 1; row < nrows-1; row++)
      {
			for (col = 1; col < ncols-1; col++)
			{
				if (m_lake[row][col]==0 && m_DAMBREAK[row][col]<=0) {
                                        
					//*******************************************/
					/* CONTINUITY EQUATION --> h(t+1)           */
					//*******************************************/
					// x direction
					// right intercell
					if (m_u1[row][col]>0 && m_u1[row][col+1]>0) {
						Fdx = m_u1[row][col]*m_h1[row][col];
					} else if (m_u1[row][col]<0 && m_u1[row][col+1]<0) {
						Fdx = m_u1[row][col+1]*m_h1[row][col+1] ;
					} else {
						u_dx = (m_u1[row][col]+m_u1[row][col+1])/2.0;
						if ( (u_dx < 0 && m_u1[row][col+1]==0) || (u_dx > 0 && m_u1[row][col]==0))
							u_dx=0;
						if (u_dx>=0) {
							h_dx = max(m_h1[row][col]+m_z[row][col]-m_z[row][col+1],0);
						} else {
							h_dx = max(m_h1[row][col+1]+m_z[row][col+1]-m_z[row][col],0);
						}
						Fdx = h_dx * u_dx;
					}

					// left intercell
					if (m_u1[row][col-1]>0 && m_u1[row][col]>0) {
						Fsx = m_u1[row][col-1] * m_h1[row][col-1];
					} else if (m_u1[row][col-1]<0 && m_u1[row][col]<0) {
						Fsx = m_u1[row][col] * m_h1[row][col];
					} else {
						u_sx = (m_u1[row][col-1]+m_u1[row][col])/2.0;
						if ( (u_sx < 0 && m_u1[row][col]==0) || (u_sx > 0 && m_u1[row][col-1]==0))
							u_sx = 0;
						if (u_sx>=0) {
							h_sx = max(m_h1[row][col-1]+m_z[row][col-1]-m_z[row][col],0);
						} else {
							h_sx = max(m_h1[row][col]+m_z[row][col]-m_z[row][col-1],0);
						}
						Fsx = h_sx * u_sx;
					}



					if(m_DAMBREAK[row][col+1]>0 && ((m_h1[row][col]+m_z[row][col]) < (m_h1[row][col+1]+m_z[row][col+1]))){
						Fdx = -m_h1[row][col+1]*velocita_breccia(method,m_h2[row][col+1]);
						if (m_h1[row][col+1]==0)
							Fdx=0.0;
					}
					if (m_DAMBREAK[row][col-1]>0 && ((m_h1[row][col]+m_z[row][col]) < (m_h1[row][col-1]+m_z[row][col-1]))){
						Fsx = m_h1[row][col-1]*velocita_breccia(method,m_h2[row][col-1]);
						if (m_h1[row][col-1]==0)
							Fsx=0.0;
					}
					F = Fdx - Fsx;

					// dGup =m_v1[row][col] * m_h1[row][col] ;irezione y
					// intercella up
					if (m_v1[row][col]>0 && m_v1[row-1][col]>0) {
						Gup = m_v1[row][col] * m_h1[row][col];
					} else if (m_v1[row][col]<0 && m_v1[row-1][col]<0) {
						Gup = m_v1[row-1][col] * m_h1[row-1][col];
					} else {
						v_up = (m_v1[row][col]+m_v1[row-1][col])/2.0;
					   if ( (v_up<0 && m_v1[row-1][col]==0) || (v_up>0 && m_v1[row][col]==0))
					   	v_up=0;
					   if (v_up>=0){
					   	h_up = max(m_h1[row][col]+m_z[row][col]-m_z[row-1][col],0);
					  	} else {
					  		h_up = max(m_h1[row-1][col]+m_z[row-1][col]-m_z[row][col],0);
					  	}
						Gup = h_up * v_up;
					}

					// intercella down
					if (m_v1[row+1][col]>0 && m_v1[row][col]>0) {
						Gdw = m_v1[row+1][col] * m_h1[row+1][col];
					} else if (m_v1[row+1][col]<0 && m_v1[row][col]<0) {
						Gdw = m_v1[row][col] * m_h1[row][col];
					} else {
						v_dw = (m_v1[row][col]+m_v1[row+1][col])/2.0;
					   if ((v_dw<0 && m_v1[row][col]==0) || (v_dw>0 && m_v1[row+1][col]==0))
					   	v_dw = 0;
					   if (v_dw>=0) {
					   	h_dw = max(m_h1[row+1][col]+m_z[row+1][col]-m_z[row][col],0);
					   } else {
					   	h_dw = max(m_h1[row][col]+m_z[row][col]-m_z[row+1][col],0);
					   }

					   Gdw = h_dw * v_dw;
					}

					if (m_DAMBREAK[row-1][col]>0 && ((m_h1[row][col]+m_z[row][col]) < (m_h1[row-1][col]+m_z[row-1][col]))){
						Gup = -m_h1[row-1][col]*velocita_breccia(method,m_h1[row-1][col]);
						if (m_h1[row-1][col]==0)
							Gup=0.0;
					}
					if (m_DAMBREAK[row+1][col]>0 && ((m_h1[row][col]+m_z[row][col]) < (m_h1[row+1][col]+m_z[row+1][col]))){
						Gup = m_h1[row-1][col]*velocita_breccia(method,m_h1[row+1][col]);
						if (m_h1[row+1][col]==0)
							Gdw=0.0;
					}
					G = Gup - Gdw;

					//equazione
					m_h2[row][col] = m_h1[row][col] - timestep / res_ew * F - timestep / res_ns * G;

					/*if ((row==20||row==21||row==22||row==23)&&(col==18||col==19)){
						printf("EQ. CONTINUITA' --> row:%d, col:%d\n)",row, col);
						printf("m_h1[row][col]:%f,m_u1[row][col]:%f,m_v1[row][col]:%f",m_h1[row][col],m_u1[row][col],m_v1[row][col]);
						printf("m_h1[row][col+1]:%f,m_h1[row][col-1]:%f,m_h1[row+1][col]:%f, m_h1[row-1][col]:%f\n",m_h1[row][col+1],m_h1[row][col-1],m_h1[row+1][col], m_h1[row-1][col]);
						printf("h_dx:%f, h_sx:%f, h_up%f, h_dw:%f\n",h_dx, h_sx, h_up, h_dw);
						printf("m_u1[row][col+1]:%f,m_u1[row][col-1]:%f,m_v1[row+1][col]:%f, m_v1[row-1][col]:%f\n",m_u1[row][col+1],m_u1[row][col-1],m_v1[row+1][col], m_v1[row-1][col]);
						printf("v_up: %f,v_dw:%f,u_dx:%f,u_sx:%f \n",v_up, v_dw, u_dx, u_sx);
						printf("Fdx: %f,Fsx: %f, F: %f, Gup:%f, Gdw:%f, G: %f\n",Fdx, Fsx, F,Gup, Gdw, G);
						printf("m_h2(row,col): %f\n \n", m_h2[row][col]);
					}*/




					if( (row==1 || row==(nrows-2) || col==1 || col==(ncols-2)) && (m_v2[1][col]>0 || m_v2[nrows-2][col]<0 || m_u1[row][1]<0 || m_u1[row][ncols-2]>0 )){
						if (warn1==0){
							G_warning("At the time %.3f the computational region is smaller than inundation",t);
							warn1=1;
						}
					}
					if (m_h2[row][col]<0){
						/*G_warning("At the time %f h is lesser than 0 h(%d,%d)=%f",t, row,col,m_h2[row][col]);
					   printf("row:%d, col:%d, H minore di zero: %.30lf)",row, col, m_h2[row][col]);
					   printf("DATI:\n");
						printf("row:%d,col%d,hmin:%g,h2:%.30lf \n ",row,col,hmin,m_h2[row][col]);
						printf("m_z[row][col]:%f\n", m_z[row][col]);
						printf("m_h1[row][col]:%.30lf\n",m_h1[row][col]);
						printf("m_u1[row][col]:%.30lf,m_v1[row][col]:%.30lf\n",m_u1[row][col], m_v1[row][col]);
						printf("m_z[row][col+1]:%f,m_z[row][col-1]:%f,m_z[row+1][col]:%f, m_z[row-1][col]:%f\n",m_z[row][col+1],m_z[row][col-1],m_z[row+1][col], m_z[row-1][col]);
						printf("m_h1[row][col+1]:%.30lf,m_h1[row][col-1]:%.30lf,m_h1[row+1][col]:%.30lf, m_h1[row-1][col]:%.30lf\n",m_h1[row][col+1],m_h1[row][col-1],m_h1[row+1][col], m_h1[row-1][col]);
						printf("h_dx:%f, h_sx:%f, h_up%f, h_dw:%f\n",h_dx, h_sx, h_up, h_dw);
						printf("m_u1[row][col+1]:%.30lf,m_u1[row][col-1]:%.30lf,m_v1[row+1][col]:%.30lf, m_v1[row-1][col]:%.30lf\n",m_u1[row][col+1],m_u1[row][col-1],m_v1[row+1][col], m_v1[row-1][col]);
						printf("timestep: %.30lf, res_ew: %.30lf, res_ns:%.30lf\n",timestep, res_ew, res_ns);
						printf("v_up: %f,v_dw:%f,u_dx:%f,u_sx:%f \n",v_up, v_dw, u_dx, u_sx);
						printf("Fdx: %.30lf,Fsx: %.30lf, F: %.30lf, Gup:%.30lf, Gdw:%.30lf, G: %.30lf\n",Fdx, Fsx, F,Gup, Gdw, G);
						printf("row: %d, col %d, m_h1(row,col): %.30lf, m_h2(row,col): %.30lf \n",row, col,m_h1[row][col], m_h2[row][col]);
						printf("m_DAMBREAk(ROW,COL):%f \n",m_DAMBREAK[row][col]);
						while(!getchar()){ }*/
						m_h2[row][col]=0;
					}

				} // fine continuita' a valle


				if (method==1 || method==2){
					//*******************************************************************
					// calcolo portata Q uscente dal lago solo nel caso di Hp stramazzo
					/* HP: method 1 or 2   */
					if (m_DAMBREAK[row][col]>0 ){
						if ((m_z[row][col]+m_h1[row][col])>(m_z[row][col+1]+m_h1[row][col+1])){
							if (t==timestep)
							Q = Q + m_h1[row][col]* velocita_breccia(method,m_h1[row][col]) * res_ns;
							m_u1[row][col]= velocita_breccia(method,m_h1[row][col]);
						} else if ((m_z[row][col]+m_h1[row][col])>(m_z[row][col-1]+m_h1[row][col-1])){
							Q = Q + m_h1[row][col]* velocita_breccia(method,m_h1[row][col]) * res_ns;
							m_u1[row][col]= -velocita_breccia(method,m_h1[row][col]);
						}
						if ((m_z[row][col]+m_h1[row][col])>(m_z[row+1][col]+m_h1[row+1][col])){
							Q = Q + m_h1[row][col]* velocita_breccia(method,m_h1[row][col]) * res_ew;
							m_v1[row][col]=velocita_breccia(method,m_h1[row][col]);
						} else if ((m_z[row][col]+m_h1[row][col])>(m_z[row-1][col]+m_h1[row-1][col])){
							Q = Q + m_h1[row][col]* velocita_breccia(method,m_h1[row][col]) * res_ew;
							m_v1[row][col]=-velocita_breccia(method,m_h1[row][col]);
						}
					}
				}
		}} //end two for cicles


		//*****************************************************************************
		// abbassamento lago (siccome c'e due volte fare poi una function)
		//*****************************************************************************
		if (method==1 || method==2){
			/* calcolo l'abbassamento sul lago*/
			if (num_cell!=0) {
				fall = (Q * timestep-vol_res) / (num_cell * res_ew * res_ns);
			} else {
				if (warn2==0){
					G_warning("At the time %.0fs no water go out from lake",t);
					warn2=1;
				}
			}
			vol_res=0.0;
			Q=0.0;

			for (row = 1; row < nrows-1; row++)
		   {
				for (col = 1; col < ncols-1; col++)
				{
					if (m_DAMBREAK[row][col]>0){
						m_h2[row][col]=m_h1[row][col]-fall;
						if (m_h2[row][col]<=0) {
							m_h2[row][col]=0.0;
							if (m_h1[row][col]>0) {
								num_break--;
								/*if (num_break==0){
									G_warning("At the time %.0fs no water go out from lake",t);
								}*/
							}
						}
					}
					if (m_lake[row][col]==1){
						m_h2[row][col]=m_h1[row][col]-fall;
						if (m_h2[row][col]<=0) {
							vol_res = vol_res-m_h2[row][col]*res_ew * res_ns;
							m_lake[row][col]=-1;
							m_h2[row][col]=0.0;
							num_cell--;
						}
					}
			}}//end two for cicles
		} //end if
                

		// DESCRIPTION OF METHOD (italian --> TRASLATE)
		//**********************************************************************************
		// terzo ciclo completo sulla matrice: applico le  -->
		// EQUAZIONI DEL MOTO IN DIREZIONE X e Y
		// e quindi calcolo u(t+1) e v(t+1)
		//
		// NOTA:
		// u(i,j) e v (i,j) sono le velocita' medie della cella i,j
		/*******************************************************************/
		for (row = 1; row < nrows-1; row++)
		{
			for (col = 1; col < ncols-1; col++)
			{
				if (m_lake[row][col]==0 && m_h2[row][col]>=hmin){

					/**********************************************************************************************************************/
					/* EQUAZIONE DEL MOTO IN DIREZIONE X */
					// right intercell
					if (m_u1[row][col]>0 && m_u1[row][col+1]>0) {
						Fdx = m_u1[row][col] * m_u1[row][col] * m_h1[row][col];
					} else if (m_u1[row][col]<0 && m_u1[row][col+1]<0) {
						Fdx = m_u1[row][col+1] * m_u1[row][col+1] * m_h1[row][col+1] ;
					} else {
						u_dx = (m_u1[row][col]+m_u1[row][col+1])/2.0;
						if ( (u_dx < 0 && m_u1[row][col+1]==0) || (u_dx > 0 && m_u1[row][col]==0))
							u_dx=0;
						if (u_dx>=0) {
							h_dx = max(m_h1[row][col]+m_z[row][col]-m_z[row][col+1],0);
						} else {
							h_dx = max(m_h1[row][col+1]+m_z[row][col+1]-m_z[row][col],0);
						}
						Fdx = h_dx * u_dx * u_dx;
					}

					// left intercell
					if (m_u1[row][col-1]>0 && m_u1[row][col]>0) {
						Fsx = m_u1[row][col-1] * m_u1[row][col-1] * m_h1[row][col-1];
					} else if (m_u1[row][col-1]<0 && m_u1[row][col]<0) {
						Fsx = m_u1[row][col] * m_u1[row][col] * m_h1[row][col];
					} else {
						u_sx = (m_u1[row][col-1]+m_u1[row][col])/2.0;
						if ( (u_sx < 0 && m_u1[row][col]==0) || (u_sx > 0 && m_u1[row][col-1]==0))
							u_sx = 0;
						if (u_sx>=0) {
							h_sx = max(m_h1[row][col-1]+m_z[row][col-1]-m_z[row][col],0);
						} else {
							h_sx = max(m_h1[row][col]+m_z[row][col]-m_z[row][col-1],0);
						}
						Fsx = h_sx * u_sx * u_sx;
					}

					if(m_DAMBREAK[row][col+1]>0 && ((m_h2[row][col]+m_z[row][col]) < (m_h2[row][col+1]+m_z[row][col+1]))){
						Fdx = m_h1[row][col+1]* pow(-velocita_breccia(method,m_h1[row][col+1]),2.0);  // -vel al quadrato perde il segno meno
						if (m_h2[row][col+1]==0)
							Fdx=0.0;
					}
					if (m_DAMBREAK[row][col-1]>0 && ((m_h2[row][col]+m_z[row][col]) < (m_h2[row][col-1]+m_z[row][col-1]))){
						Fsx = m_h1[row][col-1]*pow(velocita_breccia(method,m_h1[row][col-1]),2.0);
						if (m_h2[row][col-1]==0)
							Fsx=0.0;
					}
					F = Fdx - Fsx;

					//y
					// intercella up
					if (m_v1[row][col]>0 && m_v1[row-1][col]>0) {
						Gup = m_v1[row][col] * m_u1[row][col] * m_h1[row][col];
					} else if (m_v1[row][col]<0 && m_v1[row-1][col]<0) {
						Gup = m_v1[row-1][col] * m_u1[row-1][col] * m_h1[row-1][col];
					} else {
						v_up = (m_v1[row][col]+m_v1[row-1][col])/2.0;
					   if ( (v_up<0 && m_v1[row-1][col]==0) || (v_up>0 && m_v1[row][col]==0))
					   	v_up=0;
					   u_up = (m_u1[row][col]+m_u1[row-1][col])/2.0;
					   if (v_up>=0){
					   	h_up = max(m_h1[row][col]+m_z[row][col]-m_z[row-1][col],0);
					  	} else {
					  		h_up = max(m_h1[row-1][col]+m_z[row-1][col]-m_z[row][col],0);
					  	}
						Gup = h_up * v_up * u_up;
					}

					// intercella down
					if (m_v1[row+1][col]>0 && m_v1[row][col]>0) {
						Gdw = m_v1[row+1][col] * m_u1[row+1][col] * m_h1[row+1][col];
					} else if (m_v1[row+1][col]<0 && m_v1[row][col]<0) {
						Gdw = m_v1[row][col] * m_u1[row][col] * m_h1[row][col];
					} else {
						v_dw = (m_v1[row][col]+m_v1[row+1][col])/2.0;
					   if ((v_dw<0 && m_v1[row][col]==0) || (v_dw>0 && m_v1[row+1][col]==0))
					   	v_dw = 0;
					   u_dw = (m_u1[row][col]+m_u1[row+1][col])/2.0;
					   if (v_dw>=0) {
					   	h_dw = max(m_h1[row+1][col]+m_z[row+1][col]-m_z[row][col],0);
					   } else {
					   	h_dw = max(m_h1[row][col]+m_z[row][col]-m_z[row+1][col],0);
					   }
					   Gdw = h_dw * u_dw * v_dw;
					}


					if(m_DAMBREAK[row-1][col]>0 && ((m_h2[row][col]+m_z[row][col]) < (m_h2[row-1][col]+m_z[row-1][col]))){
						Gup = m_h1[row-1][col] * (-velocita_breccia(method,m_h1[row-1][col]) * m_u1[row-1][col]);
						if (m_h2[row-1][col]==0)
							Gup=0.0;
					}
					if (m_DAMBREAK[row+1][col]>0 && ((m_h2[row][col]+m_z[row][col]) < (m_h2[row+1][col]+m_z[row+1][col]))){
						Gdw = m_h1[row+1][col] *(velocita_breccia(method,m_h1[row+1][col]) * m_u1[row+1][col]);
						if (m_h2[row+1][col]==0)
							Gdw=0.0;
					}
					G = Gup - Gdw;


					//courant number  --> UPWIND METHOD
					if(m_u1[row][col]>0 && m_u1[row][col+1]>0 && m_u1[row][col-1]>0){
						test=1;
						dZ_dx_down = ( (m_h2[row][col+1] + m_z[row][col+1]) - (m_h2[row][col] + m_z[row][col] )) / res_ew;
						if (m_h2[row][col-1]==0 && m_z[row][col-1]>(m_h2[row][col] + m_z[row][col])){
							dZ_dx_up = 0;
						} else {
							dZ_dx_up = ( (m_h2[row][col] + m_z[row][col]) - (m_h2[row][col-1] + m_z[row][col-1]) ) / res_ew;
						}
						cr_down = (timestep / res_ew)* (fabs(m_u1[row][col+1]) + fabs(m_u1[row][col]))/2.0;
						cr_up = (timestep / res_ew)* (fabs(m_u1[row][col]) + fabs(m_u1[row][col-1]))/2.0;
						Z_piu = 0.0;
						Z_meno = 0.0;
					} else if (m_u1[row][col]<0 && m_u1[row][col-1]<0 && m_u1[row][col+1]<0){
						test=2;
						dZ_dx_down = ( (m_h2[row][col] + m_z[row][col]) - (m_h2[row][col-1] + m_z[row][col-1]) ) / res_ew;
						if (m_h2[row][col+1]==0 && m_z[row][col+1]> (m_h2[row][col] + m_z[row][col])){
							dZ_dx_up = 0;
						} else {
							dZ_dx_up = ( (m_h2[row][col+1] + m_z[row][col+1]) - (m_h2[row][col] + m_z[row][col]) ) / res_ew;
						}
						cr_down = (timestep / res_ew)* (fabs(m_u1[row][col]) + fabs(m_u1[row][col-1]))/2.0;
						cr_up = (timestep / res_ew)* (fabs(m_u1[row][col+1]) + fabs(m_u1[row][col]))/2.0;
						Z_piu = 0.0;
						Z_meno = 0.0;
					} else {
						test=3;
						if (m_h2[row][col+1]==0 && m_z[row][col+1]> (m_h2[row][col] + m_z[row][col])){
							Z_piu=(m_h2[row][col] + m_z[row][col]);
						}else{
							Z_piu = ((m_h2[row][col+1] + m_z[row][col+1])+(m_h2[row][col] + m_z[row][col]) ) / 2;
						}
						if (m_h2[row][col-1]==0 && m_z[row][col-1]>(m_h2[row][col] + m_z[row][col])){
							Z_meno = (m_h2[row][col] + m_z[row][col]);
						}else{
							Z_meno = ((m_h2[row][col-1] + m_z[row][col-1])+(m_h2[row][col] + m_z[row][col]) ) / 2;
						}
						dZ_dx_down = (Z_piu - Z_meno)/res_ew;
						dZ_dx_up = (Z_piu - Z_meno)/res_ew;
						cr_down = (timestep / res_ew)* (fabs(m_u1[row][col+1]) + 2*fabs(m_u1[row][col]) + fabs(m_u1[row][col-1]) )/4.0;
						cr_up = (timestep / res_ew)* (fabs(m_u1[row][col+1]) + 2*fabs(m_u1[row][col]) + fabs(m_u1[row][col-1]) )/4.0;
					}

						if (m_DAMBREAK[row][col+1] > 0 && m_h2[row][col+1] == 0){
							test=4;
							dZ_dx_up = 0.0;
							dZ_dx_down=( (m_h2[row][col] + m_z[row][col]) - (m_h2[row][col-1] + m_z[row][col-1]) ) / res_ew;
						}
						if (m_DAMBREAK[row][col-1] > 0 && m_h2[row][col-1] == 0){
							test=5;
							dZ_dx_up = 0.0;
							dZ_dx_down = ( (m_h2[row][col+1] + m_z[row][col+1]) - (m_h2[row][col] + m_z[row][col] )) / res_ew;
						}
						dZ_dx = (1-sqrt(cr_down)) * dZ_dx_down + sqrt(cr_up) * dZ_dx_up;
						//dZ_dx = 0.5 * dZ_dx_down + 0.5 * dZ_dx_up;

					if (m_h1[row][col]<hmin			)
						R_i=hmin;
					else
						R_i=m_h1[row][col];

					u = m_u1[row][col];
					v = m_v1[row][col];
					V = sqrt(pow(u,2.0) + pow(v,2.0));
					S = (- g * m_h2[row][col] * dZ_dx) -g*(pow(m_m[row][col],2.0) * u * V / pow(R_i,(1.0/3.0)));

			   	if (m_DAMBREAK[row][col] > 0){
			   		if ((m_z[row][col]+m_h2[row][col]) > (m_z[row][col+1]+m_h2[row][col+1]))
			   			m_u2[row][col] = velocita_breccia(method,m_h2[row][col]);  // velocita' sullo stramazzo
			   		else if ((m_z[row][col] + m_h2[row][col]) > (m_z[row][col-1] + m_h2[row][col-1]))
			   			m_u2[row][col] = - velocita_breccia(method,m_h2[row][col]);  // velocita' sullo stramazzo
			   		else
			   			m_u2[row][col] = 0.0;
			   	}else{
						m_u2[row][col] = 1.0 / m_h2[row][col] * (m_h1[row][col] * m_u1[row][col] - timestep / res_ew * F - timestep / res_ns * G + timestep * S );
						}
					// no velocita' contro la diga
	 				/*if (m_z[row][col+1]> water_elevation && m_u2[row][col]>0)
		     			m_u2[row][col]=0.0;
     				if (m_z[row][col-1] > water_elevation && m_u2[row][col]<0)
       				m_u2[row][col]=0.0;*/


       			if ((timestep/res_ew*(fabs(m_u2[row][col])+sqrt(g*m_h2[row][col])))>1.0){
				G_warning("\nATTENTION: at time: %f the Courant-Friedrich-Lewy stability condition isn't respected",t);
       				/*G_message("velocita' lungo x\n");
				G_message("row:%d, col%d \n",row,col);
				G_message("dZ_dx_down:%f, dZ_dx_up:%f,cr_up:%f, cr_down:%f\n" , dZ_dx_down,dZ_dx_up, cr_up, cr_down);
				G_message("Z_piu:%f,Z_meno:%f\n", Z_piu, Z_meno);
				G_message("dZ_dx:%f\n",dZ_dx);
				G_message("m_h1[row][col]:%f, m_h2[row][col]:%f, m_z[row][col]:%f\n",m_h1[row][col], m_h2[row][col], m_z[row][col]);
				G_message("m_h2[row][col+1]:%f, m_z[row][col+1]:%f,m_h2[row][col-1]:%f, m_z[row][col-1]:%f \n",m_h2[row][col+1], m_z[row][col+1],m_h2[row][col-1],m_z[row][col-1]);
				G_message("m_h2[row][col+1]:%f,m_h2[row][col-1]:%f,\n",m_h2[row][col+1],m_h2[row][col-1]);
				G_message("h_up: %f,h_dw:%f,h_dx:%f,h_sx:%f \n",h_up, h_dw, h_dx, h_sx);
				G_message("Fdx: %f,Fsx: %f, F: %f, Gup:%f, Gdw:%f, G: %.60lf,  S: %.60lf \n",Fdx, Fsx, F,Gup, Gdw, G, S);
				G_message("m_u1[row][col-1]:%f, m_h1[row][col-1]:%f, m_u1[row][col+1]:%f, m_h1[row][col+1]:%f\n",m_u1[row][col-1], m_h1[row][col-1], m_u1[row][col+1], m_h1[row][col+1]);
				G_message("timestep:%f, res_ew:%f\n",timestep,res_ew);
		 		G_message("m_u2[row][col]:%f,m_u1[row][col]:%f\n\n", m_u2[row][col],m_u1[row][col]);
				G_warning("   ");*/
       			}

					if (fabs(m_u2[row][col]>=1000 )){
						G_warning("At the time %f u(%d,%d)=%f", t, row,col,m_u2[row][col]);
		         }
	/*************************************************************************************************************************************************


	/*************************************************************************************************************************************************
					/* EQUAZIONE DEL MOTO IN DIREZIONE Y */
					// right intercell
					if (m_u1[row][col]>0 && m_u1[row][col+1]>0) {
						Fdx = m_u1[row][col] * m_v1[row][col] * m_h1[row][col];
					} else if (m_u1[row][col]<0 && m_u1[row][col+1]<0) {
						Fdx = m_u1[row][col+1] * m_v1[row][col+1] * m_h1[row][col+1] ;
					} else {
						u_dx = (m_u1[row][col]+m_u1[row][col+1])/2.0;
						if ( (u_dx < 0 && m_u1[row][col+1]==0) || (u_dx > 0 && m_u1[row][col]==0))
							u_dx=0;
						v_dx = (m_v1[row][col]+m_v1[row][col+1])/2.0;
						if (u_dx>=0) {
							h_dx = max(m_h1[row][col]+m_z[row][col]-m_z[row][col+1],0);
						} else {
							h_dx = max(m_h1[row][col+1]+m_z[row][col+1]-m_z[row][col],0);
						}
						Fdx = h_dx * u_dx * v_dx;
					 }

					// left intercell
					if (m_u1[row][col-1]>0 && m_u1[row][col]>0) {
						Fsx = m_u1[row][col-1] * m_v1[row][col-1] * m_h1[row][col-1];
					} else if (m_u1[row][col-1]<0 && m_u1[row][col]<0) {
						Fsx = m_u1[row][col] * m_v1[row][col] * m_h1[row][col];
					} else {
						u_sx = (m_u1[row][col-1]+m_u1[row][col])/2.0;
						if ( (u_sx < 0 && m_u1[row][col]==0) || (u_sx > 0 && m_u1[row][col-1]==0))
							u_sx = 0;
						v_sx = (m_v1[row][col-1]+m_v1[row][col])/2.0;
						if (u_sx>=0) {
							h_sx = max(m_h1[row][col-1]+m_z[row][col-1]-m_z[row][col],0);
						} else {
							h_sx = max(m_h1[row][col]+m_z[row][col]-m_z[row][col-1],0);
						}
						Fsx = h_sx * u_sx * v_sx;
					}

					if(m_DAMBREAK[row][col+1]>0.0 && ((m_h2[row][col]+m_z[row][col]) < (m_h2[row][col+1]+m_z[row][col+1]))){
						Fdx = m_h1[row][col+1]* (-velocita_breccia(method,m_h1[row][col+1])) * m_v1[row][col+1];
						if (m_h2[row][col+1]==0)
							Fdx=0.0;
					}
					if (m_DAMBREAK[row][col-1]>0.0 && ((m_h2[row][col]+m_z[row][col]) < (m_h2[row][col-1]+m_z[row][col-1]))){
						Fsx = m_h1[row][col-1] * velocita_breccia(method,m_h1[row][col-1]) * m_v1[row][col-1];
						if (m_h2[row][col-1]==0)
							Fsx=0.0;
					}
					F = Fdx - Fsx;


					//y
					// intercella up
					if (m_v1[row][col]>0 && m_v1[row-1][col]>0) {
						Gup = m_v1[row][col] * m_v1[row][col] * m_h1[row][col];
					} else if (m_v1[row][col]<0 && m_v1[row-1][col]<0) {
						Gup = m_v1[row-1][col] * m_v1[row-1][col] * m_h1[row-1][col];
					} else {
						v_up = (m_v1[row][col]+m_v1[row-1][col])/2.0;
					   if ( (v_up<0 && m_v1[row-1][col]==0) || (v_up>0 && m_v1[row][col]==0))
					   	v_up=0;
					  if (v_up>=0){
					   	h_up = max(m_h1[row][col]+m_z[row][col]-m_z[row-1][col],0);
					  	} else {
					  		h_up = max(m_h1[row-1][col]+m_z[row-1][col]-m_z[row][col],0);
					  	}
					   Gup = h_up * v_up * v_up;
					}

					// intercella down
					if (m_v1[row+1][col]>0 && m_v1[row][col]>0) {
						Gdw = m_v1[row+1][col] * m_v1[row+1][col] * m_h1[row+1][col];
					} else if (m_v1[row+1][col]<0 && m_v1[row][col]<0) {
						Gdw = m_v1[row][col] * m_v1[row][col] * m_h1[row][col];
					} else {
						v_dw = (m_v1[row][col]+m_v1[row+1][col])/2.0;
					   if ((v_dw<0 && m_v1[row][col]==0) || (v_dw>0 && m_v1[row+1][col]==0))
					   	v_dw = 0;
					   if (v_dw>=0) {
					   	h_dw = max(m_h1[row+1][col]+m_z[row+1][col]-m_z[row][col],0);
					   } else {
					   	h_dw = max(m_h1[row][col]+m_z[row][col]-m_z[row+1][col],0);
					   }
						Gdw = h_dw * v_dw * v_dw;
					}


					if(m_DAMBREAK[row-1][col]>0.0 && ((m_h2[row][col]+m_z[row][col]) < (m_h2[row-1][col]+m_z[row-1][col]))){
						Gup = m_h1[row-1][col]* pow((-velocita_breccia(method,m_h1[row-1][col])),2.0); // -0.4 al quadrato perde il segno meno
						if(m_h2[row-1][col]==0)
							Gup=0.0;
					}
					if (m_DAMBREAK[row+1][col]>0.0 && ((m_h2[row][col]+m_z[row][col]) < (m_h2[row+1][col]+m_z[row+1][col]))){
						Gdw = m_h1[row+1][col]* pow((velocita_breccia(method,m_h1[row+1][col])),2.0);
						if(m_h2[row+1][col]==0)
							Gdw=0.0;
					}
					G = Gup - Gdw;


					//courant number  --> UPWIND METHOD
					if (m_v1[row][col]>0 && m_v1[row-1][col]>0 && m_v1[row+1][col]>0){
						dZ_dy_down = ((m_h2[row-1][col] + m_z[row-1][col]) - (m_h2[row][col] + m_z[row][col]) ) / res_ns;
						if (m_h2[row+1][col]==0 &&  m_z[row+1][col]>(m_h2[row][col] + m_z[row][col])) {
							dZ_dy_up = 0;
						} else {
							dZ_dy_up = ((m_h2[row][col] + m_z[row][col]) - (m_h2[row+1][col] + m_z[row+1][col]) ) / res_ns;
						}
						cr_down = (timestep / res_ns) * fabs(m_v1[row-1][col] + m_v1[row][col])/2.0;
						cr_up = (timestep / res_ns) * (fabs(m_v1[row][col]) + fabs(m_v1[row+1][col]))/2.0;
						Z_piu= 0.0;
						Z_meno=0.0;
					} else if (m_v1[row][col]<0 && m_v1[row+1][col]<0 && m_v1[row-1][col]<0){
						dZ_dy_down = ((m_h2[row][col] + m_z[row][col]) - (m_h2[row+1][col] + m_z[row+1][col]) ) / res_ns;
						if (m_h2[row-1][col]==0 &&  m_z[row-1][col]>(m_h2[row][col] + m_z[row][col])) {
							dZ_dy_up = 0;
						} else {
							dZ_dy_up = ((m_h2[row-1][col] + m_z[row-1][col]) - (m_h2[row][col] + m_z[row][col]) ) / res_ns;
						}
						cr_down = (timestep / res_ns) * fabs(m_v1[row][col] + m_v1[row+1][col])/2.0;
						cr_up = (timestep / res_ns) * fabs(m_v1[row-1][col] + m_v1[row][col])/2.0;
						Z_piu= 0.0;
						Z_meno=0.0;
					} else {
						if (m_h2[row-1][col]==0 &&  m_z[row-1][col]>(m_h2[row][col] + m_z[row][col])) {
							Z_piu = (m_h2[row][col] + m_z[row][col]);
						} else {
							Z_piu = ((m_h2[row-1][col] + m_z[row-1][col]) + (m_h2[row][col] + m_z[row][col]) )/2.0;
					   }
					   if (m_h2[row+1][col]==0 &&  m_z[row+1][col]>(m_h2[row][col] + m_z[row][col])) {
							Z_meno = (m_h2[row][col] + m_z[row][col]);
						} else {
							Z_meno = ((m_h2[row][col] + m_z[row][col]) + (m_h2[row+1][col] + m_z[row+1][col]) ) /2.0;
						}
						dZ_dy_down = (Z_piu - Z_meno)/res_ns;
						dZ_dy_up = (Z_piu - Z_meno)/res_ns;
						cr_down = (timestep / res_ns)* (fabs(m_u1[row+1][col]) + 2*fabs(m_u1[row][col]) + fabs(m_u1[row-1][col]) )/4;
						cr_up = (timestep / res_ns)* (fabs(m_u1[row+1][col]) + 2*fabs(m_u1[row][col]) + fabs(m_u1[row-1][col]) )/4;
					}
					if (m_DAMBREAK[row-1][col]>0.0 && m_h2[row-1][col]==0.0){
						dZ_dy_up=0.0;
						dZ_dy_down = ((m_h2[row][col] + m_z[row][col]) - (m_h2[row+1][col] + m_z[row+1][col]) ) / res_ns;
					}
					if (m_DAMBREAK[row+1][col]>0.0 && m_h2[row+1][col]==0.0){
						dZ_dy_up=0.0;
						dZ_dy_down = ((m_h2[row-1][col] + m_z[row-1][col]) - (m_h2[row][col] + m_z[row][col]) ) / res_ns;
					}
						dZ_dy = (1-sqrt(cr_down)) * dZ_dy_down + sqrt(cr_up) * dZ_dy_up;
						//dZ_dy = 0.5 * dZ_dy_down + 0.5 * dZ_dy_up;


					if (m_h1[row][col]<hmin)
						R_i=hmin;
					else
						R_i=m_h1[row][col];

					u = m_u1[row][col];
					v = m_v1[row][col];
					V = sqrt(pow(u,2.0) + pow(v,2.0));
					S = (- g * m_h2[row][col] * dZ_dy) - g*( pow(m_m[row][col],2.0) * v * V / pow(R_i,(1.0/3.0)) );



					if (m_DAMBREAK[row][col] > 0.0 ){
						if ((m_z[row][col]+m_h2[row][col]) >  (m_z[row-1][col] + m_h2[row-1][col]))
						   m_v2[row][col] = velocita_breccia(method,m_h2[row][col]);  // velocita sullo stramazzo
						else if ((m_z[row][col]+m_h2[row][col]) >  (m_z[row+1][col] + m_h2[row+1][col]))
							m_v2[row][col] = -velocita_breccia(method,m_h2[row][col]);  // velocita sullo stramazzo
						else
							m_v2[row][col] = 0.0;
					}else{
						m_v2[row][col] = 1.0 / m_h2[row][col] * (m_h1[row][col] * m_v1[row][col] - timestep / res_ew * F  - timestep / res_ns * G + timestep * S);
						}

	 				// no velocita' contro la diga
					/*if (m_z[row-1][col] > water_elevation && m_v2[row][col] >0)
	 					m_v2[row][col]=0.0;
			  		if (m_z[row+1][col] > water_elevation && m_v2[row][col] < 0 )
	 					m_v2[row][col]=0.0;*/

					if ((timestep/res_ns*(abs(abs(m_v2[row][col])+sqrt(g*m_h2[row][col]))))>1){
						G_warning("\nATTENTION: at time: %f the Courant-Friedrich-Lewy stability condition isn't respected",t);
						/*G_message("EQ. MOTO DIR Y' --> row:%d, col:%d\n)",row, col);
						G_message("m_h1[row][col]:%f,m_u1[row][col]:%f,m_v1[row][col]:%f",m_h1[row][col],m_u1[row][col],m_v1[row][col]);
						G_message("m_h1[row][col+1]:%f,m_h1[row][col-1]:%f,m_h1[row+1][col]:%f, m_h1[row-1][col]:%f\n",m_h1[row][col+1],m_h1[row][col-1],m_h1[row+1][col], m_h1[row-1][col]);
						G_message("h_dx:%f, h_sx:%f, h_up%f, h_dw:%f\n",h_dx, h_sx, h_up, h_dw);
						G_message("m_u1[row][col+1]:%f,m_u1[row][col-1]:%f,m_v1[row+1][col]:%f, m_v1[row-1][col]:%f\n",m_u1[row][col+1],m_u1[row][col-1],m_v1[row+1][col], m_v1[row-1][col]);
						G_message("v_up: %f,v_dw:%f,u_dx:%f,u_sx:%f \n",v_up, v_dw, u_dx, u_sx);
						G_message("Fdx: %f,Fsx: %f, F: %f, Gup:%f, Gdw:%f, G: %f\n",Fdx, Fsx, F,Gup, Gdw, G);
						G_message("m_h2[row][col+1]:%f,m_h2[row][col-1]:%f,m_h2[row+1][col]:%f, m_h2[row-1][col]:%f\n",m_h1[row][col+1],m_h1[row][col-1],m_h1[row+1][col], m_h1[row-1][col]);
						G_message("dZ_dy_down:%f, dZ_dy_up:%f,cr_up:%f, cr_down:%f\n" , dZ_dy_down,dZ_dy_up, cr_up, cr_down);
						G_message("Z_piu:%f,Z_meno:%f\n", Z_piu, Z_meno);
						G_message("dZ_dy:%f,\n",dZ_dy);
						G_message("u:%f,v:%f,V:%f\n", u,v,V);
						G_message("R_i:%f,manning[row][col]:%f\n", R_i, m_m[row][col]);
						G_message("S=%f\n",S);
						G_message("m_v2(row,col): %f\n \n", m_v2[row][col]);
						G_warning("   ");*/
					}



					//************** stampa  ********************************************************
					//if ((t>6.8 && m_v2[row][col]!=m_v1[row][col]) && (row==87) && (col == 193)) {
					/*if (fabs(m_v2[row][col])>=1000.0){
						G_warning("At the time %f v(%d,%d)=%f", t, row,col,m_v2[row][col]);
					}*/



			 } else {
			   // tolgo h<hmin quando si svuota
			   m_u2[row][col] = 0.0;
			   m_v2[row][col] = 0.0;
			 }
	    	/* close the cicle */
			}
	}
         

//*************************************** overwriting *********************************************
    timestep_ct=0;
    if (t<TSTOP){
    /* open new cicle */
    	for (row = 1; row < nrows-1; row++)
      	{
	   	for (col = 1; col < ncols-1; col++)
			 	{
                                //********************************************************************				
				// timestep optimization using the CFL stability condition 
				if (m_h2[row][col]>=hmin ){
                         		timestep_ct_temp = max( (fabs(m_u2[row][col])+sqrt(g*m_h2[row][col]))/res_ew , (fabs(m_v2[row][col])+sqrt(g*m_h2[row][col]))/res_ns );
					if(timestep_ct_temp>timestep_ct){
						timestep_ct=timestep_ct_temp;
						//G_message("t=%f,row=%d,col=%d,timestep_ct=%f,m_u2=%f m_v2=%f m_h2=%f",t,row,col,timestep_ct,m_u2[row][col],m_v2[row][col],m_h2[row][col]);
					}
				}
				
                                //********************************************************************

			 	m_u1[row][col] =  m_u2[row][col];
				m_v1[row][col] =  m_v2[row][col];
				if (OUT_HMAX) {
					if (m_hmax[row][col]<m_h2[row][col]){
						m_hmax[row][col]=m_h2[row][col];
						m_t_hmax[row][col]=t;
						velocity=sqrt(pow(m_u1[row][col],2.0) + pow(m_v1[row][col],2.0));
						m_i_hmax[row][col]=velocity*m_h2[row][col];
					}
				}
				if (OUT_VMAX) {
					velocity=sqrt(pow(m_u1[row][col],2.0) + pow(m_v1[row][col],2.0));
					if (m_vmax[row][col]<velocity){
						m_vmax[row][col]=velocity;
						m_i_vmax[row][col]=velocity*m_h2[row][col];
						m_t_vmax[row][col]=t;
						if (flag_d->answer) {
						   if (m_u1[row][col]==0 && m_v1[row][col]==0){
								m_dir_vmax[row][col] =0;
							} else if (m_u1[row][col]>0 && m_v1[row][col]>0) {
								m_dir_vmax[row][col] = 180/PI*atan(fabs(m_v1[row][col])/fabs(m_u1[row][col]));
							} else if (m_u1[row][col]==0 && m_v1[row][col]>0) {
								m_dir_vmax[row][col] = 90.0;
							} else if (m_u1[row][col]<0 && m_v1[row][col]>0) {
								m_dir_vmax[row][col] = 90.0 + 180/PI*atan(fabs(m_u1[row][col])/fabs(m_v1[row][col]));
							} else if (m_u1[row][col]<0 && m_v1[row][col]==0) {
								m_dir_vmax[row][col] = 180.0;
							} else if (m_u1[row][col]<0 && m_v1[row][col]<0) {
								m_dir_vmax[row][col] = 180.0 + 180/PI*atan(fabs(m_v1[row][col])/fabs(m_u1[row][col]));
							} else if (m_u1[row][col]==0 && m_v1[row][col]<0) {
								m_dir_vmax[row][col] = 270.0;
							} else if (m_u1[row][col]>0 && m_v1[row][col]<0) {
								m_dir_vmax[row][col] = 270.0 + 180/PI*atan(fabs(m_u1[row][col])/fabs(m_v1[row][col]));
							} else if (m_u1[row][col]>0 && m_v1[row][col]==0) {
								m_dir_vmax[row][col] = 0.0;
							}
						}
					}
				}
				if (OUT_IMAX) {
					velocity=sqrt(pow(m_u1[row][col],2.0) + pow(m_v1[row][col],2.0));
					i=velocity*m_h2[row][col];
					if (m_imax[row][col]<i){
						m_imax[row][col]=i;
						m_t_imax[row][col]=t;
					}
				}
                                if (OUT_WAVEFRONT) {
					if (m_wavefront[row][col]==0.0 && m_h2[row][col]>hmin){
						m_wavefront[row][col]=t;
					}
				}

				m_h1[row][col] =  m_h2[row][col];
	}}}
	
	//******************************   new timestep   ******************************************
        timestep=0.1/timestep_ct;
        //******************************************************************************************

        //G_message("timestep =%f,t=%f",timestep,t);
	
	/*if (input_DELTAT->answer != NULL){
	G_message("pippo");
	}*/
   	//*****************************************************************************************
        /* if T e' piu' o meno multiplo di DELTAT allora scrivi output con istante */
     	//if ((m*DELTAT-t) <= timestep && m*DELTAT < TSTOP || (pp<=ntimes && times[pp] <= timestep)) {
     	//if ( ((m*DELTAT-t) <= timestep && m*DELTAT < TSTOP || ( ( pp<ntimes && (times[pp]-t) < timestep) && (parm.opt_t->answer != NULL) )) && (input_DELTAT->answer != NULL) ) 
	if (( ((m*DELTAT-t) <= timestep) && (m*DELTAT < TSTOP)&& (input_DELTAT->answer != NULL)) || ( (pp<ntimes && (times[pp]-t) < timestep) && (parm.opt_t->answer != NULL))) {
		//G_message("pippo");
     	 	if ((m*DELTAT-t) <= timestep && m*DELTAT < TSTOP) {
				/* devo cambiare il nome del raster e aggiungere ogni volta _timestep*/
				if (OUT_H) {
					sprintf(name1,"%s%d",OUT_H,m*DELTAT);
				}
				if (OUT_VEL){
					sprintf(name2,"%s%d",OUT_VEL,m*DELTAT);
					sprintf(name3,"%s%s%d","dir_",OUT_VEL,m*DELTAT);
				}
				G_message("Time: %d, writing the output maps",m*DELTAT);
				m++;
			} else {
				pp++;
				sprintf(name1,"%s%s%d","opt_",OUT_H,pp);
				sprintf(name2,"%s%s%d","opt_",OUT_VEL,pp);
				sprintf(name3,"%s%s%s%d","opt_","dir_",OUT_VEL,pp);
				G_message("Time: %lf, writing an optional output maps %d",t, pp);
			}
			if (OUT_H) {
				if ( (outfd_H = G_open_raster_new (name1,DCELL_TYPE)) < 0)
					G_fatal_error (_("Could not open <%s>"),name1);
			}
			if (OUT_VEL) {
				if ( (outfd_VEL = G_open_raster_new (name2,DCELL_TYPE)) < 0)
					G_fatal_error (_("Could not open <%s>"),name2);
			}
			if (flag_d->answer) {
				if ( (outfd_VEL_DIR = G_open_raster_new (name3,DCELL_TYPE)) < 0)
					G_fatal_error (_("Could not open <%s>"),name3);
			}
			/* allocate output buffer */
			if (OUT_VEL) {
				outrast_VEL = G_allocate_d_raster_buf();
			}
			if (OUT_H) {
				outrast_H = G_allocate_d_raster_buf();
			}
			if (flag_d->answer) {
				outrast_VEL_DIR = G_allocate_d_raster_buf();
			}
		 	for (row = 0; row < nrows; row++) {
			 	for (col = 0; col < ncols; col++) {
					/* copy matrix in buffer */
					if (OUT_VEL) {
						((DCELL *) outrast_VEL)[col] = sqrt(pow(m_u1[row][col],2.0) + pow(m_v1[row][col],2.0));
					}
					if (flag_d->answer) {
						if (m_u1[row][col]==0 && m_v1[row][col]==0){
							G_set_d_null_value(&outrast_VEL_DIR[col],1);
						} else if (m_u1[row][col]>0 && m_v1[row][col]>0) {
							((DCELL *) outrast_VEL_DIR)[col] = 180/PI*atan(fabs(m_v1[row][col])/fabs(m_u1[row][col]));
						} else if (m_u1[row][col]==0 && m_v1[row][col]>0) {
							((DCELL *) outrast_VEL_DIR)[col] = 90.0;
						} else if (m_u1[row][col]<0 && m_v1[row][col]>0) {
							((DCELL *) outrast_VEL_DIR)[col] = 90.0 + 180/PI*atan(fabs(m_u1[row][col])/fabs(m_v1[row][col]));
						} else if (m_u1[row][col]<0 && m_v1[row][col]==0) {
							((DCELL *) outrast_VEL_DIR)[col] = 180.0;
						} else if (m_u1[row][col]<0 && m_v1[row][col]<0) {
							((DCELL *) outrast_VEL_DIR)[col] = 180.0 + 180/PI*atan(fabs(m_v1[row][col])/fabs(m_u1[row][col]));
						} else if (m_u1[row][col]==0 && m_v1[row][col]<0) {
							((DCELL *) outrast_VEL_DIR)[col] = 270.0;
						} else if (m_u1[row][col]>0 && m_v1[row][col]<0) {
							((DCELL *) outrast_VEL_DIR)[col] = 270.0 + 180/PI*atan(fabs(m_u1[row][col])/fabs(m_v1[row][col]));
						} else if (m_u1[row][col]>0 && m_v1[row][col]==0) {
							((DCELL *) outrast_VEL_DIR)[col] = 0.0;
						}
					}
					if (OUT_H) {
						((DCELL *) outrast_H)[col] = m_h1[row][col];
					}

				} /* end_col*/
				 /*copia righe !!! */
				 if (OUT_VEL) {
				 	G_put_d_raster_row(outfd_VEL,outrast_VEL);
				 }
				 if (OUT_H) {
				 	G_put_d_raster_row(outfd_H,outrast_H);
				 }
				 if (flag_d->answer) {
				 	G_put_d_raster_row(outfd_VEL_DIR,outrast_VEL_DIR);
				 }
			  }	// end row
				if (OUT_VEL) {
					G_free(outrast_VEL);
				}
				if (OUT_H) {
					G_free(outrast_H);
				}
				if (flag_d->answer) {
					G_free(outrast_VEL_DIR);
				}
				/* chiudi i file */
				if (OUT_VEL) {
					G_close_cell (outfd_VEL);
				}
		  		if (OUT_H) {
		  			G_close_cell (outfd_H);
		  		}
		  		if (flag_d->answer) {
		  			G_close_cell (outfd_VEL_DIR);
		  		}
			} // end if
        //G_message("timestep =%f,t=%f",timestep,t);
	} // end time loop

//*******************************************************************
// write final flooding map
if(OUT_H) {
	sprintf(name1,"%s%d",OUT_H,TSTOP);
}
if(OUT_VEL) {
	sprintf(name2,"%s%d",OUT_VEL,TSTOP);
	sprintf(name3,"%s%s%d","dir_",OUT_VEL,TSTOP);
}
if(OUT_H) {
	if ( (outfd_H = G_open_raster_new (name1,DCELL_TYPE)) < 0)
		G_fatal_error (_("Could not open <%s>"),name1);
}
if(OUT_VEL) {
	if ( (outfd_VEL = G_open_raster_new (name2,DCELL_TYPE)) < 0)
		G_fatal_error (_("Could not open <%s>"),name2);
}
if (flag_d->answer) {
	if ( (outfd_VEL_DIR = G_open_raster_new (name3,DCELL_TYPE)) < 0)
		G_fatal_error (_("Could not open <%s>"),name3);
}

if (OUT_HMAX){
	sprintf(name4,"%s%s",OUT_HMAX,"_time");
	sprintf(name5,"%s%s",OUT_HMAX,"_intensity");
	if ( (outfd_HMAX = G_open_raster_new (OUT_HMAX,DCELL_TYPE)) < 0)
		G_fatal_error (_("Could not open <%s>"),OUT_HMAX);
	if ( (outfd_T_HMAX = G_open_raster_new (name4,DCELL_TYPE)) < 0)
		G_fatal_error (_("Could not open <%s>"),OUT_HMAX);
	if ( (outfd_I_HMAX = G_open_raster_new (name5,DCELL_TYPE)) < 0)
		G_fatal_error (_("Could not open <%s>"),OUT_HMAX);
}

if (OUT_VMAX){
	sprintf(name6,"%s%s",OUT_VMAX,"_time");
	sprintf(name7,"%s%s",OUT_VMAX,"_intensity");
	if (flag_d->answer) {
		sprintf(name8,"%s%s",OUT_VMAX,"_dir");
	}
	if ( (outfd_VMAX = G_open_raster_new (OUT_VMAX,DCELL_TYPE)) < 0)
		G_fatal_error (_("Could not open <%s>"),OUT_VMAX);
	if ( (outfd_T_VMAX = G_open_raster_new (name6,DCELL_TYPE)) < 0)
		G_fatal_error (_("Could not open <%s>"),OUT_VMAX);
	if ( (outfd_I_VMAX = G_open_raster_new (name7,DCELL_TYPE)) < 0)
		G_fatal_error (_("Could not open <%s>"),OUT_VMAX);
	if (flag_d->answer) {
		if ( (outfd_DIR_VMAX = G_open_raster_new (name8,DCELL_TYPE)) < 0)
			G_fatal_error (_("Could not open <%s>"),OUT_VMAX);
	}
}
if (OUT_IMAX){
	sprintf(name9,"%s%s",OUT_IMAX,"_time");
	if ( (outfd_IMAX = G_open_raster_new (OUT_IMAX,DCELL_TYPE)) < 0)
		G_fatal_error (_("Could not open <%s>"),OUT_IMAX);
	if ( (outfd_T_IMAX = G_open_raster_new (name9,DCELL_TYPE)) < 0)
		G_fatal_error (_("Could not open <%s>"),OUT_IMAX);
}
if (OUT_WAVEFRONT){
	if ( (outfd_WAVEFRONT = G_open_raster_new (OUT_WAVEFRONT,DCELL_TYPE)) < 0)
		G_fatal_error (_("Could not open <%s>"),OUT_WAVEFRONT);
}

/* allocate output buffer */
if (OUT_H) {
	outrast_H = G_allocate_d_raster_buf();
}
if (OUT_VEL) {
	outrast_VEL = G_allocate_d_raster_buf();
}
if(flag_d->answer) {
	outrast_VEL_DIR = G_allocate_d_raster_buf();
}

if (OUT_HMAX){
	outrast_HMAX = G_allocate_d_raster_buf();
	outrast_I_HMAX = G_allocate_d_raster_buf();
	outrast_T_HMAX = G_allocate_d_raster_buf();
}
if (OUT_VMAX){
	outrast_VMAX = G_allocate_d_raster_buf();
	outrast_I_VMAX = G_allocate_d_raster_buf();
	outrast_T_VMAX = G_allocate_d_raster_buf();
	if(flag_d->answer) {
		outrast_DIR_VMAX = G_allocate_d_raster_buf();
	}
}

if (OUT_IMAX){
	outrast_IMAX = G_allocate_d_raster_buf();
	outrast_T_IMAX = G_allocate_d_raster_buf();
}
if (OUT_WAVEFRONT){
	outrast_WAVEFRONT = G_allocate_d_raster_buf();
}

G_percent(nrows, nrows, 1);	/* finish it */

for (row = 0; row < nrows; row++){
   G_percent (row, nrows, 2);
	for (col = 0; col < ncols; col++) {
		/* copy matrix in buffer */
		if (OUT_H) {
			((DCELL *) outrast_H)[col] = m_h1[row][col];
		}
		if (OUT_VEL) {
			((DCELL *) outrast_VEL)[col] = sqrt(pow(m_u1[row][col],2.0) + pow(m_v1[row][col],2.0));
		}
		if (flag_d->answer) {
			if (m_u1[row][col]==0 && m_v1[row][col]==0){
				G_set_d_null_value(&outrast_VEL_DIR[col],1);
			} else if (m_u1[row][col]>0 && m_v1[row][col]>0) {
				((DCELL *) outrast_VEL_DIR)[col] = 180/PI*atan(fabs(m_v1[row][col])/fabs(m_u1[row][col]));
			} else if (m_u1[row][col]==0 && m_v1[row][col]>0) {
				((DCELL *) outrast_VEL_DIR)[col] = 90.0;
			} else if (m_u1[row][col]<0 && m_v1[row][col]>0) {
				((DCELL *) outrast_VEL_DIR)[col] = 90.0 + 180/PI*atan(fabs(m_u1[row][col])/fabs(m_v1[row][col]));
			} else if (m_u1[row][col]<0 && m_v1[row][col]==0) {
				((DCELL *) outrast_VEL_DIR)[col] = 180.0;
			} else if (m_u1[row][col]<0 && m_v1[row][col]<0) {
				((DCELL *) outrast_VEL_DIR)[col] = 180.0 + 180/PI*atan(fabs(m_v1[row][col])/fabs(m_u1[row][col]));
			} else if (m_u1[row][col]==0 && m_v1[row][col]<0) {
				((DCELL *) outrast_VEL_DIR)[col] = 270.0;
			} else if (m_u1[row][col]>0 && m_v1[row][col]<0) {
				((DCELL *) outrast_VEL_DIR)[col] = 270.0 + 180/PI*atan(fabs(m_u1[row][col])/fabs(m_v1[row][col]));
			} else if (m_u1[row][col]>0 && m_v1[row][col]==0) {
				((DCELL *) outrast_VEL_DIR)[col] = 0.0;
			}
		}

		// output HMAX
		if (OUT_HMAX){
			if(m_hmax[row][col]==0){
				G_set_d_null_value(&outrast_HMAX[col], 1);
			} else {
				((DCELL *) outrast_HMAX)[col] = m_hmax[row][col];
			}
			if(m_hmax[row][col]==0){
				G_set_d_null_value(&outrast_T_HMAX[col], 1);
			} else {
				((DCELL *) outrast_T_HMAX)[col] = m_t_hmax[row][col];
			}
			if(m_hmax[row][col]==0){
				G_set_d_null_value(&outrast_I_HMAX[col], 1);
			} else {
				((DCELL *) outrast_I_HMAX)[col] = m_i_hmax[row][col];
			}
		}
		// output VMAX
		if (OUT_VMAX){
			if(m_vmax[row][col]==0){
				G_set_d_null_value(&outrast_VMAX[col], 1);
			} else {
				((DCELL *) outrast_VMAX)[col] = m_vmax[row][col];
			}
		   if(m_vmax[row][col]==0){
				G_set_d_null_value(&outrast_I_VMAX[col], 1);
			} else {
				((DCELL *) outrast_I_VMAX)[col] = m_i_vmax[row][col];
			}
			if(m_vmax[row][col]==0){
				G_set_d_null_value(&outrast_T_VMAX[col], 1);
			} else {
				((DCELL *) outrast_T_VMAX)[col] = m_t_vmax[row][col];
			}
			if (flag_d->answer) {
				if(m_vmax[row][col]==0){
					G_set_d_null_value(&outrast_DIR_VMAX[col], 1);
				} else {
					((DCELL *) outrast_DIR_VMAX)[col] = m_dir_vmax[row][col];
				}
			}
		}
		// output IMAX
		if (OUT_IMAX){
			if(m_imax[row][col]==0){
				G_set_d_null_value(&outrast_IMAX[col], 1);
			} else {
				((DCELL *) outrast_IMAX)[col] = m_imax[row][col];
			}
		   if(m_imax[row][col]==0){
				G_set_d_null_value(&outrast_T_IMAX[col], 1);
			} else {
				((DCELL *) outrast_T_IMAX)[col] = m_t_imax[row][col];
			}
		}
		// output WAVEFRONT
		if (OUT_WAVEFRONT){
			if(m_wavefront[row][col]==0){
				G_set_d_null_value(&outrast_WAVEFRONT[col], 1);
			} else {
				((DCELL *) outrast_WAVEFRONT)[col] = m_wavefront[row][col];
			}
		}

	} /* end_col*/
	if(OUT_H) {
 		if (G_put_raster_row (outfd_H, outrast_H, TYPE_DOUBLE) < 0)
			G_fatal_error (_("Cannot write to <%s>"),name2);
	}
	if (OUT_VEL) {
		if (G_put_raster_row (outfd_VEL,outrast_VEL, TYPE_DOUBLE) < 0)
			G_fatal_error (_("Cannot write to <%s>"),name1);
 	}
	if (flag_d->answer) {
		if (G_put_raster_row (outfd_VEL_DIR, outrast_VEL_DIR, TYPE_DOUBLE) < 0)
			G_fatal_error (_("Cannot write to <%s>"),name3);
	}

	if (OUT_HMAX){
		if (G_put_raster_row (outfd_HMAX, outrast_HMAX, TYPE_DOUBLE) < 0)
			G_fatal_error (_("Cannot write to <%s>"),OUT_HMAX);
		if (G_put_raster_row (outfd_T_HMAX, outrast_T_HMAX, TYPE_DOUBLE) < 0)
			G_fatal_error (_("Cannot write to <%s>"),OUT_HMAX);
		if (G_put_raster_row (outfd_I_HMAX, outrast_I_HMAX, TYPE_DOUBLE) < 0)
			G_fatal_error (_("Cannot write to <%s>"),OUT_HMAX);
	}

	if (OUT_VMAX){
		if (G_put_raster_row (outfd_VMAX, outrast_VMAX, TYPE_DOUBLE) < 0)
			G_fatal_error (_("Cannot write to <%s>"),OUT_VMAX);
		if (G_put_raster_row (outfd_T_VMAX, outrast_T_VMAX, TYPE_DOUBLE) < 0)
			G_fatal_error (_("Cannot write to <%s>"),OUT_VMAX);
		if (G_put_raster_row (outfd_I_VMAX, outrast_I_VMAX, TYPE_DOUBLE) < 0)
			G_fatal_error (_("Cannot write to <%s>"),OUT_VMAX);
		if (flag_d->answer) {
			if (G_put_raster_row (outfd_DIR_VMAX, outrast_DIR_VMAX, TYPE_DOUBLE) < 0)
				G_fatal_error (_("Cannot write to <%s>"),OUT_VMAX);
		}
	}

	if (OUT_IMAX){
		if (G_put_raster_row (outfd_IMAX, outrast_IMAX, TYPE_DOUBLE) < 0)
			G_fatal_error (_("Cannot write to <%s>"),OUT_IMAX);
		if (G_put_raster_row (outfd_T_IMAX, outrast_T_IMAX, TYPE_DOUBLE) < 0)
			G_fatal_error (_("Cannot write to <%s>"),OUT_IMAX);
	}
        
        if (OUT_WAVEFRONT){
		if (G_put_raster_row (outfd_WAVEFRONT, outrast_WAVEFRONT, TYPE_DOUBLE) < 0)
			G_fatal_error (_("Cannot write to <%s>"),OUT_WAVEFRONT);
	}
 	//G_message("pippo, row=%d e nrows=%d", row, nrows);
}	// end row
/* chiudi i file */


  	if (OUT_H) {
  		G_message("Writing the output final map %s", name1);
  	}
  	if (OUT_VEL) {
    		G_message("Writing the output final map %s",name2);
  	}
	if (OUT_HMAX) {
  		G_message("Writing the output final map %s, corresponding time and intensity", OUT_HMAX);
  	}
	if (OUT_VMAX) {
  		G_message("Writing the output final map %s, corresponding time and intensity", OUT_VMAX);
  	}
	if (OUT_IMAX) {
  		G_message("Writing the output final map %s and corresponding time", OUT_IMAX);
  	}
	if (OUT_WAVEFRONT) {
  		G_message("Writing the output final map %s", OUT_WAVEFRONT);
  	}

G_percent(nrows, nrows, 1);	/* finish it */
if (OUT_VEL) {
	G_free(outrast_VEL);
	G_close_cell (outfd_VEL);
}
if (OUT_H) {
	G_free(outrast_H);
	G_close_cell (outfd_H);
}
if (flag_d->answer) {
	G_free(outrast_VEL_DIR);
	G_close_cell (outfd_VEL_DIR);
}

if (OUT_HMAX){
	G_free(outrast_HMAX);
	G_close_cell (outfd_HMAX);
	G_free(outrast_I_HMAX);
	G_close_cell (outfd_I_HMAX);
	G_free(outrast_T_HMAX);
	G_close_cell (outfd_T_HMAX);
}
if (OUT_VMAX){
	G_free(outrast_VMAX);
	G_close_cell (outfd_VMAX);
	G_free(outrast_I_VMAX);
	G_close_cell (outfd_I_VMAX);
	G_free(outrast_T_VMAX);
	G_close_cell (outfd_T_VMAX);
	if (flag_d->answer) {
		G_free(outrast_DIR_VMAX);
		G_close_cell (outfd_DIR_VMAX);
	}
}
if (OUT_IMAX){
	G_free(outrast_IMAX);
	G_close_cell (outfd_IMAX);
	G_free(outrast_T_IMAX);
	G_close_cell (outfd_T_IMAX);
}
if (OUT_WAVEFRONT){
	G_free(outrast_WAVEFRONT);
	G_close_cell (outfd_WAVEFRONT);
}

//************************************************************************
// da sistemare
//************************************************************************
/* add command line incantation to history file */
//G_short_history(result, "raster", &history);
//G_command_history(&history);
//G_write_history(result, &history);


/* deallocate memory matrix */
G_free_fmatrix(m_DAMBREAK);
G_free_fmatrix(m_m);
G_free_fmatrix(m_z);
G_free_dmatrix(m_h1);
G_free_dmatrix(m_h2);
G_free_dmatrix(m_u1);
G_free_dmatrix(m_u2);
G_free_dmatrix(m_v1);
G_free_dmatrix(m_v2);
G_free_imatrix(m_lake);
if( parm.opt_t->answer != NULL){
	G_free_vector(times);pp=0;
}
if (OUT_HMAX){
	G_free_fmatrix(m_hmax);
	G_free_fmatrix(m_t_hmax);
	G_free_fmatrix(m_i_hmax);
}
if (OUT_VMAX){
	G_free_fmatrix(m_vmax);
	G_free_fmatrix(m_i_vmax);
	G_free_fmatrix(m_t_vmax);
	if (flag_d->answer) {
		G_free_fmatrix(m_dir_vmax);
	}
}
if (OUT_IMAX){
	G_free_fmatrix(m_imax);
	G_free_fmatrix(m_t_imax);
}
if (OUT_WAVEFRONT){
	G_free_fmatrix(m_wavefront);
}




} //END MAIN.C







