/****************************************************************************
 *
 * MODULE:       r.dfwalk
 * AUTHOR(S):    Massimiliano Cannata - massimiliano.cannata supsi.ch
 *              
 * PURPOSE:      Calculate debris-flow hazard map
 *
 * COPYRIGHT:    (C) 2002,2007 by the GRASS Development Team
 *
 *               This program is free software under the GNU General Public
 *   	    	 License (>=v2). Read the file COPYING that comes with GRASS
 *   	    	 for details.
 *
 *********************************************************************************/
#define D_SLIM  "0.36" /* pendenza limite Grenzgefaelle - Slope limit value [tan]  */
#define D_ESP   "1.8"  /* coefficiente di espandimento - expansion coefficient [-] */
#define D_NRW   "100"  /* number of random walk for each starting cell [-]         */
#define D_VOL   "2000" /* total volume   [m3]         */
#define D_SMEAN "0.1" /* pendenza complessiva Pauschälgefalle - mean slope of conoid [tan]*/
#define D_PMD   "70"  /* perla mass-to-drag rough coefficient [m] */
#define D_PS    "0.1"  /* perla slippering rough coefficient [-] */
#define D_Sed_VL    "14"  /* maximun velocity at wich there is deposition [m/s] */
#define D_Sed_SL    "0.4"  /* maximun slope at wich there is deposition [tan] */
#define D_Sed_SHM   "1.0"  /* maximum deposition height in a single cell for each RW due to slope  [m] */
#define D_Sed_VHM   "1.0"  /* maximum deposition height in a single cell for each RW due to velocity  [m] */
#define D_Sed_SHM0  "20"  /* maximum volume deposition possible per random walk [m3] */
#define BACK 50

#define LENGHT 2000
#define START_DIM 10

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <grass/gis.h>
#include <grass/glocale.h>
#include "nrutil.h"
#include "dfFunctions.h"
#include "malloc.h"



    /* determine the inputmap type -> G_raster_map_type
     * return file descriptors -> G_open_cell_old
     * control if we can open input rasters -> G_get_cellhd */



/*
 * main function
 * it copies raster input raster map, calling the appropriate function for each
 * data type (CELL, DCELL, FCELL)
 */
int main(int argc, char *argv[])
{

    /* input raster structure as defined in dfFunctions.h */
    raster_struct S_dem, S_start, S_pdir; /* S_obj */
    /* output raster structure as defined in dfFunctions.h */
    raster_struct S_nrRW, S_vel, S_sed;
    /* arrays of structures as defined in dfFunctions.h */
    rw_point *rw_start, *tmp_start, cell_to;
    rw_cell *rw_cells, *tmp_rw;
    rw_mask mask;
    /* counters and supporting variables */
    int n_rw_cells;
    int ia, ri, n_start, n_alloc_start;
    int n_alloc_rw, n_rw, i_rw;
    int m, n, nn, r, c, count, stop, lenght;
    int row_to, col_to;
    int DirMode;

    int nrows, ncols;
    int row, col;

    /* flag descriptors */		
    int verbose, maxvel, number_rw;

    /* GRASS module for parsing arguments */
    struct GModule *module;	

    /* options */
    struct Option *input_dem, *input_start, *input_pdir, *input_slope, *input_aspect; /*input raster *input_obj */
    struct Option *input_slim, *input_esp, *input_vol, *input_rw, *input_smean; /* coefficient */
    struct Option *input_perla_m2d, *input_perla_rough; /* perla coefficients */
    struct Option *input_sed_Vlim, *input_sed_Slim,*input_sed_Shmax, *input_sed_Vhmax,*input_sed_0hmax ; /* tresholds */
    struct Option *output_nrRW, *output_vel, *output_sed; /*output raster*/

    /* input options parameters */
    int rw;    
    double slim, esp, smean, vol_tot;
    double perla_m2d, perla_rough, sed_Vlim, sed_Slim, sed_Shmax, sed_Vhmax, sed_0hmax;

    /* flags */
    struct Flag *flag1, *flag2, *flag3;		

    /* initialize GIS environment */
    G_gisinit(argv[0]);		/* reads grass env, stores program name to G_program_name() */

    /* initialize module */
    module = G_define_module();
    module->keywords = _("raster, debris-flow, hazard");
    module->description = _("Calculate debris-flow hazard map");

    /* Define the different options as defined in gis.h */
    
    /* INPUT RASTERS */
    input_dem              = G_define_standard_option(G_OPT_R_ELEV);
    input_dem->key	   = "dsm";
    input_dem->type        = TYPE_STRING;
    input_dem->required    = YES;
    input_dem->guisection  = _("Input_options");
    input_dem->description = _("Input dem raster");
    
    input_start              = G_define_standard_option(G_OPT_R_INPUT);
    input_start->key         = "start";
    input_start->type        = TYPE_STRING;
    input_start->required    = YES;
    input_start->guisection  = _("Input_options");
    input_start->description = _("Input zone file");    
        
    input_pdir              = G_define_standard_option(G_OPT_R_INPUT);
    input_pdir->key	    = "pdir";
    input_pdir->type        = TYPE_STRING;
    input_pdir->required    = NO;
    input_pdir->guisection  = _("Input_options");
    input_pdir->description = _("Input raster preferential directions");

    /* PARAMETERS */
    input_slim              = G_define_option();
    input_slim->key	    = "slim";
    input_slim->type        = TYPE_DOUBLE;
    input_slim->answer      = D_SLIM;
    input_slim->required    = NO;
    input_slim->guisection  = _("Input_options");
    input_slim->description = _("Slope limit value [%]");
	
    input_esp              = G_define_option();
    input_esp->key	   = "esp";
    input_esp->type        = TYPE_DOUBLE;
    input_esp->answer      = D_ESP;
    input_esp->required    = NO;
    input_esp->guisection  = _("Input_options");
    input_esp->description = _("expansion coefficient [-]");

    input_vol              = G_define_option();
    input_vol->key	  = "vol";
    input_vol->type        = TYPE_INTEGER;
    input_vol->answer      = D_VOL;
    input_vol->required    = NO;
    input_vol->guisection  = _("Input_options");
    input_vol->description = _("Total volumel [m3]");
    
    input_rw              = G_define_option();
    input_rw->key	  = "nrw";
    input_rw->type        = TYPE_INTEGER;
    input_rw->answer      = D_NRW;
    input_rw->required    = NO;
    input_rw->guisection  = _("Input_options");
    input_rw->description = _("number of random walk for each starting cell [-]");

    input_smean              = G_define_option();
    input_smean->key	     = "smean";
    input_smean->type        = TYPE_DOUBLE;
    input_smean->answer      = D_SMEAN;
    input_smean->required    = NO;
    input_smean->guisection  = _("Input_options");
    input_smean->description = _("mean slope of conoid [%]");

    input_perla_m2d              = G_define_option();
    input_perla_m2d->key	 = "pmd";
    input_perla_m2d->type        = TYPE_DOUBLE;
    input_perla_m2d->answer      = D_PMD;
    input_perla_m2d->required    = NO;
    input_perla_m2d->guisection  = _("Input_options");
    input_perla_m2d->description = _("perla mass-to-drag rough coefficient [m]");

    input_perla_rough              = G_define_option();
    input_perla_rough->key	   = "ps";
    input_perla_rough->type        = TYPE_DOUBLE;
    input_perla_rough->answer      = D_PS;
    input_perla_rough->required    = NO;
    input_perla_rough->guisection  = _("Input_options");
    input_perla_rough->description = _("perla slippering rough coefficient [-]");

    input_sed_Vlim		= G_define_option();
    input_sed_Vlim->key	        = "vl";
    input_sed_Vlim->type        = TYPE_DOUBLE;
    input_sed_Vlim->answer      = D_Sed_VL;
    input_sed_Vlim->required    = NO;
    input_sed_Vlim->guisection  = _("Input_options");
    input_sed_Vlim->description = _("maximun speed at wich there is deposition [m/s]");

    input_sed_Slim    		= G_define_option();
    input_sed_Slim->key         = "sl";
    input_sed_Slim->type        = TYPE_DOUBLE;
    input_sed_Slim->answer      = D_Sed_SL;
    input_sed_Slim->required    = NO;
    input_sed_Slim->guisection  = _("Input_options");
    input_sed_Slim->description = _("maximun slope at wich there is deposition [decimal value]");

    input_sed_Shmax   		 = G_define_option();
    input_sed_Shmax->key	 = "shm";
    input_sed_Shmax->type        = TYPE_DOUBLE;
    input_sed_Shmax->answer      = D_Sed_SHM;
    input_sed_Shmax->required    = NO;
    input_sed_Shmax->guisection  = _("Input_options");
    input_sed_Shmax->description = _("max deposition possible due to slope in a single cell for a RW [m]");
    
    input_sed_Vhmax              = G_define_option();
    input_sed_Vhmax->key	 = "vhm";
    input_sed_Vhmax->type        = TYPE_DOUBLE;
    input_sed_Vhmax->answer      = D_Sed_VHM;
    input_sed_Vhmax->required    = NO;
    input_sed_Vhmax->guisection  = _("Input_options");
    input_sed_Vhmax->description = _("max deposition possible due to speed in a single cell for a RW [m]");
    
    input_sed_0hmax              = G_define_option();
    input_sed_0hmax->key	 = "sh0m";
    input_sed_0hmax->type        = TYPE_DOUBLE;
    input_sed_0hmax->answer      = D_Sed_SHM0;
    input_sed_0hmax->required    = NO;
    input_sed_0hmax->guisection  = _("Input_options");
    input_sed_0hmax->description = _("max volume deposition possible per random walk [m3]");

    /* OUTPUT RASTERS */
    output_nrRW              = G_define_standard_option(G_OPT_R_OUTPUT);
    output_nrRW->key	     = "nrrw";
    output_nrRW->type        = TYPE_STRING;
    output_nrRW->required    = YES;
    output_nrRW->guisection  = _("Output_options");
    output_nrRW->description = _("map random walks");
    
    output_vel              = G_define_standard_option(G_OPT_R_OUTPUT);
    output_vel->key	    = "vel";
    output_vel->type        = TYPE_STRING;
    output_vel->required    = YES;
    output_vel->guisection  = _("Output_options");
    output_vel->description = _("raster max velocities");

    output_sed              = G_define_standard_option(G_OPT_R_OUTPUT);
    output_sed->key	    = "hgt";
    output_sed->type        = TYPE_STRING;
    output_sed->required    = YES;
    output_sed->guisection  = _("Output_options");
    output_sed->description = _("raster sedimentation");

    /* Define the different flags */
    flag1 = G_define_flag();
    flag1->key = 'q';
    flag1->description = _("Quiet");
    
    flag2 = G_define_flag();
    flag2->key = 'm';
    flag2->description = _("Calculate max velocity instead of mean");

    flag3 = G_define_flag();
    flag3->key = 'n';
    flag3->description = _("Use a number of random walk instead of volume");


    /* options and flags parser */
    if (G_parser(argc, argv))
	exit(EXIT_FAILURE);

	
	//back = atof(BACK);
    /**************************************/
    /* stores options and flags to variables   */
    /**************************************/ 
    /* input raster maps */
    S_dem.name   = input_dem->answer;             
    S_start.name = input_start->answer;
    S_pdir.name  = input_pdir->answer;
     
    /* treshold parameters */
    slim = atof(input_slim->answer);
    esp = atof(input_esp->answer);
    smean = atof(input_smean->answer);
    rw = atoi(input_rw->answer);
    vol_tot = atof(input_vol->answer);
    
    
    /* perla model parameters */
    perla_m2d = atof(input_perla_m2d->answer);
    perla_rough = atof(input_perla_rough->answer);
    //sscanf(input_perla_m2d->answer, "%lf", &(perla_m2d));
    //sscanf(input_perla_rough->answer, "%lf", &(perla_rough));
      
    /* sedimentation thresholds */
    sed_Vlim = atof(input_sed_Vlim->answer);
    sed_Slim = atof(input_sed_Slim->answer);
    sed_Vhmax = atof(input_sed_Vhmax->answer);
    sed_Shmax = atof(input_sed_Shmax->answer);
    sed_0hmax = atof(input_sed_0hmax->answer);
       
    /* output raster maps */
    S_nrRW.name = output_nrRW->answer;
    S_vel.name  = output_vel->answer;
    S_sed.name  = output_sed->answer;

    /* input flags */
    verbose = (!flag1->answer);
    maxvel = (flag2->answer);
    number_rw = (flag3->answer);
    
    /* set input map type */       
    SetMapset(&S_dem);
    SetMapset(&S_start);
    if(S_pdir.name != NULL) { SetMapset(&S_pdir); }
    /* open raster files */
    SetOpen(&S_dem);
    SetOpen(&S_start);
    if(S_pdir.name != NULL) { SetOpen(&S_pdir); }

    G_debug(3, "number of rows %d", S_dem.cellhd.rows);

    /* Allocate Input buffer */
    AllocateInBuf(&S_dem);
    AllocateInBuf(&S_start);
    if(S_pdir.name != NULL) { AllocateInBuf(&S_pdir); }
    
    /* get current region dimension */
    nrows = G_window_rows();
    ncols = G_window_cols();
  
    /* Allocate Output buffer and open file descriptor using specific data_type */
    AllocateOpenOut(&S_nrRW,CELL_TYPE);
    AllocateOpenOut(&S_vel,FCELL_TYPE);
    AllocateOpenOut(&S_sed,FCELL_TYPE);
    
    /* allocate matrix to store raster in memory */
    AllocateMatrix (&S_dem, nrows, ncols,0);
    if(S_pdir.name != NULL) { 
        AllocateMatrix (&S_pdir, nrows, ncols,0); 
    } else {
        S_pdir.data_type=FCELL_TYPE;
        AllocateMatrix (&S_pdir, nrows, ncols,1);
    }
    
    /* allocate array to store start points */
    rw_start = (rw_point *)G_malloc(sizeof(rw_point) * START_DIM );
    n_start = 0; //count found starting cells
    ia = 1; //count reallocation
    
    /* allocate matrix to store output raster in memory */
    AllocateMatrix(&S_nrRW,nrows, ncols,0);
    AllocateMatrix(&S_vel,nrows, ncols,0);
    AllocateMatrix(&S_sed,nrows, ncols,0);
    
    /* processing data message */
    if (verbose)
        G_message (_("Scanning input data:... "));

    /* go trough raster maps */
    for (row = 0; row < nrows; row++) {
    	if (verbose)
	       G_percent(row, nrows, 2);
	/* read input row */
        GetInRow(&S_dem, row);
        GetInRow(&S_start, row);
        if(S_pdir.name != NULL) { GetInRow (&S_pdir, row); }
	/* process the data */
	for (col = 0; col < ncols; col++) {
            /* copy raster value to matrix */
            GetCellValue2Matrix(&S_dem, row, col);
            if(S_pdir.name != NULL) { GetCellValue2Matrix(&S_pdir, row, col); }
            /* set new starting point in array and reallocate space if needed*/
            if( ((CELL *) S_start.in_buf)[col] != 0 ) {
                if ( (n_start!=0) && (n_start % START_DIM)==0 ) {
                    ia++;
                    tmp_start = (rw_point *)G_realloc(rw_start, sizeof(rw_point) * START_DIM *ia );
                    if (tmp_start != NULL ) { rw_start = tmp_start; }
                }
                rw_start[n_start].r = row;
                rw_start[n_start].c = col;
                n_start++;
            }   
        } //end row
    } // end col
     
     /* indicative number of random walks for each cell */
     i_rw = vol_tot/(sed_0hmax*n_start);
     
     
      
     //debug message	    
     /*G_message("Reading data completed, The total number of starting cell is n_start=%d\n",n_start);
     while(!getchar()){ }*/

    if (verbose)
        G_percent(row, nrows, 2);

   /****************************************************************/
   /* At this point we have:                                	*/
   /* - S_dem, S_slope, S_aspect & S_pdir  with filled matrices    */ 
   /* - rw_start with starting (row,col)                    	*/
   /****************************************************************/
   /*             START THE DEBRIS FLOW PROCESSING          	*/
   /****************************************************************/
   /*                                                      	*/
   /*                                                       	*/
   /*                                                       	*/

   srand(time(0));
   /* process path variables */
   rw_cell *cell_path;                    /* array to store paths */
   int pk;                                /* packages counter (number of random walk sent) */
   int step=0, nrw=0, rw_stop=0;                            /* step counter (number of cells passed during diffusion)*/
   int cr,cc;                             /* current row and current column values */
   int idx, idx_from, idx_next;           /* index of the selected cell within the mask */
   int next_r, next_c;
   double alpha, teta;                    /* parameters of the Perla equation */
   double step_lenght;                    /* distance between current and next cell */
   double delta_slope;                    /* slope difference between current and next cell */
   double sed_veloc, sed_slope;  	  /* sedimantation due to velocity and slope principles */
   double sed_res, sed_tot;               /* sedimentation residual volume for each package */
   double sed_next;
   double Vol ;                        /* sedimentation of the next cell */
   double vel_next, v_star;               /* velocity of the next cell */
   int j,bk,s;                            /* counters */
   int directions[9];                     /* directions in case 2 */
   int low_nr;
   int delta_r,delta_c;
   double dem_v,sed_v,vel_v,nrRW_v,vel_max,sed_hmax;
   double slope_n,slope_c,DHc, DN2c, DE2c,DHn,DN2n,DE2n,DN2,DE2,DL,DH,DZ, slope_s,sed_s;
   int stopbk;
   double beta_prev, beta;

   lenght = LENGHT / S_dem.cellhd.ew_res; /* approximated number of cells for single path */
  
   /* simulate process message */
   if (verbose) {
        G_percent_reset();
        G_message (_("Running random walks packages...\n"));
   }
 
   /* loop trough each starting cell (rw_start[ ]) */
   //while(getchar() != 'y') { printf("Sending %d packages.\n",n_start+1*rw); }         
nrw=0;
	
	while ( Vol < vol_tot && rw_stop!=1  ) {
		nrw++;
        
        if (number_rw) {
		if (nrw==rw) { rw_stop=1; }
	} else {
		if (nrw == (5*i_rw)) { rw_stop=1; }	
	}

   	if (number_rw) {
		if (verbose) { G_percent(nrw, rw, 2); }
	} else {
		if (verbose) { G_percent(Vol, vol_tot, 2); }	
	}
	if (nrw==5*i_rw) {
		G_warning ("The number of random walk is very high (%d).\n The sedimentation volume is: %f.\n The calculation will be stopped",nrw,Vol); 
	}
	
	for(m=0;m<n_start;m++){
        
        
        //if (verbose) { G_percent(m, n_start, 2); }
	
	
	
	//while(getchar() != 'y') { printf("\nStarting cell %d of %d.\n",m+1,n_start); }                  
        //printf("m=%d",m);
        //while(!getchar()){ }
        /* send rw packages */
	//for(pk=0;pk<rw;pk++){


//	while(Vol < (vol_tot/n_start)){

	    
	    //while(getchar() != 'y') { printf("\nSending package %d of %d.\n",pk+1,rw); }
            /* initialize package variables */
            stop = FALSE;          
            step = 0; /* nn = index of the cell within the random walk */
            //sed_res=sed_0hmax;
            ia = 1; // count reallocations
            /* allocate array to store cells path and values (cell_path[ ]) */
            cell_path = (rw_cell *)G_malloc(sizeof(rw_cell) * lenght);
            /* calculate the singles cells path */        
            while (stop != TRUE) {            
                /* if it is full then realloc cells path array */
                if ( (step!=0) && (step+1 % lenght)==0 ) {
                    ia++;
		    while(getchar() != 'y') { printf("\nREALLOC"); }
                    tmp_rw = (rw_cell  *)G_realloc(cell_path, sizeof(rw_cell ) * lenght * ia);
                    if (tmp_rw != NULL) { cell_path = tmp_rw; }
                }
                 
                // set current row and col
                if (step==0) {
                    cell_path[step].point.r = rw_start[m].r;
                    cell_path[step].point.c = rw_start[m].c;                   
                    cell_path[step].Vel = 0;
                    cell_path[step].Sed = 0;
                    cell_path[step].Pend = 0;
                    // ERRORE non erano definiti cr e cc
                    /*cr = cell_path[step].point.r;
                    cc = cell_path[step].point.c;
                    cell_path[step].h = (double) GetValue(&S_dem,cr,cc);*/
                    //forse basta la seguente riga 
                     cell_path[step].h = (double) GetValue(&S_dem,rw_start[m].r,rw_start[m].c);
                    //printf("Zero step: %d \n",step);
                }
                cr = cell_path[step].point.r;
                cc = cell_path[step].point.c;
                /* get the 8 neighbours cells of the current cell if not on the region border */
                if(cr==0 || cr==nrows-1 || cc==0 || cc==ncols-1) {
                	G_warning (_("Package at the bounduary of the current region: please enlarge your region!"));
			vel_next=0;
		} else {
		     	/* set the current 3x3 mask */
		    	SetMask(&mask, &S_dem, &S_pdir, cr, cc);
			/* get the previous cell index in current mask */
			if(step>0) {
	                    delta_r = cell_path[step-1].point.r - cell_path[step].point.r;
	                    delta_c = cell_path[step-1].point.c - cell_path[step].point.c;
	                    /*if(delta_r==0 && delta_c==0){
				printf("c'e' qualcosa che non va");
				while(!getchar()){ }
			    }*/  
		        } else {
		            delta_r = 0;
		            delta_c = 0;
		        }
						
			/* get the index of the previous cell */
			idx_from = GetIdx(delta_r, delta_c);
				
		        /* check diffusion mode */
			//while(getchar() != 'y') { printf("slim=%lf",slim);	}			
		        DirMode = SetDirectionMode(&mask, slim, idx_from);
		        cell_path[step].mode = DirMode;
				
			//printf("Mode: %d\n",DirMode);
			int k;
			//for(k=0;k<9;k++){
			//	printf("mask.h[%d]: %lf\n",k,mask.h[k]);
			//	printf("mask.s[%d]: %lf\n",k,mask.s[k]);
			//}

			/* DirMode = 1 -> direzione fissa SFD */
			/* DirMode = 2 -> random walk MFD */
			/* DirMode = 3 -> pianura MFD */
			/* DirMode = 4 -> sink */




			/* calculate values for the next cell */
			switch (DirMode) {
			case 1: /* deterministic path: the package follows the max slope */

				/* get the index of the next cell (maximum slope) */
				idx = maximum(mask.w);
				//printf("next cell idx=%d",idx);					

				/* set next cell row and column */
			        cell_path[step+1].point.r = GetNewRowIdx(idx,cr);
			        cell_path[step+1].point.c = GetNewColIdx(idx,cc);
				
			        /* set slope */
			        cell_path[step+1].Pend = mask.s[idx];
			        if (cell_path[step+1].Pend<=0) {
			        	printf("ecco");
			        	while(!getchar()){ }
			        }
				
			        /* set height */
			        cell_path[step+1].h = mask.h[idx];      
			        
			        /* estimate velocity with Perla model */
				beta = atan( cell_path[step+1].Pend );
				beta_prev = atan( cell_path[step].Pend );
			        alpha = 9.81 * ( sin(beta) - perla_rough * cos(beta) );
			        //alpha = 9.81 * sin(cell_path[step+1].Pend) - perla_rough * cos(cell_path[step+1].Pend);
			        if (idx == 0 || idx == 2 || idx == 6 || idx == 8) {
			            step_lenght = sqrt(2) * S_dem.cellhd.ns_res;
			        } else {  
			            step_lenght = S_dem.cellhd.ns_res;
			        }
			        /*if( cell_path[step+1].Pend > cell_path[step].Pend ) {
			            delta_slope = cell_path[step+1].Pend - cell_path[step].Pend; 
			        } else {
			            delta_slope = 0;
			        }*/
			        teta = (-2 * step_lenght ) /  perla_m2d;

				if(cell_path[step+1].Pend <= cell_path[step].Pend) {
					v_star = cell_path[step].Vel * cos(beta_prev-beta);
					//v_star = pow(cell_path[step].Vel,2) * cos(beta_prev-beta);
				} else {
					v_star = cell_path[step].Vel;
				}
			
				vel_next = alpha * perla_m2d * (1 - exp(teta)) + pow(v_star,2) * exp(teta);
				//vel_next = alpha * perla_m2d * (1 - exp(teta)) + v_star * exp(teta);
				//while(getchar() != 'y') { 
				//printf("Calculated values:\nalpha=%lf\nstep_lenght=%lf\ndelta_slope=%lf\nteta=%lf\n",alpha, step_lenght,delta_slope,teta);
				//printf("(1 - exp(teta)=%lf\npow(cell_path[step].Vel,2)=%lf\ncos(delta_slope)=%lf\nexp(teta)=%lf\npow(v_star,2)=%lf\n",(1 - exp(teta)), pow(cell_path[step].Vel,2),cos(delta_slope),exp(teta),pow(v_star,2));
				//	printf("\n SED\n cell_path[step+1].Vel=%lf - sed_Vlim=%lf",cell_path[step+1].Vel,sed_Vlim);
				//	printf("\n SED\n cell_path[step+1].Sed=%lf - sed_Slim=%lf",cell_path[step+1].Pend,sed_Slim);
				//}
						
				if(vel_next>0) {
					vel_next = sqrt(vel_next);
				} else {
					vel_next = 0;
				}
					
			        /* estimate sedimentation */
			        //if ( cell_path[step+1].Vel  < sed_Vlim ) {
				if ( vel_next  < sed_Vlim ) {
	                    		sed_veloc = sed_Vhmax * ( 1 - (vel_next / sed_Vlim ));
					//sed_veloc = min_among_2 (sed_veloc,sed_Vlim);
					sed_veloc = min_among_2 (sed_veloc,sed_Vhmax);
	                	} else {
	                    		sed_veloc = 0;
	                	}

			        if ( cell_path[step+1].Pend  < sed_Slim ) {
			            sed_slope = sed_Shmax * ( 1 - (cell_path[step+1].Pend / sed_Slim ));
				    //sed_slope = min_among_2 (sed_slope,sed_Slim);
				    sed_slope = min_among_2 (sed_slope,sed_Shmax);
			        } else {
			            sed_slope = 0;
			        }
					
				if(sed_veloc < sed_slope) {
					sed_next = sed_veloc;
				} else {
					sed_next = sed_slope;
				}
					
	                	//sed_next = min_among_2 (sed_veloc,sed_slope);
				if (sed_next>0){
					printf("vel_next=%f, sed_next=%f\n",vel_next, sed_next);
					while(!getchar()){ }
				}
				//while(getchar() != 'y') { 
				//printf("Calculated values 1: vel=%lf - Sed=%lf\n",vel_next, sed_next);						
				//printf("Calculated values 1: sed_v=%lf - Sed_s=%lf\n",sed_veloc, sed_slope);							
				//}						
		    	break; 
				    
			case 2: /* stocastic path: the package follows a statistica extraction */
				
				/* choose the next cell */                        
				if(step>0) {
				    delta_r = cell_path[step-1].point.r - cell_path[step].point.r;
				    delta_c = cell_path[step-1].point.c - cell_path[step].point.c;
				} else {
				    delta_r = 0;
				    delta_c = 0;
				}
				
				/* get the index of the next cell (maximum slope) */    
				idx = get_cell_mode2(&mask, slim, esp, delta_r, delta_c); 
				
				/* set next cell row and column */
				cell_path[step+1].point.r  = GetNewRowIdx(idx,cr);
				cell_path[step+1].point.c = GetNewColIdx(idx,cc);
				
				/* set slope */
				cell_path[step+1].Pend = mask.s[idx];
				if (cell_path[step+1].Pend<=0) {
			        	printf("ecco2");
			        	while(!getchar()){ }
			        }
				
				/* set height */
				cell_path[step+1].h = mask.h[idx];
				
				/* estimate velocity with Perla model */
				beta = atan( cell_path[step+1].Pend );
			    	beta_prev = atan( cell_path[step].Pend );
				alpha = 9.81 * ( sin(beta) - perla_rough * cos(beta) );
				//alpha = 9.81 * sin(cell_path[step+1].Pend) - perla_rough * cos(cell_path[step+1].Pend);
				if (idx == 0 || idx == 2 || idx == 6 || idx == 8) {
				    step_lenght = sqrt(2) * S_dem.cellhd.ns_res;
				} else {  
				    step_lenght = S_dem.cellhd.ns_res;
				}
				/*if( cell_path[step+1].Pend > cell_path[step].Pend ) {
				    delta_slope = cell_path[step+1].Pend - cell_path[step].Pend; 
				} else {
				    delta_slope = 0;
				}*/
		        	teta = (-2 * step_lenght ) /  perla_m2d;
				
				if(cell_path[step+1].Pend <= cell_path[step].Pend) {
					v_star = cell_path[step].Vel * cos(beta_prev-beta);
					//v_star = pow(cell_path[step].Vel,2) * cos(beta_prev-beta);
				} else {
					v_star = cell_path[step].Vel;
				}
		
				vel_next = alpha * perla_m2d * (1 - exp(teta)) + pow(v_star,2) * exp(teta);
				//vel_next = alpha * perla_m2d * (1 - exp(teta)) + v_star * exp(teta);

				if(vel_next>0) {
					vel_next = sqrt(vel_next);
				} else {
					vel_next = 0;
				}
				
				/* estimate sedimentation */
				if ( vel_next  < sed_Vlim ) {
				    sed_veloc = sed_Vhmax * ( 1 - (vel_next / sed_Vlim ));
				    //sed_veloc = min_among_2 (sed_veloc,sed_Vlim);
				    sed_veloc = min_among_2 (sed_veloc,sed_Vhmax);
				} else {
				    sed_veloc = 0;
				}
				if ( cell_path[step+1].Pend  < sed_Slim ) {
				    sed_slope = sed_Shmax * ( 1 - (cell_path[step+1].Pend / sed_Slim ));
				    //sed_slope = min_among_2 (sed_slope,sed_Slim);
				    sed_slope = min_among_2 (sed_slope,sed_Shmax);
				} else {
				    sed_slope = 0;
				}
				
				if(sed_veloc < sed_slope) {
					sed_next = sed_veloc;
				} else {
					sed_next = sed_slope;
				}
				
				//sed_next = min_among_2 (sed_veloc,sed_slope);
				
				//while(getchar() != 'y') { 
				//printf("Calculated values 1: vel=%lf - Sed=%lf\n",vel_next, sed_next);						
				//printf("Calculated values 1: sed_v=%lf - Sed_s=%lf\n",sed_veloc, sed_slope);							
				//}	                       
	            	break;               
	            
	            	case 3: /* horizontal plane diffusion: stocastic approach */
	                
				/* choose the next cell */                   
			        if(step>0) {
			            delta_r = cell_path[step-1].point.r - cell_path[step].point.r;
			            delta_c = cell_path[step-1].point.c - cell_path[step].point.c;
			        } else {
			            delta_r = 0;
			            delta_c = 0;
			        }
		      
	                	/* get the index of the next cell (maximum slope) */ 
				//while(getchar() != 'y') { printf("pre get cell 3"); }						
	                	idx = get_cell_mode3(&mask, slim, esp, delta_r, delta_c);                       
				//while(getchar() != 'y') { printf("post get cell 3"); }
	                
				/* set next cell row and column */
			        cell_path[step+1].point.r  = GetNewRowIdx(idx,cr);
			        cell_path[step+1].point.c = GetNewColIdx(idx,cc);
			        /* set slope */
			        cell_path[step+1].Pend = mask.s[idx];
			        if (cell_path[step+1].Pend<0) {
			        	printf("ecco3");
			        	while(!getchar()){ }
			        }
			        /* set height */
			        cell_path[step+1].h = mask.h[idx]; 
		                
				/* estimate velocity with Perla model */                        
				beta = atan( cell_path[step+1].Pend );
				beta_prev = atan( cell_path[step].Pend );
				alpha = 9.81 * ( sin(beta) - perla_rough * cos(beta) );
				//alpha = 9.81 * sin(cell_path[step+1].Pend) - perla_rough * cos(cell_path[step+1].Pend);
				if (idx == 0 || idx == 2 || idx == 6 || idx == 8) {
				    step_lenght = sqrt(2) * S_dem.cellhd.ns_res;
				} else {  
				    step_lenght = S_dem.cellhd.ns_res;
				}
			        /*if( cell_path[step+1].Pend > cell_path[step].Pend ) {
			            delta_slope = beta -beta_prev;
			        } else {
			            delta_slope = 0;
			        }*/
			        teta = (-2 * step_lenght ) /  perla_m2d;
					
				if(cell_path[step+1].Pend <= cell_path[step].Pend) {
					v_star = cell_path[step].Vel * cos(beta_prev-beta);
					//v_star = pow(cell_path[step].Vel,2) * cos(beta_prev-beta);
				} else {
					v_star = cell_path[step].Vel;
				}
			
				vel_next = alpha * perla_m2d * (1 - exp(teta)) + pow(v_star,2) * exp(teta);
				//vel_next = alpha * perla_m2d * (1 - exp(teta)) + v_star * exp(teta);
			
				if(vel_next>0) {
					vel_next = sqrt(vel_next);
				} else {
					vel_next = 0;
				}
					
			        /* estimate sedimentation */
			        if ( vel_next  < sed_Vlim ) {
			            sed_veloc = sed_Vhmax * ( 1 - (vel_next / sed_Vlim ));
			            sed_veloc = min_among_2(sed_slope,sed_Vhmax);
			        } else {
			            sed_veloc = 0;
			        }
			        if ( cell_path[step+1].Pend  < sed_Slim ) {
			            sed_slope = sed_Shmax * ( 1 - (cell_path[step+1].Pend / sed_Slim ));
			            sed_slope = min_among_2(sed_slope,sed_Shmax);
			        } else {
			            sed_slope = 0;
			        }
					
				if(sed_veloc < sed_slope) {
					sed_next = sed_veloc;
				} else {
					sed_next = sed_slope;
				}
				//printf("vel=%f, Pend=%f, sed_next=%f\n",vel_next,cell_path[step+1].Pend, sed_next);
					
	                	//sed_next = min_among_2 (sed_veloc,sed_slope);
					
				//printf("Calculated values 3: vel=%lf - Sed=%lf\n",vel_next, sed_next);							
	                
	                /* calculate values (velocity, sedimentation, number of pass) */               
	            	break; 

			case 4: /* sink: go back and fill with sediment the hole */
  						
				//for(k=0;k<9;k++){
				//	printf("mask.h[%d]: %lf - ",k,mask.h[k]);
				//	printf("mask.s[%d]: %lf\n",k,mask.s[k]);
				//}
			        //choose the next cell 
				if(step>0) {
			            delta_r = cell_path[step-1].point.r - cell_path[step].point.r;
			            delta_c = cell_path[step-1].point.c - cell_path[step].point.c; 
			        } else {
			            delta_r = 0;
			            delta_c = 0;
			        }
						
				// get the index of the previous cell 
				idx_from = GetIdx(delta_r, delta_c);
	                	
					
				/* get lower cell in the mask except the last one*/
				idx_next = GetLower(&mask,idx_from);
				//printf("next_idx=%d",idx_next);				
					
				/* get row and col of the next cell */
				next_r = GetNewRowIdx(idx_next,cr);
	                	next_c = GetNewColIdx(idx_next,cc);
					
				/* go back in path by BACK meters or to the starting cell and get distance */
				bk=0; 
				//DL = 0;
				DL= sqrt(pow((next_r - cell_path[step].point.r)*S_dem.cellhd.ns_res,2) + pow((next_c - cell_path[step].point.c)*S_dem.cellhd.ns_res,2));
				while( DL < BACK && bk<step ){
					bk++;
					DN2 = pow((cell_path[step-bk].point.r - cell_path[step-bk+1].point.r)*S_dem.cellhd.ns_res,2);
					DE2 = pow((cell_path[step-bk].point.c - cell_path[step-bk+1].point.c)*S_dem.cellhd.ew_res,2);
					DL += sqrt( DN2 + DE2  );								
				}
				//while(getchar() != 'y') { 
				//	printf("back of %d steps of %d\n",bk,step);
				//	printf("current[%d][%d] -- backcell[%d][%d]\n",cr,cc,cell_path[step-bk].point.r,cell_path[step-bk].point.c);

				//}

				/* estimate the mean slope among 50m back cell and the current cell */
				slope_s = (cell_path[step-bk].h - mask.h[idx_next]) / DL;
				//slope_s = (cell_path[step-bk].h - cell_path[idx_next].h) / DL;
					
				//printf("estimated slope (backcell-next) = %lf\n",slope_s);						
				
				/* calculate sedimentation in the sink cell */				
				DL= sqrt(pow((next_r - cell_path[step].point.r)*S_dem.cellhd.ns_res,2) + pow((next_c - cell_path[step].point.c)*S_dem.cellhd.ns_res,2));
				sed_s = slope_s * DL  + mask.h[idx_next] - cell_path[step].h;
				if (sed_s > 0) {
					cell_path[step].Sed = sed_s;
				}
				
				/* calculate sedimentation in the cells above sink */
				for(s=1;s<=bk;s++){
					DN2 = pow((cell_path[step-s].point.r - cell_path[step-s+1].point.r)*S_dem.cellhd.ew_res,2);
					DE2 = pow((cell_path[step-s].point.c - cell_path[step-s+1].point.c)*S_dem.cellhd.ew_res,2);
					DL += sqrt(DN2+DE2);
					sed_s = (slope_s) * DL + mask.h[idx_next] - cell_path[step-s].h;
					if (sed_s > 0) {
						cell_path[step-s].Sed = sed_s;
					}
				//printf("estimated sedimantation (%d) = %lf\n",step-s, sed_s);							
				}
					
				/* stop the path */
				vel_next=0;
				sed_next=0;

					//printf("Calculated values 4: vel=%lf - Sedcur=%lf\n",vel_next, sed_s);
					//while(getchar() != 'y') { printf("SINK FOUND"); } 						
						
		                 	/* calculate values (velocity, sedimentation, number of pass) */               
		            	break;
	            	




		        	}  /* END SWITCH - found next cell in path */  
				
			} /* END IF OUT OF REGION */
			/* set calculated values in path array if not stopped */
			if(vel_next<=0) {
			    stop=TRUE;
			} else {
			    //cell_path[step+1].Vel = vel_next;
			    //cell_path[step+1].Sed = sed_next;
			    cell_path[step+1].Vel = vel_next;
			    cell_path[step+1].Sed = sed_next;
			}   
				
			//printf("Calculated values: step=%d vel=%lf - Sed=%lf\n",step, vel_next, sed_next);				
			step = step + 1;			
		} /* END WHILE - path of package done (stop=TURE)*/
		//if(step<5) {
		//printf("package %d of %d stopped. Executed %d steps\n",pk+1,rw,step); 
		//}
		/* calculate total sedimented height */
		sed_tot = 0;
		for(j=1;j<step;j++) {
			sed_tot += cell_path[j].Sed;
			//printf("step=%d, row=%d, col=%d ",j ,cell_path[j].point.r, cell_path[j].point.c);
			//printf("mode: %d ",cell_path[j].mode);
			//printf("quota: %lf - ",cell_path[j].h);
			//printf("Vel: %lf - ",cell_path[j].Vel);
			//printf("Sed: %lf - ",cell_path[j].Sed);
			//printf("Slope: %lf\n",cell_path[j].Pend);			   
		}

		//while(getchar() != 'y') { printf("...DONE!\n"); }

		/* store path in matrices */
		for(j=1;j<step;j++) {

			//printf("\n\nStoring cell values in matrices");            
			       
			/* store path values in matrices */
			row = cell_path[j].point.r;
			col = cell_path[j].point.c;      
			//printf("row=%d, col=%d sed =%f - ",row,col,cell_path[j].Sed);
					   
			/* velocity sum or max */
		   	vel_v = GetValue(&S_vel,row,col);
			if (maxvel) {
				vel_v = max_among_2(vel_v,cell_path[j].Vel);
			} else {
				vel_v += cell_path[j].Vel;
			}        
			SetValue(&S_vel,row,col,vel_v);

			/* sedimentation sum */
			sed_v = GetValue(&S_sed,row,col);
			//printf("sed[%d][%d]=%f\n",row,col,sed_v);
			//printf("\nstored value Sed: %lf - calculated value: %lf - sum of path: %lf",sed_v,cell_path[j].Sed,sed_tot);			   
			sed_hmax = sed_0hmax/(S_dem.cellhd.ew_res*S_dem.cellhd.ns_res);
			//printf("\n max value Sed: %lf - ",sed_hmax);				   
			if( sed_tot > sed_hmax ) {
				sed_v += (sed_hmax * cell_path[j].Sed) / sed_tot;
				//printf("\nsed excess, adding value: %lf - ",(sed_hmax * cell_path[j].Sed) / sed_tot);					   
			} else {
				sed_v += cell_path[j].Sed;
				//printf("\nsed no-excess, adding value: %lf - ",cell_path[j].Sed);					   				   
			}
			//printf("\nsetted value Sed: %lf - ",sed_v);			   			   
			SetValue(&S_sed,row,col,sed_v);
			//printf("sed[%d][%d]=%f\n",row,col,sed_v);
	       		
	       		       		
	       		/* RW number counter */
			nrRW_v = GetValue(&S_nrRW,row,col);
			if(cell_path[j].Vel > 0 ) {
				if(nrRW_v >= rw*n_start) {
				   nrRW_v = rw;
				} else {
				   nrRW_v++;
				}
			}
			SetValue(&S_nrRW,row,col,nrRW_v);			   
				   
			/* dem update */
			dem_v = GetValue(&S_dem,row,col);
			SetValue(&S_dem,row,col,dem_v+sed_v);
					   
			//while(getchar() != 'y') { printf("\n- Sedimentation stored"); }                
			//printf("step=%d, row=%d, col=%d ",j ,row, col);
			//printf("mode: %d ",cell_path[j].mode);
			//printf("nrRW: %d - ",nrRW_v);
			//printf("Vel: %lf - ",vel_v);
			//printf("Sed: %lf - ",sed_v);
			//printf("Dem: %lf\n",dem_v+sed_v);

		}
		// G_free(cell_path);
		free(cell_path);
		//while(getchar() != 'y') { printf("package %d done\n",pk); }		   
        //} /* END FOR - all packages sent */
        /* total volume update*/
		Vol += min_among_2(sed_hmax,sed_tot)*(S_dem.cellhd.ew_res*S_dem.cellhd.ns_res);
		//printf("Vol=%f --> ",Vol);
		/*if (min_among_2(sed_hmax,sed_tot)==0) {
			Vol= vol_tot/n_start+1;
		}*/
 	// } /* END WHILE on the total volume */
        } /* END FOR - all start cells processed */
        m=0;
    //printf("m=%d, Vol=%f --> ",m,Vol);
    //Vol=0;
    //printf("Vol=%f\n",Vol);
    //while(!getchar()){ } 
    printf("All packages sent!\n\n");       
   // } /* END FOR - all start cells processed */
   } /* END WHILE on the total volume */

   if (number_rw) {
	G_message ("The total sedimentation volume is %f m3\n", Vol);
   } else {
	G_message ("The number of random walk carried out is:%d\n", nrw);	
   }

/************************ FINE CALCOLO ***************************************************/

	printf("Writing rasters");   

    /* write output matrix to output raster map */
    for (row = 0; row < nrows; row++) {
        for (col = 0; col < ncols; col++) {
            ((CELL *) S_nrRW.out_buf)[col] = S_nrRW.cell_matrix[row][col];
            ((FCELL *) S_sed.out_buf)[col]  = S_sed.fcell_matrix[row][col];
            if (S_nrRW.cell_matrix[row][col] >0) {
                ((FCELL *) S_vel.out_buf)[col]   = S_vel.fcell_matrix[row][col] / S_nrRW.cell_matrix[row][col];
            } else {
                ((FCELL *) S_vel.out_buf)[col]   = 0;
            }
      }

        
	     if (G_put_raster_row(S_nrRW.fd, S_nrRW.out_buf, CELL_TYPE) < 0)
	         G_fatal_error(_("Cannot write to <%s>"), S_nrRW.name);
	     if (G_put_raster_row(S_vel.fd, S_vel.out_buf, FCELL_TYPE) < 0)
	         G_fatal_error(_("Cannot write to <%s>"), S_vel.name);
	     if (G_put_raster_row(S_sed.fd, S_sed.out_buf, FCELL_TYPE) < 0)
	         G_fatal_error(_("Cannot write to <%s>"), S_sed.name);
	}

printf("...DONE!\n");

struct History history_nrRW, history_vel, history_sed; 

printf("Writing history");

    /* add command line incantation to history file */
    G_short_history(S_nrRW.name, "raster", &history_nrRW);
    G_command_history(&history_nrRW);
    G_write_history(S_nrRW.name, &history_nrRW);
    G_short_history(S_vel.name, "raster", &history_vel);
    G_command_history(&history_vel);
    G_write_history(S_vel.name, &history_vel);
    G_short_history(S_sed.name, "raster", &history_sed);
    G_command_history(&history_sed);
    G_write_history(S_sed.name, &history_sed);

printf("...DONE!\n");


printf("Cleaning memory");
//while(!getchar()){ }
    /* memory cleanup */
//while(getchar() != 'y') { printf("...1 DONE!\n"); }        
    free_close_rast_struct(&S_dem, nrows, ncols);
//while(getchar() != 'y') { printf("...2 DONE!\n"); }
    /*  NON SO MA NON VA IL FREE CELL */
    free(rw_start);
    //G_free(rw_start);    
    // free_close_rast_struct(&S_start, nrows, ncols);
//while(getchar() != 'y') { printf("...3 DONE!\n"); }     
    //if(obj != NULL)
    //    free_close_rast_struct(&S_obj, nrows, ncols);
    if(S_pdir.name != NULL)
        free_close_rast_struct(&S_pdir, nrows, ncols);
//while(getchar() != 'y') { printf("...4 DONE!\n"); }          
    free_close_rast_struct(&S_nrRW, nrows, ncols);
//while(getchar() != 'y') { printf("...5 DONE!\n"); }      
    free_close_rast_struct(&S_vel, nrows, ncols);
//while(getchar() != 'y') { printf("...6 DONE!\n"); }      
    free_close_rast_struct(&S_sed, nrows, ncols);    

//while(getchar() != 'y') { 
	printf("...ALL DONE!\n"); 
//	}

    exit(EXIT_SUCCESS);
}
