/*****************************************************************************
*
* MODULE:	r.rockcone
*
* AUTHOR:	Massimiliano Cannata - massimiliano.cannata[]supsi.ch (2005)
*
* PURPOSE:	Estimate the area potentially affected by gravitative falls by
*               using a cone as energy surface.
*
* COPYRIGHT:	(C) 2006 by Istituto Scienze della Terra
*
*		This program is free software under the GNU General Public
*		Licence (>=2). Read the file COPYING that cames with GRASS
*		for details.
*
/***************************************************************************/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <grass/gis.h>
#include <math.h>
#include <grass/glocale.h>

#define min(a,b) ((a) < (b) ? (a) : (b))
#define max(a,b) ((a) > (b) ? (a) : (b))

#define ALLOC_DIM 10000
#define PI 3.14159265

int main(int argc, char *argv[]){

	typedef struct
    	{
       	   double sx, sy, sz, dir;
    	} startpt;

	startpt *a_start,*tmp_start;

	struct Cell_head cellhd;
	struct Cell_head window;
	struct History history;
	
	/* input-output raster files */
	int infd_ELEV,infd_START,infd_ASP;
	int outfd_F,outfd_V,outfd_P,outfd_C;
	
	/* mapset name locator */
	char *mapset_ELEV,*mapset_START,*mapset_ASP;

	/* buffer for input-output raster */
	DCELL *inrast_ELEV;			/* elevation map [m]*/
	DCELL *inrast_ASP;			/* aspect map [m]*/
	CELL  *inrast_START;			/* start map [not null cells]*/
	CELL  *outrast_F;			/* count frequency [-] */ 
	DCELL *outrast_V;		        /* velocity [m/s] */ 
	DCELL *outrast_P;	                /* potential energy */
	DCELL *outrast_C;			/* cinetic energy */

	/* cell values */
	DCELL d_ELEV;
	DCELL d_ASP;
	CELL  d_START;
	DCELL d_V,d_P,d_C;
        CELL c_F;
	
	/* cell counters */
	int nrows, ncols;
	int row, col;

	int start_count, k, ia;
	int i_freq;
        double i_max,i_min,i_sum,zi;
	double ee,nn;
	double Dz,Dzi,Dhi,Di;
        double north, east;
	
	/* variables to handle user inputs */
	char *ELEV, *ASP, *START, *OUT_F, *OUT_V, *OUT_P, *OUT_C, *METHOD;
        double ANG, FVEL, MASS, DIR;
	
	/* GRASS structure */
	struct GModule *module;
	struct Option *input_ELEV, *input_ASP, *input_START, *input_ANG, *input_DIR, *input_FVEL, *input_MASS, *input_METHOD;
	struct Option *output_F,*output_V,*output_P,*output_C;

	G_gisinit(argv[0]);
	
	module = G_define_module();
	module->keywords = _("raster, rockfall");
	module->description = _("Estimate the area potentially affected by gravitative falls");
	
	/* Define different options */
	input_ELEV = G_define_option();
	input_ELEV->key	= "elev";
	input_ELEV->type = TYPE_STRING;
	input_ELEV->required = YES;
	input_ELEV->gisprompt = "old,cell,raster";
	input_ELEV->description = _("Name of elevation raster map");

	input_ASP = G_define_option();
	input_ASP->key	= "aspect";
	input_ASP->type = TYPE_STRING;
	input_ASP->required = YES;
	input_ASP->gisprompt = "old,cell,raster";
	input_ASP->description = _("Name of aspect raster map");
	
	input_START = G_define_option();
	input_START->key = "start";
	input_START->type = TYPE_STRING;
	input_START->required = YES;
	input_START->gisprompt = "old,cell,raster";
	input_START->description = _("Name of input starting cells raster map");
		
	input_ANG = G_define_option();
	input_ANG->key	= "deg";
	input_ANG->type = TYPE_DOUBLE;
	input_ANG->required = YES;
	input_ANG->multiple   = NO;
	input_ANG->answer     = "35";
	input_ANG->gisprompt = _("0-90");
	input_ANG->description = _("Degree of cone inclination from the horizont [0<x<90]");

	input_DIR = G_define_option();
	input_DIR->key	= "dir";
	input_DIR->type = TYPE_INTEGER;
	input_DIR->required = NO;
	input_DIR->multiple   = NO;
	input_DIR->answer     = "90";
//	input_DIR->gisprompt = _("0-360");
	input_DIR->description = _("Degree of direction angle range [0<x<360]");

	input_FVEL = G_define_option();
	input_FVEL->key	= "fv";
	input_FVEL->type = TYPE_STRING;
	input_FVEL->required = NO;
	input_FVEL->multiple   = NO;
	input_FVEL->answer     = "1.0";
//	input_FVEL->gisprompt = _("0-1");
	input_FVEL->description = _("Velocity factor reduction [0<x<1]");

	input_MASS = G_define_option();
	input_MASS->key	= "mass";
	input_MASS->type = TYPE_STRING;
	input_MASS->required = NO;
	input_MASS->multiple   = NO;
//	input_MASS->answer     = "100";
	input_MASS->description = _("Mean block mass [Kg]");

	output_F = G_define_option();
	output_F ->key = "freq";
	output_F ->type = TYPE_STRING;
	output_F ->required = NO;
	output_F ->gisprompt = "new,cell,raster";
	output_F ->description = _("Output frequency raster map");
	output_F ->guisection  = _("Output_options");

	output_V = G_define_option();
	output_V ->key = "vel";
	output_V ->type = TYPE_STRING;
	output_V ->required = NO;
	output_V ->gisprompt = "new,cell,raster";
	output_V ->description = _("Output velocity raster map");
	output_V ->guisection  = _("Output_options");

/*M
	output_P = G_define_option();
	output_P ->key = "epot";
	output_P ->type = TYPE_STRING;
	output_P ->required = NO;
	output_P ->gisprompt = "new,cell,raster";
	output_P ->description = _("Output potential energy raster map");
	output_P ->guisection  = _("Output_options");
*/
	output_C = G_define_option();
	output_C ->key = "ecin";
	output_C ->type = TYPE_STRING;
	output_C ->required = NO;
	output_C ->gisprompt = "new,cell,raster";
	output_C ->description = _("Output kinetic energy raster map");
	output_C ->guisection  = _("Output_options");

	input_METHOD = G_define_option();
	input_METHOD->key	= "method";
	input_METHOD->key_desc   = "max|ave|min" ;
	input_METHOD->type = TYPE_STRING;
	input_METHOD->required = NO;
	input_METHOD->multiple   = NO;
	input_METHOD->answer     = "max";
	input_METHOD->options    = "max,ave,min";
	input_METHOD->description = _("Method for cones velocity and energy aggregation");

	if (G_parser(argc, argv))
		exit(EXIT_FAILURE);
	
	/* get entered parameters */
	ELEV=input_ELEV->answer;
	ASP=input_ASP->answer;
	START=input_START->answer;
	OUT_F=output_F->answer;
	OUT_V=output_V->answer;
//	OUT_P=output_P->answer;
	OUT_C=output_C->answer;
	METHOD=input_METHOD->answer;
	ANG = atof(input_ANG->answer);
	DIR = atof(input_DIR->answer);
	FVEL = atof(input_FVEL->answer);
	if(input_MASS->answer)
	    MASS = atof(input_MASS->answer);
        //sscanf(input_ANG->answer, "%f", &ANG);
        //sscanf(input_FVEL->answer, "%f", &FVEL);
        //sscanf(input_MASS->answer, "%f", &MASS);
		
	/* find maps in mapset */
	mapset_ELEV = G_find_cell2 (ELEV, "");
	if (mapset_ELEV == NULL)
	        G_fatal_error (_("cell file [%s] not found"), ELEV);

	mapset_START = G_find_cell2 (START, "");
	if (mapset_START == NULL)
	        G_fatal_error (_("cell file [%s] not found"), START);

	mapset_ASP = G_find_cell2 (ASP, "");
	if (mapset_ASP == NULL)
	        G_fatal_error (_("cell file [%s] not found"), ASP);
        

	/* check legal output name */
	if (OUT_F) {
		if (G_legal_filename (OUT_F) < 0)
			G_fatal_error (_("[%s] is an illegal name"), OUT_F);
	}
	if (OUT_V) {
		if (G_legal_filename (OUT_V) < 0)
			G_fatal_error (_("[%s] is an illegal name"), OUT_V);
	}
/*M
	if (OUT_P) {
		if (G_legal_filename (OUT_P) < 0)
			G_fatal_error (_("[%s] is an illegal name"), OUT_P);
	}
*/
	if (OUT_C) {
		if (G_legal_filename (OUT_C) < 0)
			G_fatal_error (_("[%s] is an illegal name"), OUT_C);
	}

	/* check at least one output is defined */
//	if( !(OUT_F || OUT_V || OUT_P || OUT_C) )
	if( !(OUT_F || OUT_V || OUT_C) )
		G_fatal_error (_("you must define at least one output map"));

	if( OUT_C && !MASS )
		G_fatal_error (_("for energy map calculation the \"mass\" parameters is required"));

/****************** collecting start cells **************************

	/* open input raster files */
	if ( (infd_ELEV = G_open_cell_old (ELEV, mapset_ELEV)) < 0)
		G_fatal_error (_("Cannot open cell file [%s]"), ELEV);
	if ( (infd_ASP = G_open_cell_old (ASP, mapset_ASP)) < 0)
		G_fatal_error (_("Cannot open cell file [%s]"), ASP);
	if ( (infd_START = G_open_cell_old (START, mapset_START)) < 0)
		G_fatal_error (_("Cannot open cell file [%s]"),START);	

	if (G_get_cellhd (ELEV, mapset_ELEV, &cellhd) < 0)
		G_fatal_error (_("Cannot read file header of [%s]"), ELEV);
	if (G_get_cellhd (ASP, mapset_ASP, &cellhd) < 0)
		G_fatal_error (_("Cannot read file header of [%s]"), ASP);
	if (G_get_cellhd (START, mapset_START, &cellhd) < 0)
		G_fatal_error (_("Cannot read file header of [%s]"), START);
		
	/* Allocate input buffer */
	inrast_ELEV	= G_allocate_d_raster_buf();
	inrast_ASP	= G_allocate_d_raster_buf();
	inrast_START	= G_allocate_c_raster_buf();	

	/* get windows rows & cols */
	nrows	= G_window_rows();
	ncols	= G_window_cols();
	G_get_window(&window); 

	/* alloc start array */
	a_start = (startpt *)G_malloc(sizeof(startpt) * ALLOC_DIM );
	start_count = 0;
    	ia = 1; //count reallocation

    	G_debug(3, "starting source detection...\n");


	for (row = 0; row < nrows; row++)
	{
		G_percent (row, nrows, 2);
				
		/* read a line input maps into buffers*/	
		if (G_get_d_raster_row (infd_ELEV, inrast_ELEV, row) < 0)
			G_fatal_error (_("Could not read from <%s>"),ELEV);
		if (G_get_d_raster_row (infd_ASP, inrast_ASP, row) < 0)
			G_fatal_error (_("Could not read from <%s>"),ASP);
		if (G_get_c_raster_row (infd_START, inrast_START, row) < 0)
			G_fatal_error (_("Could not read from <%s>"),START);
	
		/* read every cell in the line buffers */
		for (col = 0; col < ncols; col++)
		{
		    d_ELEV = ((DCELL *) inrast_ELEV)[col];
		    d_ASP = ((DCELL *) inrast_ASP)[col];
		    d_START= ((CELL *) inrast_START)[col];

  		    if (!G_is_c_null_value(&d_START)) {
		        /* set new starting point in the array and reallocate space if needed*/
                        if ( (start_count!=0) && (start_count % ALLOC_DIM)==0 ) {
                            ia++;
                            tmp_start = (startpt *)G_realloc(a_start, sizeof(startpt) * ALLOC_DIM *ia );
                            if (tmp_start != NULL ) 
                                a_start = tmp_start;
                        }
                        a_start[start_count].sy = G_row_to_northing(row+0.5,&window);
                        a_start[start_count].sx = G_col_to_easting(col+0.5,&window);
                        a_start[start_count].sz = d_ELEV;
                        a_start[start_count].dir = d_ASP;
                        start_count++;
            	    }
		}
	}
	G_free(inrast_ELEV);
	G_free(inrast_ASP);
	G_free(inrast_START);
	G_close_cell (infd_ELEV);
	G_close_cell (infd_ASP);
	G_close_cell (infd_START);

    	G_debug(3, "source detection completed: selected %d source cells!\n", start_count);

//printf("start count: %i, z: %f", start_count,a_start[start_count-1].sz);
//while(getchar() != 'y') { }

/******************* running operation **************************

	/* open input raster files */
	if ( (infd_ELEV = G_open_cell_old (ELEV, mapset_ELEV)) < 0)
		G_fatal_error (_("Cannot open cell file [%s]"), ELEV);
	if (G_get_cellhd (ELEV, mapset_ELEV, &cellhd) < 0)
		G_fatal_error (_("Cannot read file header of [%s]"), ELEV);		
	/* Allocate input buffer */
	inrast_ELEV	= G_allocate_d_raster_buf();

	/* open output raster */
	if (OUT_F) {
		if ( (outfd_F = G_open_raster_new (OUT_F,CELL_TYPE)) < 0)
			G_fatal_error (_("Could not open <%s>"),OUT_F);
	}
	if (OUT_V) {
		if ( (outfd_V = G_open_raster_new (OUT_V,DCELL_TYPE)) < 0)
			G_fatal_error (_("Could not open <%s>"),OUT_V);
	}
/*M
	if (OUT_P) {
		if ( (outfd_P = G_open_raster_new (OUT_P,DCELL_TYPE)) < 0)
			G_fatal_error (_("Could not open <%s>"),OUT_P);
	}
*/
	if (OUT_C) {
		if ( (outfd_C = G_open_raster_new (OUT_C,DCELL_TYPE)) < 0)
			G_fatal_error (_("Could not open <%s>"),OUT_C);
	}

	/* Allocate output buffer */
	if (OUT_F)
		outrast_F = G_allocate_c_raster_buf();
	if (OUT_V)
		outrast_V = G_allocate_d_raster_buf();
/*M
	if (OUT_P)
		outrast_P = G_allocate_d_raster_buf();
*/
	if (OUT_C) 
		outrast_C = G_allocate_d_raster_buf();

	/* get windows rows & cols */
	nrows	= G_window_rows();
	ncols	= G_window_cols();
	G_get_window(&window); 

    	G_debug(3, "calculating shadows!\n");

	for (row = 0; row < nrows; row++)
	{
		G_percent (row, nrows, 2);

/*M
printf("nrows: %i - ncols: %i", nrows,ncols);
while(getchar() != 'y') { }
*/

		/* read a line input maps into buffers*/	
		if (G_get_d_raster_row (infd_ELEV, inrast_ELEV, row) < 0)
			G_fatal_error (_("Could not read from <%s>"),ELEV);

		/* read every cell in the line buffers */
		for (col = 0; col < ncols; col++)
		{
		    i_freq = 0;
		    i_max = -1;
		    i_min = 1000000;
		    i_sum = 0;
		    d_ELEV = ((DCELL *) inrast_ELEV)[col];
		    if (!G_is_d_null_value(&d_ELEV)) {
		        north = G_row_to_northing(row+0.5,&window);
                        east  = G_col_to_easting(col+0.5,&window);

		        for ( k = 0; k < start_count; k++) {
			    ee = pow(a_start[k].sx-east,2);
			    nn = pow(a_start[k].sy-north,2);
			    Dz = sqrt( ee + nn ) * tan(ANG*PI/180);
			    zi = a_start[k].sz - Dz;

			    Di = atan2(north-a_start[k].sy, east-a_start[k].sx) * 180 / PI;
			    if(a_start[k].dir>180) { a_start[k].dir = a_start[k].dir - 360; }
			    Dhi = fabs( a_start[k].dir - Di );
			    if (Dhi>180) { Dhi = 360-Dhi;}

//printf("Dhi: %f - asp: %f - Di: %f\n", Dhi, a_start[k].dir, Di);


/*M
printf("xS: %f - yS: %f - zS: %f\n",a_start[k].sx,a_start[k].sy,a_start[k].sz);
printf("xi: %f - yi: %f - zi: %f\n",east,north,d_ELEV);
printf("row: %i - col: %i - Ds: %f - zi: %f - d_ELEV: %f - Dz: %f - a_start[k].sz: %f\n", row, col, sqrt( ee + nn ), zi,d_ELEV,Dz,a_start[k].sz);
printf("Dhi: %f - asp: %f - dir: %f\n", Dhi, a_start[k].dir, atan2(north-a_start[k].sy, east-a_start[k].sx) * 180 / PI);
while(getchar() != 'y') { }
*/

//			    if(zi > d_ELEV && Di<a_start[k].dir-(DIR/2) && Di>a_start[k].dir+(DIR/2) ) {
			    if(zi > d_ELEV && Dhi < DIR/2) {
			        Dzi = zi - d_ELEV;
			        i_freq++;
			        i_max = max(i_max,Dzi);
			        if(i_min > Dzi){ i_min = Dzi; };
			        i_sum = i_sum + Dzi;
/*M
printf("FATTO - row: %i - col: %i - freq: %i - i_max: %f - i_min: %f - i_sum: %f Dzi: %f Dhi: %f\n", row, col, i_freq,max(i_max,Dzi),min(i_min,Dzi),i_sum,Dzi,Dhi);
while(getchar() != 'y') { }
*/	
  			    }
		        }
		    }

		    /* write output values */
		    if (OUT_F) {
		         ((CELL *) outrast_F)[col] = i_freq;
		    }

		    if (OUT_V) { /* Vel = fv sqrt(2 g ∆h) */
		        if (strcmp (METHOD, "max") == 0) {
		            if (i_freq == 0) {
				G_set_d_null_value(&outrast_V[col], 1);
			    } else {
			        ((DCELL *) outrast_V)[col] = sqrt(i_max*2*9.80665)*FVEL;
			    } 
		        }
			if (strcmp (METHOD, "min") == 0) {
			    if (i_freq == 0) {
				G_set_d_null_value(&outrast_V[col], 1);
			    } else {
			        ((DCELL *) outrast_V)[col] = sqrt(i_min*2*9.80665)*FVEL;
			    } 
		        }
			if (strcmp (METHOD, "ave") == 0) {
			    if (i_freq == 0) {
				G_set_d_null_value(&outrast_V[col], 1);
			    } else {
			        ((DCELL *) outrast_V)[col] = sqrt((i_sum/i_freq)*2*9.80665)*FVEL;
			    } 
		        }
		    }
/*M
		    if (OUT_P) { // Ep = m g ∆h 
		        if (strcmp (METHOD, "max") == 0) {
			    if (i_freq == 0) {
				G_set_d_null_value(&outrast_P[col], 1);
			    } else {
			        ((DCELL *) outrast_P)[col] = i_max*9.80665*MASS;
			    } 
		        }
			if (strcmp (METHOD, "min") == 0) {
			    if (i_freq == 0) {
				G_set_d_null_value(&outrast_P[col], 1);;
			    } else {
			       ((DCELL *) outrast_P)[col] = i_min*9.80665*MASS;
			    } 
		        }
			if (strcmp (METHOD, "ave") == 0) {
			    if (i_freq == 0) {
				G_set_d_null_value(&outrast_P[col], 1);
			    } else {
			        ((DCELL *) outrast_P)[col] = (i_sum/i_freq)*9.80665*MASS;
			    } 
		        }
		    }
*/
/*M
		    if (OUT_C) { /* Ec = 1/2 m v *
		        if (strcmp (METHOD, "max") == 0) {
			    if (i_freq == 0) {
				G_set_d_null_value(&outrast_C[col], 1);
			    } else {
			        ((DCELL *) outrast_C)[col] = 0.5 * MASS * sqrt(i_max*2*9.80665)*FVEL;
			    } 
		        }
			if (strcmp (METHOD, "min") == 0) {
			    if (i_freq == 0) {
				G_set_d_null_value(&outrast_C[col], 1);
			    } else {
			        ((DCELL *) outrast_C)[col] =  0.5 * MASS * sqrt(i_min*2*9.80665)*FVEL;
			    } 
		        }
			if (strcmp (METHOD, "ave") == 0) {
			    if (i_freq == 0) {
				G_set_d_null_value(&outrast_C[col], 1);
			    } else {
			        ((DCELL *) outrast_C)[col] =  0.5 * MASS * sqrt((i_sum/i_freq)*2*9.80665)*FVEL;
			    } 
		        }
*/
		    if (OUT_C) { /* Ec = m g ∆h*/
		        if (strcmp (METHOD, "max") == 0) {
			    if (i_freq == 0) {
				G_set_d_null_value(&outrast_C[col], 1);
			    } else {
			        ((DCELL *) outrast_C)[col] = MASS * 9.80665 * i_max * pow(FVEL,2);
			    } 
		        }
			if (strcmp (METHOD, "min") == 0) {
			    if (i_freq == 0) {
				G_set_d_null_value(&outrast_C[col], 1);
			    } else {
			        ((DCELL *) outrast_C)[col] = MASS * 9.80665 * i_min * pow(FVEL,2);
			    } 
		        }
			if (strcmp (METHOD, "ave") == 0) {
			    if (i_freq == 0) {
				G_set_d_null_value(&outrast_C[col], 1);
			    } else {
			        ((DCELL *) outrast_C)[col] =  MASS * 9.80665 * i_min * (i_sum/i_freq)  * pow(FVEL,2);
			    } 
		        }

		    }

		} //end col

    		G_debug(3, "writing row %d !\n",row);

		if (OUT_F) {
		    if (G_put_c_raster_row (outfd_F, outrast_F) < 0)
		        G_fatal_error (_("Cannot write to <%s>"),OUT_F);
		}

		if (OUT_V) {
		    if (G_put_d_raster_row (outfd_V, outrast_V) < 0)
		        G_fatal_error (_("Cannot write to <%s>"),OUT_V);
		}

/*M
		if (OUT_P) {
		    if (G_put_d_raster_row (outfd_P, outrast_P) < 0)
		        G_fatal_error (_("Cannot write to <%s>"),OUT_P);
		}
*/
		if (OUT_C) {
		    if (G_put_d_raster_row (outfd_C, outrast_C) < 0)
		        G_fatal_error (_("Cannot write to <%s>"),OUT_C);
		}


	} //end row

    	G_debug(3, "command executed! Bye!\n");

G_free(inrast_ELEV);
G_free(outrast_F);
G_free(outrast_V);
//G_free(outrast_P);
G_free(outrast_C);

G_close_cell (infd_ELEV);
G_close_cell (outfd_F);
G_close_cell (outfd_V);
//G_close_cell (outfd_P);
G_close_cell (outfd_C);

G_free(a_start);

if (OUT_F) {
    G_short_history(OUT_F, "raster", &history);
    G_command_history(&history);
    G_write_history(OUT_F, &history);
}

if (OUT_V) {
    G_short_history(OUT_V, "raster", &history);
    G_command_history(&history);
    G_write_history(OUT_V, &history);
}

/*M
if (OUT_P) {
    G_short_history(OUT_P, "raster", &history);
    G_command_history(&history);
    G_write_history(OUT_P, &history);
}
*/
if (OUT_C) {
    G_short_history(OUT_C, "raster", &history);
    G_command_history(&history);
    G_write_history(OUT_C, &history);
}

exit(EXIT_SUCCESS);

//return (0);
}
