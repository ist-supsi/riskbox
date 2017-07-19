#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <grass/gis.h>
#include <grass/glocale.h>
#include <math.h>
#include "nrutil.h"
#include "dfFunctions.h"
#include "malloc.h"

void SetMapset(raster_struct *rst){
	rst->mapset_name = G_find_cell2(rst->name, "");
	if (rst->mapset_name == NULL)
		G_fatal_error(_("cell file [%s] not found"), rst->name);
	return;
}

void SetOpen(raster_struct *rst){
	rst->data_type = G_raster_map_type(rst->name, rst->mapset_name);
	if ((rst->fd = G_open_cell_old(rst->name, rst->mapset_name)) < 0)
	G_fatal_error(_("Cannot open raster map [%s]"), rst->name);
	if (G_get_cellhd(rst->name, rst->mapset_name, &rst->cellhd) < 0)
	G_fatal_error(_("Cannot read file header of [%s]"), rst->name);
	return;
}

void AllocateInBuf(raster_struct *rst){
	rst->in_buf   = G_allocate_raster_buf(rst->data_type);
	return;
}

void AllocateOpenOut(raster_struct *rst,RASTER_MAP_TYPE data_type){
	rst->data_type=data_type;
	if (G_legal_filename(rst->name) < 0)
	 G_fatal_error(_("[%s] is an illegal name"), rst->name);
	if (data_type==CELL_TYPE) {
		rst->out_buf = G_allocate_raster_buf(CELL_TYPE);
		if ((rst->fd = G_open_raster_new(rst->name, CELL_TYPE)) < 0)
			 G_fatal_error(_("Could not open <%s>"), rst->name);
	} else if (data_type==FCELL_TYPE) {
		rst->out_buf = G_allocate_raster_buf(FCELL_TYPE);
		if ((rst->fd = G_open_raster_new(rst->name, FCELL_TYPE)) < 0)
			 G_fatal_error(_("Could not open <%s>"), rst->name);
	} else if (data_type==DCELL_TYPE) {
		rst->out_buf = G_allocate_raster_buf(DCELL_TYPE);
		if ((rst->fd = G_open_raster_new(rst->name, DCELL_TYPE)) < 0)
			 G_fatal_error(_("Could not open <%s>"), rst->name);
	} else {
		G_fatal_error(_("Could not detect <%s> data type"), rst->name);
	}
	return;
}

void AllocateMatrix (raster_struct *rst,int nrows, int ncols,double init){
	int r,c;
	rst->rows = nrows;
	rst->cols = ncols;
	if(rst->data_type==CELL_TYPE) {
		rst->cell_matrix = Cmatrix(0, nrows, 0, ncols);
		for (r=0;r<nrows;r++){
			for (c=0;c<ncols;c++){
				rst->cell_matrix[r][c]=(CELL) init;
			}
		}
	} else if (rst->data_type==FCELL_TYPE){
		rst->fcell_matrix = Fmatrix(0, nrows, 0, ncols);
		for (r=0;r<nrows;r++){
			for (c=0;c<ncols;c++){
				rst->fcell_matrix[r][c]=(FCELL) init;
			}
		}
	} else if (rst->data_type==DCELL_TYPE){
		rst->dcell_matrix = Dmatrix(0, nrows, 0, ncols);
		for (r=0;r<nrows;r++){
			for (c=0;c<ncols;c++){
				rst->dcell_matrix[r][c]=(DCELL) init;
			}
		}
	} else {
		G_fatal_error(_("Could not detect <%s> data type and allocate matrix"), rst->name);
	}
	return;
}

void GetInRow(raster_struct *rst, int row) {
	if (G_get_raster_row(rst->fd, rst->in_buf, row, rst->data_type) < 0)
	    G_fatal_error(_("Could not read from <%s>"), rst->name);
    return;
}

void free_close_rast_struct(raster_struct *rst, int nrows, int ncols) {    
	free(rst->in_buf);
	//G_free(rst->in_buf);
	//while(getchar() != 'y') { printf("...in buf cleaned!\n"); }   	     
	free(rst->out_buf);
	//G_free(rst->out_buf);
//while(getchar() != 'y') { printf("...out buf cleaned!\n"); }   	     
	switch (rst->data_type) {
        case CELL_TYPE:
//while(getchar() != 'y') { printf("...CELL RASTER!\n"); }               
	        free_Cmatrix(rst->cell_matrix,0, nrows, 0, ncols);
//while(getchar() != 'y') { printf("...3 DONE!\n"); }
        break;
        case FCELL_TYPE:
//while(getchar() != 'y') { printf("...FCELL RASTER!\n"); }               
		    free_Fmatrix(rst->fcell_matrix,0, nrows, 0, ncols);
//while(getchar() != 'y') { printf("...3 DONE!\n"); }		            
        break;
        case DCELL_TYPE:
//while(getchar() != 'y') { printf("...DCELL RASTER!\n"); }               
		    free_Dmatrix(rst->dcell_matrix,0, nrows, 0, ncols);
//while(getchar() != 'y') { printf("...3 DONE!\n"); }		            
        break;
    }
    G_close_cell(rst->fd);
    return;
}
    
    
void GetCellValue2Matrix(raster_struct *rst, int row, int col) {
	/* use different function for each data type */
    switch (rst->data_type) {
        case CELL_TYPE:
            rst->cell_matrix[row][col] = ((CELL *) rst->in_buf)[col];
        break;
        case FCELL_TYPE:
            rst->fcell_matrix[row][col] = ((FCELL *) rst->in_buf)[col];
        break;
        case DCELL_TYPE:
            rst->dcell_matrix[row][col] = ((DCELL *) rst->in_buf)[col];
        break;
    }
    return;
}

void SetMask(rw_mask *m, raster_struct *dem,raster_struct *pdir, int r,int c) {
    int i;
	// problema sui bordi if(cr!=0 && cc!=0)

	/* set mask height */
	if(dem->data_type==CELL_TYPE) {
		m->h[0] = (double) dem->cell_matrix[r-1][c-1];     /* N-W */
		m->h[1] = (double) dem->cell_matrix[r-1][c];       /* N */
		m->h[2] = (double) dem->cell_matrix[r-1][c+1];   /* N-E */
		m->h[3] = (double) dem->cell_matrix[r][c-1];       /* W */
		m->h[4] = (double) dem->cell_matrix[r][c];          /* center cell */
		m->h[5] = (double) dem->cell_matrix[r][c+1];     /* E */
		m->h[6] = (double) dem->cell_matrix[r+1][c-1];  /* S-W */
		m->h[7] = (double) dem->cell_matrix[r+1][c];     /* S */
		m->h[8] = (double) dem->cell_matrix[r+1][c+1]; /* S-E */
	} else if (dem->data_type==FCELL_TYPE) {
		m->h[0] = (double) dem->fcell_matrix[r-1][c-1];
		m->h[1] = (double) dem->fcell_matrix[r-1][c];
		m->h[2] = (double) dem->fcell_matrix[r-1][c+1];
		m->h[3] = (double) dem->fcell_matrix[r][c-1];
		m->h[4] = (double) dem->fcell_matrix[r][c];
		m->h[5] = (double) dem->fcell_matrix[r][c+1];
		m->h[6] = (double) dem->fcell_matrix[r+1][c-1];
		m->h[7] = (double) dem->fcell_matrix[r+1][c];
		m->h[8] = (double) dem->fcell_matrix[r+1][c+1];
	} else {
		m->h[0] = (double) dem->dcell_matrix[r-1][c-1];
		m->h[1] = (double) dem->dcell_matrix[r-1][c];
		m->h[2] = (double) dem->dcell_matrix[r-1][c+1];
		m->h[3] = (double) dem->dcell_matrix[r][c-1];
		m->h[4] = (double) dem->dcell_matrix[r][c];
		m->h[5] = (double) dem->dcell_matrix[r][c+1];
		m->h[6] = (double) dem->dcell_matrix[r+1][c-1];
		m->h[7] = (double) dem->dcell_matrix[r+1][c];
		m->h[8] = (double) dem->dcell_matrix[r+1][c+1];
	}

	/* set mask slope */
	for (i=0; i<9; i++) {
		if(i==0 || i==2 || i==6 || i==8) {
			//m->s[i] = atan( (m->h[4] - m->h[i]) / (sqrt(2)*(dem->cellhd.ew_res)) );
			m->s[i] = (m->h[4] - m->h[i]) / (sqrt(2)*(dem->cellhd.ew_res)) ;
		} else {
			//m->s[i] = atan( (m->h[4] - m->h[i]) / (dem->cellhd.ew_res) );
			m->s[i] = (m->h[4] - m->h[i]) / (dem->cellhd.ew_res);
		}
	}

	//print debug 
	/*for (i=0; i<9; i++) {
      		printf("m->h[%d]=%lf - ",i,m->h[i]);
      		printf("m->s[%d]=%lf - ",i,m->s[i]);
      		//printf("res=%f - ",dem->cellhd.ew_res);
       		//printf("PDIR=%lf - ", pdir->fcell_matrix[r-1][c-1]);
   	}
   	printf("\n");*/
       

	/* set mask weighted slope */
	if(pdir->data_type==CELL_TYPE) {
		m->w[0] = (double) pdir->cell_matrix[r-1][c-1] * m->s[0];
		m->w[1] = (double) pdir->cell_matrix[r-1][c] * m->s[1];
		m->w[2] = (double) pdir->cell_matrix[r-1][c+1] * m->s[2];
		m->w[3] = (double) pdir->cell_matrix[r][c-1] * m->s[3];
		m->w[4] = (double) pdir->cell_matrix[r][c] * m->s[4];
		m->w[5] = (double) pdir->cell_matrix[r][c+1] * m->s[5];
		m->w[6] = (double) pdir->cell_matrix[r+1][c-1] * m->s[6];
		m->w[7] = (double) pdir->cell_matrix[r+1][c] * m->s[7];
		m->w[8] = (double) pdir->cell_matrix[r+1][c+1] * m->s[8];
	} else if (pdir->data_type==FCELL_TYPE) {
		m->w[0] = (double) pdir->fcell_matrix[r-1][c-1] * m->s[0];
		m->w[1] = (double) pdir->fcell_matrix[r-1][c] * m->s[1];
		m->w[2] = (double) pdir->fcell_matrix[r-1][c+1] * m->s[2];
		m->w[3] = (double) pdir->fcell_matrix[r][c-1] * m->s[3];
		m->w[4] = (double) pdir->fcell_matrix[r][c] * m->s[4];
		m->w[5] = (double) pdir->fcell_matrix[r][c+1] * m->s[5];
		m->w[6] = (double) pdir->fcell_matrix[r+1][c-1] * m->s[6];
		m->w[7] = (double) pdir->fcell_matrix[r+1][c] * m->s[7];
		m->w[8] = (double) pdir->fcell_matrix[r+1][c+1] * m->s[8];
	} else {
		m->w[0] = (double) pdir->dcell_matrix[r-1][c-1] * m->s[0];
		m->w[1] = (double) pdir->dcell_matrix[r-1][c] * m->s[1];
		m->w[2] = (double) pdir->dcell_matrix[r-1][c+1] * m->s[2];
		m->w[3] = (double) pdir->dcell_matrix[r][c-1] * m->s[3];
		m->w[4] = (double) pdir->dcell_matrix[r][c] * m->s[4];
		m->w[5] = (double) pdir->dcell_matrix[r][c+1] * m->s[5];
		m->w[6] = (double) pdir->dcell_matrix[r+1][c-1] * m->s[6];
		m->w[7] = (double) pdir->dcell_matrix[r+1][c] * m->s[7];
		m->w[8] = (double) pdir->dcell_matrix[r+1][c+1] * m->s[8];
	}
	   
	/* set mask weight */
	if(pdir->data_type==CELL_TYPE) {
		m->p[0] = (double) pdir->cell_matrix[r-1][c-1];
		m->p[1] = (double) pdir->cell_matrix[r-1][c];
		m->p[2] = (double) pdir->cell_matrix[r-1][c+1];
		m->p[3] = (double) pdir->cell_matrix[r][c-1];
		m->p[4] = (double) pdir->cell_matrix[r][c];
		m->p[5] = (double) pdir->cell_matrix[r][c+1];
		m->p[6] = (double) pdir->cell_matrix[r+1][c-1];
		m->p[7] = (double) pdir->cell_matrix[r+1][c];
		m->p[8] = (double) pdir->cell_matrix[r+1][c+1];
	} else if (pdir->data_type==FCELL_TYPE) {
		m->p[0] = (double) pdir->fcell_matrix[r-1][c-1];
		m->p[1] = (double) pdir->fcell_matrix[r-1][c];
		m->p[2] = (double) pdir->fcell_matrix[r-1][c+1];
		m->p[3] = (double) pdir->fcell_matrix[r][c-1];
		m->p[4] = (double) pdir->fcell_matrix[r][c];
		m->p[5] = (double) pdir->fcell_matrix[r][c+1];
		m->p[6] = (double) pdir->fcell_matrix[r+1][c-1];
		m->p[7] = (double) pdir->fcell_matrix[r+1][c];
		m->p[8] = (double) pdir->fcell_matrix[r+1][c+1];
	} else {
		m->p[0] = (double) pdir->dcell_matrix[r-1][c-1];
		m->p[1] = (double) pdir->dcell_matrix[r-1][c];
		m->p[2] = (double) pdir->dcell_matrix[r-1][c+1];
		m->p[3] = (double) pdir->dcell_matrix[r][c-1];
		m->p[4] = (double) pdir->dcell_matrix[r][c];
		m->p[5] = (double) pdir->dcell_matrix[r][c+1];
		m->p[6] = (double) pdir->dcell_matrix[r+1][c-1];
		m->p[7] = (double) pdir->dcell_matrix[r+1][c];
		m->p[8] = (double) pdir->dcell_matrix[r+1][c+1];
	}
	return;
} // end SetMask


//??????????????????????
void SetMaskWeight(rw_mask *m, raster_struct *pdir,int r,int c) {
	int i;
	// problema sui bordi if(cr!=0 && cc!=0)
	if(pdir->data_type==CELL_TYPE) {
		m->w[0] = pdir->cell_matrix[r-1][c-1] * m->s[0];
		m->w[1] = pdir->cell_matrix[r-1][c] * m->s[1];
		m->w[2] = pdir->cell_matrix[r-1][c+1] * m->s[2];
		m->w[3] = pdir->cell_matrix[r][c-1] * m->s[3];
		m->w[4] = pdir->cell_matrix[r][c] * m->s[4];
		m->w[5] = pdir->cell_matrix[r][c+1] * m->s[5];
		m->w[6] = pdir->cell_matrix[r+1][c-1] * m->s[6];
		m->w[7] = pdir->cell_matrix[r+1][c] * m->s[7];
		m->w[8] = pdir->cell_matrix[r+1][c+1] * m->s[8];
	} else if (pdir->data_type==FCELL_TYPE) {
		m->w[0] = pdir->fcell_matrix[r-1][c-1] * m->s[0];
		m->w[1] = pdir->fcell_matrix[r-1][c] * m->s[1];
		m->w[2] = pdir->fcell_matrix[r-1][c+1] * m->s[2];
		m->w[3] = pdir->fcell_matrix[r][c-1] * m->s[3];
		m->w[4] = pdir->fcell_matrix[r][c] * m->s[4];
		m->w[5] = pdir->fcell_matrix[r][c+1] * m->s[5];
		m->w[6] = pdir->fcell_matrix[r+1][c-1] * m->s[6];
		m->w[7] = pdir->fcell_matrix[r+1][c] * m->s[7];
		m->w[8] = pdir->fcell_matrix[r+1][c+1] * m->s[8];
	} else {
		m->w[0] = pdir->dcell_matrix[r-1][c-1] * m->s[0];
		m->w[1] = pdir->dcell_matrix[r-1][c] * m->s[1];
		m->w[2] = pdir->dcell_matrix[r-1][c+1] * m->s[2];
		m->w[3] = pdir->dcell_matrix[r][c-1] * m->s[3];
		m->w[4] = pdir->dcell_matrix[r][c] * m->s[4];
		m->w[5] = pdir->dcell_matrix[r][c+1] * m->s[5];
		m->w[6] = pdir->dcell_matrix[r+1][c-1] * m->s[6];
		m->w[7] = pdir->dcell_matrix[r+1][c] * m->s[7];
		m->w[8] = pdir->dcell_matrix[r+1][c+1] * m->s[8];
	}
	return;
}

double GetValue(raster_struct *rst,int r,int c) {
	double value;
	if(rst->data_type==CELL_TYPE) {
		value = (double) rst->cell_matrix[r][c];
	} else if (rst->data_type==FCELL_TYPE) {
		value = (double) rst->fcell_matrix[r][c];
	} else if (rst->data_type==DCELL_TYPE) {
		value = (double) rst->dcell_matrix[r][c];
	} else {
		G_fatal_error(_("Could not detect <%s> data type in GetValue"), rst->name);
	}
	return value;
}

void SetValue(raster_struct *rst,int r,int c,double value) {
	if(rst->data_type==CELL_TYPE) {
		rst->cell_matrix[r][c] = (CELL) value;
	} else if (rst->data_type==FCELL_TYPE) {
		rst->fcell_matrix[r][c] = (FCELL) value;
	} else if (rst->data_type==DCELL_TYPE) {
		rst->dcell_matrix[r][c] = (DCELL) value;
	} else {
		G_fatal_error(_("Could not detect <%s> data type in GetValue"), rst->name);
	}
	return;
}


int SetDirectionMode(rw_mask *m, double lim, int idx_from){
	int mode, i;
	/* mode = 1 -> direzione fissa SFD */
	/* mode = 2 -> random walk MFD */
	/* mode = 3 -> pianura MFD */
	/* mode = 4 -> sink */
	for (i=0; i<9; i++) {
		//while(getchar() != 'y') { printf("pendenza-slim: %lf - %lf!\n",m->s[i],lim); } 
		//if(fabs(m->w[i]) > lim && i!=4 && i!=idx_from) {
		if(m->w[i] >= lim && i!=4 && i!=idx_from) {
			//printf("mode1\n");
			mode = 1;
			return mode;
		}
	}
	for (i=0; i<9; i++) {
		if(m->h[i] < m->h[4] && i!=4 && i!=idx_from) {
			//printf("mode2\n");
			mode = 2;
			return mode;
		}
	}
	for (i=0; i<9; i++) {
		if(m->h[i] == m->h[4] && i!=4  && i!=idx_from) {
			//printf("mode3\n");
			mode = 3;             
			return mode;
		}
	}
	for (i=0; i<9; i++) {
		if(m->h[i] > m->h[4] && i!=4  && i!=idx_from) {
			//printf("mode3\n");
			mode = 4;             
			return mode;
		}
	}
	printf("c'e' qualcosa che non va!\n");
	while(!getchar()){ }
	return mode;
}

int  maximum( double values[9] ) {
	int  i, index;
	double max_value;
	max_value = values[4];
	index=4;
	for( i = 0; i < 9; i++ ) {
		if( values[i] > max_value ) {
			max_value = values[i];
			index = i;
		}
	}
	return index;
}

void  lower ( double* h, int* index ) {
    int i;		
	for (i=0; i < 9; i++ ) {
		index[i]=-1;
	}
	for( i = 0; i < 9; i++ ) {
		if( h[i] < h[4] ) {
			index[i] = i;
		}
	}
	return;
}


int get_cell_mode2(rw_mask *m, double slim, double a, int delta_r, int delta_c) {
	int extracted, idx, to_idx, from_idx;
	int from; /* index of the cell the flow came from */
	int i, j, count=0, to_r, to_c;
	double gamma_i[9];
	double gamma_max=0, prob_tot=0;
	double *probability, *upper_bound;
	int *selected_idx;
	double max_slope=-99;
	
	
	
	/* calculate gamma i and max */
	for( i = 0; i < 9; i++ ) {
		gamma_i[i] = m->s[i] / slim;
		max_slope = max_among_2(m->s[i],max_slope);
		//printf("gamma_i[%d]=%lf - ", i, gamma_i[i]);
		//printf("m->s[%d]=%lf - ", i, m->s[i]);
	}
	gamma_max = max_slope/ slim;
	//print debug 
	//printf("\n- calculated gamma max = %lf --> ", gamma_max);
	//printf("calculated gamma_max^alfa = %lf", pow(gamma_max,a) );
	from_idx = GetIdx(delta_r , delta_c);
	
	/* count selectable cells */
	for( i = 0; i < 9; i++ ) {
		if( gamma_i[i] >= pow(gamma_max,a) && i!=4 && gamma_max>0 && gamma_max<1 && i!=from_idx ) {
	    		count++;
		}
		
		//debug: aggiunta condizione gamma_max>=1)
		/*if( gamma_i[i] >= pow(gamma_max,a) && i!=4 && gamma_max<1) {
	    		count++;
	    	} else if (gamma_i[i]==gamma_max && i!=4 && gamma_max>=1){
	    		count++;
		}*/
	    	
	 }
	 
	 /*for (i=0; i<9; i++) {
      		printf("m->h[%d]=%lf - ",i,m->h[i]);
      		printf("m->s[%d]=%lf - ",i,m->s[i]);
      		//printf("res=%f - ",dem->cellhd.ew_res);
       		//printf("PDIR=%lf - ", pdir->fcell_matrix[r-1][c-1]);
   	}
   	printf("\n");*/
	
	//printf("\n --> count selectable cells = %d\n",count);

	/* allocate array */
	selected_idx = (int*)malloc(sizeof(int)*(count+1));
	probability = (double*)malloc(sizeof(double)*(count+1));
	upper_bound = (double*)malloc(sizeof(double)*(count+1));
	upper_bound[0] = -0.01;
	upper_bound[count] = 100;
	selected_idx[0] = -1;
	probability[0] = 0;
	//printf ("pippo2\n");
	/* calculate selected cells index and probability adding 0.5 to preferred direction*/
	j=1; 
	//printf("\n-  delta_r=%d, delta_c=%d",delta_r,delta_c);
	to_r = delta_r * -1;
	to_c = delta_c * -1;
	//printf("\n-  to_r=%d, to_c=%d\n",to_r,to_c);
	to_idx = GetIdx(to_r,to_c);
	//to_idx = GetIdx(delta_r * -1,delta_c * -1);
	//printf("-  delta_r=%d, delta_c=%d, to_idx=%d",delta_r,delta_c,to_idx);
	for( i = 0; i < 9; i++ ) {
		//if( (gamma_i[i] >= pow(gamma_max,a) || gamma_i[i] == gamma_max)  && i!=4 ) {
		if( gamma_i[i] >= pow(gamma_max,a)  && i!=4) {
			selected_idx[j] = i;
			if (i==to_idx) {
		        	probability[j] = m->s[i] * 1.5 * m->w[i];
	    	    	} else {
	        		probability[j] = m->s[i] * m->w[i];
	    	    	}
	    	    	//printf("\nprobability[%d]=%lf - ", j, probability[j]);
		    	//printf("\nselected_idx[%d]=%d - ", j, selected_idx[j]);
	    	    	j++; 
		}
    	}    
    
    	/* calculate standardization value */
    	for( i = 1; i < count+1; i++ ) {
        	prob_tot += probability[i];
        	//printf("i=%d,selected_idx[i]=%d,probability[i]=%f,prob_tot)%f\n",i,selected_idx[i],probability[i],prob_tot);
    	}
    	
	//printf("\n-  calculate standardization value: %lf\n",prob_tot);
	/*	
	printf("gamma_max=%lf\n",gamma_max);
     	/* calculate bounduary */
    	for( i = 1; i < count+1; i++ ) {
        	upper_bound[i] = upper_bound[i-1] + probability[i]*100.01/prob_tot;
        	//printf("upper_bound[%d]=%lf - ", i, upper_bound[i]);
        	//printf("selected_idx[%d]=%d - ", i, selected_idx[i]);
        	//printf("m->s[%d]=%lf - ",selected_idx[i],m->s[selected_idx[i]]);
        	//printf("gamma_i[%d]=%lf\n", selected_idx[i],gamma_i[selected_idx[i]]);
    	}
	

    /* extract random value between 0-100 */
    //extracted = rand() / (RAND_MAX / 101 + 1);
    extracted = rand()%100;
    //printf("pippo4\n");
    //printf("\n-  extracted random=%d",extracted);    
    /* set selected cell index */
    for( i = 1; i < count+1; i++ ) {
	//while(getchar() != 'y') { 		
        //printf("\n-  selected_idx[%d]=%d \n", i, selected_idx[i]);
	//printf("\n-  upper_bound[%d]=%lf \n", i, upper_bound[i]);
	//printf("\n-  upper_bound[%d-1]=%lf \n", i-1, upper_bound[i-1]);
	//}
	//printf("i=%d,upper_bound[i]=%f\n",i,upper_bound[i]);
        if(extracted <= upper_bound[i] && extracted > upper_bound[i-1]) {
            idx = selected_idx[i];
            //printf("-set selected cell index=%d, extracted=%d\n",idx, extracted);
        }
    }
    //while(!getchar()){ }
    //if(idx==8 || idx==6) { while(getchar() != 'y') {} }

    //G_free(selected_idx);
    free(selected_idx);
    //G_free(probability);
    free(probability);
    //G_free(upper_bound);
    free(upper_bound);
    if(idx==4){
    	printf("no!!!!!!!");
    }
    return idx;
} 

int get_cell_mode3(rw_mask *m, double slim, double a, int delta_r, int delta_c) {
    int extracted, idx, to_idx;
    int from_idx,idx_cw,idx_ccw; /* index of the cell the flow came from */
    int i, j, count=0;
    double prob_tot=0;
    double *probability, *upper_bound;
	int *selected_idx;
    	//printf("\n -  delta_r=%d, delta_c=%d\n",delta_r,delta_c);
	from_idx = GetIdx(delta_r , delta_c);
    /* count selectable cells */
    for( i = 0; i < 9; i++ ) {
        if( m->h[4] == m->h[i] && i!=4 && i!=from_idx ) {
            count++;
        }
    }

    //printf("\n- count selectable cells = %d\n",count);
    
    /* allocate array */
    selected_idx = (int*)malloc(sizeof(int)*(count+1));
    probability = (double*)malloc(sizeof(double)*(count+1));
    upper_bound = (double*)malloc(sizeof(double)*(count+1));
    
    upper_bound[0] = -0.01;
    upper_bound[count] = 100;
    selected_idx[0] = -1;
    probability[0] = 0;

    /* calculate selected cells index and probability adding 0.5 to preferred direction*/
    j=1; 
    to_idx = GetIdx(delta_r * -1 , delta_c * -1 );
    idx_cw = GetIdxCW(to_idx);
    idx_ccw = GetIdxCCW(to_idx);
    
	//printf("\n-  delta_r=%d, delta_c=%d, to_idx=%d\n",delta_r,delta_c,to_idx);
	// modifica Roberto da testare
	for( i = 0; i < 9; i++ ) {
		if( m->h[4] == m->h[i] && i!=4 && i!=from_idx ) {
            		selected_idx[j] = i;
			if(i==to_idx) {
				probability[j] = pow(count,-1) * 1.5* m->p[i];
			} else if (i==idx_cw) {
				probability[j] = pow(count,-1) * 0.5 * m->p[i];
			} else if (i==idx_ccw) {
				probability[j] = pow(count,-1) * 0.5 * m->p[i];
			} else {
				/* forse bisognerebbe impedire che finisca qua se anche ci fosse un muro prima si ferma qua e sedimenta, 
				/ al passo dopo o alcuni andrebbe avanti in un altra direzione */
				probability[j] = pow(count,-1) * 0.1 * m->p[i];
			}
			j++; 
        	}
        	/*if( m->h[4] == m->h[i] && i!=4 && i!=from_idx ) {
            		selected_idx[j] = i;
			if(i==to_idx) {
				probability[j] = 70 * m->p[i];
			} else if (i==idx_cw) {
				probability[j] = 10 * (100/count) * m->p[i];
			} else if (i==idx_ccw) {
				probability[j] = 10 * (100/count) * m->p[i];
			} else {
				probability[j] = m->p[i];
			}
		
			//printf("\nweight[%d]=%lf - ", i, m->p[i]);
            		//printf("\nprobability[%d]=%lf - ", j, probability[j]);
			//printf("\nselected_idx[%d]=%d - ", j, selected_idx[j]);
            		j++; 
        	}*/
    	}
    	    
    /* calculate standardization value */
    for( i = 1; i < count+1; i++ ) {
        prob_tot = prob_tot + probability[i];
        //printf("i=%d,selected_idx[i]=%d,probability[i]=%f,prob_tot)%f\n",i,selected_idx[i],probability[i],prob_tot);
    }
     
     
    //printf("\n-  calculate standardization value: %lf\n",prob_tot);
	
     /* calculate bounduary */
    for( i = 1; i < count+1; i++ ) {
        upper_bound[i] = upper_bound[i-1] + (probability[i]*100.01)/prob_tot;
        //printf("upper_bound[%d]=%lf \n ", i, upper_bound[i]);
        //printf("selected_idx[%d]=%d \n", i, selected_idx[i]);
    }
    //printf("ecco4 RAND_MAX=%d\n",RAND_MAX);
    /* extract random value between 0-100 */
    //extracted = rand() / (RAND_MAX / 101 + 1);
    extracted = rand()%100;
    //printf("\n-  extracted random=%d",extracted);
    //printf("ecco5\n");   
    /* set selected cell index */
    for( i = 1; i < count+1; i++ ) {
        //printf("\n-  selected_idx[%d]=%d \n", i, selected_idx[i]);
	//printf("\n-  upper_bound[%d]=%lf \n", i, upper_bound[i]);
	//printf("\n-  upper_bound[%d-1]=%lf \n", i-1, upper_bound[i-1]);
	if(extracted <= upper_bound[i] && extracted > upper_bound[i-1]) {
            idx = selected_idx[i];
            //printf("-set selected cell index=%d, extracted=%d\n",idx, extracted);
            //printf("\n-  set selected cell index=%d",idx);
        } 
    }
   //while(!getchar()){ }
//if(idx==8 || idx==6) { while(getchar() != 'y') {} }

    //G_free(selected_idx);
    free(selected_idx);
    //G_free(probability);
    free(probability);
    //G_free(upper_bound);
    free(upper_bound);
    return idx;
} 


int get_cell_mode4(rw_mask *m, double slim, double a, int delta_r, int delta_c) {
    int extracted, idx, to_idx;
    int from; /* index of the cell the flow came from */
    int i, j, count=0;
    double prob_tot=0;
    double *probability, *upper_bound;
	int *selected_idx;
    
    /* count selectable cells */
    for( i = 0; i < 9; i++ ) {
        if( m->h[4] == m->h[i] && i!=4 ) {
            count++;
        }
    }

//printf("\n- count selectable cells = %d",count);
    
    /* allocate array */
    selected_idx = (int*)malloc(sizeof(int)*(count+1));
    probability = (double*)malloc(sizeof(double)*(count+1));
    upper_bound = (double*)malloc(sizeof(double)*(count+1));
    
    upper_bound[0] = -0.01;
    upper_bound[count] = 100;
	selected_idx[0] = -1;
	probability[0] = 0;

    /* calculate selected cells index and probability adding 0.5 to preferred direction*/
    j=1; 
    to_idx = GetIdx(delta_r * -1 , delta_c * -1 );
	//printf("\n-  delta_r=%d, delta_c=%d, to_idx=%d",delta_r,delta_c,to_idx);
    
	for( i = 0; i < 9; i++ ) {
        if( m->h[4] == m->h[i] && i!=4 ) {
            selected_idx[j] = i;
            if (i==to_idx) {
                probability[j] = (100/count) * (m->p[i]*10);
            } else {
                probability[j] = (100/count) * m->p[i];
            }
			//printf("\nweight[%d]=%lf - ", i, m->p[i]);
            //printf("\nprobability[%d]=%lf - ", j, probability[j]);
			//printf("\nselected_idx[%d]=%d - ", j, selected_idx[j]);
            j++; 
        }
    }    
    
    /* calculate standardization value */
    for( i = 1; i < count+1; i++ ) {
        prob_tot = prob_tot + probability[i];
    }
	
//printf("\n-  calculate standardization value: %lf\n",prob_tot);
	
     /* calculate bounduary */
    for( i = 1; i < count+1; i++ ) {
        upper_bound[i] = upper_bound[i-1] + (probability[i]*100.01)/prob_tot;
        //printf("upper_bound[%d]=%lf - ", i, upper_bound[i]);
        //printf("selected_idx[%d]=%d \n", i, selected_idx[i]);
    }
 
    /* extract random value between 0-100 */
    extracted = rand() / (RAND_MAX / 101 + 1);
	//printf("\n-  extracted random=%d",extracted);    
    /* set selected cell index */
    for( i = 1; i < count+1; i++ ) {
        //printf("\n-  selected_idx[%d]=%d \n", i, selected_idx[i]);
		//printf("\n-  upper_bound[%d]=%lf \n", i, upper_bound[i]);
		//printf("\n-  upper_bound[%d-1]=%lf \n", i-1, upper_bound[i-1]);

        if(extracted <= upper_bound[i] && extracted > upper_bound[i-1]) {
            idx = selected_idx[i];
			//printf("\n-  set selected cell index=%d",idx);
        }
    }

//if(idx==8 || idx==6) { while(getchar() != 'y') {} }

    //G_free(selected_idx);
    free(selected_idx);
    //G_free(probability);
    free(probability);
    //G_free(upper_bound);
    free(upper_bound); 
    return idx;
} 


int GetNewRowIdx(int idx,int r) {
    int row;
    if(idx==0 || idx==1 || idx==2)
        row = r - 1;
    if(idx==3 || idx==4 || idx==5)
        row = r ;
    if(idx==6 || idx==7 || idx==8)
        row = r + 1;
    return row;
}

int GetNewColIdx(int idx,int c) {
      int col;
      if(idx==0 || idx==3 || idx==6)
          col = c - 1;
      if(idx==1 || idx==4 || idx==7)
          col = c ;
      if(idx==2 || idx==5 || idx==8)
          col = c + 1;
      return col;    
}

int GetIdx(int rfrom, int cfrom) {
    int idx;
    if (rfrom==-1 && cfrom ==-1) {
        idx=0;
    } else if (rfrom==-1 && cfrom ==0) {
        idx =1;        
    } else if (rfrom==-1 && cfrom ==1) {
        idx =2;
    } else if (rfrom==0 && cfrom ==-1) {
        idx =3;
    } else if (rfrom==0 && cfrom ==0) {
        idx =4;
    } else if (rfrom==0 && cfrom ==1) {
        idx =5;
    } else if (rfrom==1 && cfrom ==-1) {
        idx =6;
    } else if (rfrom==1 && cfrom ==0) {
        idx =7;
    } else if (rfrom==1 && cfrom ==1) {
        idx =8;
    } else {
        G_fatal_error(_("Could not get index in mask, values must be between -1 and 1"));
    }
    return idx;
}

int GetLower(rw_mask *m,int idx_from) {
	int i, index;
	double min=100000;
	for( i = 0; i < 9; i++ ) {
		if(m->h[i]<min && i!=4 && i!=idx_from) {
			min = m->h[i];
			index=i;
		}
	}
	return index;
}

int GetIdxCW(int idx_to) {
	int index;
	switch (idx_to) {
		case 0:
			index=1;
		break;
		case 1:
			index=2;
		break;
		case 2:
			index=5;
		break;
		case 5:
			index=8;
		break;
		case 8:
			index=7;
		break;
		case 7:
			index=6;
		break;
		case 6:
			index=3;
		break;
		case 3:
			index=0;
		break;
	}
	return index;
}

int GetIdxCCW(int idx_to) {
	int index;
	switch (idx_to) {
		case 0:
			index=3;
		break;
		case 1:
			index=0;
		break;
		case 2:
			index=1;
		break;
		case 5:
			index=2;
		break;
		case 8:
			index=5;
		break;
		case 7:
			index=8;
		break;
		case 6:
			index=7;
		break;
		case 3:
			index=6;
		break;
	}
	return index;
}
