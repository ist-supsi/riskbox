#ifndef _DF_FUNCTIONS_H_
#define _DF_FUNCTIONS_H_

#define min_among_2(a,b) ((a) <= (b) ? (a) : (b))
#define max_among_2(a,b) ((a) >= (b) ? (a) : (b))

    typedef struct
    {
       char *name;
       char *mapset_name;
       int fd;
       RASTER_MAP_TYPE data_type;      /* which raster type CELL_TYPE, FCELL_TYPE, DCELL_TYPE */
       int rows, cols;
       struct History hystory;
       struct Cell_head cellhd;
       CELL **cell_matrix;             /*The data is stored in an one dimensional array internally */
       FCELL **fcell_matrix;           /*The data is stored in an one dimensional array internally */
       DCELL **dcell_matrix;           /*The data is stored in an one dimensional array internally */
       void *in_buf;
       unsigned char *out_buf;
    } raster_struct;

    typedef struct
    {
       int r, c;
    } rw_point;

    typedef struct
    {
       rw_point point;
       float Vel, Sed, Pend, h;
	   int mode;
    } rw_cell;

    typedef struct
    {

       double h[9], w[9], s[9], p[9];
    } rw_mask;


/* ANSI functions prototype */
void SetRastName(raster_struct rst,struct Option opt);
void SetMapset(raster_struct *rst);
void SetOpen(raster_struct *rst);
void AllocateInBuf(raster_struct *rst);
void AllocateOpenOut(raster_struct *rst,RASTER_MAP_TYPE data_type);
void AllocateMatrix (raster_struct *rst,int nrows, int ncols,double init);
void GetInRow(raster_struct *rst, int row);
void GetCellValue2Matrix(raster_struct *rst, int row, int col) ;
void SetMask(rw_mask *m, raster_struct *dem,raster_struct *pdir, int r,int c);
void SetMaskWeight(rw_mask *m, raster_struct *pdir,int r,int c);
int SetDirectionMode(rw_mask *m, double lim, int idx_from);
int maximum( double [] );
int GetNewRowIdx(int idx,int r);
int GetNewColIdx(int idx,int c);
int GetIdx(int rfrom, int cfrom);
int GetLower(rw_mask *m,int idx_from);
int GetIdxCCW(int idx_to);
int GetIdxCW(int idx_to);
void free_close_rast_struct(raster_struct *rst, int nrows, int ncols) ;
double GetValue(raster_struct *rst,int r,int c);
void SetValue(raster_struct *rst,int r,int c,double value);
void  lower ( double* h, int* index );
int  count_lower (double* values );
int get_cell_mode2(rw_mask *m, double slim, double a, int delta_c, int delta_r);
int get_cell_mode3(rw_mask *m, double slim, double a, int delta_r, int delta_c);


#endif /* _DF_FUNCTIONS_H_ */
