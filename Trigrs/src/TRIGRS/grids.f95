! By Rex L. Baum, 1 April 2004
module grids
	integer,allocatable:: pf2(:),indx(:),nxt(:),nv(:),nvu(:)
	integer,allocatable:: dsctr(:),dsc(:),zo(:) 
	real,allocatable:: rikzero(:)
	real,allocatable::  rik(:),rik1(:),ri(:),rizero(:),pf1(:)
	real,allocatable:: temp(:),ro(:),wf(:),ir(:),tfg(:)
	real,allocatable:: zmax(:),slo(:),depth(:)
	real,allocatable:: zfmin(:),fsmin(:),pmin(:) 
end module grids