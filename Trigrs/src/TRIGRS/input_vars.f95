! By Rex L. Baum, 1 April 2004
module input_vars
	logical:: ans,outp(8),rodoc,lskip,lany,llus,lps0,unsat0
	logical,allocatable:: unsat(:) 
 	integer:: imax,row,col,nwf,tx,nmax
	integer:: flag,nper
  	integer:: nzs,mmax
	integer:: nzon,nout
	integer,allocatable:: ksav(:)
	real:: uww,zmin,t,dep,czmax,crizero,slomin
	real,allocatable:: ths(:),thr(:),alp(:),dif(:),c(:),phi(:)
	real,allocatable:: ks(:),uws(:),capt(:),cri(:),tsav(:)
	character (len=5):: flowdir
end module input_vars
