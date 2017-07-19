! By Rex L. Baum, 1 April 2004
module model_vars
	integer,parameter:: double=kind(1d0)
	integer:: nts,kper,nmax1,nmax2,nmn,nmin
	integer,allocatable:: jsav(:)
	real:: test1,dg2rad
	real,allocatable:: q(:),qb(:)
	real (double):: eps,tmin,tmax,ts,tinc,qt,tns,beta,qmax !
	real (double):: test,nodat,sumex,dusz,dcf,vf0 !(:)
	real (double):: celsiz,param(6),parami(6),ti,tis,pi,smt,lard
	real (double),allocatable:: p(:),ptran(:),pzero(:),bline(:)
	real (double),allocatable:: r(:),fc(:),fw(:),thz(:),kz(:),tcap(:)
	real (double),allocatable:: trz(:),uwsp(:),gs(:),qtime(:)
end module model_vars

