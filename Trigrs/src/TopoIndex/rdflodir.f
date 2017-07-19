c !  
c !  subroutine to determine flow direction from grid cells
c !  
c !  Rex L. Baum, USGS, 8 January 2000
c !  latest revsion 27 May 2008  
c !
	subroutine rdflodir(rc,nrow,ncol,dir,direc,u,nodat,mnd,u1)
	implicit none
	integer rc,u,u1
	integer j,nrow,ncol,m,mnd
	integer dir(rc)
	double precision param(6),celsiz,nodat,cns,cew
	double precision east,west,north,south
	character*14 header(6)
	character*255 direc
c !  
	open(u,file=direc,status='old',err=23)
	do 100, m=1,6
	read(u,*) header(m),param(m)
  100	continue
c ! set default value of nodat & celsiz for use with GRASS GIS ascii files   
  	nodat=-9999.d0
  	celsiz=-10.d0
	do 110, m=1,6
	if (trim(header(m)).eq.'ncols') ncol=int(param(m))
	if (trim(header(m)).eq.'nrows') nrow=int(param(m))
	if (trim(header(m)).eq.'cellsize') celsiz=param(m)
	if (trim(header(m)).eq.'cols:') ncol=int(param(m))
	if (trim(header(m)).eq.'rows:') nrow=int(param(m))
	if (trim(header(m)).eq.'east:') east=param(m)
	if (trim(header(m)).eq.'west:') west=param(m)
	if (trim(header(m)).eq.'north:') north=param(m)
	if (trim(header(m)).eq.'south:') south=param(m)
	if (trim(header(m)).eq.'NODATA_value') then
	  nodat=param(m)
	  mnd=m
	end if
	if (trim(header(m)).eq.'nodata_value') then
	  nodat=param(m)
	  mnd=m
	end if
  110	continue
  	if (celsiz.le.0) then
  	  cew=abs(east-west)/ncol
  	  cns=abs(north-south)/nrow
  	  if (cew.eq.cns) then
  	    celsiz=cew
  	  else
  	    celsiz=sqrt(cew*cns)
  	    write(*,*) 'Rectangular cells ',cew, ' X ', cns
  	    write(u1,*) 'Rectangular cells ',cew, ' X ', cns
  	  end if
  	end if
  	if (ncol*nrow .gt. rc) then
  	  write(*,*) 'Grid file exceeds array size'
	  write (*,*) '--> ',trim(direc)
  	  write(*,*) 'Check intialization file row and column values.'
  	  write(u1,*) 'Grid file exceeds array size'
	  write (u1,*) '--> ',trim(direc)
  	  write(u1,*) 'Check intialization file row and column values.'
	  pause 'Press RETURN to exit'
	  close(u)
	  close(u1)
	  stop
  	end if
	do 120, m=1,nrow
  	read(u,*,end=130) (dir(j+(m-1)*ncol), j=1,ncol)
  120	continue
  	close(u)
  	return
  130	continue
  	write (*,*) 'Error reading direction grid file.'
  	write (u1,*) 'Error reading direction grid file.'
	pause 'Press return/enter to quit program'
	stop
   23	continue
   	write (*,*) 'Error opening file < ',direc,' >'
	write (*,*) 'Check file name and location'
   	write (u1,*) 'Error opening file < ',direc,' >'
	write (u1,*) 'Check file name and location'
	pause 'Press RETURN to continue'
	stop
	end
