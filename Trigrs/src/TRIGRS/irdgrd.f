c !  subroutine to read an ascii grid file of integer values
c !  and store in a 1-d array
c ! Rex L. Baum, USGS, spring 2001
c !  Dec 2003, added code to handle GRASS GIS ASCII grid files
c ! 
	subroutine irdgrd(grd,pth,ncol,nrow,celsiz,inodat,mnd,
     +	y,y1,ctr,imax,itemp,u,infil,param,header,u1)
	implicit none
	integer grd,pth,i,m,ncol,nrow,ctr,imax,u,inodat,u1
	integer y(imax),y1(grd),itemp(pth),mnd
	double precision param(6),celsiz,cns,cew,nodat
	double precision east,west,north,south
	character*14 header(6)
	character*255 infil
c !  	
	open(u,file=infil,status='old',err=10)
	do 200, m=1,6
	read(u,*) header(m),param(m)
  200	continue
c ! set default value of nodat & celsiz for use with GRASS GIS ascii files   
  	inodat=-9999
  	celsiz=-10.
	do 210, m=1,6
	if (trim(header(m)).eq.'ncols') ncol=int(param(m))
	if (trim(header(m)).eq.'nrows') nrow=int(param(m))
	if (trim(header(m)).eq.'cellsize') celsiz=param(m)
	if (trim(header(m)).eq.'cols:') ncol=int(param(m))
	if (trim(header(m)).eq.'rows:') nrow=int(param(m))
	if (trim(header(m)).eq.'east:') east=param(m)
	if (trim(header(m)).eq.'west:') west=param(m)
	if (trim(header(m)).eq.'north:') north=param(m)
	if (trim(header(m)).eq.'south:') south=param(m)
c ! 7/19/06 added pointer for nodata postion in param array	
	if (trim(header(m)).eq.'NODATA_value') then
	  nodat=int(param(m))
	  mnd=m
	end if
	if (trim(header(m)).eq.'nodata_value') then
	  nodat=int(param(m))
	  mnd=m 
	end if
	inodat=int(nodat)
  210	continue
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
  	if (ncol*nrow .gt. grd) then
  	  write(*,*) 'Grid file exceeds array size'
	  write (*,*) '--> ',trim(infil)
  	  write(*,*) 'Check intialization file row and column values.'
  	  write(u1,*) 'Grid file exceeds array size'
	  write (u1,*) '--> ',trim(infil)
  	  write(u1,*) 'Check intialization file row and column values.'
	  pause 'Press RETURN to exit'
	  close(u)
	  close(u1)
	  stop '-11'
  	end if
	ctr=0
	do 120, m=1,nrow
c !  next sequence of lines read data in but skips no_data values
c !  count maintained by ctr should coincide with node numbers from GIS
c !  y1() keeps track of positions of nodata values so that results
c !  can be written out in grid format.
  	read(u,*,end=125) (itemp(i), i=1,ncol)
	do 250, i=1,ncol
	  y1(i+(m-1)*ncol)=itemp(i)
	  if(itemp(i).ne.inodat) then
	  ctr=ctr+1
	  if (ctr>imax) then
  	    write(*,*) 'Number of data cells exceeds array size'
	    write (*,*) '--> ',trim(infil)
  	    write(*,*) 'Check imax value in intialization file.'
  	    write(u1,*) 'Number of data cells exceeds array size'
	    write (u1,*) '--> ',trim(infil)
  	    write(u1,*) 'Check imax value in intialization file.'
	    pause 'Press RETURN to exit'
  	    close(u)
	    close(u1)
	    stop '-12'
	  end if
	  y(ctr)=itemp(i)
	  end if
  250	continue
  120	continue
  125	close(u)
	return
   10	continue
   	write (*,*) '*** Error opening input file ***'
	write (*,*) '--> ',infil
	write (*,*) 'Check file name and location'
   	write (u1,*) '*** Error opening input file ***'
	write (u1,*) '--> ',infil
	write (u1,*) 'Check file name and location'
	pause 'Press RETURN to exit'
	close(u)
	close(u1)
	stop '-10'
	end