	subroutine trini(ulog,uini,dg2rad,myinit)
! reads initialization file, R.L. Baum, USGS	
	use input_file_defs; use input_vars
	implicit none
	integer:: i,j,iz,linct
	integer, intent(in):: ulog,uini
	real, intent(in):: dg2rad
	character (len=31):: scratch
!-----modified
	character (len=255) :: myinit
	init=trim(myinit)//'tr_in.txt'
	open (uini,file=trim(init),status='old',err=201)
	write (*,*) 'Opening default initialization file'
!-----------
	!init='tr_in.txt'
	!inquire (file=trim(init),exist=ans)
	!if(ans) then
	!  open (uini,file=trim(init),status='old',err=201)
	!  write (*,*) 'Opening default initialization file'
	!else
	!  write (*,*) 'Cannot locate default initialization file, <tr_in.txt>'
	  !-------------------
	  !call getarg( 1, argv )
	  !open (uini,file=trim(argv),status='old',err=201)
	  !write (*,*) 'letto'
	  !-------------------
	  !write (*,*) 'Type name of initialization file and'
	  !write (*,*) 'press RETURN to continue'
	  !read (*,'(a)') init
	  !open (uini,file=trim(init),status='old',err=201)
	!end if
	write (ulog,*) 'initialization file -->',init
	write (ulog,*) '-- LISTING OF INITIALIZATION FILE --'	
! write copy of data to log file
	linct=1	
  	read (uini,'(a)',err=420) heading; linct=linct+1
	read (uini,'(a)',err=420) title; linct=linct+1
  	write (ulog,*) heading
  	write (ulog,*) title
  	write (*,*) title
  	read (uini,'(a)',err=420) heading; linct=linct+1
  	read (uini,*,err=420) imax,row,col,nwf,tx,nmax; linct=linct+1
  	write (ulog,*) heading
  	if(nmax<2) nmax=2 ! set minimum value for nmax
  	write (ulog,*) imax,row,col,nwf,tx,nmax
  	read (uini,'(a)',err=420) heading; linct=linct+1
	read (uini,*,err=420) nzs,mmax,nper,zmin,uww,t,nzon
	linct=linct+1
  	write (ulog,*) heading
  	write (ulog,*) nzs,mmax,nper,zmin,uww,t,nzon
  	read (uini,'(a)',err=420) heading; linct=linct+1
        read (uini,*,err=420) czmax,dep,crizero,slomin; linct=linct+1
  	write (ulog,*) heading
  	write (ulog,*) czmax,dep,crizero,slomin
  	slomin=slomin*dg2rad ! convert minimum slope angle to radians
! allocate & read arrays for zone properties and initial conditions  	
  	allocate (c(nzon),phi(nzon),uws(nzon),dif(nzon),&
	  & ks(nzon),ths(nzon),thr(nzon),alp(nzon),unsat(nzon))
	c=0;phi=0
	uws=0;dif=0;ks=0;ths=0;thr=0;alp=0
	unsat=.true.;unsat0=.false.
  	do i=1,nzon
  	  read (uini,*,err=420) scratch,iz ! property zone number
  	  linct=linct+1
  	  read (uini,'(a)',err=420) heading; linct=linct+1
	  read (uini,*,err=420) c(iz),phi(iz),&
	  & uws(iz),dif(iz),ks(iz),ths(iz),thr(iz),alp(iz) ! zone parameters
	  linct=linct+1
  	  write (ulog,*) trim(scratch),': ',iz
  	  write (ulog,*) heading
  	  write (ulog,*) c(iz),phi(iz),&
	  & uws(iz),dif(iz),ks(iz),ths(iz),thr(iz),alp(iz) 
	  if(ths(iz)<thr(iz)) then 
	    write(*,*)'Error, Theta-resid. > Theta-sat. for property zone',iz
	    write(*,*)'Saturated infiltration model will be used for cells in zone', iz,'.'
	    write(ulog,*)'Error, Theta-resid. > Theta-sat. for property zone',iz
	    write(ulog,*)'Saturated infiltration model will be used for cells in zone', iz,'.'
	    unsat(iz)=.false.
	  end if
	  if(alp(iz)<=0) then
	    unsat(iz)=.false.
	    write(*,*)'Negative or zero value of Alpha for property zone',iz
	    write(*,*)'Saturated infiltration model will be used for cells in zone', iz,'.'
	    write(ulog,*)'Negative or zero value of Alpha for property zone',iz
	    write(ulog,*)'Saturated infiltration model will be used for cells in zone', iz,'.'
	  end if
	  if(unsat(iz)) then
	    unsat0=.true. ! tracks whether any property zones are unsat.
	    write(*,*)'Unsaturated infiltration model will be used for cells in zone', iz,'.'
	    write(ulog,*)'Unsaturated infiltration model will be used for cells in zone', iz,'.'
	  end if
  	end do
  	phi=phi*dg2rad ! convert phi angles to radians
! Allocate & read arrays for storm period data	
	allocate (cri(nper),capt(nper+2),rifil(nper)) 
  	read (uini,'(a)',err=420) heading; linct=linct+1
	read (uini,*,err=420) (cri(j), j=1,nper) ! List of rainfall rates
	linct=linct+1
  	write (ulog,*) heading
  	write (ulog,*) (cri(j), j=1,nper)
  	read (uini,'(a)',err=420) heading; linct=linct+1
	read (uini,*,err=420) (capt(j), j=1,nper+1) ! List of times corresponding to change of rate
	capt(nper+2)=t ! for cases where t>capt(nper+1)
	linct=linct+1
  	write (ulog,*) heading
  	write (ulog,*) (capt(j), j=1,nper+1)
!  path names of input files
  	read (uini,'(a)',err=420) heading; linct=linct+1
	read (uini,'(a)',err=420) slofil ! File name of slope angle grid (slofil)
	linct=linct+1
  	write (ulog,*) heading
  	write (ulog,*) slofil
  	read (uini,'(a)',err=420) heading; linct=linct+1
	read (uini,'(a)',err=420) zonfil ! File name of property zone grid (zonfil)
	linct=linct+1
  	write (ulog,*) heading
  	write (ulog,*) zonfil
  	read (uini,'(a)',err=420) heading; linct=linct+1
	read (uini,'(a)',err=420) zfil ! File name of depth grid (zfil)
	linct=linct+1
  	write (ulog,*) heading
  	write (ulog,*) zfil
  	read (uini,'(a)',err=420) heading; linct=linct+1
	read (uini,'(a)',err=420) depfil ! File name of initial depth of water table grid (depfil)
	linct=linct+1
  	write (ulog,*) heading
  	write (ulog,*) depfil
  	read (uini,'(a)',err=420) heading; linct=linct+1
	read (uini,'(a)',err=420) rizerofil ! File name of initial infiltration rate grid (rizerofil)
	linct=linct+1
  	write (ulog,*) heading
  	write (ulog,*) rizerofil
  	read (uini,'(a)',err=420) heading; linct=linct+1
  	write (ulog,*) heading
	do j=1,nper
	  read (uini,'(a)',err=421) rifil(j) ! List of file names of rainfall intensity for each period, (rifil())
	  linct=linct+1	
  	  write (ulog,*) rifil(j)
  	end do
  	read (uini,'(a)',err=421) heading; linct=linct+1
	read (uini,'(a)',err=421) nxtfil ! File name of grid of D8 runoff receptor cell numbers (nxtfil)
	  linct=linct+1	
  	write (ulog,*) heading
  	write (ulog,*) nxtfil
  	read (uini,'(a)',err=421) heading; linct=linct+1
	read (uini,'(a)',err=421) ndxfil ! File name of list of defining runoff computation order (ndxfil)
	  linct=linct+1	
  	write (ulog,*) heading
  	write (ulog,*) ndxfil
  	read (uini,'(a)',err=421) heading; linct=linct+1
	  linct=linct+1	
	read (uini,'(a)',err=421) dscfil ! File name of list of all runoff receptor cells  (dscfil)
  	write (ulog,*) heading
  	write (ulog,*) dscfil
  	read (uini,'(a)',err=421) heading; linct=linct+1
	read (uini,'(a)',err=421) wffil ! File name of list of runoff weighting factors  (wffil)
	  linct=linct+1	
  	write (ulog,*) heading
  	write (ulog,*) wffil
!  location of output files	
  	read (uini,'(a)',err=421) heading; linct=linct+1
	read (uini,'(a)',err=421) folder ! Folder where output grid files will be stored  (folder)
	  linct=linct+1	
  	write (ulog,*) heading
  	write (ulog,*) folder
!  output-file ID code
  	read (uini,'(a)',err=421) heading; linct=linct+1
	read (uini,'(a)',err=421) suffix ! Identification code to be added to names of output files (suffix)
	  linct=linct+1	
  	write (ulog,*) heading
  	write (ulog,*) suffix
!  output file selections	
  	read (uini,'(a)',err=421) heading; linct=linct+1
	read (uini,*,err=421) rodoc ! Save grid files of runoff?
	  linct=linct+1	
  	write (ulog,*) heading
  	write (ulog,*) rodoc
  	read (uini,'(a)',err=421) heading; linct=linct+1
	read (uini,*,err=421) outp(3) ! Save grid of minimum factor of safety?
  	write (ulog,*) heading
  	write (ulog,*) outp(3)
  	read (uini,'(a)',err=421) heading; linct=linct+1
	read (uini,*,err=421) outp(4) ! Save grid of depth of minimum factor of safety?
	  linct=linct+1	
  	write (ulog,*) heading
  	write (ulog,*) outp(4)
  	read (uini,'(a)',err=421) heading; linct=linct+1
	read (uini,*,err=421) outp(5) ! Save grid of pore pressure at depth of minimum factor of safety?
	  linct=linct+1	
  	write (ulog,*) heading
  	write (ulog,*) outp(5)
  	read (uini,'(a)',err=421) heading; linct=linct+1
	read (uini,*,err=421) outp(6) ! Save grid files of actual infiltration rate?
	  linct=linct+1	
  	write (ulog,*) heading
  	write (ulog,*) outp(6)
  	read (uini,'(a)',err=421) heading; linct=linct+1
	read (uini,*,err=421) outp(7) ! Save grid files of unsaturated zone basal flux?
	  linct=linct+1	
  	write (ulog,*) heading
  	write (ulog,*) outp(7)
  	read (uini,'(a)',err=421) heading; linct=linct+1
  	read (uini,*,err=421) flag ! Save listing of pressure head and factor of safety ("flag")?
	  linct=linct+1	
  	write (ulog,*) heading
  	write (ulog,*) flag
  	read (uini,'(a)',err=421) heading; linct=linct+1
  	read (uini,*,err=421) nout ! Number of times to save output grids
	  linct=linct+1	
  	if(nout<1) nout=1 ! must save at least one time; negative number not allowed
  	write (ulog,*) heading
  	write (ulog,*) nout 
  	allocate (tsav(nout),ksav(nout))
  	tsav=0.;ksav=0
  	read (uini,'(a)',err=420) heading; linct=linct+1 
  	read (uini,*,err=420) (tsav(j), j=1,nout) ! Times of output grids
	  linct=linct+1	
  	write (ulog,*) heading
  	write (ulog,*) (tsav(j), j=1,nout)
! user options  	
  	read (uini,'(a)',err=420) heading; linct=linct+1
	read (uini,*,err=420) lskip ! Skip other timesteps?
	  linct=linct+1	
  	write (ulog,*) heading
  	write (ulog,*) lskip
  	read (uini,'(a)',err=420) heading; linct=linct+1
	read (uini,*,err=420) lany ! Use analytic solution for fillable porosity?
	  linct=linct+1	
  	write (ulog,*) heading
  	write (ulog,*) lany
  	read (uini,'(a)',err=420) heading; linct=linct+1
	read (uini,*,err=420) llus ! Estimate positive pressure head in rising water table zone
	  linct=linct+1	
  	write (ulog,*) heading
  	write (ulog,*) llus
  	read (uini,'(a)',err=420) heading; linct=linct+1
	read (uini,*,err=420) lps0 ! Use psi0=-1/alpha?
	  linct=linct+1	
  	write (ulog,*) heading
  	write (ulog,*) lps0
  	read (uini,'(a)',err=420) heading; linct=linct+1
	read (uini,*,err=420) outp(8) ! Log mass balance results?
	  linct=linct+1	
  	write (ulog,*) heading
  	write (ulog,*) outp(8)
 	read (uini,'(a)',err=420) heading; linct=linct+1
	read (uini,*,err=420) flowdir ! Specify flow direction
  	write (ulog,*) heading
  	write (ulog,*) flowdir
	close(uini)
	write (ulog,*) '-- END OF INITIALIZATION DATA --'	
	write (ulog,*) ''
	write (ulog,*) title
	write (ulog,*) ''
  	do iz=1,nzon
	  if(alp(iz)>=0 .and. unsat(iz)) then
	    write(*,*) '******** Zone ',iz, ' *********'
	    write(*,*)'Using unsaturated infiltration model.'
	    write(*,*) 'Cells where water table is shallower than '
	    write(*,*) '           ', 1./alp(iz)
	    write(*,*) 'treated as tension saturated--Saturated infiltration model used.'
	    write(ulog,*) '******** Zone ',iz, ' *********'
	    write(ulog,*)'Using unsaturated infiltration model.'
	    write(ulog,*) 'Cells where water table is shallower than '
	    write(ulog,*) '           ', 1./alp(iz)
	    write(ulog,*) 'treated as tension saturated--Saturated infiltration model used.'
	  end if
	end do
	return
  201	continue
	write (*,*) '*** Error opening intialization file ***'
	write (*,*) '--> ',trim(init)
	write (*,*) 'Check file location and name'
     	write (ulog,*) '*** Error opening intialization file ***'
     	write (ulog,*) '--> ',trim(init)
     	write (ulog,*) 'Check file location and name'
	pause 'Press RETURN to exit'
	stop '201'
  420	continue
  	write (*,*) 'Error reading initialization file'
	write (*,*) '--> ',trim(init), 'at line ',linct
	write (*,*) 'Check file contents and organization'
  	write (ulog,*) 'Error reading initialization file'
	write (ulog,*) '--> ',trim(init), 'at line ',linct
	write (ulog,*) 'Check file contents and organization'
	pause 'Press RETURN to exit'
	stop '420'
  421	continue
  	write (*,*) 'Error reading initialization file'
	write (*,*) '--> ',trim(init), 'at line ',linct
	write (*,*) 'Check file contents and organization'
	write (*,*) 'Number of file names/place holders for rainfall data'
	write (*,*) 'must equal nper.  List each on a separate line.'
  	write (ulog,*) 'Error reading initialization file'
	write (ulog,*) '--> ',trim(init), 'at line ',linct
	write (ulog,*) 'Check file contents and organization'
	write (ulog,*) 'Number of file names/place holders for rainfall data'
	write (ulog,*) 'must equal nper.  List each on a separate line.'
	pause 'Press RETURN to exit'
	stop '421'
	end subroutine trini
