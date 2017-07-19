! Program to (1) index an elevation grid (DEM) from highest point to lowest
! and to (2) compute weighting factors for dispersing flow among neighboring
! downslope cells in preparation for infilitration/runoff analysis.
! Indexing allows the cells to be analyzed from highest to lowest
! in a simple mass balance model for infiltration and runoff.
! indexing is checked and corrected against the ordering of subjacent
! cells along the steepest downslope path
!  
! by Rex L. Baum, USGS, 31 Jan 2002
! 
	program topoindex
 	implicit none
	integer:: rc,row,col,prm,aif
	integer,parameter:: double=kind(1d0)
	integer:: u(10),nrow,ncol,nodata,itmax,cctr,dsctr
	integer:: a,b,t1,t2,m,mnd,clcnt,imax,i,itctr,j,ia
	integer,allocatable:: cels(:),rndx(:),ordr(:)
	integer,allocatable:: cel(:),indx(:),cell(:)
	integer,allocatable:: list(:),lkup(:),dir(:)
	logical:: ans,op(5)
	real:: pwr
	real,allocatable:: z(:),z1(:),temp(:)
	real (double):: param(6),celsiz, nodat,nodat1,ti,aa
	character (len=255):: init,demfil,dirfil,heading
	character (len=18):: opthfil='TIdsneiList_' 
	character (len=18):: celfil='TIdscelGrid_'
 	character (len=18):: ndxfil='TIcelindxGrid_'
	character (len=18):: ndxlst='TIcelindxList_'
 	character (len=18):: drofil='TIflodirGrid_'
	character (len=18):: dscfil='TIdscelList_'
	character (len=18):: wffil='TIwfactorList_'
 	character (len=224):: folder
 	character (len=255):: outfil,outfil1,title
	character (len=31):: scratch
	character (len=14):: header(6),vdate
	character (len=10):: tm
	character (len=8):: dt,suffix,vrsn
	character (len=1):: sp
!------modified
	character argv*255
	call getarg( 1, argv )
!-------
	data u/11,12,13,14,15,16,17,18,19,20/ 	
! First executable statement.....	
     	call date_and_time(dt,tm)
! current version number and date     	
     	vdate='27 May 2008'; vrsn='1.0.03'
	sp=char(32)
!------modified
	outfil=trim(argv)//'TopoIndexLog.txt'
!------
	open(u(9),file=trim(outfil),err=25)
     	write (u(9),*) ''
	write (u(9),*) 'Starting TopoIndex'
	write (u(9),*) 'Version: ',trim(vrsn),', ',vdate
     	write (u(9),*) 'Date: ',dt(5:6),'/',dt(7:8),'/',dt(1:4)
	write (u(9),*) 'Time: ',tm(1:2),':',tm(3:4),':',tm(5:6)
	write (*,*) '  TopoIndex: Topographic Indexing and'
	write (*,*) ' flow distribution factors for routing'
	write (*,*) ' runoff through Digital Elevation Models'
	write (*,*) '            By Rex L. Baum'
	write (*,*) '       U.S. Geological Survey'
	write (*,*) '     Wersion ',trim(vrsn),', ',trim(vdate)
	write (*,*) '-----------------------------------------'
	write (*,*) 'Portions of this program include material'
	write (*,*) '           copyrighted (C) by'
	write (*,*) '      Absoft Corporation 1988-2008.'
	write (*,*) ''
!  open and read initialization file
!------modified
	init=trim(argv)//'tpx_in.txt' 
!------
	inquire (file=trim(init),exist=ans)
	if(ans)then
 	  open(u(1),file=trim(init),status='old',err=21)
	else
	  write (*,*) 'Cannot locate default initialization file'
	  write (*,*) 'Type name of initialization file and'
	  write (*,*) 'press RETURN to continue'
	  read (*,'(a)') init
	  open (u(1),file=trim(init),status='old',err=21)
	end if
	write(u(9),*) ''
	write(u(9),*) '-- LISTING OF INITIALIZATION DATA --'
	read(u(1),'(a)',err=30) heading
	read(u(1),'(a)',err=30) title
	write(u(9),*) heading
	write(u(9),*) title
	read(u(1),'(a)',err=30) heading
	read(u(1),*,err=30) row,col,aif
	write(u(9),*) heading
	write(u(9),*) row,col,aif
	read(u(1),'(a)',err=30) heading
	read(u(1),*,err=30) pwr,itmax
	write(u(9),*) heading
	write(u(9),*) pwr,itmax
	read(u(1),'(a)',err=30) heading
	read(u(1),'(a)',err=30) demfil
	write(u(9),*) heading
	write(u(9),*) demfil
	read(u(1),'(a)',err=30) heading
	read(u(1),'(a)',err=30) dirfil
	write(u(9),*) heading
	write(u(9),*) dirfil
	read(u(1),'(a)',err=30) heading
	read(u(1),*,err=30) op(1)
	write(u(9),*) heading
	write(u(9),*) op(1)
	read(u(1),'(a)',err=30) heading
	read(u(1),*,err=30) op(2)
	write(u(9),*) heading
	write(u(9),*) op(2)
	read(u(1),'(a)',err=30) heading
	read(u(1),*,err=30) op(3)
	write(u(9),*) heading
	write(u(9),*) op(3)
	read(u(1),'(a)',err=30) heading
	read(u(1),*,err=30) op(4)
	write(u(9),*) heading
	write(u(9),*) op(4)
	read(u(1),'(a)',err=30) heading
	read(u(1),*,err=30) op(5)
	write(u(9),*) heading
	write(u(9),*) op(5)
	read(u(1),'(a)',err=30) heading
	read(u(1),'(a)',err=30) folder
	write(u(9),*) heading
	write(u(9),*) folder
	read(u(1),'(a)',err=30) heading
	read(u(1),'(a)',err=30) suffix
	write(u(9),*) heading
	write(u(9),*) suffix
	write(u(9),*) '-- END OF INITIALIZATION DATA --'
	write(u(9),*) ''
	write(u(9),*) title
	write(u(9),*) ''
	close(u(1))
	write(*,*) title
! allocate arrays
	rc=row*col
	prm=2*col
	allocate (cels(rc),rndx(rc),ordr(rc))
	allocate (cel(rc),indx(rc),cell(rc))
	allocate (list(prm),lkup(rc),dir(rc))
	allocate (z(rc),z1(rc),temp(col))
! initialize arrays and counter
	 cels=0
	 rndx=0
	 ordr=0
	 cel=0
	 indx=0
	 cell=0
	 list=0
	 lkup=0
	 dir=0
	 z=0.
	 z1=0.
	 temp=0.
	clcnt=0
	imax=rc
! open and read elevation grid file
	write(*,*) 'Reading elevation grid data'
	call srdgrd1(rc,col,ncol,nrow,celsiz,nodat,&
     &	z,z1,clcnt,imax,temp,u(2),demfil,param,header,cel,u(9))
	write (u(9),*) 'Parameters for file--> ', demfil
	write (u(9),*) 'Data cells, Rows, Columns'
	write (u(9),*) clcnt,nrow,ncol
! index elevation grid from highest point to lowest
	if (clcnt==1) indx(1)=1
 	if (clcnt>1) then
 	  call sindex(clcnt,z,indx)
	  do i=1,clcnt
	   t1=indx(i)
	   lkup(t1)=i
  	  end do
	end if
	write(*,*) 'Initial elevation indexing completed'
! Check indexing against list of subjacent cells, 
!  Each cell should have a lower "order" or rank than its subjacent cell
!  so that it will be processed sooner in the infiltration/ruonoff routing model.
!  Find subjacent cell from a grid of flow directions at each cell
!    
! Import flow direction grid 
	  inquire (file=trim(dirfil),exist=ans)
	  if (ans) then
     	   write (u(9),*) 'Reading flow-direction data'
  	   call rdflodir(rc,nrow,ncol,dir,dirfil,u(7),nodat1,mnd,u(9))
	   nodata=int(nodat1) 
	   if (aif==1) then
! transform flow direction grid from ARC/GRID numbering to TRIGRS numbering scheme
	     write (*,*) 'Converting directional data'
	     call mpfldr(rc,dir,nodata)
	     if (op(5)) then
	     	outfil=trim(folder)//trim(drofil)//trim(suffix)//'.asc'
	      open (u(4),file=trim(outfil),status='unknown',err=25)
	      param(mnd)=nodat1
	      do m=1,6
	       aa=param(m)
	       ia=int(param(m))
	       if(abs(aa-ia)<=tiny(aa)) then
 	         write(scratch,'(i16)') ia
	       else
 	         write(scratch,*) aa
	       end if
   	       scratch=adjustl(scratch)
	       write(u(4),'(t1,a,t15,a)') trim(header(m)),trim(scratch)
  	      end do
	      do i=1,nrow
	        do j=1,ncol
	          write(scratch,'(i14)') dir(j+(i-1)*ncol)
		  scratch=adjustl(scratch)
		  if (j/=ncol) then
	            write(u(4),'(a,a1,$)') trim(scratch),sp
		  else
	            write(u(4),'(a)') trim(scratch)
		  end if
	        end do
	      end do
	      close(u(4))
	     end if
	   end if
	  else
! directional data not available, so save cell count and quit 	  
     	   write (u(9),*) ''
	   write (u(9),*) 'Parameters for file--> ', demfil
	   write (u(9),*) 'Data cells, Rows, Columns'
	   write (u(9),*) clcnt,nrow,ncol
     	   write (u(9),*) ''
     	   write (u(9),*) 'Flow-direction data not available'
     	   write (*,*) 'Flow-direction data not available'
     	   call date_and_time(dt,tm)
     	   write (u(9),*) 'Date: ',dt(5:6),'/',dt(7:8),'/',dt(1:4)
	   write (u(9),*) 'Time: ',tm(1:2),':',tm(3:4),':',tm(5:6)
	   close (u(9))
	   stop
	  end if
!  Analyze flow direction grid to find downslope neighbor in D8 flow direction
     	   write (*,*) 'Finding D8 neighbor cells'
	  if(op(1)) then
	    outfil=trim(folder)//trim(opthfil)//trim(suffix)//'.txt'
 	    open (u(8),file=trim(outfil),err=25)
 	    call nxtcel(nrow,ncol,rc,prm,u(8),nodat,nodata,&
     &	     z1,dir,cell,list,cels,u(9))
     	    close(u(8))
	  else
	    u(8)=-1
 	    call nxtcel(nrow,ncol,rc,prm,u(8),nodat,nodata,&
     &	     z1,dir,cell,list,cels,u(9))
	   end if
! Correct initial ordering using D8 neighbor cells
     	write (*,*) 'Correcting cell index numbers'
	do itctr=1,itmax
	  write (u(9),*) 'iteration',itctr
  	  cctr=0
	  do i=1,imax
     	   a=lkup(indx(i))	
	   b=lkup(cels(indx(i)))
	   if (b>a) then
	    cctr=cctr+1
	    t1=indx(a)
	    t2=indx(b)
	    lkup(cels(indx(i)))=a
     	    lkup(indx(i))=b
	    indx(a)=t2
	    indx(b)=t1
	   end if
	   rndx(i)=indx(clcnt+1-i)
	   a=rndx(i)
	   ordr(a)=i
  	  end do
  	  write (u(9),*) 'iteration ',itctr, ' corrections ',cctr
	  if (cctr==0) exit
	end do
  	if (cctr>0) then
	  write (u(9),*) 'Corrections did not converge after',&
     &    itmax,' iterations'
          close (u(9))
     	  stop
	end if
! compute weighting factors
     	write (*,*) 'Computing weighting factors'
	outfil=trim(folder)//trim(dscfil)//trim(suffix)//'.txt'
	outfil1=trim(folder)//trim(wffil)//trim(suffix)//'.txt'
 	call slofac(nrow,ncol,rc,z1,celsiz,nodat,nodata,cel,&
     &	u(3),outfil,pwr,cels,ordr,u(10),outfil1,dsctr,dir)
! write grid of D8 downslope neighbor cell numbers
!  Shows where runoff goes from each cell
!  if it follows the steepest downslope path
     	write (*,*) 'Saving results to disk'
	ti=tiny(param(1))
	if (op(2)) then
	outfil=trim(folder)//trim(celfil)//trim(suffix)//'.asc'
 	call isvgrd(cels,rc,z1,nrow,ncol,u(4),nodat,nodat1,mnd,param,&
     &	u(9),outfil,ti,header) ! nodat1 is the nodata value for integer grids
	end if
! write grid of computation order
	if (op(3)) then
	outfil=trim(folder)//trim(ndxfil)//trim(suffix)//'.asc'
   	call isvgrd(ordr,rc,z1,nrow,ncol,u(5),nodat,nodat1,mnd,param,&
     &	u(9),outfil,ti,header)
	end if
! write list of cell number and cell index
	if (op(4)) then
	outfil=trim(folder)//trim(ndxlst)//trim(suffix)//'.txt'
	open(u(6),file=trim(outfil),status='unknown',err=25)
	do i=1,clcnt
	 write(u(6),*) i,rndx(i)
  	end do
	close(u(6))
	end if
	write (u(9),*) 'Parameters for file--> ', demfil
	write (u(9),*) 'Exponent ',pwr
	write (u(9),*) 'Data cells, Rows, Columns, Downslope cells'
	write (u(9),*) clcnt,nrow,ncol,dsctr
     	write (u(9),*) ''
     	write (u(9),*) 'TopoIndex finished normally'
     	call date_and_time(dt,tm)
     	write (u(9),*) 'Date: ',dt(5:6),'/',dt(7:8),'/',dt(1:4)
	write (u(9),*) 'Time: ',tm(1:2),':',tm(3:4),':',tm(5:6)
	close(u(9)) 
	stop
! Report errors when opening files	
   21	continue
	write (*,*) '*** Error opening intialization file ***'
	write (*,*) '--> ',init
	write (*,*) 'Check file name and location'
	pause 'Press RETURN to exit'
	stop
   25	continue
  	write (*,*) 'Error opening output file'
	write (*,*) '--> ',outfil
	write (*,*) 'Check file path and status'
	pause 'Press RETURN to exit'
	stop
   30	continue
	write (*,*) '*** Error reading file ***'
	write (*,*) '--> ',init
	write (*,*) 'Check file format and data'
	pause 'Press RETURN to exit'
	stop
	end program topoindex
