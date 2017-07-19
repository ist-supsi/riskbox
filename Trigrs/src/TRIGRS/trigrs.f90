!  Program to compute pore-pressure response and 
!  factor of safety for saturated and unsaturated infiltration.
!  by Rex L. Baum and W.Z. Savage, USGS
!
	program trigrs
	use input_file_defs; use input_vars
	use grids; use model_vars
	implicit none
 	integer:: grd
	integer:: i,j,k,imx1,m,mnd
  	integer:: nodata,sctr
	integer:: ncol,nrow,u(25),maxzo,ncc,nccs
	real::x1
	character (len=1):: tb
	character (len=255):: outfil 
  	character (len=14):: fminfil='TRfs_min_'
  	character (len=14):: zfminfil='TRz_at_fs_min_'
  	character (len=14):: pminfil='TRp_at_fs_min_'
  	character (len=18):: profil='TRlist_z_p_fs_'
	character (len=14):: header(6)
  	character (len=13):: ncvfil='TRnon_convrg_'  	
	character (len=8):: date
	character (len=10):: time
	character (len=4):: stp
	character (len=31):: scratch,irfil
	character (len=6):: vrsn
	character (len=11):: bldate
	character (len=224)::mypath
	logical :: lwarn

!------modified
	character argv*255
	call getarg( 1, argv )
!-------

! first executable statement ............	
     	call date_and_time(date,time)
	test=-9999.D0; test1=-9999.
	u=(/11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,&
     &	28,29,30,31,32,33,34,35/)
	pi=3.141592653589793
	dg2rad=pi/180.D0
	vrsn='2.0.00'; bldate='27 May 2008'
	write (*,*) ''
	write (*,*) 'TRIGRS: Transient Rainfall Infiltration'
	write (*,*) 'and Grid-based Regional Slope-Stability'
	write (*,*) '               Analysis'
	write (*,*) '       Version ', vrsn,', ',bldate
	write (*,*) '  By Rex L. Baum and William Z. Savage'
	write (*,*) '       U.S. Geological Survey'
	write (*,*) '-----------------------------------------'
	write (*,*) 'Portions of this program include material'
	write (*,*) '           copyrighted (C) by'
	write (*,*) '      Absoft Corporation 1988-2008.'
	write (*,*) ''
 	tb=char(9)
!  open log file
!----modified
	outfil = trim(argv)//'TrigrsLog.txt'
	open (u(19),file=trim(outfil),status='unknown',err=410)
!------------
	write (u(19),*) ''
     	write (u(19),*) 'Starting TRIGRS ', vrsn,' ',bldate
     	write (u(19),*) 'Date ',date(5:6),'/',date(7:8),'/',date(1:4)
	write (u(19),*) 'Time ',time(1:2),':',time(3:4),':',time(5:6)
!  read initialization file
	!call trini(u(19),u(7),dg2rad)
!---------modified
	call trini(u(19),u(7),dg2rad,argv)
!----------
! Allocate & initialize arrays needed for runoff routing
	grd=row*col
	imx1=imax
	allocate (pf2(grd),indx(imax),nxt(imax))
	allocate (dsctr(imax+1),dsc(nwf),slo(imax))
	allocate (pf1(grd),rizero(imax))
	allocate (ri(imax),rik(imax*nper),ro(imax))
	allocate (wf(nwf),rikzero(imax),temp(col))
	allocate (depth(imax),zmax(imax))
	allocate (zo(imax),ir(imax),tfg(imax))
	 pf2=0
	 indx=0
	 nxt=0
	 dsctr=0
	 dsc=0
	 zo=0
	 pf1=0.
	 slo=0.
	 rizero=0.
	 ri=0.
	 rik=0.
	 ro=0.
	 wf=0.
	 rikzero=0.
	 temp=0.
	 ir=0.
	 depth=0.
	 zmax=0.
! *****************************************************************
!  read gridded data from GIS
	write (*,*) 'Reading input grids'
     	write(u(19),*) 'Input file name,            Cell count'
!  read slope angles
	call srdgrd(grd,col,ncol,nrow,celsiz,nodat,&
     &	slo,pf1,sctr,imax,temp,u(1),slofil,param,header,u(19))
     	  write(u(19),*) 'Slope angle grid'
    	  write(u(19),*) trim(slofil),sctr,' data cells'
	  if (imax>sctr) imx1=sctr
	  write(u(19),*) 'Set number of cells for pore-pressure&
	  & factor of safety computations to ',imx1
	  if(imx1<imax) then
	    write(*,*) 'Check initialization file: '
	    write(*,*) 'Value of imax exceeds number of grid cells.'
	    write(u(19),*) 'Check initialization file: '
	    write(u(19),*) 'Value of imax > number of grid cells.'
	  end if
	  slo=slo*dg2rad ! convert slope angles to radians
!  read property zone numbers, zo
	if(nzon==1) then
	  zo=1 ! if only one zone, all values of zone grid equal 1.
	  write(*,*) 'One property zone, no grid required!'
	  write(u(19),*) 'One property zone, no grid required!'
	else
	  call irdgrd(grd,col,ncol,nrow,celsiz,nodata,mnd,&
          &zo,pf2,sctr,imax,temp,u(15),zonfil,parami,header,u(19))
    	  write(u(19),*) 'Property zone grid'
      	  write(u(19),*) trim(zonfil),sctr,' data cells'
     	  if(sctr/=imx1) then
     	    write (*,*) 'Grid mismatch ',trim(zonfil)
     	    write (*,*) 'Correct property-zone grid and/or initializtion file.' 
     	    write (u(19),*) 'Grid mismatch ',trim(zonfil)
     	    write (u(19),*) 'Correct property-zone grid and/or initializtion file.' 
     	    close(u(19))
     	    pause 'Press return/enter to quit'
     	    stop '-1'
     	  end if
     	  maxzo=maxval(zo)
     	  if (maxzo/=nzon) then
     	    write (*,*) 'Maximum zone number does not equal number of property zones!'
     	    write (*,*) 'Correct property-zone grid and/or initializtion file.' 
     	    write (u(19),*) 'Maximum zone number does not equal number of property zones!'
     	    write (u(19),*) 'Correct property-zone grid and/or initializtion file.' 
     	    close(u(19))
     	    pause 'Press return/enter to quit'
     	    stop '-1'
     	  end if
	end if
! *********************     	
!  read background infiltration rate, Isub0 
	if (crizero.lt.0) then 
	  call srdgrd(grd,col,ncol,nrow,celsiz,nodat,&
          &rizero,pf1,sctr,imax,temp,u(16),rizerofil,param,header,u(19))
    	  write(u(19),*) 'Background infiltration rate grid'
     	  write(u(19),*) trim(rizerofil),sctr,' data cells'
     	  if(sctr/=imx1) write (u(19),*) 'Grid mismatch ',trim(rizerofil)
  	else
	  rizero=crizero
  	end if
!  read initial depth to water table, 
	if (dep.lt.0) then  
	  call srdgrd(grd,col,ncol,nrow,celsiz,nodat,&
          &depth,pf1,sctr,imax,temp,u(10),depfil,param,header,u(19))
    	  write(u(19),*) 'Initial water-table depth grid'
     	  write(u(19),*) trim(depfil),sctr,' data cells'
     	  if(sctr/=imx1) write (u(19),*) 'Grid mismatch ',trim(depfil)
  	else
	  depth=dep
  	end if
!  read depth to base of potential slide, zmax
	if (czmax.lt.0) then 
	call srdgrd(grd,col,ncol,nrow,celsiz,nodat,&
     &	zmax,pf1,sctr,imax,temp,u(11),zfil,param,header,u(19))
    	write(u(19),*) 'Maximum depth grid'
     	write(u(19),*) trim(zfil),sctr,' data cells'
     	if(sctr/=imx1) write (u(19),*) 'Grid mismatch ',trim(zfil)
  	else
	zmax=czmax
  	end if
! *****************************************************************
     	write(u(19),*) '---------------******---------------'
! test and adjust (if necessary) steady background infiltration rates
	call steady(sumex,u(19),imx1)
! conduct runoff routing and adjust transient infiltration rates
     	call rnoff(grd,sumex,imx1,celsiz,param,parami,nodat,&
        & nodata,mnd,sctr,ncol,nrow,header,test1,u,25)
! Deallocate arrays that are no longer needed
	deallocate (ri,ro,wf)
	deallocate (pf2,indx,nxt,dsctr,dsc)
     	write(u(19),*) '---------------******---------------'
     	write(u(19),*) 'Input file name,          Cell count'
! *****************************************************************
! compute pore pressure distributions for either fully saturated or
! partially saturated conditions.
! Partially saturated zone overlies saturated zone
! Allocate and initialize new arrays
	  allocate (fsmin(imax*nout),pmin(imax*nout),zfmin(imax*nout))
	  allocate (p(nzs+1),ptran(nzs+1),pzero(nzs+1),bline(nzs+1))
	  allocate (fc(nzs+1),fw(nzs+1),thz(nzs+1),kz(nzs+1),trz(nzs+1))
	  allocate (nvu(imax),nv(imax),uwsp(nzs+1),gs(nzon))
	  fsmin=0.
	  zfmin=0.
	  pmin=0.
	  p=0.
	  ptran=0.
	  pzero=0.
	  bline=0.
	  fc=0.
	  fw=0.
	  nv=0
	  nvu=0
! determine number of time steps needed	  
	  kper=nper
	  if (t>capt(nper+1)) then
	    kper=nper+1 
	  else
    	    do k=1,nper ! find the period that contains t
    	      if(t>=capt(k) .and. t<=capt(k+1)) kper=k
            end do
	  end if 
	  if (tx<1) tx=1
	  nts=kper*tx ! number of time-steps from 0 to t
 	  tns=float(nts)
 	  tmin=0.
 	  tmax=t 
	  tinc=(tmax-tmin)/tns
! compute output pointers	
	  allocate(jsav(nts+1))  
	  jsav=0
	  write (u(19),*) '******** Output times ********'
	  write (u(19),*) 'number, timestep #,  time'
	  write (*,*) '******** Output times ********'
	  write (*,*) 'number, timestep #,  time'
	  lwarn=.false.
	  do k=1,nout
	    ts=tmin
	    do j=1,nts
	      if(tsav(k)>=ts .and. tsav(k)<ts+tinc) then
	        if(tsav(k)/=ts) lwarn=.true.
	        jsav(j)=k
	        ksav(k)=j
	        tsav(k)=ts
	        exit
	      else if(tsav(k)>=tmax) then
	        jsav(nts+1)=k
	        ksav(k)=nts+1
	        tsav(k)=tmax
	        exit
	      end if
	      ts=ts+tinc
	    end do
	    if(lwarn) then
	    write(u(19),*) 'One or more specified output times unavailable, '
	    write(u(19),*) 'Nearest available time substitued.'
	    write(*,*) 'One or more specified output times unavailable, '
	    write(*,*) 'Nearest available time substitued.'
	    end if 
	    write(u(19),*) k,ksav(k),tsav(k)
	    write(*,*) k,ksav(k),tsav(k)
	  end do
! allocate and initialize additonal model arrays	  
	  allocate (r(nmax),q(kper),qtime(2*nts+1),qb(nts),tcap(nts+2))
	  if(outp(7)) then
	    allocate(rik1(imax*nts))
	    rik1=0.
	  end if
	  eps=1.0e-18
	  nmax2=0
	  nmin=1+nmax
	  ncc=0;nccs=0
	  tis=tiny(x1)
	  write(*,*) 'Starting computations of pressure head &
          &and factor of safety'	
	if(unsat0) then
 	  smt=0.05d0; lard=12.d0 ! test values for early-time
	  do j=1,nzon 
	    gs(j)=((uws(j)/uww)-ths(j))/(1-ths(j))
	  end do
	  if(mmax.lt.0) then ! infinite depth model
	    mmax=20
	    write(*,*) 'Calling unsaturated infinite-depth model'
	    call unsinf(imx1,u(19),u(2),ncc,nccs,profil)
     	  else ! finite depth model
	    write(*,*) 'Calling unsaturated finite-depth model'
	    call unsfin(imx1,u(19),u(2),ncc,nccs,profil)
	  end if
	else ! Saturated zone extends to ground surface
	  write(*,*) 'Ignoring unsaturated zone'
	  if(tx==1 .and. nout==1) then
	    deallocate (rizero,ks,ir)
!  compute pore-pressure distributions and factor of safety
	    outfil=trim(folder)//trim(profil)//trim(suffix)//'.txt'
	    nccs=0
	    if(mmax.lt.0) then ! infinite depth model
	      mmax=20
	      write(*,*) 'Calling saturated infinite-depth model'
	      call iverson(imx1,u(2),outfil,u(19))
     	    else ! finite depth model
	      write(*,*) 'Calling saturated finite-depth model'
	      call savage(imx1,u(2),outfil,u(19),nccs)
	    end if
	  else
	    if(mmax.gt.0) then	
	      write(*,*) 'Calling multistep saturated finite-depth model'
	      call satfin(imx1,u(19),u(2),nccs,profil)
	    else
	      write(*,*) 'Calling multistep saturated infinite-depth model'
	      call satinf(imx1,u(19),u(2),nccs,profil)
	    end if
	  end if
	end if
! *****************************************************************
!  write output grid files 
  101	continue
	write(*,*) 'Saving results'
	ti=tiny(param(m))
	do j=1,nout
	  write(stp,'(i4)') j
	  stp=adjustl(stp)
	  if (outp(3)) then ! minimum factor of safety
	    tfg=0.
	    do i=1,imx1
	      tfg(i)=fsmin(i+(j-1)*imax)
	    end do
	    outfil=trim(folder)//trim(fminfil)//trim(suffix)//'_'//trim(stp)//'.asc'
   	    call ssvgrd(tfg,imax,pf1,nrow,ncol,u(4),test1,param,u(19),&
            &outfil,ti,header)
	  end if
	  if (outp(4)) then ! depth of minimum factor of safety
	    tfg=0.
	    do i=1,imx1
	      tfg(i)=zfmin(i+(j-1)*imax)
	    end do
	    outfil=trim(folder)//trim(zfminfil)//trim(suffix)//'_'//trim(stp)//'.asc'
   	    call ssvgrd(tfg,imax,pf1,nrow,ncol,u(5),test1,param,u(19),&
            &outfil,ti,header)
	  end if
 	  if (outp(5)) then ! pressure head at depth of minimum factor of safety
	    tfg=0.
	    do i=1,imx1
	      tfg(i)=pmin(i+(j-1)*imax)
	    end do
 	    outfil=trim(folder)//trim(pminfil)//trim(suffix)//'_'//trim(stp)//'.asc'
   	    call ssvgrd(tfg,imax,pf1,nrow,ncol,u(6),test1,param,u(19),&
            &outfil,ti,header)
	  end if
	end do
	if (outp(7) .and. unsat0) then ! incremental basal flux, unsaturated zone
	    do j=1,nts
	      ir=0.
	      do i=1,imx1
	        ir(i)=rik1(i+(j-1)*imax)
	      end do
	      ti=tiny(param(m))
	      write(scratch,'(i6)') j
	      scratch=adjustl(scratch)
	      irfil='TRunszfluxTS'//trim(scratch)//trim(suffix)//'.asc'
	      outfil=trim(folder)//trim(irfil)
   	      call ssvgrd(ir,imax,pf1,nrow,ncol,u(6),test1,&
     &	      param,u(19),outfil,ti,header)
	    end do
	  end if	  
 	if (ncc>0) then ! non-convergent cells, unsaturated zone
 	  outfil=trim(folder)//ncvfil//'UZ_'//trim(suffix)//'.asc'
   	  call isvgrd(nvu,imax,pf1,nrow,ncol,u(7),test,test,mnd,&
          & parami,u(19),outfil,ti,header)
	end if
 	if (nccs>0) then ! non-convergent cells, saturated zone
 	  outfil=trim(folder)//ncvfil//'SZ_'//trim(suffix)//'.asc'
   	  call isvgrd(nv,imax,pf1,nrow,ncol,u(7),test,test,mnd,&
          & parami,u(19),outfil,ti,header)
	end if
   	write (*,*) 'TRIGRS finished!'
     	write (u(19),*) 'TRIGRS finished normally'
     	call date_and_time(date,time)
     	write (u(19),*) 'Date ',date(5:6),'/',date(7:8),'/',date(1:4)
	write (u(19),*) 'Time ',time(1:2),':',time(3:4),':',time(5:6)
	close (u(19))
	stop '0'
! Error reporting 	
  410	continue
  	write (*,*) 'Error opening output file'
	write (*,*) '--> ',outfil
	write (*,*) 'Check file path and status'
  	write (u(19),*) 'Error opening output file'
	write (u(19),*) '--> ',outfil
	write (u(19),*) 'Check file path and status'
	pause 'Press RETURN to exit'
	stop '410'
   	end program trigrs

