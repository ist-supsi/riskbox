!  Implementation of Iverson's (2000) method of computing pore pressure
!  for rain infiltration.
!  by W.Z. Savage, spring 2001, with modifications by R.L. Baum
	subroutine iverson(imx1,u1,profil,ulog)
	use grids; use input_vars
	use model_vars
 	implicit none
	integer:: j,i,u1,ulog,imx1,n 
	real:: finf 
	real (double) :: derfc,a1,b1,ff,zns,zinc,z
	real (double) :: zstar,tstar,x1,x2,x3,x4
	real (double) :: rf1,rf2,rf3,rf4,rfa,rfb,rf
	real (double) :: fs,rslo,rphi,fmn,ptest
	real (double) :: newdep,captstar1,captstar2,tdif1,tdif2
	character (len=255) profil
	write(ulog,*) 'Starting saturated-zone'
	write(ulog,*) 'computations for infinite-depth'
	write(*,*) 'Starting saturated-zone'
	write(*,*) 'computations for infinite-depth'
	pi=3.141592653589793
	dg2rad=pi/180.D0
!  maximum value of Factor of Safety	
	finf=10.
!  open file for storing complete data on each cell
	if (flag==-1 .or. flag==-2) then
	  open (u1,file=profil,err=70)
	  write(u1,'(a)') 'TRIGRS depth profiles at each cell'
	  write(u1,'(a)') 'Infinite depth no-flow boundary'
	  write(u1,*) 'Time:',t
	  write(u1,'(a)') 'Cell Number, Slope angle'
	end if
  	if (flag==-1) write(u1,'(a)') 'Z         P         FS'
  	if (flag==-2) write(u1,'(a)') &
         & 'Z         P      Pzero     Ptran      Pbeta       FS'
!  loop steps through all grid cells
	write(*,*) 'Cells completed: '
	grid_loop: do i=1,imx1
	  rslo=slo(i) 
	  if(rslo<slomin) then
	    fsmin(i)=finf+1.
	    zfmin(i)=zmax(i)
	    pmin(i)=0.
	    if (mod(i,2000)==0) write (*,*) i
	    cycle grid_loop
	  end if
	  if (flag==-1 .or. flag==-2) write (u1,'(i12,f6.1)')i,rslo/dg2rad
	  rphi=phi(zo(i))
	  a1=sin(rslo)
	  b1=cos(rslo)
	  select case (flowdir) ! set value of beta (Iverson's beta line)
	    case ('slope')
	    beta=b1*b1
	    case ('hydro')
	    beta=1.d0
	    case default
	    beta=b1*(b1-rikzero(i))
	  end select
	  if(abs(b1-rikzero(i))<1.e-6) beta=0.d0
	  if (abs(rslo)>1.e-5) then
	    ff=tan(rphi)/tan(rslo)
	  else
!  set factor of safety to fixed value for flat slopes 	  
	    ff=finf
	  end if
	  zns=float(nzs)
	  zinc=(zmax(i)-zmin)/zns
	  z=zmin
	  fmn=1.e25
	  rf=0.0
	  z_loop: do j=1,nzs+1
	    if (abs(a1)>1.e-5) then
	      fc(j)=c(zo(i))/(uws(zo(i))*z*a1*b1)
	    else
	      fc(j)=0.d0
	    end if
	    pzero(j)=beta*(z-depth(i))
	    zstar=z**2/(4.*dif(zo(i))/(b1*b1))
	    tstar=t/zstar
	    rf=0.0
	    temporal_loop: do n=1,nper
	      captstar1=capt(n)/zstar
	      tdif1=tstar-captstar1
	      if(tdif1 > 0.0) then 
	        x1=1./tdif1
	        x2=1./(sqrt(tdif1))
 	        rf1=sqrt(1./(x1*pi))*exp(-x1)
  	        rf2=derfc(x2)
	        rfa=rf1-rf2
	      else
	        rfa=0.0
	      end if
	      captstar2=capt(n+1)/zstar
	      tdif2=tstar-captstar2
	      if(tdif2 > 0.0) then	
	        x3=1./tdif2
	        x4=1./(sqrt(tdif2))
	        rf3=sqrt(1./(x3*pi))*exp(-x3)
	        rf4=derfc(x4)
	        rfb=rf3-rf4
	      else
	        rfb=0.0
 	      end if
	      rf=rf+rik(i+(n-1)*imax)*(rfa-rfb)
	    end do temporal_loop
  	    ptran(j)=z*rf
 	    p(j)=pzero(j)+ptran(j)
	    bline(j)=z*beta
	    ptest=p(j)-bline(j)
	    if(ptest > 0.0) then
	      p(j)=bline(j)
	    end if
	    if (abs(a1)>1.e-5) then
	      fw(j)=-(p(j)*uww*tan(rphi))/(uws(zo(i))*z*a1*b1)
	    else
	      fw(j)=0.d0
	    end if
	    z=z+zinc
	  end do z_loop
! find new height of rising water table in zones of upward seepage   
	  if(rikzero(i)<0.0) then
 	    zinc=(zmax(i)-zmin)/zns
	    z=zmin
	    newdep=0.0
   	    z_loop_a: do j=1,nzs+1
	      if(p(j)<0.0) newdep=z
	      z=z+zinc
   	    end do z_loop_a
! adjust presures 
   	    z=zmin
	    z_loop_b: do j=1,nzs+1
	      if(p(j)>0.0 .and. z<newdep) p(j)=0.d0
	      if(p(j)>=0.0 .and. z>=newdep) p(j)=beta*(z-newdep)
	      z=z+zinc	    
	    end do z_loop_b
   	    end if
   	    z=zmin
	  fs_loop: do j=1,nzs+1
	    fs=ff+fw(j)+fc(j)
! frictional strength cannot be less than zero 	
	    if ((ff+fw(j))<0.) fs=fc(j)
 	    if (fs>finf) fs=finf
 	    if (z<=1.e-02) fs=finf 
  	    if (flag==-1) write(u1,'(6(g12.5,1x):)') z,p(j),fs
  	    if (flag==-2) write(u1,'(6(g12.5,1x):)') z,p(j),pzero(j),ptran(j),&
             & bline(j),fs
	    if (fs<fmn) then
	      fmn=fs
	      zfmin(i)=z
	      pmin(i)=p(j)
	    end if
	    z=z+zinc
	  end do fs_loop
!  next statement assumes that computations begin at surface and work downward   
	  fsmin(i)=fmn
	  if (mod(i,2000)==0) write (*,*) i
   	end do grid_loop
   	write(*,*) imax, ' cells completed' 
   	write(ulog,*) imax, ' cells completed'
	if (flag==-1 .or. flag==-2) close(u1)
	return
  70	continue
  	write (*,*) 'Error opening output file'
	write (*,*) '--> ',profil
	write (*,*) 'Check file path and status'
	pause 'Press RETURN to exit'
	stop '-71'
	end subroutine iverson
