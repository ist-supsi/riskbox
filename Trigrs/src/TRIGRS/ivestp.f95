!  Implementation of Iverson's (2000) method of computing pore pressure
!  for rain infiltration.
!  by W.Z. Savage, spring 2001, with modifications by R.L. Baum
!  (both) USGS
!  
	subroutine ivestp(u1,rikf,profil,&
        & ulog,i,rf,nccs,lcvs,nmx)
	use grids; use input_vars
	use model_vars
	implicit none
	integer:: i,j,jf,u1,ulog,nccs,nmx 
	integer:: n,nn
	real:: rikf(nts+1),finf 
	real (double):: derfc,a1,b1,ff,zns,zinc,z,t0,znew
	real (double):: zstar,tstar,x1,x2,x3,x4
	real (double):: rf1,rf2,rf3,rf4,rfa,rfb,rf(nzs+1)
	real (double):: fs,rslo,rphi,fmn,ptest
	real (double):: captstar1,captstar2,tdif1,tdif2
	real (double):: dusz1,ddg2rad,newdep
	logical:: lcvs 
	character(len=255):: profil
	pi=3.141592653589793
	ddg2rad=pi/180.D0
!  maximum value of Factor of Safety	
 	finf=10.
 	nmx=0
 	nmn=1+mmax
 	rslo=slo(i) 
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
	else !  set factor of safety to fixed value for flat slopes
	  ff=finf
	end if
	zns=float(nzs)
	zinc=(zmax(i)-zmin)/zns
	dusz1=0.
	temporal_loop: do n=1,nts+1
	  t0=tcap(n)
	  fmn=1.e25
	  jf=jsav(n)
	  if (flag==-1 .or. flag==-2) then
	    write (u1,'(i12,f6.1,i12,2x,g14.8)')i,rslo/ddg2rad,n,t0
	  end if
	  z=zmin
	  Z_loop: do j=1,nzs+1
	    znew=z
	    if(znew < 1.0e-30) znew =1.0e-30
	    if (abs(a1)>1.e-5) then
	      fc(j)=c(zo(i))/(uws(zo(i))*znew*a1*b1)
	    else
	      fc(j)=0.d0
	    end if
	    pzero(j)=beta*(z-depth(i))
	    if(z<dusz1 .or. rikf(n)==0.) then
	      rf(j)=0.0
	    else
	      rf(j)=0.0
	      zstar=z**2/(4.*dif(zo(i))/(b1*b1))
	      tstar=t0/zstar
	      temporal_loop_1: do nn=1,nper
	        captstar1=capt(nn)/zstar 
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
	        captstar2=capt(nn+1)/zstar
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
 	        rf(j)=rf(j)+rik(i+(nn-1)*imax)*(rfa-rfb)
  	      end do temporal_loop_1
	    end if
	    bline(j)=z*beta
  	    ptran(j)=z*rf(j)
 	    p(j)=pzero(j)+ptran(j)
	    ptest=p(j)-bline(j)
	    if(ptest > 0.0) then
	      p(j)=bline(j)
	    end if
	    z=z+zinc
   	  end do   Z_loop
! find new height of rising water table in zones of upward seepage   
	  if(rikzero(i)<0.0) then
 	    zinc=(zmax(i)-zmin)/zns
	    z=zmin
	    newdep=0.0
   	    do j=1,nzs+1
	      if(p(j)<0.0) newdep=z
	      z=z+zinc
            end do
! adjust presures 
   	    z=zmin
	    do j=1,nzs+1
	      if(p(j)>0.0 .and. z<newdep) p(j)=0.d0
	      if(p(j)>=0.0 .and. z>=newdep) p(j)=beta*(z-newdep)
	      z=z+zinc	    
   	    end do
   	  end if
! Compute factor of safety & save results  	  
   	  z=zmin
	  Z_FS_loop: do j=1,nzs+1 
	    if (abs(a1)>1.e-5 .and. z>1.e-30) then
	      fw(j)=-(p(j)*uww*tan(rphi))/(uws(zo(i))*z*a1*b1)
	    else
	      fw(j)=0.d0
	    end if
	    fs=ff+fw(j)+fc(j)
! frictional strength cannot be less than zero 	
	    if ((ff+fw(j))<0.) fs=fc(j)
 	    if (fs>finf) fs=finf
 	    if (z<=1.e-02) fs=finf 
  	    if (flag==-1) write(u1,'(6(g12.5,1x):)') z,p(j),fs
  	    if (flag==-2) write(u1,'(6(g12.5,1x):)') z,p(j)&
             & ,pzero(j),ptran(j),bline(j),fs
	    if(jf>0) then
	      if (fs<fmn) then
	        fmn=fs
	        zfmin(i+(jf-1)*imax)=z
	        pmin(i+(jf-1)*imax)=p(j)
	      end if
	    end if
	    z=z+zinc
   	  end do Z_FS_loop
	  if (jf>0) fsmin(i+(jf-1)*imax)=fmn
  	end do temporal_loop
	return
	end subroutine ivestp
