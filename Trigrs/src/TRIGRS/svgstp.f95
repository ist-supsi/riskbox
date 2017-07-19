!  finite depth diffusion model for rain infiltration.
!  by W.Z. Savage, Spring 2001, with modifications by R.L. Baum
!  (both) USGS	
	subroutine svgstp(u1,rikf,profil,&
        & ulog,i,rf,nccs,lcvs,nmx)
	use grids; use input_vars
	use model_vars
	implicit none
	integer:: i,j,u1,nmx,ulog,n1,nccs,n,m,mm,jf
	real:: rikf(nts+1) 
	real (double):: finf,t1,t2,term1,term2,rn,znew,t0
	real (double):: fierfc1,fierfc2,fierfc3,fierfc4
	real (double):: a1,b1,zns,zinc,tdif1,tdif2,z,newdep
	real (double):: ar1,ar2,ar3,ar4,fs,rf(nzs+1)
	real (double):: rfa,rfb,derfc,ff,rslo,rphi,fmn,ptest
	real (double):: tol,delt1,delt2,t1old,t2old
	real (double):: ddg2rad,dusz1,dlz
	logical:: lcvs 
	character (len=255)::  profil
	pi=3.141592653589793
	ddg2rad=pi/180.D0
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
	dlz=zmax(i)-dusz1
	temporal_loop: do m=1,nts+1 
	  t0=tcap(m) 
	  fmn=1.e25
	  jf=jsav(m) 
	  if (flag==-1 .or. flag==-2) then 
	    write (u1,'(i12,f6.1,i12,2x,g14.8)')i,rslo/ddg2rad,m,t0
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
	    if(z<dusz1 .or. rikf(m)==0.)then
	      rf(j)=0.0 
	    else
	      rf(j)=0.0 
	      	temporal_loop_1: do mm=1,nper
	        tdif1=t0-capt(mm)
	        if(tdif1 > 0.0) then
!  corrected diffusivity term in next line (divide by b1*b1) 	      
                  t1=sqrt(tdif1*dif(zo(i))/(b1*b1))
	          if (t1<1.0e-29) t1=1.0e-29
	          term1=0.0
	          series_a: do n=1,mmax
	            rn=float(n)
	            ar1=((2.*rn-1.)*dlz-(zmax(i)-z))/(2.*t1)
	            ar2=((2.*rn-1.)*dlz+(zmax(i)-z))/(2.*t1)
	            fierfc1=exp(-ar1**2)/sqrt(pi)-ar1*derfc(ar1)
	            fierfc2=exp(-ar2**2)/sqrt(pi)-ar2*derfc(ar2)
! test for convergence of series to within 1/10000 of previous value	      
	            t1old=term1
	            tol=term1/10000.
   	            term1=term1+fierfc1+fierfc2
 	            delt1=abs(term1-t1old)
 	            n1=n
 	            if(delt1<=tol) exit
   	          end do series_a
                  if(lcvs .and. delt1>tol) then
                    nccs=nccs+1
                    nv(i)=1
                    lcvs=.false.
                  end if
  	          if(n1>nmx) nmx=n1
   	          if(n1<nmn) nmn=n1
	          rfa=2.*t1*term1
	        else
	          rfa=0.0
	        end if
	        tdif2=t0-capt(mm+1)
	        if(tdif2 > 0.0) then
!  corrected diffusivity term in next line (divide by b1*b1) 	      
                  t2=sqrt(tdif2*dif(zo(i))/(b1*b1))
	          if (t2<1.0e-29) t2=1.0e-29
	          term2=0.0 
	          series_b: do n=1,mmax
	            rn=float(n)
	            ar3=((2.*rn-1.)*dlz-(zmax(i)-z))/(2.*t2)
	            ar4=((2.*rn-1.)*dlz+(zmax(i)-z))/(2.*t2)
	            fierfc3=exp(-ar3**2)/sqrt(pi)-ar3*derfc(ar3)
	            fierfc4=exp(-ar4**2)/sqrt(pi)-ar4*derfc(ar4)
! test for convergence of series to within 1/10000 of previous value	      
	            t2old=term2
	            tol=term2/10000.
   	            term2=term2+fierfc3+fierfc4
 	            delt2=abs(term2-t2old)
 	            n1=n
 	            if(delt2<=tol) exit
   	          end do series_b
   	          if(lcvs .and. delt2>tol) then
   	            nccs=nccs+1
                    nv(i)=1
   	            lcvs=.false.
   	          end if
   	          if(n1>nmx) nmx=n1
   	          if(n1<nmn) nmn=n1
	          rfb=2.*t2*term2
	        else
	          rfb=0.0
	        end if
   	        rf(j)=rf(j)+rik(i+(mm-1)*imax)*(rfa-rfb)
   	      end do temporal_loop_1 
	    end if
	    bline(j)=z*beta
            ptran(j)=rf(j)
   	    p(j)=pzero(j)+ptran(j)
	    ptest=p(j)-bline(j)
	    if(ptest > 0.0) then
	      p(j)=bline(j)
	    end if
	    z=z+zinc
   	  end do Z_loop  
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
  	    if (flag==-2) write(u1,'(6(g12.5,1x):)') z,p(j),&
  	    & pzero(j),ptran(j),bline(j),fs
	    if (jf>0) then
	      if (fs<fmn) then
	        fmn=fs 
	        zfmin(i+(jf-1)*imax)=z
	        pmin(i+(jf-1)*imax)=p(j)
	      end if
	    end if
	    z=z+zinc
   	  end do Z_FS_loop
!  next statement assumes that computations begin at surface and work downward   
	  if (jf>0) fsmin(i+(jf-1)*imax)=fmn
   	    end do temporal_loop 
	return
	end subroutine svgstp
