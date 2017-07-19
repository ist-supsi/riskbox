	subroutine flux(ic,iper,t0,j1,&
     &lcv,ncc,nv0,lwt)
	use grids; use input_vars
	use model_vars
	implicit none
	integer:: ic
	integer:: j1,iper,i,i1,k,k1,ncc,nv0
	real (double):: al,qa 
	real (double):: t0,tf,q1a,q1b,q2a,q2b,qta,qtb,bot,qtop
	real (double):: q2old,delq2,tdif1,tdif2,tstar
	real (double):: q0a,q0b 
	real (double):: z,psih,qlb,ck,th,b1 
	logical:: lcv,lwt
	b1=cos(slo(ic))
! coordinate transformation, corrected 9/5/06	
	al=alp(zo(ic))*b1*b1
	qa=rizero(ic)
	tf=al*ks(zo(ic))/(ths(zo(ic))-thr(zo(ic)))
	qt=0.0; qta=0.0; qtb=0.0 
	nmax1=0
	nmn=nmax+1
	do  i=1,iper
	  i1=1+(i-1)*tx
      	  tdif1=t0-capt(i)
	  if(tdif1 > 0.0) then
	    tstar=tf*tdif1
	    if(tstar < smt) then
! Early-time solution (ETS) ...............
		z=dusz
		qlb=ks(zo(ic))
		call smallt(tstar,al,qa,q(i),dusz,&
                 &qlb,ths,thr,qta,psih,ck,th,z)  
	        nmn=1
	        write(*,*) 'time, t*', t0,tstar, ' using ETS for basal flux'
	    else    
! later-time solution ...............	    
   	      q2a=0.0
	      do k=1,nmax
	        bot=1.+al*dusz/2.+2.*al*dusz*r(k)**2
 	        qtop=r(k)*sin(r(k)*al*dusz)*&
                &exp(-(r(k)**2)*tstar)
     	        q2old=q2a
	        q2a=q2a+qtop/bot
	        delq2=abs((q2a-q2old)/q2a)
	        k1=k
	        if(delq2<=1.0e-06) exit 
	        if(abs(q2a)<=eps .and. k>3) exit 
            end do
              if(lcv) then
                if((delq2>1.0e-06) .and. (abs(q2a)>eps .and. k>3)) then
                  ncc=ncc+1
                  nv0=1
                  lcv=.false.
                end if
              end if
 	      if(k1>nmax1) nmax1=k1
	      if(k1<nmn) nmn=k1
	      q0a=q(i)
	      q1a=4.d0*(q(i)-qa)*exp(al*dusz/2.)
	      qta=q0a-q1a*q2a*exp(-tstar/4.)
  	    end if
	  else 
	    qta=0.0
	  end if
	  tdif2=t0-capt(i+1)
	  if(tdif2 > 0.0) then
	    tstar=tf*tdif2
	    if(tstar < smt) then
! Early-time solution ...............	    
	      call smallt(tstar,al,qa,q(i),dusz,&
               &qlb,ths,thr,qtb,psih,ck,th,z)  
	      nmn=1
	        write(*,*) 'time, t*', t0,tstar, ' using ETS for basal flux'
	    else
! later-time solution ...............	    
              q2b=0.0
	      do  k=1,nmax
                bot=1.+al*dusz/2.+2.*al*dusz*r(k)**2
	        qtop=r(k)*sin(r(k)*al*dusz)*&
                &exp(-(r(k)**2)*tstar)
     	        q2old=q2b
	        q2b=q2b+qtop/bot
	        delq2=abs((q2b-q2old)/q2b)
	        k1=k
	        if(delq2<=1.0e-06) exit 
	        if(abs(q2b)<=eps .and. k>3) exit 
   	      end do
              if(lcv) then
                if((delq2>1.0e-06) .and. (abs(q2b)>eps .and. k>3)) then
                  ncc=ncc+1
                  nv0=1
                  lcv=.false.
                end if
              end if
	      if(k1>nmax1) nmax1=k1
	      if(k1<nmn) nmn=k1
	      q0b=q(i)
	      q1b=4.d0*(q(i)-qa)*exp(al*dusz/2.)
	      qtb=q0b-q1b*q2b*exp(-tstar/4.)
   	    end if
      	  else
            qtb=0.0
      	  end if
      	    qt=qt+qta-qtb
	if(qmax<qt .or. qt<0) then 
	  write(*,*) 'Error computing basal flux, output exceeds input!'
    	  write(*,*) 'Cell, Time, Max. Input flux, Basal flux: ',ic,t0,qmax,qt
    	  write(*,*) ''
	end if
    	end do
	return
	end subroutine flux


