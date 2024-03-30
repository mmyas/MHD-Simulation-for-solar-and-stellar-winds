subroutine avmnmx(idx)
  use nrtype
  USE vars
  implicit none

  integer(I4B),intent(in) :: idx
  integer(I4B) :: ip
  real(DP) :: Tmaxpr,rvpr,vpr

  if(idx==0)then  ! initial
     vsum = 0.d0
     vmin = 1.d10
     vmax = 0.d0
     rvsum = 0.d0
     rvmin = 1.d10 
     rvmax = 0.d0
     Tsum = 0.d0
     Tmin = 1.d10 
     Tmax = 0.d0
     tmsum = 0.d0
  else if(idx==2)then ! final
     vsum = vsum / tmsum
     rvsum = rvsum / tmsum
     Tsum = Tsum / tmsum

     vsum = vsum * akm
     vmax = vmax * akm
     vmin = vmin * akm
     Tsum = Tsum * spetkl
     Tmax = Tmax * spetkl
     Tmin = Tmin * spetkl
     rvsum = rvsum * akm * rdkm*rdkm * rhocgs * 4.d15 * pi * 1.586d-26
     rvmax = rvmax * akm * rdkm*rdkm * rhocgs * 4.d15 * pi * 1.586d-26
     rvmin = rvmin * akm * rdkm*rdkm * rhocgs * 4.d15 * pi * 1.586d-26
     open(unit=52,file='avmnmx2.dat')
     write(52,*)rvsum,rvmax,rvmin,vsum,vmax,vmin,Tsum,Tmax,Tmin
     close(52)
  else ! average
     tmsum = tmsum + dt
     Tmaxpr = 0.d0
     do ip=1,np
        if(speene(ip) > Tmaxpr)Tmaxpr = speene(ip)
     end do
     Tsum = Tsum + Tmaxpr*dt
     if(Tmax<Tmaxpr)Tmax = Tmaxpr
     if(Tmin>Tmaxpr)Tmin = Tmaxpr
     rvpr = 0.d0
     vpr = 0.d0
     do ip=np-200,np-101
        rvpr = rvpr + densit(ip)*dxdt(ip)*radius(ip)*radius(ip)
        vpr = vpr + dxdt(ip)
     end do
     rvpr = 0.01d0 * rvpr
     vpr = 0.01d0 * vpr
     vsum = vsum + vpr*dt
     rvsum = rvsum + rvpr*dt
     if(vmax<vpr)vmax = vpr
     if(vmin>vpr)vmin = vpr
     if(rvmax<rvpr)rvmax = rvpr
     if(rvmin>rvpr)rvmin = rvpr
  end if

end subroutine avmnmx
!*******************************************************************
!subroutine engflx(ef)
subroutine engflx(index)
  use nrtype
  USE vars
  implicit none       

!  ------------------------------------------------------------------
!  real*8,intent(out) :: ef(9,nd)
  integer,intent(in) :: index
  real(DP) :: radint,Aflow,pflux,drad,Temp,fltb,dradflx,mu,xe
  integer(I4B) :: ip,j

  if(index==0)then
     ef(1:10,1:np) = 0.d0
     tefsum = 0.d0
     return
  end if

!--- print out energy flux ----------------------------
  radint = 0.d0
  do j=1,5
     radems(j) = 0.d0
  enddo
  tefsum = tefsum + dt
  do ip=1,np-1
     Aflow = radius(ip) * radius(ip) * fltb(radius(ip))
     ef(1,ip) = ef(1,ip) - dt * Aflow * fc2*Bparal(ip)* (Btrans(ip)*vtrans(ip) &
          & + Btrant(ip)*vtrant(ip) )
!     pflux = Aflow * densit(ip) * dxdt(ip)
     pflux = dt * Aflow * densit(ip) * dxdt(ip)
     ef(10,ip) = ef(10,ip) + pflux
     ef(2,ip) = ef(2,ip) + pflux * .5d0 * dxdt(ip)*dxdt(ip)
     ef(3,ip) = ef(3,ip) + pflux * .5d0 * (vtrans(ip)*vtrans(ip) &
          & + vtrant(ip)*vtrant(ip) )
     ef(4,ip) = ef(4,ip) + pflux * (speene(ip) + prsgas(ip) / densit(ip) )
     ef(5,ip) = ef(5,ip) - pflux * grav / radius(ip) 
     ef(6,ip) = ef(6,ip) + dt * dxdt(ip)*Aflow &
          & * fc2*(Btrans(ip)*Btrans(ip)+Btrant(ip)*Btrant(ip)) 
     if(ip==1)then
        ef(7,ip) = ef(7,ip) - dt * Aflow * cond*speene(ip)**2.5d0 &
             & * (speene(ip+1)-speene(ip-1))/(radius(ip+1)-radius(ip))
     else 
        ef(7,ip) = ef(7,ip) - dt * Aflow * cond*speene(ip)**2.5d0 &
             & * (speene(ip+1)-speene(ip-1))/(radius(ip+1)-radius(ip-1))
     endif
     if(ip <= 1)then
!     if(ip.le.50)then !for ?00km for Sun (sw657 model)
!         if(ip.le.82)then ! for 700km of Sun (sw657 model)
!         if(ip.le.106)then ! for 1000km of Sun (sw657 model)
        radint = 0.d0
        drad = 0.d0
     else 
!            drad = radcol(ip)*Aflow * .5*(radius(ip+1)-radius(ip-1)) 
                       ! along flux tube 
        drad = radcol(ip)* .5d0 *(radius(ip+1)-radius(ip-1)) 
        dradflx = drad * Aflow
        drad = drad * radius(ip)*radius(ip) ! for spherical 
!        radint = radint + drad
        radint = radint + dradflx
     endif
     ef(8,ip) = ef(8,ip) + dt * radint 
     ef(9,ip) = ef(1,ip) + ef(2,ip) + ef(3,ip) + ef(4,ip) + ef(5,ip) &
          & + ef(6,ip) + ef(7,ip) + ef(8,ip)

     call tp_mu(speene(ip),Temp,mu)
     if(t>t_e2T_rss.and.Temp<2.d4)call e2T(ip,Temp,mu,xe)
!     Temp = speene(ip)*spetkl            
!     if(Temp < 3.d3)then
!        Temp = 2.*Temp
!     else if(Temp < 1.d4)then
!        Temp = .5714*Temp + 4285.7
!     endif

     if(Temp < 7.d3)then
!         if(Temp.lt.1.e4)then
        radems(1) = radems(1) + drad
     else if(Temp < 2.d4)then
        radems(2) = radems(2) + drad
     else if(Temp < 5.d5)then
        radems(3) = radems(3) + drad
     else if(Temp < 1.d6)then
        radems(4) = radems(4) + drad
     else
        radems(5) = radems(5) + drad
     endif
!        write(*,*)ip,radcol(ip)
!        write(13,'(1p10e13.5)')radius(ip),trflux,qkinpl,qkintr,entha
!     &        ,grv,btr,cdflx,radnow,ten
  enddo

!  write(*,*)grav,1000,radius(1000),ef(5,1000)/ef(10,1000)*radius(1000) &
!       &,ef(5,1000),ef(10,1000),tefsum
!       & ,5000,radius(5000),ef(5,5000)/ef(10,5000)*radius(5000)

end subroutine engflx

!*******************************************************************
subroutine conduction
  use nrtype
  USE vars
  implicit none       
  
! This subroutine calculate (Lagrange) conduction & radiation cooling 

  real(DP), parameter :: epsln = 1.d-4
  real(DP),allocatable :: a(:),b(:),c(:),r(:),dsp(:),dspsum(:),spetmp(:)
  real(DP),allocatable :: colles(:),drm(:),drp(:)
  real(DP) :: r1,r2,r3,A1,A2,A3,spav,dtcd,spt,Temp,umenf,rhonow,dnsuse &
       & ,qlogT,radvol,drdvde,dTemp,qlogTpl,radvolpl,err,errmax,dttot &
       & ,fltb,rdctn,radcoltab,xe,condmu !,rdesc
  integer(I4B) :: ip,idtdt,icd
  
  allocate(a(nd),b(nd),c(nd),r(nd),dsp(nd),dspsum(nd),spetmp(nd) &
       & ,colles(nd),drm(nd),drp(nd))
!      goto 121 ! skip conduction

  ip=1
  r1 = radius(ip)  ! we have no mesh at 0
  r2 = radius(ip)
  r3 = radius(ip+1)
  A1 = r1*r1*fltb(r1)
  A2 = r2*r2*fltb(r2)
  A3 = r3*r3*fltb(r3)
  drm(ip) = (dxplac0+dradiu(ip))*dradiu(ip)*A2/(A1+A2)
  drp(ip) = (dradiu(ip)+dradiu(ip+1))*dradiu(ip)*A2/(A3+A2) 
                                ! 0.5 is cancelled both in nemer.& denom.
  do ip =2,np
     r1 = radius(ip-1)
     r2 = radius(ip)
     r3 = radius(ip+1)
     A1 = r1*r1*fltb(r1)
     A2 = r2*r2*fltb(r2)
     A3 = r3*r3*fltb(r3)
     drm(ip) = (dradiu(ip-1)+dradiu(ip))*dradiu(ip)*A2/(A1+A2)
     drp(ip) = (dradiu(ip)+dradiu(ip+1))*dradiu(ip)*A2/(A3+A2)
                                ! 0.5 is cancelled both in nemer.& denom.
  enddo

!----outflow boundary-----------------------
  speene(np) = speene(np-1)*(radius(np-1)/radius(np))**.285d0
  speene(np+1) = speene(np)*(radius(np)/radius(np+1))**.285d0
!------------------------------------------
  do ip=2,np-1
     if(speene(ip)<0.d0)then
        spav = .5d0*(speene(ip-1) + speene(ip+1))
        if(spav>0.d0)then
           speene(ip) = spav
        else
           speene(ip) = .9d0*max(speene(ip-1),speene(ip+1))
        endif
        write(*,*)'T1 < 0', ip,speene(ip),spav
     endif
  end do

!      cdfc = 1.

  dtcd = min(dt,1.d-9)
!      dtcd = min(dt,2.e-7)
  dttot = 0.d0 ! for division of conduction

 152    continue

  idtdt = 1 ! which is used for dtcd
  spt = 0.d0
  do ip = 1,np
     dspsum(ip) = 0. ! initialize total improvement
     spt = abs(speene(ip) ) + spt 
  enddo
  spt = spt / dble(np-2)

  do ip=1,np+1
!     if(ip<1000.and.speene(ip)>1000.d0)write(*,*)'bef',t,ip,radius(ip),speene(ip)*spetkl
     spetmp(ip) = speene(ip)
  enddo

!      call clscnd(dtcd,spetmp,colles)

!      do 109 icd = 1,ncd
! 151    continue                  ! return here if dtcd is reduced 
  do icd = 1,50
     ip = 1
! conduction part
!      sp1 = (.5*(speene0 + spetmp(ip)))
!      sp2 = (.5*(spetmp(ip) + spetmp(ip+1)))
!      sptt1 = sp1**1.5
!      sptt2 = sp2**1.5
! radiation part
!     Temp = spetmp(ip)*spetkl
!     if(Temp<3.d3)then
!        Temp = 2.d0*Temp
!     else if(Temp<1.d4)then
!!         umenf = dmin1(2.d0,(1.d4 - Temp)/7.d3 + 1.d0) ! 7.e3 = 1.e4-3.e3
!!         Temp = Temp * umenf
!        Temp = .5714d0*Temp + 4285.7d0
!        umenf = 1.6d0
!     else 
!        umenf = 1.d0
!     endif
     call tp_mu(speene(ip),Temp,umenf)
     if(t>t_e2T_rss.and.Temp<2.d4)call e2T(ip,Temp,umenf,xe)
     rhonow = rhocgs*densit(ip)
     dnsuse = min(rhonow , rhocrt) / rhocgs
!      dnsuse = dmin1(rhocgs*densit(ip) , rhocrt) / rhocgs
!      if(Temp.le.1.e4)then
!         radvol = 0.
!         drdvde = 0.
!      else
!         qlogT = log10(Temp)
!         radvol = 10.**rdesc(qlogT) * radfac
!         dTemp = 0.01 * Temp
!         qlogTpl = log10(Temp + dTemp)
!         radvolpl = 10.**rdesc(qlogTpl) * radfac
!         drdvde = (radvolpl - radvol)*spetkl/dTemp 
!      endif
!      if(Temp.le.5.e3)then
!c      if(Temp.le.3.e3)then
!         radvol = 0.
!         drdvde = 0.
!      else if(Temp.le.2.e4 .and. rhonow.gt.rhocrt)then 
     if(Temp<=Tcrt .and. rhonow>rhocrt)then 
!     if(Temp<2.d4 .and. rhonow>rhocrt)then 
                                           !     chromospheric cooling
!        qlogT = 4.302d0
!        radvol = 10.d0**rdesc(qlogT) * radfac ! /umenf/umenf
                                ! no need for umenf since cool ~ rho       
        radvol = lambdacrt * radfac
!        if(Temp < 4.d3)radvol = radvol*exp(.2d0*Temp - 8.d2)
        if(Temp < Tclso)radvol = radvol*exp(.2d0*(Temp - Tclso))
                                !     gradual switch-off; tau_T=1K
        radcol(ip) = radvol * densit(ip) * dnsuse 
                                !     preserve for flux after switch-off
!         if(densit(ip).gt.1.e-6)radvol = radvol*1.e-6/densit(ip)
        drdvde = 0.d0
     else
        qlogT = log10(Temp)
!        radvol = 10.d0**rdesc(qlogT) * radfac
        radvol = 10.d0**radcoltab(qlogT) * radfac
        dTemp = 0.01d0 * Temp
        qlogTpl = log10(Temp + dTemp)
!        radvolpl = 10.d0**rdesc(qlogTpl) * radfac
        radvolpl = 10.d0**radcoltab(qlogTpl) * radfac
        drdvde = (radvolpl - radvol)*spetkl/dTemp 
        radcol(ip) = radvol * densit(ip) * dnsuse ! preserve for flux
     endif
! calculation of matrix elements
!      a(ip) = colles(ip)*cond * sptt1 * (.75*spetmp(ip) - 1.75*speene0) 
!     &     / drm(ip)
!      b(ip) = densit(ip)/dtcd + colles(ip)*cond 
!     &     * (sptt1*(1.75*spetmp(ip)-.75*speene0) / drm(ip)
!     &     + sptt2*(1.75*spetmp(ip)-.75*spetmp(ip+1)) / drp(ip) )
!     &     + drdvde * densit(ip) * dnsuse * cdfc
!      c(ip) = - colles(ip)*cond * sptt2 
!     &     * (1.75*spetmp(ip+1)-.75*spetmp(ip)) / drp(ip)
!      r(ip) = colles(ip)*cond 
!     &     * ( sp2**2.5 * ( spetmp(ip+1) - spetmp(ip) ) / drp(ip)
!     &     - sp1**2.5 * ( spetmp(ip) - speene0 ) / drm(ip) )
!     &     - densit(ip)/dtcd * dspsum(ip)
!     &     - radvol * densit(ip) * dnsuse * cdfc
!     if(Temp>=1.d0)then  ! should be 1e4!!
!        a(ip) = - cond *  speene0**2.5d0 / drm(ip)
!        b(ip) = densit(ip)/dtcd + cond * spetmp(ip)**2.5d0 &
!             & * (1.d0/drm(ip) + 1.d0/drp(ip)) + drdvde * densit(ip) * dnsuse 
!        c(ip) = - cond  * spetmp(ip+1)**2.5d0 / drp(ip)
!        r(ip) = cond  * .2857d0 * ( ( spetmp(ip+1)**3.5d0 - spetmp(ip)**3.5d0)&
!             & / drp(ip) - ( spetmp(ip)**3.5d0 - speene0**3.5d0 ) / drm(ip)) &
!             & - densit(ip)/dtcd * dspsum(ip) - radvol * densit(ip) * dnsuse 
!     else

     condmu = cond*(umenf/mufull)**3.5d0 ! mu is normalized for mufull
     a(ip) = - condmu *  speene0**2.5d0 / drm(ip)
     b(ip) = densit(ip)/dtcd + condmu * spetmp(ip)**2.5 &
          & * (1.d0/drm(ip) + 1.d0/drp(ip)) + drdvde * densit(ip) * dnsuse 
     c(ip) = - condmu  * spetmp(ip+1)**2.5d0 / drp(ip)
     r(ip) = condmu  * .2857d0 * ( ( spetmp(ip+1)**3.5d0 &
          & - spetmp(ip)**3.5d0 ) / drp(ip) - ( spetmp(ip)**3.5d0 &
          & - speene0**3.5d0 ) / drm(ip)) - densit(ip)/dtcd * dspsum(ip) &
          & - radvol * densit(ip) * dnsuse
!     endif

!     if(colles(ip).lt..9999)then
!         abcnd = colles(ip)*cond 
!     &        * (sptt2*sp2*(spetmp(ip+1)-spetmp(ip))/drp
!     &        - sptt1*sp1*(spetmp(ip)-speene0)/drm)
!         cacnd = abcnd/(spetmp(ip+1)-speene0)
!         b(ip) = b(ip) + abcnd/spetmp(ip)
!         a(ip) = a(ip) - cacnd
!         c(ip) = c(ip) + cacnd
!      endif
!      write(*,*)ip,r(ip)

     do ip = 2,np
! conduction part
!         sp1 = (.5*(spetmp(ip-1) + spetmp(ip)))
!         sp2 = (.5*(spetmp(ip) + spetmp(ip+1)))
!         sptt1 = sp1**1.5
!         sptt2 = sp2**1.5
! radiation part
!        Temp = spetmp(ip)*spetkl
!        if(Temp<3.d3)then
!           Temp = 2.d0*Temp
!        else if(Temp<1.d4)then
!!            umenf = dmin1(2.d0,(1.d4 - Temp)/7.d3 + 1.d0) ! 7.e3 = 1.e4-3.e3
!!            Temp = Temp * umenf
!           Temp = .5714d0*Temp + 4285.7d0
!           umenf = 1.6d0
!!            write(*,*)ip,Temp,spetmp(ip),umenf
!        else
!           umenf = 1.d0
!        endif
        call tp_mu(speene(ip),Temp,umenf)
        if(t>t_e2T_rss.and.Temp<2.d4)call e2T(ip,Temp,umenf,xe)
        rhonow = rhocgs*densit(ip)
        dnsuse = min(rhonow , rhocrt) / rhocgs
!         dnsuse = dmin1(rhocgs*densit(ip) , rhocrt) / rhocgs
!         if(Temp.le.1.e4)then
!            radvol = 0.
!            drdvde = 0.
!         else
!            qlogT = log10(Temp)
!            radvol = 10.**rdesc(qlogT) * radfac
!            dTemp = 0.01 * Temp
!            qlogTpl = log10(Temp + dTemp)
!            radvolpl = 10.**rdesc(qlogTpl) * radfac
!            drdvde = (radvolpl - radvol)*spetkl/dTemp 
!         endif
!         if(Temp.le.5.e3)then
!c         if(Temp.le.3.e3)then
!            radvol = 0.
!            drdvde = 0.
!         else if(Temp.le.2.e4 .and. rhonow.gt.rhocrt)then   
        if(Temp<=Tcrt .and. rhonow > rhocrt)then
!        if(Temp<=2.d4 .and. rhonow>rhocrt)then
                                   ! chromospheric cooling
!           qlogT = 4.302d0 ! T=2.e4
!           radvol = 10.d0**rdesc(qlogT) * radfac ! /umenf/umenf
                 ! no need for umenf since cool ~ rho
           radvol = lambdacrt * radfac
!           if(Temp<4.d3)radvol = radvol*exp(.2d0*Temp - 8.d2)
           if(Temp<Tclso)radvol = radvol*exp(.2d0*(Temp - Tclso))
                                !     gradual switch-off; tau_T=1K
           radcol(ip) = radvol * densit(ip) * dnsuse 
                                !     preserve for flux after switch-off
!            if(densit(ip).gt.1.e-6)radvol = radvol*1.e-6/densit(ip)
           drdvde = 0.d0
        else
           qlogT = log10(Temp)
!           radvol = 10.d0**rdesc(qlogT) * radfac
           radvol = 10.d0**radcoltab(qlogT) * radfac
!           if(ip==500.and.qlogT>4.d0)write(*,*)qlogT,rdesc(qlogT)
           ! ,radcoltab(qlogT)
           dTemp = 0.01d0 * Temp
           qlogTpl = log10(Temp + dTemp)
!           radvolpl = 10.d0**rdesc(qlogTpl) * radfac
           radvolpl = 10.d0**radcoltab(qlogTpl) * radfac
           drdvde = (radvolpl - radvol)*spetkl/dTemp 
           
!           if(Temp<4.d3)then ! Low T & Low density
           if(Temp<Tclso)then
!              radvol = radvol*exp(.2d0*Temp - 8.d2)
!              drdvde = drdvde*exp(.2d0*Temp - 8.d2)
              radvol = radvol*exp(.2d0*(Temp - Tclso))
              drdvde = drdvde*exp(.2d0*(Temp - Tclso))
           endif
           radcol(ip) = radvol * densit(ip) * dnsuse 
                                ! preserve for flux after switch-off
!           if(ip==1000)write(*,*)radius(ip),radcol(ip),radvol,qlogT,densit(ip)&
!                & ,dnsuse,radcoltab(qlogT),radfac          
        endif
!        if(t>1.19416d-2.and.ip<=5592.and.ip>=5590)write(*,*)radvol,Temp
! calculation of matrix elements
!         a(ip) = colles(ip)*cond * sptt1 * (.75*spetmp(ip) 
!     &        - 1.75*spetmp(ip-1)) / drm(ip)
!         b(ip) = densit(ip)/dtcd + colles(ip)*cond 
!     &        * ( sptt1*(1.75*spetmp(ip)-.75*spetmp(ip-1)) / drm(ip)
!     &        + sptt2*(1.75*spetmp(ip)-.75*spetmp(ip+1)) / drp(ip) )
!     &        + drdvde * densit(ip) * dnsuse *cdfc
!         c(ip) = - colles(ip)*cond * sptt2 * (1.75*spetmp(ip+1) 
!     &        - .75*spetmp(ip)) / drp(ip)
!         r(ip) = colles(ip)*cond 
!     &        * (sp2**2.5 * ( spetmp(ip+1) - spetmp(ip) ) / drp(ip)
!     &        - sp1**2.5 * ( spetmp(ip) - spetmp(ip-1) ) / drm(ip) )
!     &        - densit(ip)/dtcd * dspsum(ip)
!     &        - radvol * densit(ip) * dnsuse * cdfc
!        if(Temp>=1.d0)then  ! should be 1.e4!!
!           a(ip) = - cond * spetmp(ip-1)**2.5d0 / drm(ip)
!           b(ip) = densit(ip)/dtcd + cond * spetmp(ip)**2.5d0 &
!                & * (1.d0/drm(ip) + 1.d0/drp(ip)) + drdvde * densit(ip)*dnsuse 
!           c(ip) = - cond  * spetmp(ip+1)**2.5d0 / drp(ip)
!           r(ip) = cond * .2857d0 * ( ( spetmp(ip+1)**3.5d0 &
!                & - spetmp(ip)**3.5d0 ) / drp(ip) - ( spetmp(ip)**3.5d0 &
!                & - spetmp(ip-1)**3.5d0 ) / drm(ip)) - densit(ip)/dtcd &
!                & * dspsum(ip) - radvol * densit(ip) * dnsuse
!!           if(ip==150)write(*,*)spetmp(ip-1:ip+1)
!        else
        condmu = cond*(umenf/mufull)**3.5d0 ! mu is normalized for mufull
        a(ip) = - condmu * spetmp(ip-1)**2.5d0 / drm(ip)
        b(ip) = densit(ip)/dtcd + condmu * spetmp(ip)**2.5d0 &
             & * (1.d0/drm(ip) + 1.d0/drp(ip)) + drdvde * densit(ip)*dnsuse 
        c(ip) = - condmu * spetmp(ip+1)**2.5d0 / drp(ip)
        r(ip) = condmu * .2857d0 * ( ( spetmp(ip+1)**3.5d0 &
             & - spetmp(ip)**3.5d0 ) / drp(ip) - ( spetmp(ip)**3.5d0 &
             & - spetmp(ip-1)**3.5d0 ) / drm(ip)) - densit(ip)/dtcd &
             & * dspsum(ip) - radvol * densit(ip) * dnsuse 
!        endif

!        if(t>1.19416d-2.and.ip==5591)write(*,'(1p9e12.4)')r(ip),Temp &
!             & ,spetmp(ip),drp(ip),drm(ip),dspsum(ip),radvol,dnsuse,dtcd
       
!         if(colles(ip).lt..9999)then
!            abcnd = colles(ip)*cond 
!     &           * (sptt2*sp2*(spetmp(ip+1)-spetmp(ip))/drp
!     &           - sptt1*sp1*(spetmp(ip)-spetmp(ip-1))/drm)
!            cacnd = abcnd/(spetmp(ip+1)-spetmp(ip-1))
!            b(ip) = b(ip) + abcnd/spetmp(ip)
!            a(ip) = a(ip) - cacnd
!            c(ip) = c(ip) + cacnd
!         endif
!     if(Temp.ge.2.e4 .and. densit(ip).ge.1.e-4)write(*,*)
!     &        'rad in dense region',ip,Temp,qlogT,rdesc(qlogT)
!     &        ,densit(ip),r(ip)/b(ip),speene(ip),rdesc(qlogT)
!     &        ,10.**rdesc(qlogT),-radvol*densit(ip)*densit(ip)
!     &        ,cond*sp1**2.5 * ( speene(ip+1) - speene(ip) ) / drm
!         if(ip.eq.697)write(*,*)ip,dtcd,densit(ip),b(ip),dradiu(ip-1)
!     &        ,dradiu(ip),dradiu(ip+1),speene(ip-1),speene(ip)
!     &        ,speene(ip+1)  
!    if(ip.ge.200)write(*,*)ip,Temp,- radvol*densit(ip)*densit(ip)
!     &        ,sp2**2.5 * ( speene(ip) - speene(ip-1) ) / drm
     end do

!      write(*,*)ip,np,r(ip)

     call tridag(a,b,c,r,dsp,np)
!     if(t>1.19416d-2)write(*,*)t,r(5591)
!        do ip = 1,np
!           write(*,'(i4,1p5e12.4)')ip,a(ip),b(ip),c(ip),r(ip),dsp(ip)
!        end do
!     end if
     rdctn = 1.d0
     
     do ip = 2,np-1 
        if(abs(dsp(ip))> .5d0*abs(spetmp(ip)))then 
           rdctn = .5d0*abs(spetmp(ip)/dsp(ip) ) 
           dsp(ip) = dsp(ip)*rdctn
        endif
!         spetmp(ip) =  speene(ip) + dsp(ip)
     end do
     do ip=2,np-1
        spetmp(ip) = spetmp(ip) + dsp(ip)
     end do
!--- outflow boundary -------------------------------------------
!     spetmp(np) = spetmp(np-1)*(radius(np-1)/radius(np))**.285
!     dsp(np) =  spetmp(np) - speene(np)
!----------------------------------------------------------------
!---- sum of improvment & error estimate ----------------------------
     errmax = 0.d0
     err = 0.d0
     do ip = 2,np - 1
        dspsum(ip) = dsp(ip) + dspsum(ip)
!         err = err + abs(dsp(ip) ) 
        err = abs(dsp(ip)/speene(ip))
!         if(err.gt.errmax)errmax = err
        if(err>errmax)then
!           if(icd>30)write(*,*)ip,err,spetmp(ip),dsp(ip),2.d0*spetmp(ip)*spetkl
           errmax = err
        endif
!         speene(ip) = spetmp(ip)
        if(errmax==0.d0)write(*,*)'dsp&err',t,ip,spetmp(ip),dsp(ip),speene(ip),a(ip),b(ip),c(ip)
     end do                  ! end of 1 step of conduction part
!      err = err/dfloat(np-2)
!      do 124 ip = 2,np 
!         speene(ip) = spetmp(ip)
! 124  continue
!      write(*,*)icd,errmax
!      if (err.le. rdctn*epsln*spt) goto 121 ! iteration converged
     if(errmax<=epsln)goto 121 ! iteration converged

!      call clscnd(dtcs,spetmp,colles)
  end do
!      cdfc = .1*cdfc ! reduce cond & rad
  dtcd = .1d0*dtcd
  write(*,*)'No convergence in conduction; reduce',dtcd,dt,dttot
  goto 152

121 continue
  do ip=2,np-1
     speene(ip) = spetmp(ip)
!      if(ip<1000.and.speene(ip)>1000.d0)write(*,*)'aft',ip,radius(ip),speene(ip)*spetkl
  enddo
!-------------outgoing boundary ------------------------
  speene(np) = speene(np-1)*(radius(np-1)/radius(np))**.285d0
  speene(np+1) = speene(np)*(radius(np)/radius(np+1))**.285d0
!------------------------------------------
  dttot = dtcd + dttot
!  if(t>1.19416d-2)write(*,*)'conduction iteration',t,icd,dt,dttot,errmax

  if(dttot<.3d0*dt)then
!      dtcd = min(dt - dttot,2.d-7)
     dtcd = min(dt - dttot,dtcd0)
     goto 152
  else if(dttot<.99d0*dt)then
!     dtcd = min(dt - dttot,2.d-7)
     dtcd = min(dt - dttot,dtcd0)
!         write(*,*)'conduction division',dt,dtcd,dttot
!         call clscnd(dtcd,spetmp,colles)
     goto 152
  endif
   
!  write(*,*)t,radius(1000),radcol(1000),speene(1000),densit(1000)

  deallocate(a,b,c,r,dsp,dspsum,spetmp,colles,drm,drp)

end subroutine conduction
 
!-----------------------------------------------------
real(dp) function radcoltab(log10T)
  use vars
  implicit none

  real(DP),intent(in) :: log10T  
  integer(I4B) :: it

  if(log10T<=Ttab(1))then
     radcoltab = radtab(1)
  else if(log10T>=Ttab(nradcol))then
     radcoltab = radtab(nradcol)
  else
     it = int((log10T-Ttab0)/dTtab)+1

     radcoltab = radtab(it) + (radtab(it+1)-radtab(it)) &
          & / (Ttab(it+1)-Ttab(it)) * (log10T-Ttab(it)) 
!     write(*,*)'radcol',it,log10T,Ttab(it),Ttab(it+1)
  end if

end function radcoltab
!------------------------------------------------------
subroutine calcTcr(Tcrtrad)!(zmetal) !(rhocr,zmetal)
  use nrtype
  use vars
  implicit none

!  real(DP),intent(in) :: zmetal
  real(DP),intent(out) :: Tcrtrad
  real(DP) :: th
  integer(I4B) :: i, iztab
  integer(I4B),parameter :: nztab=4
  real(DP) :: ztab(nztab)
  data Ztab/1.d-10,1.d-2,1.d-1,1.d0/
  real(DP),allocatable :: radtab1(:),radtab2(:)

  if(zmetal==1.0d0)then
!     open(unit=17,file='SD93_cie_z1e0.dat',status='old')
     open(unit=17,file='SD93_cie_z7e-1.dat',status='old') ! from Asplund+ 2009
     do i = 1,nradcol
        read(17,*)Ttab(i),radtab(i)
     end do
     close(17)
  else if(zmetal==0.0d0)then
     open(unit=17,file='SD93_cie_z0e0.dat',status='old')
     do i = 1,nradcol
        read(17,*)Ttab(i),radtab(i)
     end do
     close(17)
  else if(zmetal==1.0d-1)then
!     open(unit=17,file='SD93_cie_z1e-1.dat',status='old')
     open(unit=17,file='SD93_cie_z7e-2.dat',status='old')
     do i = 1,nradcol
        read(17,*)Ttab(i),radtab(i)
     end do
     close(17)
  else if(zmetal==1.0d-2)then
!     open(unit=17,file='SD93_cie_z1e-2.dat',status='old')
     open(unit=17,file='SD93_cie_z7e-3.dat',status='old')
     do i = 1,nradcol
        read(17,*)Ttab(i),radtab(i)
     end do
     close(17)
  else if(zmetal<ztab(1).or. zmetal>ztab(nztab))then
     write(*,*)'no table for zmetal'
     stop
  else
     do iztab = 1,nztab
        if((zmetal-ztab(iztab))*(zmetal-ztab(iztab+1)) < 0.d0)exit
     end do
     select case(iztab)
     case(1)
        open(unit=171,file='SD93_cie_z0e0.dat',status='old')
        open(unit=172,file='SD93_cie_z7e-3.dat',status='old')
     case(2)
        open(unit=171,file='SD93_cie_z7e-3.dat',status='old')
        open(unit=172,file='SD93_cie_z7e-2.dat',status='old')
     case(3)
        open(unit=171,file='SD93_cie_z7e-2.dat',status='old')
        open(unit=172,file='SD93_cie_z7e-1.dat',status='old')
     case default
        write(*,*)'no table for zmetal interpolation'
        stop
     end select
     allocate(radtab1(nradcol),radtab2(nradcol))
     do i = 1,nradcol
        read(171,*)Ttab(i),radtab1(i)
     end do
     do i = 1,nradcol
        read(172,*)Ttab(i),radtab2(i)
     end do
!     write(*,*)iztab
     do i = 1,nradcol
        radtab(i) = radtab1(i) + (radtab2(i) - radtab1(i)) &
             & / (log10(ztab(iztab+1)) - log10(ztab(iztab))) *(log10(zmetal)-log10(ztab(iztab)))
!        write(*,*)i,radtab(i),radtab1(i),radtab2(i)
     end do
     deallocate(radtab1,radtab2)
     close(171); close(172)
  end if
  
  Ttab0 = Ttab(1)
  dTtab = Ttab(2) - Ttab(1)  
  if(dTtab/=(Ttab(3)-Ttab(2)))then
     write(*,*)'dTtab is not constant'
  end if
!  write(*,*)Ttab(1:50)

  lambdacrt = log10(4.5d-39/rhocrt*(0.2 + 0.8*zmetal))
  do i = 1,nradcol-1
     th = (lambdacrt - radtab(i)) * (lambdacrt - radtab(i+1))
     if(th<=0.d0)exit
  end do     
  Tcrtrad = Ttab(i) + (lambdacrt - radtab(i)) / (radtab(i+1) - radtab(i)) &
       & * (Ttab(i+1) - Ttab(i))

  lambdacrt = 10.d0**lambdacrt

!  write(*,*)'Tcr',i,lambdacrt,radtab(i),radtab(i+1),Ttab(i),Ttab(i+1),Tcrtrad
!  write(*,*)'Tcr',lambdacrt,Ttab0,dTtab,Tcrtrad

end subroutine calcTcr

!*******************************************************************
! FUNCTION rdesc(T1)
!   use nrtype 
!   USE vars
!   implicit none 
!   real(DP) :: rdesc
!   real(DP),intent(in) :: T1
!! This functon returs log10(Lambda (erg cm^3 s^-1)) for input log10(T(K))
!
!   rdesc = cofrad(1) + cofrad(2)*T1 + cofrad(3)*T1*T1 + cofrad(4)*T1*T1*T1    
!
! end FUNCTION rdesc
!*******************************************************************
 SUBROUTINE tridag(a,b,c,r,u,n)
   use nrtype
!   use vars
   implicit none
   INTEGER(I4B),intent(in) :: n
   integer(I4B),parameter :: NMAX=200000
   REAL(DP),intent(in) :: a(n),b(n),c(n),r(n)
   real(DP),intent(out) :: u(n)
   INTEGER(I4B) :: j
   REAL(DP) :: bet,gam(NMAX)

   if(b(1)==0.d0)then
      write(*,*) 'tridag: rewrite equations'
      stop
   end if
   bet=b(1)
   u(1)=r(1)/bet
   do j=2,n
      gam(j)=c(j-1)/bet
      bet=b(j)-a(j)*gam(j)
!        if(j.eq.608)write(*,*)j,bet,b(j),a(j),gam(j),c(j-1)
      if(bet==0.d0)then
         write(*,*)'tridag failed',j,b(j),a(j),gam(j)
!         pause
         stop
      endif
      u(j)=(r(j)-a(j)*u(j-1))/bet
!      if(t>1.19416d-2)write(*,'(i4,1p4e12.4)')j,u(j),r(j),a(j),bet
   end do
   do j=n-1,1,-1
      u(j)=u(j)-gam(j+1)*u(j+1)
   end do
 
END SUBROUTINE tridag
!*******************************************************************
function fltb(r)
  use nrtype
  USE vars
  implicit none   

  real(DP) :: fltb
  real(DP),intent(in) :: r
  real(DP) :: ps,f1,ex,fltbc

! This subroutine calculate area expansion factor, f, on r.

! In th corona => currently comment out

!  ps = (r-re1)/sgm
!  if(ps>10.d0)then
!     fltb=fmax
!     goto 11
!  endif
!  f1 = 1.d0 - (fmax-1.d0)*exp((1.d0-re1)/sgm)
!  ex = exp(ps)
!  
!  fltb=(fmax*ex+f1)/(ex+1.)
!11 continue

! In the chromosphere

  ps=(r-re1c)/sgmc
  if(ps>10.d0)then
     fltbc=fmaxc
  else
     f1 = 1.d0 - (fmaxc-1.d0)*exp((1.d0-re1c)/sgmc)
     ex = exp(ps)
     fltbc = (fmaxc*ex+f1) / (ex+1.d0)
  endif

!  fltb = fltb * fltbc
  fltb = fltbc
  
end function fltb

!**************************************************************************
FUNCTION rtnewt(funcd,x1,x2,xacc)
  use nrtype
  implicit none
  INTEGER(I4B),parameter :: JMAX=20
  REAL(DP) :: rtnewt
  real(DP),intent(in) :: x1,x2,xacc
  EXTERNAL :: funcd
  INTEGER(I4B) :: j
  REAL(DP) :: df,dx,f
  rtnewt=.5d0*(x1+x2)
  do j=1,JMAX
     call funcd(rtnewt,f,df)
     dx=f/df
     rtnewt=rtnewt-dx
     if((x1-rtnewt)*(rtnewt-x2)<0.d0)then
        write(*,*) 'rtnewt jumped out of brackets'
        stop
     end if
     if(abs(dx)<xacc) return
  end do
!  pause 'rtnewt exceeded maximum iterations'
  write(*,*) 'rtnewt exceeded maximum iterations'
  stop
END FUNCTION rtnewt

!***********************************************************************
!subroutine centmass(r,f,df)
!  implicit none
!
!  common/centm/dns,ddndr,drh
!
!      r2 = r*r
!      r3 = r2*r
!      
!      f = .33333 * dns * r3 + .25 * ddndr * r3*r 
!     &     - .25 * ddndr * drh*drh*drh*drh
!
!      df = dns * r2 + ddndr * r3
!
!c      write(*,*)r,drh,dns,ddndr,f,df
!
!      return
!      end

!***********************************************************************
subroutine outgobdsc
  use nrtype
  USE vars
  implicit none      

  real(DP) :: err,Btrtmp,Btttmp,dnstmp,pgtmp,pttmp,vxtmp,vtrtmp,vtttmp,sqdnsi &
       & ,sound,Alfvn,sp1,sp2,sp3,sp4,fast,slow,alf2,fst2,slw2 &
       & ,chr1,chr2,chr3,chr4,chr5,chr6,chr7,rad0,rad1,dr0&
       & ,f0,f1,geo1,geo2,alcmpy,alcmpz,alelmy,alelmz,div,trwvel &
       & ,qL1,qL2,qL3,qL4,qL5,qL6,qL7,qLcf,qLcs,qLca,qLmf,qLms,qLma&
       & ,elm1,elm2,elm3,elm4,elm5,elm6,elm7,elm8,dPtdt,dvxdt,Bymag,Bzmag &
       & ,Ptot,dPgdt,dBtrdt,dBttdt,dvtrdt,dvttdt,dPmdt,dPmdtdf,ddndt &
       & ,dBtrdt1,dBttdt1,dvtrdt1,dvttdt1,dPmdt1,dPmdtdf1,ddndt1,dPgdt1,dvxdt1 &
       & ,Bprs,fltb 
  integer(I4B) :: iter      

  do iter = 1,100
     err = 0.d0
     if(iter == 1)then
        Btrans(np) = Btrbd2
        Vtrans(np) = vtrbd2
        Btrant(np) = Bttbd2
        Vtrant(np) = vttbd2
        Btrtmp = Btrbd2
        Btttmp = Bttbd2
        dnstmp = densit(np)
        pgtmp = prsgas(np)
        pttmp = pressu(np)
        vxtmp = dxdt(np)
        vtrtmp = vtrbd2
        vtttmp = vttbd2
     endif
     sqdnsi = fc1 / sqrt(dnstmp)
         
     sound = sqrt (gamma*pgtmp / dnstmp)
     Alfvn = Bparal(np) * sqdnsi
     sp1 = sound*sound + fc2 * (Bparal(np)*Bparal(np) + Btrtmp*Btrtmp &
          & + Btttmp*Btttmp) / dnstmp
     sp2 = 2.*sound * Alfvn       
     sp3 = sqrt(sp1 + sp2) 
     sp4 = sqrt(sp1 - sp2)
     
     fast = .5 * (sp3 + sp4)
     slow = .5 * (sp3 - sp4)

     chr1 = vxtmp + fast
     chr2 = vxtmp + Alfvn
     chr3 = vxtmp + slow
     chr4 = vxtmp
     chr5 = vxtmp - slow
     chr6 = vxtmp - Alfvn
     chr7 = vxtmp - fast 

     alf2 = Alfvn*Alfvn
     fst2 = fast*fast
     slw2 = slow*slow
     elm1 = fst2 - alf2
     elm2 = slw2 - alf2
!      write(*,*)fast,Alfven(np),slow,sounds(np)

!---------- elements for spherical terms-----------------------
     rad0 = radius(np)
     dr0 = .1*dradiu(np)
     rad1 = radius(np) + dr0
     f0 = fltb(rad0)
     f1 = fltb(rad1)
     geo1 = (rad1*sqrt(f1)/(rad0*sqrt(f0)) - 1.d0) / dr0
     geo2 = (rad1*rad1*f1 / (rad0*rad0*f0) - 1.d0) / dr0
     alcmpy = ( vtrtmp * vxtmp - fc2 * Bparal(np) * Btrtmp / dnstmp ) * geo1 
     alcmpz = ( vtttmp * vxtmp - fc2 * Bparal(np) * Btttmp / dnstmp ) * geo1 
     alelmy = ( Btrtmp * vxtmp - Bparal(np) * vtrtmp ) * geo1
     alelmz = ( Btttmp * vxtmp - Bparal(np) * vtttmp ) * geo1
     div = dnstmp * vxtmp * geo2 * sound * sound 
     trwvel =  ( - dnstmp * ( vtrtmp*vtrtmp + vtttmp*vtttmp ) &
          & + fc2 * (Btrtmp*Btrtmp + Btttmp*Btttmp ) ) * geo1
!----------------------------------------------------------------

     if(chr1 > 0.d0)then
        qL1 = chr1*( dnstmp * elm1 * fast * (vxtmp - dxdt(np-1)) - fc2 &
             & * Bparal(np) * fast * ( Btrtmp * (vtrtmp - vtrans(np-1) ) &
             & + Btttmp  * (vtttmp - vtrant(np-1) ) ) + elm1 &
             & * (pgtmp - prsgas(np-1)) + fc3 * fst2 * (Btrtmp*Btrtmp &
             & + Btttmp*Btttmp - Btrans(np-1)*Btrans(np-1) &
             & - Btrant(np-1)*Btrant(np-1) ) ) / dradiu(np-1) + elm1 * fast &
             & * dnstmp * grav / (radius(np)*radius(np-1))
!     &           + trwvel * fast * elm1 - fc2 * Bparal(np) * fast 
!     &           * ( alcmpy * Btrtmp + alcmpz * Btttmp )
!     &           + fst2 * fc2 * (alelmy * Btrtmp + alelmz * Btttmp)
!     &           + div * elm1
!            qL1 = chr1*( dnstmp * elm1 * fast * (vxtmp - dxdt(np-1)) 
!     &           - fc2 * Bparal(np) * fast * ( Btrtmp  
!     &           * (vtrtmp - vtrans(np-1) ) + Btttmp  * (vtttmp 
!     &           - vtrant(np-1) ) ) + elm1 * (pttmp - pressu(np-1))  
!     &           + fc3 * alf2 * (Btrtmp*Btrtmp 
!     &           + Btttmp*Btttmp - Btrans(np-1)*Btrans(np-1) 
!     &           -Btrant(np-1)*Btrant(np-1) ) ) / dradiu(np-1)
!            qL1 = chr1 * ( densit(np) * elm1 * fast * (vxout2 - vxout1) 
!     &           - fc2 * Bparal(np) * fast * ( .5 * (Btrbd1 
!     &           + Btrbd2 ) * (vtrbd2 - vtrbd1 ) 
!     &           + .5 * ( Bttbd2 + Bttbd2 ) * (vttbd2 
!     &           - vttbd1 ) ) + elm1 * (ptout2 - ptout1)  
!     &           + fc3 * alf2 * (Btrbd2*Btrbd2 
!     &           + Bttbd2*Bttbd2 - Btrbd1*Btrbd1 
!     &           -Bttbd1*Bttbd1 ) ) / dradiu(np-1)
!         write(*,*)'qL1',iter,qL1,elm1 
!,dnstmp,vxtmp,Btrtmp,vtrtmp,pttmp
!     &           ,dnstmp * elm1 * fast * (vxtmp - dxdt(np-1))
!     &           ,- fc2 * Bparal(np) * fast * ( Btrtmp  
!     &           * (vtrtmp - vtrans(np-1) ) + Btttmp  * (vtttmp 
!     &           - vtrant(np-1) ) ), elm1 * (pttmp - pressu(np-1)) 
!     &           ,fc3 * alf2 * (Btrtmp*Btrtmp 
!     &           + Btttmp*Btttmp - Btrans(np-1)*Btrans(np-1) 
!     &           -Btrant(np-1)*Btrant(np-1) ) 
!         write(*,*)dxdt(np-2),vxout1,dxdt(np-1),vxout2,dxdt(np)
!     &        ,pressu(np-2),ptout1,pressu(np-1),ptout2,pressu(np)
     else
        qL1 = 0.d0 
     endif
     
     if(chr2 > 0.d0)then
        qL2 = chr2*( Btttmp  * (vtrtmp - vtrans(np-1) ) - Btrtmp &
             & * (vtttmp - vtrant(np-1) ) - sqdnsi * ( Btttmp  &
             & * (Btrtmp - Btrans(np-1) ) - Btrtmp * (Btttmp - Btrant(np-1))))& 
             & / dradiu(np-1)
!     &           + alcmpy * Btttmp - alcmpz * Btrtmp 
!     &           - (alelmy * Btttmp - alelmz * Btrtmp ) * sqdnsi 
!            qL2 = chr2 * ( .5 * ( Bttbd1 + Bttbd2 ) 
!     &           * (vtrbd2 - vtrbd1 ) - .5 * (Btrbd2 
!     &           + Btrbd2 ) * (vttbd2 - vttbd1 ) 
!     &           - sqdnsi * ( .5 * ( Bttbd1 + Bttbd2 ) 
!     &           * (Btrbd2 - Btrbd1 ) - .5 * ( Btrbd1 
!     &           + Btrbd2 )  * (Bttbd2 - Bttbd1) ) ) 
!     &           / dradiu(np-1) 
     else
        qL2 = 0.d0
     endif

     if(chr3 > 0.d0)then
        qL3 = chr3*( dnstmp * elm2 * slow * (vxtmp - dxdt(np-1) ) - fc2 &
             & * Bparal(np) * slow * (  Btrtmp * (vtrtmp - vtrans(np-1) ) &
             & + Btttmp * (vtttmp - vtrant(np-1) ) ) + elm2 &
             & * (pgtmp - prsgas(np-1) ) + fc3 * slw2 * (Btrtmp*Btrtmp &
             & + Btttmp*Btttmp - Btrans(np-1)*Btrans(np-1) &
             &  - Btrant(np-1)*Btrant(np-1) ) ) / dradiu(np-1) &
             & + elm2 * slow * dnstmp * grav / (radius(np)*radius(np-1))
!     &           + trwvel * slow * elm2 - fc2 * Bparal(np) * slow 
!     &           * ( alcmpy * Btrtmp + alcmpz * Btttmp )
!     &           + slw2 * fc2 * (alelmy * Btrtmp + alelmz * Btttmp)
!     &           + div * elm2
!            qL3 = chr3*( dnstmp * elm2 * slow * (vxtmp - dxdt(np-1) ) 
!     &           - fc2 * Bparal(np) * slow * (  Btrtmp 
!     &           * (vtrtmp - vtrans(np-1) ) + Btttmp * (vtttmp 
!     &           - vtrant(np-1) ) ) + elm2 * (pttmp - pressu(np-1) )  
!     &           + fc3 * alf2 * (Btrtmp*Btrtmp 
!     &           + Btttmp*Btttmp - Btrans(np-1)*Btrans(np-1) 
!     &           -Btrant(np-1)*Btrant(np-1) ) ) / dradiu(np-1)
!            qL3 = chr3 * ( densit(np) * elm2 * slow * (vxout2 - vxout1) 
!     &           - fc2 * Bparal(np) * slow * ( .5 * ( Btrbd1 
!     &           + Btrbd2 ) * (vtrbd2 - vtrbd1 ) 
!     &           + .5 * ( Bttbd1 + Bttbd2 ) * (vttbd2 
!     &           - vttbd1 ) ) + elm2 * (ptout2 - ptout1)  
!     &           + fc3 * alf2 * (Btrbd2*Btrbd2 
!     &           + Bttbd2*Bttbd2 - Btrbd1*Btrbd1 
!     &           -Bttbd1*Bttbd1 ) ) / dradiu(np-1)
!     write(*,*)'qL3',qL3,densit(np) * elm2 * slow *(vxout2 - vxout1)
!     &        ,- fc2 * Bparal(np) * slow * .5 * ( Btrbd1 + Btrbd2 ) 
!     &        * (vtrbd2 - vtrbd1 ), elm2 * (ptout2 - ptout1)
!     &        ,fc3 * alf2 * (Btrbd2*Btrbd2 + Bttbd2*Bttbd2 
!     &        - Btrbd1*Btrbd1 - Bttbd1*Bttbd1 ), chr3 
     else
        qL3 = 0.d0 
     endif
         
     if(chr4 > 0.d0) then
        qL4 = chr4*( (pgtmp-prsgas(np-1)) - sound*sound &
             & * (dnstmp - densit(np-1) ) ) / dradiu(np-1)
!            qL4 = chr4 * ( (pgout2-pgout1) - sounds(np)*sounds(np) 
!     &           * (dnout2 - dnout1) ) / dradiu(np-1)
     else
        qL4 = 0.
     endif
     
     if(chr5 > 0.d0)then
        qL5 = chr5*( dnstmp * elm2 * slow * (vxtmp - dxdt(np-1) ) - fc2 &
             & * Bparal(np) * slow * (  Btrtmp * (vtrtmp - vtrans(np-1) ) &
             & + Btttmp * (vtttmp - vtrant(np-1) ) ) - elm2 &
             & * (pgtmp - prsgas(np-1) ) - fc3 * slw2 * (Btrtmp*Btrtmp &
             & + Btttmp*Btttmp - Btrans(np-1)*Btrans(np-1) &
             & - Btrant(np-1)*Btrant(np-1) ) ) / dradiu(np-1) + elm2 &
             & * slow * dnstmp * grav / (radius(np)*radius(np-1))
!     &           + trwvel * slow * elm2 - fc2 * Bparal(np) * slow 
!     &           * ( alcmpy * Btrtmp + alcmpz * Btttmp )
!     &           - slw2 * fc2 * (alelmy * Btrtmp + alelmz * Btttmp)
!     &           - div * elm2
!            qL5 = chr5*( dnstmp * elm2 * slow * (vxtmp - dxdt(np-1) ) 
!     &           - fc2 * Bparal(np) * slow * (  Btrtmp 
!     &           * (vtrtmp - vtrans(np-1) ) + Btttmp * (vtttmp 
!     &           - vtrant(np-1) ) ) - elm2 * (pttmp - pressu(np-1) )  
!     &           - fc3 * alf2 * (Btrtmp*Btrtmp 
!     &           + Btttmp*Btttmp - Btrans(np-1)*Btrans(np-1) 
!     &           -Btrant(np-1)*Btrant(np-1) ) ) / dradiu(np-1)
!            qL5 = chr5 * ( densit(np) * elm2 * slow * (vxout2 - vxout1) 
!     &           - fc2 * Bparal(np) * slow * ( .5 * ( Btrbd1 + 
!     &           Btrbd2 ) * (vtrbd2 - vtrbd1 ) 
!     &           + .5 * ( Bttbd1 + Bttbd2 ) * (vttbd2 
!     &           - vttbd1 ) ) - elm2 * (ptout2 - ptout1)  
!     &           - fc3 * alf2 * (Btrbd2*Btrbd2 
!     &           + Bttbd2*Bttbd2 - Btrbd1*Btrbd1 
!     &           -Bttbd1*Bttbd1 ) ) / dradiu(np-1)
!         write(*,*)'qL5',qL5,densit(np) * elm2 * slow *(vxout2 - vxout1)
!     &        ,- fc2 * Bparal(np) * slow * .5 * ( Btrbd1 + Btrbd2 ) 
!     &        * (vtrbd2 - vtrbd1 ), -elm2 * (ptout2 - ptout1)
!     &        ,fc3 * alf2 * (Btrbd2*Btrbd2 + Bttbd2*Bttbd2 
!     &        - Btrbd1*Btrbd1 - Bttbd1*Bttbd1 ),chr5  
     else
        qL5 = 0.d0 
     endif

     if(chr6 > 0.d0)then
        qL6 = chr6*( Btttmp  * (vtrtmp - vtrans(np-1) ) - Btrtmp &
             & * (vtttmp - vtrant(np-1) ) + sqdnsi * ( Btttmp &
             & * (Btrtmp - Btrans(np-1) ) - Btrtmp * (Btttmp - Btrant(np-1))))& 
             & / dradiu(np-1) 
!     &           + alcmpy * Btttmp - alcmpz * Btrtmp 
!     &           + (alelmy * Btttmp - alelmz * Btrtmp ) * sqdnsi 
!            qL6 = chr6 * ( .5 * ( Bttbd1 + Bttbd2 ) 
!     &           * (vtrbd2 - vtrbd1 ) - .5 * ( Btrbd1 
!     &           + Btrbd2 ) * (vttbd2 - vttbd1 ) 
!     &           + sqdnsi * ( .5 * ( Bttbd1 + Bttbd2 ) 
!     &           * (Btrbd2 - Btrbd1) - .5 * ( Btrbd1 
!     &           + Btrbd2 ) * (Bttbd2 - Bttbd1) ) ) 
!     &           / dradiu(np-1) 
     else
        qL6 = 0.d0
     endif

     if(chr7 > 0.d0)then
        qL7 = chr7*( dnstmp * elm1 * fast * (vxtmp - dxdt(np-1)) - fc2 &
             & * Bparal(np) * fast * ( Btrtmp * (vtrtmp - vtrans(np-1) ) &
             & + Btttmp  * (vtttmp - vtrant(np-1) ) ) - elm1 &
             & * (pgtmp - prsgas(np-1)) - fc3 * fst2 * (Btrtmp*Btrtmp &
             & + Btttmp*Btttmp - Btrans(np-1)*Btrans(np-1) &
             & - Btrant(np-1)*Btrant(np-1) ) ) / dradiu(np-1) + elm1 &
             & * fast * dnstmp * grav / (radius(np)*radius(np-1))
!     &           + trwvel * fast * elm1 - fc2 * Bparal(np) * fast 
!     &           * ( alcmpy * Btrtmp + alcmpz * Btttmp )
!     &           - fst2 * fc2 * (alelmy * Btrtmp + alelmz * Btttmp)
!     &           - div * elm1
!            qL7 = chr7*( dnstmp * elm1 * fast * (vxtmp - dxdt(np-1)) 
!     &           - fc2 * Bparal(np) * fast * ( Btrtmp  
!     &           * (vtrtmp - vtrans(np-1) ) + Btttmp  * (vtttmp 
!     &           - vtrant(np-1) ) ) - elm1 * (pttmp - pressu(np-1))  
!     &           - fc3 * alf2 * (Btrtmp*Btrtmp 
!     &           + Btttmp*Btttmp - Btrans(np-1)*Btrans(np-1) 
!     &           -Btrant(np-1)*Btrant(np-1) ) ) / dradiu(np-1)
!            qL7 = chr7 * ( densit(np) * elm1 * fast * (vxout2 - vxout1) 
!     &           - fc2 * Bparal(np) * fast * ( .5 * ( Btrbd1 
!     &           + Btrbd2 ) * (vtrbd2 - vtrbd1 ) 
!     &           + .5 * ( Bttbd1 + Bttbd2 ) * (vttbd2 
!     &           - vttbd1 ) ) - elm1 * (ptout2 - ptout1)  
!     &           - fc3 * alf2 * (Btrbd2*Btrbd2 
!     &           + Bttbd2*Bttbd2 - Btrbd1*Btrbd1 
!     &           - Bttbd1*Bttbd1 ) ) / dradiu(np-1)
     else
        qL7 = 0.d0 
     endif
         
!      write(*,*)qL1,qL2,qL3,qL4,qL5,qL6,qL7,pgout1,pgout2,vxout1,vxout2
!     &     ,Btrbnd(np-2),Btrbnd(np-1),vtrbnd(np-2),vtrbnd(np-1)
!     &     ,dxdt(np-2),dxdt(np-1),dxdt(np)

!         rad0 = radius(np)
!         dr0 = .1*dradiu(np)
!         rad1 = radius(np) + dr0
!         f0 = fltb(rad0)
!         f1 = fltb(rad1)
!         write(*,*)f0,f1,vtrtmp,Btrtmp
!         geo1 = (rad1*sqrt(f1)/(rad0*sqrt(f0)) - 1.) / dr0
!         geo2 = (rad1*rad1*f1 / (rad0*rad0*f0) - 1.) / dr0
!         alcmpy = ( vtrtmp * vxtmp - fc2 * Bparal(np) * Btrtmp  
!     &        / dnstmp ) / rad0
!         alcmpz = ( vtttmp * vxtmp - fc2 * Bparal(np) * Btttmp  
!     &        / dnstmp ) / rad0      
!         alelmy = ( Btrtmp * vxtmp - Bparal(np) * vtrtmp ) 
!     &        * geo1
!         alelmz = ( Btttmp * vxtmp - Bparal(np) * vtttmp )
!     &        * geo1
!         div = dnstmp * vxtmp * geo2 * sound * sound 
!         trwvel =  ( - dnstmp * ( vtrtmp*vtrtmp + vtttmp*vtttmp )
!     &        + fc2 * (Btrtmp*Btrtmp + Btttmp*Btttmp ) ) / rad0

     qLcf = .5d0 * (qL1 + qL7) !+ elm1 * fast * ( trwvel
!     &        + dnstmp * grav / (rad0*rad0) !grav term is included in qL
!     &        )- fc2 * Bparal(np) * fast 
!     &        * ( Btrtmp * alcmpy + Btttmp * alcmpz )
     qLcs = .5d0 * (qL3 + qL5) !+ elm2 * slow * ( trwvel
!     &        + dnstmp * grav / (rad0*rad0) 
!     &        )- fc2 * Bparal(np) * slow 
!     &        * ( Btrtmp * alcmpy + Btttmp * alcmpz )
!      write(*,*)qLcf,qL1,.5*dnstmp*elm1*fast*grav/(rad0*rad0)
!                                ! ,rad0,-densit(np) 
!     &     * ( vtrans(np)*vtrans(np) + vtrant(np)*vtrant(np) )
!     &     , fc2 * ( Btrans(np)*Btrans(np) + Btrant(np)*Btrant(np) )
!     &     ,grav / (rad0*rad0) , - fc2 * Bparal(np) * slow 
!     &     * ( Btrans(np) * alcmpy + Btrant(np) * alcmpz )
     qLca = .5d0 * (qL2 + qL6) !+ Btttmp * alcmpy 
!     &     + Btrtmp * alcmpz
     qLmf = .5d0 * (qL1 - qL7) !+ fst2 * ( Btrtmp * alelmy 
!     &     + Btttmp * alelmz ) + div * elm1
     qLms = .5d0 * (qL3 - qL5) !+ slw2 * ( Btrtmp * alelmy 
!     &     + Btttmp * alelmz ) + div * elm2
     qLma = .5 * (qL2 - qL6) !+ sqdnsi * ( - Btttmp * alelmy 
!     &     + Btrtmp * alelmz )
      
     elm3 = fst2 - slw2 
     elm4 = Alfvn*sound*elm3
!      elm5 = Btrans(np)*Btrans(np) + Btrant(np)*Btrant(np)
     elm5 = Btrtmp*Btrtmp + Btttmp*Btttmp
     elm6 = elm4 * elm5 * Bparal(np)
     elm7 = elm3 * alf2
     elm8 = elm7 * elm5

     dPtdt = (-qLmf + qLms) / elm3  
!         if(iter.gt.1)then
!            if(dPtdt1 .ne.0.)then
!               err = err + abs(dPtdt/dPtdt1 - 1.)
!            endif
!         endif   
     if(Alfvn<=0.d0)then
        dvxdt = - qLcf / (dnstmp * fst2 * fast)
        Bymag = fc3 * Btrtmp * Btrtmp
        Bzmag = fc3 * Btttmp * Btttmp
        Ptot = pgtmp + Bymag + Bzmag
        dPgdt = dPtdt * pgtmp / Ptot
        if(Btrtmp/=0.d0)then
           dBtrdt = dPtdt * Bymag / Ptot / (Btrtmp * fc2) 
           if(iter>1)then
              if(dBtrdt1 /=0.d0)then
                 err = err + abs(dBtrdt/dBtrdt1 - 1.d0)
              endif
           endif
        else
           dBtrdt = 0.d0
        endif
        if(Btttmp/=0.d0)then
           dBttdt = dPtdt * Bzmag / Ptot / (Btttmp * fc2) 
           if(iter>1)then
              if(dBttdt1 /=0.d0)then
                 err = err + abs(dBttdt/dBttdt1 - 1.d0)
              endif
           endif
        else
           dBttdt = 0.d0
        endif
        dvtrdt = 0.d0
        dvttdt = 0.d0
!     write(*,*)Ptot,prsgas(np),dPtdt,dPgdt,elm3,qL1,qLmf
     else
        dvxdt = (-slow * qLcf + fast * qLcs) / (elm4 * dnstmp )
        dPgdt = ( - slw2*qLmf + fst2*qLms ) / elm7
        if(iter>1)then
           if(dvxdt1 /=0.d0)then
              err = err + abs(dvxdt/dvxdt1 - 1.d0)
           endif
!               if(dPgdt1 .ne.0.)then
!                  err = err + abs(dPgdt/dPgdt1 - 1.)
!               endif
        endif
        if(elm5<=1.d-15)then
           dvtrdt = 0.d0 
           dvttdt = 0.d0
           dBtrdt = 0.d0
           dBttdt = 0.d0
        else
!     dvtrdt = ( Btrans(np) * ( - slow * elm2 * qLcf 
!     &           + fast * elm1 * qLcs) - elm4 * fc2 * Bparal(np) 
!     &           * Btrant(np) * qLca ) / (fc2 * elm6 )
           dvtrdt = ( Btrtmp * ( - slow * elm2 * qLcf + fast * elm1 * qLcs) &
                & - elm4 * fc2 * Bparal(np) * Btttmp * qLca ) / (fc2 * elm6 )
           dvttdt = ( Btttmp * ( - slow * elm2 * qLcf + fast * elm1 * qLcs) &
                & + elm4 * fc2 * Bparal(np) * Btrtmp * qLca ) / (fc2 * elm6 )
               
!     dBtrdt = (  Btrans(np) * sqdnsi * ( elm2*qLmf - elm1*qLms) 
!     &           + elm7 * fc2 * Btrant(np) * qLma) 
!     &           / (fc2 * elm8 * sqdnsi)
           dBtrdt = (  Btrtmp * sqdnsi * ( elm2*qLmf - elm1*qLms) &
                & + elm7 * fc2 * Btttmp * qLma) / (fc2 * elm8 * sqdnsi)
               
           dBttdt = ( Btttmp * sqdnsi * ( elm2*qLmf - elm1*qLms) &
                & - elm7 * fc2 * Btrtmp * qLma) / (fc2 * elm8 * sqdnsi)      
!               write(*,*)iter,dBtrdt,dvtrdt,Btrtmp,qLmf,qLms,fast,slow
!     &              ,fast * elm1 * qLcs, elm6,- elm4 * fc2 * Bparal(np) 
!     &              * Bttbd2 * qLca 
!               write(*,*)dBttdt,Btttmp,qLma
           if(iter > 1)then
              if(dvtrdt1 /=0.)then
                 err = err + abs(dvtrdt/dvtrdt1 - 1.d0)
!                     write(*,*)'dvy',err,dvtrdt,dvtrdt1
              endif
              if(dvttdt1 /=0.d0)then
                 err = err + abs(dvttdt/dvttdt1 - 1.d0)
              endif
              if(dBtrdt1 /=0.d0)then
                 err = err + abs(dBtrdt/dBtrdt1 - 1.d0)
              endif
              if(dBttdt1 /=0.d0)then
                 err = err + abs(dBttdt/dBttdt1 - 1.d0)
              endif
           endif
        endif
     endif
     
     dPmdt = fc2*(Btrtmp*dBtrdt + Btttmp*dBttdt)       
     dPmdtdf = dPtdt - dPgdt
!         dPgdtdf = dPtdt - dPmdt
!         if(dPmdt .ne. 0.)then
!            dBtrdt = dBtrdt * dPmdtdf / dPmdt
!            dBttdt = dBttdt * dPmdtdf / dPmdt
!         endif
!         dPgdt = dPtdt - fc2 * (Btrans(np)*dBtrdt + Btrant(np)*dBttdt)
     ddndt = (dPgdt + qL4) / (sound*sound)
     if(iter >1)then
        if(ddndt1 /= 0.)then
           err = err + abs(ddndt/ddndt1 - 1.d0)
        endif
     endif

!      write(*,*)iter,dBtrdt,dBttdt,dvxdt,dvtrdt,dvttdt,dPgdt,dPtdt,ddndt
!     &        ,qL1,qL3,qL5,qL7,qLcf,qLcs
!     &        ,Btrtmp,vtrtmp,pgtmp,dnstmp,vxtmp
!     &     ,qL1,qL2,qL3,qL4,qL5,qL6,qL7
!      write(*,*)dvxdt,qL1,qLcf,qL3,qL5,qLcs,slow,fast,densit(np),elm4
!     &     ,-slow * qLcf, + fast * qLcs
!      write(*,*)fast,slow,Alfven(np),sounds(np),qL1,qL2,qL3,qL4
!     &     ,-slow*elm2*qLcf,fast*elm1*qLcs,dvtrdt

!     if (err < 1.d-6 .and. iter > 1) goto 102
     if (err < 1.d-6 .and. iter > 1) exit

     if(iter > 1)then
        dBtrdt = .5d0*(dBtrdt + dBtrdt1)
        dBttdt = .5d0*(dBttdt + dBttdt1)
        ddndt = .5d0*(ddndt + ddndt1)
        dPgdt = .5d0*(dPgdt + dPgdt1)
        dvxdt = .5d0*(dvxdt + dvxdt1)           
!            dPtdt = .5*(dPtdt + dPtdt1)           
        dvtrdt = .5d0*(dvtrdt + dvtrdt1)
        dvttdt = .5d0*(dvttdt + dvttdt1)
     endif

     Btrtmp = Btrans(np) + dt*dBtrdt 
     Btttmp = Btrant(np) + dt*dBttdt
     dnstmp = densit(np) + dt*ddndt
     pgtmp = prsgas(np) + dt*dPgdt
     vxtmp = dxdt(np) + dt*dvxdt
!         pttmp = pressu(np) + dt*dPtdt
     vtrtmp = vtrans(np) + dt*dvtrdt
     vtttmp = vtrant(np) + dt*dvttdt

!         write(*,*)iter,dBtrdt,Btrtmp,vtrtmp,dnstmp,pgtmp,vxtmp,pttmp

!         dPtdt1 = dPtdt
     dvxdt1 = dvxdt
     dPgdt1 = dPgdt
     dvtrdt1 = dvtrdt
     dvttdt1 = dvttdt
     dBtrdt1 = dBtrdt
     dBttdt1 = dBttdt
     ddndt1 = ddndt

  end do

!102 continue
  dxdt(np) = dxdt(np) + dvxdt*dt
  densit(np) = densit(np) + ddndt*dt
  prsgas(np) = prsgas(np) + dPgdt*dt
!      pressu(np) = pressu(np) + dPtdt*dt
  Btrans(np) = Btrans(np) + dBtrdt*dt
  vtrans(np) = vtrans(np) + dvtrdt*dt
  Btrant(np) = Btrant(np) + dBttdt*dt
  vtrant(np) = vtrant(np) + dvttdt*dt      


  Bprs = fc3 * (Btrans(np)*Btrans(np) + Btrant(np)*Btrant(np))
!      Bprdf = pressu(np) - prsgas(np)
!      fr1 = Bprs / pressu(np)
!      fr2 = Bprdf / pressu(np)
!      if(fr1.lt.10.3 .and. fr2.lt.10.3)goto 103
!      rdfc = sqrt(max(Bprdf,1.e-20)/max(Bprs,1.e-20))
!      write(*,*)Bprs,Bprdf,dpgdt,dPtdt,dBtrdt
!      Btrans(np) = Btrans(np) * rdfc
!      Btrant(np) = Btrant(np) * rdfc
!      prsgas(np) = pressu(np) 
!     &     - fc3 * (Btrans(np)*Btrans(np) + Btrant(np)*Btrant(np))
! 103  continue
  if(prsgas(np) < 0.d0)then
     write(*,*)'Pg < 0 in outgoing bd' 
     prsgas(np) = 1.d-10
  endif
!      if(pressu(np).lt.0.)then
!         write(*,*)'Pt < 0 in outgoing bd' 
!         pressu(np) = prsgas(np) + Bprdf 
!      endif
  pressu(np) = prsgas(np) + Bprs
  speene(np) = prsgas(np) / (gammi1 * densit(np))
  energy(np) = speene(np) + .5*(dxdt(np)*dxdt(np) + vtrans(np)*vtrans(np) &
       & + vtrant(np)*vtrant(np) ) + Bprs /densit(np) - grav / radius(np)
  sounds(np) = sqrt( gamma*gammi1*speene(np) )
  Alfven(np) = fc1*Bparal(np)/sqrt(densit(np))

!      write(*,*)prsgas(np),pressu(np),densit(np),Alfven(np),sounds(np)
!     &     ,speene(np),energy(np),Btrans(np),vtrans(np),Btrant(np)
      
end subroutine outgobdsc

!------------------------------------------------------------------
subroutine wveng
  use nrtype
  USE vars
  implicit none

! This subroutine calculates energy flux of various waves
 
  real(DP),parameter :: dncgs=-7.d0, vkm=6.4d0
  real(DP) :: vlav(nd),dnav(nd)
  real(DP) :: fcb,r1,r2,dnt,vlt,wmd,spmin,smmax,velabs,sl,wdm,dnln,ddvtr &
       & ,dvtin,dvtout,dvlg,ddns,ddvlg,dvlout,dvlin,geo,fltb &
       & ,ddvtr3,dvtout3,dvtin3
  integer(I4B) :: i,ip

  fcb = sqrt(4.*pi*10.d0**dncgs)*1.d5*vkm
  
  do ip = 1,np
     r1 = radius(ip) - 2.d0
     r2 = radius(ip) + 2.d0
     dnt = 0.d0
     vlt = 0.d0
     wmd=0.d0
     do i=1,np
        if(radius(i) > r1 .and. radius(i) < r2)then
           spmin = min(speene(i),speene(ip))
           if(spmin > 50.d0)then
              smmax=2.d0
           else
              smmax=(.04d0*spmin)**4.
           endif
           smmax = min(smmax,max(.5d0*(radius(ip)-1.2d0),1.d-4))
           velabs = (Alfven(i) + abs(dxdt(i)) )*vkm
           sl = min(velabs*1800.d0/7.d5,smmax) ! smoothing length
           wdm = exp(-((radius(ip) - radius(i))/sl)**2)
           dnt = dnt+(log10(densit(i))+dncgs)*wdm
           vlt = vlt+dxdt(i)*vkm*wdm
           wmd = wmd+wdm
        endif
     enddo
     dnav(ip) = dnt/wmd
     vlav(ip) = vlt/wmd  ! outflow speed &
  enddo
  
  do ip=1,np
     dnln = 10.d0**dnav(ip)
     ddvtr = - fcb*Btrans(ip) / sqrt(4.d0*pi*dnln) * 1.d-5
     dvtout = .5d0*(vtrans(ip)*vkm + ddvtr)
     dvtin = vtrans(ip)*vkm - dvtout
     
     ddvtr3 = - fcb*Btrant(ip) / sqrt(4.d0*pi*dnln) * 1.d-5
     dvtout3 = .5d0*(vtrant(ip)*vkm + ddvtr3)
     dvtin3 = vtrant(ip)*vkm - dvtout3

     dvlg = dxdt(ip)*vkm - vlav(ip)
     ddns = densit(ip)*10.d0**dncgs - dnln
     ddvlg = ddns/dnln *sounds(ip)*vkm 
     dvlout = .5d0*(dvlg + ddvlg)
     dvlin = dvlg -dvlout
     
     geo = radius(ip)*radius(ip)*fltb(radius(ip))
     
     Falout(ip) = dnln*dvtout*dvtout * (Alfven(ip)*vkm + 1.5d0*vlav(ip)) &
          & * 1.d15*geo !/15.d0
     Falin(ip) = dnln*dvtin*dvtin * (Alfven(ip)*vkm - 1.5d0*vlav(ip)) &
          & * 1.e15*geo !/15.d0
     Falout3(ip) = dnln*dvtout3*dvtout3 * (Alfven(ip)*vkm + 1.5d0*vlav(ip)) &
          & * 1.d15*geo !/15.d0
     Falin3(ip) = dnln*dvtin3*dvtin3 * (Alfven(ip)*vkm - 1.5d0*vlav(ip)) &
          & * 1.e15*geo !/15.d0
     Fcsout(ip) = dnln*dvlout*dvlout * (sounds(ip)*vkm + (1.5d0+0.5d0*gamma) &
          & * vlav(ip) )* 1.d15*geo !/15.d0
  enddo

  
end subroutine wveng

!-----------------------------------------------------------------------
subroutine printvars(itime,index)
  use nrtype
  USE vars
  implicit none   

  INTEGER(I4B),intent(in) :: itime,index
  INTEGER(I4B) :: ip
  REAL(DP) :: Temp,mu,xe
  
  if(index==0)goto 99


  dxdtav(1:np) = dxdtav(1:np) / densav(1:np)
  Btraav(1:np) = sqrt(Btraav(1:np)/dtavtt)
  vtraav(1:np) = sqrt(vtraav(1:np)/densav(1:np))
  Btrtav(1:np) = sqrt(Btrtav(1:np)/dtavtt)
  vtrtav(1:np) = sqrt(vtrtav(1:np)/densav(1:np))
  spenav(1:np) = spenav(1:np) / densav(1:np)
  densav(1:np) = densav(1:np) / dtavtt
  dvbaav(1:np) = dvbaav(1:np) / (dtavtt*sqrt(densav(1:np)))
  dvbtav(1:np) = dvbtav(1:np) / (dtavtt*sqrt(densav(1:np)))

  do ip = 1,np
     call tp_mu(speene(ip),Temp,mu)
     if(t>t_e2T_rss.and.Temp<2.d4)call e2T(ip,Temp,mu,xe)
     Tpprt(ip) = Temp
     muprt(ip) = real(mu)
  end do
  
!  write(12,rec=itime)t,radius,dradiu,densit,dxdt,vtrans,vtrant &
!       & ,Btrans,Btrant,speene(1:nd),dxdtbd,real(Alfven),etaprt,eadprt,xeprt &
!       & ,real(densav),real(dxdtav),real(vtraav),real(vtrtav),real(Btraav) &
!       & ,real(Btrtav),real(spenav),real(dvbaav),real(dvbtav)
  write(12,rec=itime)t,radius,dradiu,densit,dxdt,vtrans,vtrant &
       & ,Btrans,Btrant,speene(1:nd),Tpprt,real(Alfven),etaprt,eadprt,muprt &
       & ,real(densav),real(dxdtav),real(vtraav),real(vtrtav),real(Btraav) &
       & ,real(Btrtav),real(spenav),real(dvbaav),real(dvbtav)
  call flush(12)

99 continue 
! initialization of averaged values
  dtavtt=0.d0
  dxdtav(1:np)=0.d0
  densav(1:np)=0.d0
  enrgav(1:np)=0.d0
  prgsav(1:np)=0.d0
  spenav(1:np)=0.d0
  Btraav(1:np)=0.d0
  vtraav(1:np)=0.d0
  Btrtav(1:np)=0.d0
  vtrtav(1:np)=0.d0
  dvbaav(1:np)=0.d0
  dvbtav(1:np)=0.d0
  zplus(1:np)=0.d0
  zmins(1:np)=0.d0

end subroutine printvars

!-----------------------------------------------------------------------
subroutine tavevars
  use nrtype
  USE vars
  implicit none   

  INTEGER(I4B) :: ip
  REAL(DP) :: sqdav

  dtavtt = dtavtt + dt

  do ip = 1,np
!     dxdtav(ip) = dxdtav(ip) + dt*dxdt(ip)
     dxdtav(ip) = dxdtav(ip) + dt*dxdt(ip)*densit(ip)
     densav(ip) = densav(ip) + dt*densit(ip)
!     enrgav(ip) = enrgav(ip) + dt*energy(ip)
     enrgav(ip) = enrgav(ip) + dt*energy(ip)*densit(ip)
     prgsav(ip) = prgsav(ip) + dt*prsgas(ip)
!     spenav(ip) = spenav(ip) + dt*speene(ip)
     spenav(ip) = spenav(ip) + dt*speene(ip)*densit(ip)
     Btraav(ip) = Btraav(ip) + dt*Btrans(ip)*Btrans(ip)
!     vtraav(ip) = vtraav(ip) + dt*vtrans(ip)*vtrans(ip)
     vtraav(ip) = vtraav(ip) + dt*vtrans(ip)*vtrans(ip)*densit(ip)
     Btrtav(ip) = Btrtav(ip) + dt*Btrant(ip)*Btrant(ip)
!     vtrtav(ip) = vtrtav(ip) + dt*vtrant(ip)*vtrant(ip)
     vtrtav(ip) = vtrtav(ip) + dt*vtrant(ip)*vtrant(ip)*densit(ip)
     sqdav = sqrt(densit(ip))
     dvbaav(ip) = dvbaav(ip) + dt*vtrans(ip)*btrans(ip)*sqdav
     dvbtav(ip) = dvbtav(ip) + dt*vtrant(ip)*btrant(ip)*sqdav
     zplus(ip) = zplus(ip) + dt*(vtrans(ip) - Btrans(ip)/sqdav) 
     zmins(ip) = zmins(ip) + dt*(vtrans(ip) + Btrans(ip)/sqdav) 
  end do

end subroutine tavevars
        
!----------------------------------------------------------------
subroutine resistivity
  use nrtype
  USE vars
  implicit none   

  real(DP) :: Tp,mu,xe,fltb
  integer(I4B) :: ip

  do ip = 1,np
     call tp_mu(speene(ip),Tp,mu)
     if(t>t_e2T_rss.and.Tp<2.d4)then      
        call e2T(ip,Tp,mu,xe)
     else
        xe = max(1.d0/(mu*XHyd) - (1.d0 + 0.25d0*YoX + 0.05d0*zmetal*Zsun/XHyd) &
             & ,1.d-10)
     end if
     
!     eta(ip) = 234.d0 * sqrt(Tp) / (xe *rdkm*akm) * 1.d-10
     eta(ip) = min(etafct * sqrt(Tp) / xe, etamax)
!     if(mod(ip,50)==1)write(*,*)ip,real(radius(ip)),real(Tp),real(xe) &
!          & ,real(c_h),real(eta(ip))
     etaprt(ip) = real(etafct * sqrt(Tp) / xe) ! print out
!     xeprt(ip) = real(xe)
     if(irss==1)then
        ead(ip) = 0.d0
        eadprt(ip) = 0.0
     else        
        ead(ip) = eadfct / (densit(ip)*densit(ip)*xe) &
             & *(max(0.d0,(1.d0 - xe)))**2.d0 !Khomenko & Collados 2012 eq.(20)
!             & * 0.5d0*(tanh(1.d2*(0.1d0-xe))+1.d0) ! ead => 0 for xe=>1
!        if(ip==62)write(*,*)ip,Tp,xe
        if(irss==2)then
           eadprt(ip) = real(eadfct / (densit(ip)*densit(ip)*xe) &
                & *(Btrans(ip)*Btrans(ip) + Btrant(ip)*Btrant(ip) ) )
        else
           eadprt(ip) = real(eadfct / (densit(ip)*densit(ip)*xe) &
                & *(Btrans(ip)*Btrans(ip) + Btrant(ip)*Btrant(ip) &
                & + Bplint*Bplint/(radius(ip)**4*fltb(radius(ip))**2) ))
        end if
        eadprt(ip) = eadprt(ip) * (max(0.0,(1.0 - xe)))**2.0 ! print out
     end if
!     if(ip==1000)write(*,*)"xe",radius(ip),Tp,xe
  end do
!  stop
end subroutine resistivity

!--------------------------------------------------------------------
subroutine tp_mu(spvar,Temperature,mu)
  use vars
  implicit none

  real(dp),intent(in) :: spvar
  real(dp),intent(out) :: Temperature,mu

  Temperature = spvar*spetkl
!  write(*,*)Temperature
!  if(Temperature < 2.d3)then
  if(Temperature < Tlow*mufrac)then
!     Temperature = 2.d0*Temperature
     Temperature = Temperature / mufrac
!  else if(Temperature < 1.d4)then
  else if(Temperature < Thigh)then
!     Temperature = .75d0*Temperature + 2500.d0
     Temperature = coefT1*Temperature + coefT0
  endif
  
  mu = Temperature*mufull / (spvar*spetkl)
!  write(*,*)spvar,mufrac,coefT1,coefT0,Temperature
!  stop
  
end subroutine tp_mu
!------------------------------------------------------------------
real(DP) function Temperature(spvar) 
  use nrtype
  use vars
  implicit none

  real(DP),intent(in) :: spvar
  
  Temperature = spvar*spetkl
!  if(Temperature < 3.d3)then
!  if(Temperature < 1.d3)then
!  if(Temperature < 1.de3)then
!  if(Temperature < 2.d3)then
  if(Temperature < Tlow*mufrac)then
!     Temperature = 2.d0*Temperature
     Temperature = Temperature / mufrac
!  else if(Temperature < 1.d4)then
  else if(Temperature < Thigh)then
!     Temperature = .5714d0*Temperature + 4285.7d0
!     Temperature = .8889d0*Temperature + 1111.d0
!     Temperature = .8235d0*Temperature + 1765.d0
!     Temperature = .75d0*Temperature + 2500.d0
     Temperature = coefT1*Temperature + coefT0
  endif

end function Temperature

!--------------------------------------------------------------------
subroutine e2T(ip,Tp1,mu1,xe)
  use nrtype
  use vars
  implicit none

! derive Tp (& mu) from given speene by calculating ionization

  integer(I4B),intent(in) :: ip
  real(DP),intent(out) :: Tp1,mu1,xe
  integer(I4B),parameter :: niter = 100
  integer(I4B) :: iter 
  real(DP),parameter :: epsilon = 1.d-6
  real(DP) :: spvar,Temperature,Tp,Tp0,sp0,sp1,dfsp,mu0


  spvar = speene(ip)
!  write(*,*)radius(ip),densit(ip),spvar
  
  Tp0 = Temperature(spvar) ! initial guess of Tp  
  call T2e(ip,Tp0,sp0,mu0,xe) ! derive sp0 & mu0 from Tp0
!  if(ip==5)write(*,*)'0',ip,sp0,spvar,Tp0,mu0
  if(abs((sp0-spvar)/spvar) < epsilon)then
     write(*,*)"No iteration needed",ip
     Tp1 = Tp0; mu1 = mu0
     return
  end if
  Tp1 = Tp0*1.01d0
  do iter = 1,niter
     call T2e(ip,Tp1,sp1,mu1,xe) ! derive sp1 & mu1 from Tp1
!     if(ip==5)write(*,*)iter,ip,sp1,spvar,Tp1,mu1
     dfsp = abs((sp1-spvar)/spvar)
     if(dfsp < epsilon)exit      
     Tp = Tp0 + (Tp1-Tp0) * (spvar - sp0) / (sp1 - sp0)  
     Tp0 = Tp1
     Tp1 = Tp
     sp0 = sp1
  end do
  if(iter==niter .and. dfsp>epsilon)then
     write(*,*)"No convergence of iteration"
  else
!     write(*,*)sp1,mu1,Tp1
  end if
  
end subroutine e2T

!----------------------------------------------------------------
subroutine T2e(ip,Tp,spetmp,mu,xe)
  use nrtype
  USE vars
  implicit none   

  integer(I4B),intent(in) :: ip
  real(DP),intent(in) :: Tp
  real(DP),intent(out) :: spetmp,mu,xe
!  real(DP),parameter :: bk  = 1.38066d-16, hh  = 6.626076d-27 &
!       & ,me = 9.109390d-28, c_pannekoek = 0.5d0, 
  real(DP),parameter :: wgal =1.d-14, tgal =7500.d0, hvfc = 1.207354d15 
         !hvfc = c_pannekoek*(2.d0*pi_d*me*bk/hh**2)**(1.5d0)
  real(DP),parameter :: epsrs = 1.0d-3
  real(DP),parameter :: xemin = 1.0d-15 !, Tpmin = 3.d3
  integer,parameter :: numele=11,niter = 100
!     C    O   Na    Mg     Al      Si     S     K    Ca    Cr    Fe    
  real(DP):: Ax(1:numele),IK(1:numele) !,Iev(1:numele)
! Ax : relative to H (revised to Asplund+2009)
!  data Ax/3.0d-4,7.0d-4,2.0d-6,4.4d-5,3.0d-6,3.6d-5,1.5d-5,1.2d-7,2.1d-6,5.5d-7,3.5d-5/
  data Ax/2.7d-4,4.9d-4,1.9d-6,3.4d-5,2.7d-6,3.2d-5,1.4d-5,1.2d-7,2.0d-6,4.4d-7,2.8d-5/
  data IK/1.3667d5,1.5805d5,5.9647d4,8.8658d4,6.9511d4,9.4577d4,1.2022d5 &
       & ,5.0364d4,7.0904d4,7.8685d4,9.1328d4/
!  data Iev/11.26d0,13.62d0,5.14d0,7.64d0,5.99d0,8.15d0,10.36d0,4.34d0,6.11d0,6.78d0,7.87d0/
     !IK(1:numele)  =  11604.515d0*IeV(1:numele) in K unit
  real(DP) :: amet,tbalmer,sm,vturb,hcol,wz,fte,tau21,n1nk &
       & ,c_h,c_hei,c_heii,npos,Temperature,fltb
  real(DP) :: trad,hvfctp,hvfctp_star,hvfctp_gal,ane,axii,xii,chv
  integer(I4B) :: j,i,it
  
!  sig     = 1.4d+00                ! mean mass per hydrogen nuclei 
  amet    = 0.d0 ! 5.d-04      ! contribution from photoionized metals
  tbalmer = 3000.d0 !max(3000.d0, 0.85d0*Teff)  
  trad = 0.8d0*Teff 

  sm = 0.d0
  do j = ip+1,np
     sm = sm + densit(j)*dradiu(j)
  end do
  
  vturb   = 12.85d5 * SQRT(Tp/1.d+04)    ! cm s-1
!     hcol    = rho_scale*nh  ; estimate overlying column density cm-2
!     hcol = sm * rdkm * 1.d5 * rhocgs / (sig*prtnms)
  hcol = sm * clmfct
  wz      = 0.5d0*(1.d0 -SQRT(1.d0 -1.d0/max(radius(ip), 1.0001d0)**2))
  fte     = 6.96d-8 * EXP(1.184d5/Tp + 3.946d4/tbalmer) &
       & /Tbalmer/SQRT(Tp)
!     tau21   = 0.0265d0 * 1215.67d-8 * 0.4164d0 * hcol / SQRT(pi_d) / vturb
  tau21   = 7.568d-8 * hcol  / vturb
  n1nk    = 6.265d8 * fte / (tau21*wz)
  c_h   = 10.d0**(-EXP((4.1d0-LOG10(Tp))/9.d-2))
!     c_hei = 0.1d0*10.d0**(-EXP((4.35d0-LOG10(Tp))/9.d-2)) 
!     c_heii= 0.2d0*10.d0**(-EXP((4.75d0-LOG10(Tp))/9.d-2))
  c_hei = 0.25d0*YoX*10.d0**(-EXP((4.35d0-LOG10(Tp))/9.d-2)) 
  c_heii= 0.5d0*YoX*10.d0**(-EXP((4.75d0-LOG10(Tp))/9.d-2))
  npos  = max(1.d0/(1.d0 + n1nk), c_h) + max(c_hei, c_heii)
!  if(ip==3)write(*,*)npos,n1nk,sm,hcol
  
     ! for heavy element
  If(zmetal==0.d0)then
     xe = npos
  else
     hvfctp = sqrt(Tp) * hvfc  
     hvfctp_star = hvfctp*wz*trad
     hvfctp_gal  = hvfctp*wgal*tgal
!     xe = 1.d-4 + npos              ! initial guess for ne/nH
     xe = 1.d-4 * zmetal * (1.d0 - exp(-Tp/3000.d0)) + npos 
     do it=1, niter 
        ane = xe * densit(ip) * dnnH1
!        axii=0.d0
        axii = npos 
        do i=1, numele
           chv = hvfctp_star*EXP(-IK(i)/trad) + hvfctp_gal*EXP(-IK(i)/tgal)
           xii = zmetal*Ax(i)*chv/(chv + ane) ! dependence on metallicity
           axii = axii + xii
!           ioni(i) = xii
        end do
!        err = ABS((xe-axii)/xe) 
!        if(ip==100)write(*,*)it,xplace(ip),Tp,xe,axii,npos,c_h,c_hei,c_heii &
!             & ,n1nk,vturb,hcol,tau21
!        if(ip==5)write(*,*)it,xe,axii
        if (abs((xe-axii)/xe) < epsrs)exit
        xe = sqrt(axii*xe)
!        xe = axii
!     end if
     end do
     if(it>niter)then
        write(*,*)'it>niter',ip,xe
        stop
     end if
     xe = axii
!     xe = amet + npos
  end If
  xe = max(xe,xemin)
!  write(*,*)xe
  
  mu = 1.d0 / ((1.d0 - YHe - zmetal*Zsun) * (1.d0+xe) &
       & + 0.25d0*YHe + 0.1d0*zmetal*Zsun) ! A_heavyelem = 10

  spetmp = mufull / (spetkl * mu) * Tp

!  write(*,*)ip,sm,xe,mu,Tp
!  stop
  
end subroutine T2e

!------------------------------------------------------------------------------
subroutine Teff2akm
  use nrtype
  USE vars
  implicit none   

  real(DP) :: mu
  
!  mu = 1.d0 / ((1.d0 - YHe - zmetal*Zsun) * (1.d0+xe) &
!       & + 0.25*YHe + 0.1d0*zmetal*Zsun) ! A_heavyelem = 10

  if(Teff<=6.d3)then
     mu = 1.2d0
  else
     write(*,*)"Teff > 6.d3"
     stop
  end if
  
  akm = sqrt(Teff * Rcon / mu) * 1.d-5
  write(*,*)Teff,akm,mu !,xe,npos,c_h,c_hei, c_heii,wz,YHe
  
end subroutine Teff2akm

!----------------------------------------------------------------
subroutine Teff2akm2
  use nrtype
  USE vars
  implicit none   

  real(DP),parameter :: wgal =1.d-14, tgal =7500.d0, hvfc = 1.207354d15 
  real(DP),parameter :: epsrs = 1.0d-3
  real(DP),parameter :: xemin = 1.0d-15 !, Tpmin = 3.d3
  integer,parameter :: numele=11,niter = 100
!     C    O   Na    Mg     Al      Si     S     K    Ca    Cr    Fe    
  real(DP):: Ax(1:numele),IK(1:numele) !,Iev(1:numele)
! Ax : relative to H (revised to Asplund+2009)
  data Ax/2.7d-4,4.9d-4,1.9d-6,3.4d-5,2.7d-6,3.2d-5,1.4d-5,1.2d-7,2.0d-6,4.4d-7,2.8d-5/
  data IK/1.3667d5,1.5805d5,5.9647d4,8.8658d4,6.9511d4,9.4577d4,1.2022d5 &
       & ,5.0364d4,7.0904d4,7.8685d4,9.1328d4/
  real(DP) :: amet,tbalmer,sm,vturb,hcol,wz,fte,tau21,n1nk &
       & ,c_h,c_hei,c_heii,npos,Temperature,fltb
  real(DP) :: trad,hvfctp,hvfctp_star,hvfctp_gal,ane,axii,xii,chv
  real(DP) :: mu,xe
  integer(I4B) :: i,it
  
  amet    = 0.d0 ! 5.d-04      ! contribution from photoionized metals
  tbalmer = 3000.d0 !max(3000.d0, 0.85d0*Teff) 
  trad = 0.8d0*Teff

  sm = 1.d0/3600.d0 * Rstar/Mstar * (Teff/5780.d0)
  
  vturb   = 12.85d5 * SQRT(Teff/1.d+04)    ! cm s-1
  hcol = sm * clmfct
  wz      = 0.5d0*(1.d0 -SQRT(1.d0 -1.d0/1.0001d0**2))
  fte     = 6.96d-8 * EXP(1.184d5/Teff + 3.946d4/tbalmer) &
       & /Tbalmer/SQRT(Teff)
  tau21   = 7.568d-8 * hcol  / vturb
  n1nk    = 6.265d8 * fte / (tau21*wz)
  c_h   = 10.d0**(-EXP((4.1d0-LOG10(Teff))/9.d-2))
  c_hei = 0.25d0*YoX*10.d0**(-EXP((4.35d0-LOG10(Teff))/9.d-2)) 
  c_heii= 0.5d0*YoX*10.d0**(-EXP((4.75d0-LOG10(Teff))/9.d-2))
  npos  = max(1.d0/(1.d0 + n1nk), c_h) + max(c_hei, c_heii)
!  write(*,*)npos,n1nk,c_h,sm,Teff,hcol
  
     ! for heavy element
  If(zmetal==0.d0)then
     xe = npos
  else
     hvfctp = sqrt(Teff) * hvfc  
     hvfctp_star = hvfctp*wz*trad
     hvfctp_gal  = hvfctp*wgal*tgal
     xe = 1.d-4 * zmetal * (1.d0 - exp(-Teff/3000.d0)) + npos 
     do it=1, niter 
        ane = xe * dnnH1 ! densit=1 at r=1
        axii = npos 
        do i=1, numele
           chv = hvfctp_star*EXP(-IK(i)/trad) + hvfctp_gal*EXP(-IK(i)/tgal)
           xii = zmetal*Ax(i)*chv/(chv + ane) ! dependence on metallicity
           axii = axii + xii
        end do
        if (abs((xe-axii)/xe) < epsrs)exit
        xe = sqrt(axii*xe)
     end do
     if(it>niter)then
        write(*,*)'it>niter',xe
        stop
     end if
     xe = axii
  end If
  xe = max(xe,xemin)
  
  mu = 1.d0 / ((1.d0 - YHe - zmetal*Zsun) * (1.d0+xe) &
       & + 0.25d0*YHe + 0.1d0*zmetal*Zsun) ! A_heavyelem = 10

  akm = sqrt(Rcon / mu* Teff)*1.d-5
  write(*,*)"mu=",mu,"akm=",akm
  
end subroutine Teff2akm2
  
!-------------------------------------------------------------------------------------
subroutine Teff2akm_sun
  use nrtype
  USE vars
  implicit none   

  real(DP),parameter :: wgal =1.d-14, tgal =7500.d0, hvfc = 1.207354d15 
  real(DP),parameter :: epsrs = 1.0d-3
  real(DP),parameter :: xemin = 1.0d-15 !, Tpmin = 3.d3
  integer,parameter :: numele=11,niter = 100
!     C    O   Na    Mg     Al      Si     S     K    Ca    Cr    Fe    
  real(DP):: Ax(1:numele),IK(1:numele) !,Iev(1:numele)
! Ax : relative to H (revised to Asplund+2009)
  data Ax/2.7d-4,4.9d-4,1.9d-6,3.4d-5,2.7d-6,3.2d-5,1.4d-5,1.2d-7,2.0d-6,4.4d-7,2.8d-5/
  data IK/1.3667d5,1.5805d5,5.9647d4,8.8658d4,6.9511d4,9.4577d4,1.2022d5 &
       & ,5.0364d4,7.0904d4,7.8685d4,9.1328d4/
  real(DP) :: amet,tbalmer,sm,vturb,hcol,wz,fte,tau21,n1nk &
       & ,c_h,c_hei,c_heii,npos,Temperature,fltb
  real(DP) :: trad,hvfctp,hvfctp_star,hvfctp_gal,ane,axii,xii,chv
  real(DP) :: mu,xe
  integer(I4B) :: i,it
  
  amet    = 0.d0 ! 5.d-04      ! contribution from photoionized metals
  tbalmer = 3000.d0 !max(3000.d0, 0.85d0*Teffsun) 
  trad = 0.8d0*Teffsun

  sm = 1.d0/3600.d0 !* Rstar/Mstar * (Teff/5780.d0)
  
  vturb   = 12.85d5 * SQRT(Teffsun/1.d+04)    ! cm s-1
  hcol = sm * clmfct
  wz      = 0.5d0*(1.d0 -SQRT(1.d0 -1.d0/1.0001d0**2))
  fte     = 6.96d-8 * EXP(1.184d5/Teffsun + 3.946d4/tbalmer) &
       & /Tbalmer/SQRT(Teffsun)
  tau21   = 7.568d-8 * hcol  / vturb
  n1nk    = 6.265d8 * fte / (tau21*wz)
  c_h   = 10.d0**(-EXP((4.1d0-LOG10(Teffsun))/9.d-2))
!  c_hei = 0.25d0*YoX*10.d0**(-EXP((4.35d0-LOG10(Teffsun))/9.d-2)) ! =0 for Teffsun 
!  c_heii= 0.5d0*YoX*10.d0**(-EXP((4.75d0-LOG10(Teffsun))/9.d-2))
  npos  = max(1.d0/(1.d0 + n1nk), c_h) !+ max(c_hei, c_heii)

  
     ! for heavy element
  If(zmetal==0.d0)then
     xe = npos
  else
     hvfctp = sqrt(Teffsun) * hvfc  
     hvfctp_star = hvfctp*wz*trad
     hvfctp_gal  = hvfctp*wgal*tgal
     xe = 1.d-4 * zmetal * (1.d0 - exp(-Teffsun/3000.d0)) + npos 
     do it=1, niter 
        ane = xe * dnnH1 ! densit=1 at r=1
        axii = npos 
        do i=1, numele
           chv = hvfctp_star*EXP(-IK(i)/trad) + hvfctp_gal*EXP(-IK(i)/tgal)
           xii = zmetal*Ax(i)*chv/(chv + ane) ! dependence on metallicity
           axii = axii + xii
        end do
        if (abs((xe-axii)/xe) < epsrs)exit
        xe = sqrt(axii*xe)
     end do
     if(it>niter)then
        write(*,*)'it>niter',xe
        stop
     end if
     xe = axii
  end If
  xe = max(xe,xemin)
  
  mu = 1.d0 / ((1.d0 - YHe - zmetal*Zsun) * (1.d0+xe) &
       & + 0.25d0*YHe + 0.1d0*zmetal*Zsun) ! A_heavyelem = 10

  akmsun = sqrt(Rcon / mu* Teffsun)*1.d-5
!  write(*,*)mu,akmsun
  
end subroutine Teff2akm_Sun
  
