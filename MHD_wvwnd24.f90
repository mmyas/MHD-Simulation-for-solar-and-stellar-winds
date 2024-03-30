program MUSCL
  use nrtype
  USE vars
  implicit none

  call setup
!      IF (PGOPEN('?') .LE. 0) STOP

  select case(icont)
  case(0)
     call initia
  case(1)
     call initiaread
  case(2)
     call initiaread_rdexpnd
  end select

   write(*,*) "dvelo=",dvelo,"dvtrit=",dvtrit,"dvttit=",dvttit
  
  call evolve

!      CALL PGCLOS
!  ------------------------------------------------------------------
  close(unit= 8)
  close(unit=11)
  close(unit=13)
  close(unit=16)
  close(unit=43)
  close(unit=61)
  
end program MUSCL
!--------------------------------------------------------------------------------------------------------
subroutine setup
  use nrtype
  use vars
  implicit none

  real(DP) :: Temperature,Tcrtrad, taucsc
  real(DP) :: Tmfrac,mu,xe,spetmp
  INTEGER(I4B) :: i,irec
  
! ------------------------------------------------------------------
  open(unit= 8,file='InputParameters_MHDstwind')
!  open(unit= 8,file='Input_Parameters_MHDstwind_Betelgeuse' )
  open(unit= 9,file='sp_f-1_100rnm1.dat',status='old')
  open(unit= 10,file='sp_f-1_100rnm2.dat',status='old')
  open(unit= 18,file='sp_f-1_100rnm3.dat',status='old')
!      open(unit= 9,file='spinp_10_2.dat')
!      open(unit= 9,file='randpert1e4.dat')
  open(unit=11,file='Initial_Data2')

  open(unit=13,file='energyflow.dat')
!  open(unit=16,file='SWout816_30.dat')
  !      open(unit=16,file='plebm_sw710_14_r1e-24_v3e9_p0.1_t0.05.dat')
  
  open(unit=43,file='radflx2.dat')
  open(unit=61,file='waveenergy.dat') 
!   ----------------------------
!   fc1 = 0.5d0/sqrt(pi)
  fc1 = 1.d0
  fc2 = fc1*fc1
!      fc3 = .125/pi
  fc3 = .5d0*fc2

  read(8,*)  icont,irss,itrb 
  read(8,*)  np             ! Number of Grid points
!   read(8,*)  ql             ! L1 & L2 (Length to the Boundary)
  read(8,*)  Bplint,Btrint
  read(8,*)  dvelo,dBtrit,dBttit     ! v perturbation  (km/s)
      ! Underlied field & Wave Amplitude
  read(8,*)  cflnum,tfinal,gamma ! Time Interval
  read(8,*)  facdep,critit  ! facdep,critit
  read(8,*)  Mstar,Teff,rhocgs,Rstar,rhocrt,zmetal,Tclso
  read(8,*)  fmaxc
  read(8,*)  vrot0,dslmbd0 !,tefave
  read(8,*)  drin, drout, einit, rhoturn, drfac, drcnst, dtcd0
  gammi1 = gamma - 1.0d0 
!  ----------------------------[ for Godunov ]
  critif = 1.0d0 - 1.0d-6
  rdkm = Rstar * Rsun
  rtgamm = sqrt( gamma )
!      gammp1 = 0.5*( gamma + 1.0 )
!      gammm1 = 0.5*( gamma - 1.0 )
  gammp1 = gamma + 1.0d0 
  gammx1 = gammi1/gamma
  gammx2 = gammi1/rtgamm
  gammp3 = gamma + 3.0d0
  gammm3 = gamma - 3.0d0
  gammm2 = gamma - 2.0d0

  dinit = 1.d0
  pinit = 1.d0
  vinit = 0.d0
  iwv = 0

  Tclso = Tclso*Teff
  
  YHe = Ysun + dYdZ*(zmetal - 1.d0)*Zsun
  XHyd = 1.d0 - YHe - zmetal*Zsun
  YoX = YHe / XHyd
  mufull = 1.d0 / (2.d0*XHyd + 0.75d0*YHe + 0.45d0*zmetal*Zsun)
  muzero = 1.d0 / (XHyd + 0.25d0*YHe + 0.05d0*zmetal*Zsun)
  mufrac = mufull / muzero
  Tmfrac = Thigh - Tlow*mufrac
  coefT1 = (Thigh - Tlow) / Tmfrac
  coefT0 = Thigh*Tlow*(1.d0 - mufrac) / Tmfrac

  dnnH1 = rhocgs*XHyd/ prtnms
  clmfct = rdkm * 1.d5 * dnnH1 
  
!  call Teff2akm
  call Teff2akm2
!   
  dvelo = sqrt(2.d0)*dvelo/akm
  dBtrit = sqrt(2.d0)*dBtrit/akm
  dBttit = sqrt(2.d0)*dBttit/akm

  call Teff2akm_sun ! Solar normalization
!  taucsc = Rstar/Mstar * akm*akm / (akmsun*akmsun) ! normalization for wave period
  taucsc = Rstar/Mstar * Teff / Teffsun ! normalization for wave period
  sgmc = sgmcsun * taucsc 
  re1c = 1.d0 + 2.d0*sgmc ! H/R + R

  grav = Gslm*Mstar * 1.d-15/(akm*akm*rdkm)
    ! Gslm : Msun*G
!      cond = cdcon / (Rcon*2.)**3.5 *akm**4 * 1.e15/(rdkm*rhocgs)*.7
                                ! 2. = 1/mean.mol.weight
!      radfac = 1.e-10 * rdkm/(akm*akm*akm) * rhocgs 
!     &     / (.25*prtnms*prtnms)  * 5.e-3
!      radfac0 = 1.e-10 * rdkm/(akm*akm*akm) * rhocgs 
!     &     / (.25*prtnms*prtnms) 
!      spetkl = 1.e10* akm*akm * .5/Rcon

!-----------for mu=mufull---------------------------------
  cond = cdcon / (Rcon/mufull)**3.5d0 *akm**4.d0 * 1.d15/(rdkm*rhocgs) &
       &     *gammi1**3.5d0
                                ! 2. = 1/mean.mol.weight
  radfac = 1.d-10 * rdkm/(akm*akm*akm) * rhocgs &
       &     / (mufull*mufull*prtnms*prtnms)  !* 5.d-3
!  radfac0 = 1.d-10 * rdkm/(akm*akm*akm) * rhocgs &
!       &     / (.36d0*prtnms*prtnms) 
  spetkl = 1.d10* akm*akm * mufull/Rcon * gammi1
!  Teff = Temperature(1.d0/gammi1)
  etafct = 2.34d-8 / (rdkm*akm)
  eadfct = 2.1d-16 * akm / (rhocgs*rdkm)
!--------------------------------------------------------
!      rhocrt = 1.e-15 ! critical density for optically thin condition 
!      rhocrt = 1.23e-17 ! for Tcr=4.e4
!  rhocrt = 4.9d-17 ! for Tcr = 2e4
!      rhocrt = 8.7e-17 ! for Tcr = 1.5e4
!      rhocrt = 1.25e-16 ! for Tcr = 1.4e4
! rhocrt
  call calcTcr(Tcrtrad)
  Tcrt = 10.d0**Tcrtrad!(zmetal) ! in Kelvin
  
!  write(**,*)grav,cond,radfac,spetkl,akm,Teff,Tcrt,YHe,YoX,mufull,muzero,coefT1,coefT0 &
!       & ,taucsc,Tclso,re1c,sgmc
  write(*,*)"Tcrt=",Tcrt,"YHe=",YHe,"nufull=",mufull,"muzero=",muzero,"taucsc=",taucsc&
       & ,"Tclso=",Tclso,"re1c=",re1c,"sgmc=",sgmc
  
  nd = np + 1

  if(icont==2)then
     np_ori = np
     nd_ori = nd
!     np = int(1.5*real(np))
     np = np_ori+500
     nd = np+1
  end if

  allocate(radius(nd),dradiu(nd),xplace(nd),dxplac(nd),ddxdtt(nd) &
       & ,densit(nd),volume(nd),pressu(nd),speene(0:nd),energy(nd) &
       & ,sounds(nd),dxdt(nd),dvoldt(nd),denedt(nd) &
       & ,gradd(nd), gradp(nd), gradv(nd), grade(nd) &
       & ,xlagra(nd),dxlagr(nd), prsgas(nd),Tpprt(nd) &
       & ,Bparal(nd),Alfven(nd),Btense(nd),Bcenrf(nd),Btrans(nd),vtrans(nd)&
       & ,dBtran(nd),dvtran(nd),dtdBtr(nd),dtvrtr(nd),fBtrbd(nd),vftrbd(nd)&
       & ,Btrant(nd),vtrant(nd),dBtrat(nd),dvtrat(nd),dtdBtt(nd),dtvrtt(nd) &
       & ,fBttbd(nd),vfttbd(nd),clocal(nd) &
       & ,radbnd(nd),dxdtbd(nd),rsqfbnd(nd),Btrrbd(nd),Btrlbd(nd) &
       & ,Bttrbd(nd),Bttlbd(nd) &
       & ,xplori(nd),dxplor(nd),radori(nd),drador(nd),dxlg(nd) &
       & ,fBtreb(nd),fBtteb(nd),vftreb(nd),vftteb(nd),dxdteu(nd) &
       & ,dxdtav(nd),densav(nd),enrgav(nd),prgsav(nd),spenav(nd) &
       & ,Btraav(nd),vtraav(nd),Btrtav(nd),vtrtav(nd),zplus(nd),zmins(nd) &
       & ,dvbaav(nd),dvbtav(nd),radcol(nd) &
       & ,Falout(nd),Falin(nd),Fcsout(nd),Faloav(nd),Faliav(nd),Fcsoav(nd) &
       & ,Falout3(nd),Falin3(nd),Faloav3(nd),Faliav3(nd) &
       & ,ef(10,nd),Btrprv(nd),Bttprv(nd),engtot(nd) &
       & ,eta(nd),ead(nd),ettrb(nd)) !,efav(9,nd)
  allocate(etaprt(nd),eadprt(nd),muprt(nd)) !,xeprt(nd)
  etaprt(1:nd) = 0.0; eadprt(1:nd) = 0.0; muprt(1:nd) = 1.2 !xeprt(1:nd) = 1.0d0 ! initial value
  Faloav(1:nd) = 0.0d0 ; Faliav(1:nd) = 0.0d0; Fcsoav(1:nd) = 0.0d0
  Faloav3(1:nd) = 0.0d0; Faliav3(1:nd) = 0.0d0

  inquire(iolength=irec)t,radius,dradiu,densit,dxdt,vtrans,vtrant &
       & ,Btrans,Btrant,speene(1:nd),dxdtbd,real(Alfven),etaprt,eadprt,muprt &
       & ,real(densav),real(dxdtav),real(vtraav),real(vtrtav),real(Btraav) &
       & ,real(Btrtav),real(spenav),real(dvbaav),real(dvbtav)
  open(unit=12,file='printvars.dat',access='direct',form='unformatted' &
       & ,recl = irec)
!       & ,recl=icmpl*(2+(2*nmvar+nmave)*nd))

!  ------------------------------------------------------------------
! read random perturbation
  do i = 1,nrand
!     read(9,*)trd(i),random(i)
!     trd(i) = taucsc * trd(i)
!     read(10,*)trdlg(i),randlg(i)
!     trdlg(i) = taucsc * trdlg(i)
!     read(18,*)trd3c(i),rand3c(i)
!     trd3c(i) = taucsc * trd3c(i)
     trd(i) = dble(i-1) * trandsun * taucsc 
     read(9,*)random(i)
     read(10,*)randlg(i)
     read(18,*)rand3c(i)
  enddo

  close(9);close(10);close(18)
  
end subroutine setup
!======================================================================*
!                     Prepare an Initial State                         *
!======================================================================*
subroutine initia
  use nrtype
  USE vars
  implicit none

  integer(I4B),parameter :: iran1=12869,iran2=6925,iran3=32768
  integer(I4B) :: ip,ivtr,iddns
  real(DP) :: dr,volum,press,uinit,sound !,rhoturn ,einit,drout
  real(DP) :: rnow,Aflow,Bprs,velmax,dtmin,drad,rad,spdm,fac,fltb,dttmp &
       & ,r1,r2
  integer(I4B) :: i,np_buf

  t = 0.d0
  irstr = 0

! initialization of averaged values
!   dtavtt=0.d0
!   do i = 1,nd
!      dxdtav(i)=0.d0
!      densav(i)=0.d0
!      enrgav(i)=0.d0
!      prgsav(i)=0.d0
!      spenav(i)=0.d0
!      Btraav(i)=0.d0
!      vtraav(i)=0.d0
!      Btrtav(i)=0.d0
!      vtrtav(i)=0.d0
!      zplus(i)=0.d0
!      zmins(i)=0.d0
!   enddo


!  ------------------------------------------------------------------
!  dr = 2.d-3
!  drout = 2.d-3
!  rhoturn = 2.d-6
!  einit = 25.0d0
!  dr = 4.d-5
  dr = drin
!  -------------------------------------
  volum = 1.0d0/dinit
!      press = pinit + Bpre 
  press = pinit
  uinit = pinit*volum/gammi1
  sound = sqrt( gamma*pinit*volum )
  
!  dvelo = ddens*sound/(dinit + 2.d0*ddens) ! dvelo is input
  ddens = dvelo * dinit / max((sound - 2.d0*dvelo),0.1d0) 
  dpres = gamma*ddens*press/dinit ! p perturbation
  
!  write(*,*) dvelo,ddens,dpres

  dvtrit = fc1*dBtrit/sqrt(dinit)
  
  dvttit = fc1*dBttit/sqrt(dinit)
  
  ip = 1
!  -------------------
  dradiu(ip) = dr
!      dradiu(ip) = (.5 * ( 2.e-3 + 1.e-5 ) + ( 2.e-3 - 1.e-5 ) 
!     &     / pi * atan( 3000.*(1. - 1.01) ) ) 
!     &     * 1.
  radius(ip) = 1.d0
  rnow = radius(ip)  
  Aflow = rnow*rnow*fltb(rnow)
  r1 = radius(ip) - 0.5d0*dradiu(ip)
  r2 = radius(ip) + 0.5d0*dradiu(ip)
!  dxplac(ip) = dradiu(ip) * Aflow
  dxplac(ip) = 0.5d0 * dradiu(ip) *(r1*r1*fltb(r1) + r2*r2*fltb(r2))
!  dxplac(ip) = onovth*(r2*r2*r2*fltb(r2) - r1*r1*r1*fltb(r1))
!      xplace(ip) = .333333
  xplace(ip) = 1.d0
  radori(ip) = radius(ip)
  drador(ip) = dradiu(ip)
  xplori(ip) = xplace(ip)
  dxplor(ip) = dxplac(ip) 
  radbnd(ip) = radius(ip) + .5d0*dradiu(ip)
  rsqfbnd(ip) = radbnd(ip)*radbnd(ip) * fltb(radbnd(ip))
!     xposit(ip) = xcoord
  dxdt(ip) = vinit
  ddxdtt(ip) = 0.d0
  densit(ip) = dinit
  volume(ip) = 1.d0/densit(ip)
  prsgas(ip) = press
  speene(ip) = prsgas(ip)*volume(ip)/gammi1
  sounds(ip) = sqrt(gamma*prsgas(ip)*volume(ip))
  clocal(ip) = dtmin
  dxlagr(ip) = dxplac(ip) * dinit 
!      xlagra(ip) = .333333
  xlagra(ip) = 1.d0
!  -------------------
  Bparal(ip) = Bplint / Aflow
  Alfven(ip) = fc1*Bparal(ip)/sqrt(dinit)
  Btrans(ip) = 0.d0
  vtrans(ip) = 0.d0 
  vtrant(ip) = dvttit + vrot0   ! for circular polarized wave
  Btrant(ip) = -dBttit       ! for circular polarized wave
  ettrb(ip) = 0.025d0*rdkm / dslmbd0 / sqrt(Aflow)
!----------------------------------------------
  Bprs = fc3 * (Btrans(ip)*Btrans(ip) + Btrant(ip)*Btrant(ip))
  pressu(ip) = prsgas(ip) + Bprs
  energy(ip) = speene(ip) + .5d0*(dxdt(ip)*dxdt(ip) &
       &  + vtrans(ip)*vtrans(ip) + vtrant(ip)*vtrant(ip) )  &
       &  + Bprs / densit(ip) - grav / radius(ip)
  velmax = sqrt ( sounds(ip)*sounds(ip) + fc2/densit(ip)  &
       & * (Bparal(ip)*Bparal(ip) + Btrans(ip)*Btrans(ip) &
       & + Btrant(ip)*Btrans(ip) ) ) + abs(dxdt(ip))
  dtmin = dradiu(ip)/velmax

! ------------ ghost mesh at ip = 0----------------------
  drad = dr
  rad = radius(ip) - drad
  Aflow = rad * rad *fltb(rad)
  r1 = radius(ip) - 1.5d0 * drad
  r2 = radius(ip) - 0.5d0 * drad
!  dxplac0 = drad * Aflow 
  dxplac0 = 0.5d0 * dradiu(ip) *(r1*r1*fltb(r1) + r2*r2*fltb(r2))
!  dxplac0 = onovth * (r2*r2*r2*fltb(r2) - r1*r1*r1*fltb(r1))  
  densit0 =  dinit * exp(grav*(1.d0/rad - 1.d0) )
  pressu0 = press / dinit * densit0
  dxdt0 = 0.d0
  speene0 = pinit / (densit0*gammi1)
  Btrans0 = 0.d0
  vtrans0 = 0.d0
  Btrant0 = 0.d0
  vtrant0 = 0.d0
  energy0 = speene0 + .5d0*( dxdt0*dxdt0 + vtrans0*vtrans0 &
       & + vtrant0*vtrant0 ) + fc3*( Btrans0*Btrans0 &
       & + Btrant0*Btrant0 ) / densit0 - grav / rad
!--------------------------------------------------------      
   
  ivtr = 1
  iddns = 1

!  np_buf = int(real(np)/1.5)
  np_buf = np-300

  do ip = 2,np + 1
!         xcoord = dx*float(idummy)
!        -------------------
!         dradiu(ip) = dr * radius(ip-1)**20.
!         if(radius(ip-1) .lt. 1.0099)then
!            dradiu(ip) = (.5 * ( 2.e-3 + 1.e-5 ) + ( 2.e-3 - 1.e-5 ) 
!     &           / pi * atan( 3000.*(radius(ip-1) - 1.01) ) ) 
!     &           * radius(ip-1)
!         else
!            dradiu(ip) = (.5 * ( 2.e-3 + 1.e-5 ) + ( 2.e-3 - 1.e-5 ) 
!     &           / pi * atan( 60.*(radius(ip-1) - 1.01) ) ) 
!     &           * radius(ip-1)
!          endif
!         dradiu(ip) = min(max(dr * sqrt(Alfven(ip-1)/Alfven(1))
!     &        ,dradiu(ip-1)) * radius(ip-1)**.5,5.e-3)
!     if(ip > (np-500))then
     if(ip > (np-200))then
        dradiu(ip) = 1.02d0*dradiu(ip-1)
     else if(ip > np_buf)then
        dradiu(ip) = dradiu(np_buf) * (1.0d0 + 0.001d0 &
             & * (dble(ip-np_buf))**2.d0 / dble(np-np_buf))
        !     else if(densit(ip-1)>rhoturn)then
!!            dradiu(ip) = min(max(dr * (Alfven(ip-1)/Alfven(1))**.5
!!     &           ,dradiu(ip-1)) * radius(ip-1)**.5,5.e-3)
!        dradiu(ip) = min(max(dr * (Alfven(ip-1)/Alfven(1))**.3d0 &
!             & ,dradiu(ip-1)),drout)
     else if(radius(ip-1)<drcnst)then
        dradiu(ip) = drin
     else
!        dradiu(ip) = min(dradiu(ip-1)*radius(ip-1),drout)
        dradiu(ip) = min(dradiu(ip-1)+drfac*drin*radius(ip-1),drout)
!        dradiu(ip) = min(dradiu(ip-1)*1.005d0,drout)
     endif
!         if(ip.eq.2)dradiu(ip) = .01*dradiu(ip)
     radius(ip) = radius(ip-1)+.5d0*( dradiu(ip-1) + dradiu(ip) )
     rnow = radius(ip)  
     Aflow = rnow*rnow*fltb(rnow)
     r1 = radius(ip) - 0.5d0*dradiu(ip)
     r2 = radius(ip) + 0.5d0*dradiu(ip)
!     dxplac(ip) = dradiu(ip) * Aflow
!     write(*,*)Aflow,0.5d0*(r1*r1*fltb(r1)+r2*r2*fltb(r2))
     dxplac(ip) = dradiu(ip) * 0.5d0*(r1*r1*fltb(r1)+r2*r2*fltb(r2))
!     dxplac(ip) = onovth*(r2*r2*r2*fltb(r2) - r1*r1*r1*fltb(r1))

     xplace(ip) = xplace(ip-1) + .5d0*(dxplac(ip-1)+dxplac(ip))
     radori(ip) = radius(ip)
     drador(ip) = dradiu(ip)
     xplori(ip) = xplace(ip)
     dxplor(ip) = dxplac(ip) 
     radbnd(ip) = radius(ip) + .5d0*dradiu(ip)
     rsqfbnd(ip) = radbnd(ip)*radbnd(ip) * fltb(radbnd(ip))
!        xposit(ip) = xcoord
     dxdt(ip) = vinit !* radius(ip)**2
     ddxdtt(ip) = 0.d0
!         densit(ip) = max(1.e-10*dinit,
!     &        dinit * dexp(grav*(1./radius(ip) - 1.) ) )
!         if(prsgas(ip-1).lt.2.e-1)then
     if(densit(ip-1) < rhoturn)then
!        fac = 1.d0 + 2.d0/spdm
        fac = 2.d0
        spdm =  min(fac*spdm,einit)
!        spdm = einit
        prsgas(ip) = prsgas(ip-1) * exp(grav/(spdm*gammi1) &
             & *(1.d0/radius(ip) - 1.d0/radius(ip-1) ))
        speene(ip) = 1.d0/gammi1
        densit(ip) = prsgas(ip)/(spdm*gammi1)
     else
        speene(ip) = 1.d0/gammi1
!            densit(ip) = dinit * exp(grav*(1./radius(ip) - 1.) )
        densit(ip) = dinit * exp(1.d0*grav*(1.d0/radius(ip) - 1.d0) )
        prsgas(ip) = press / dinit * densit(ip)
        spdm = speene(ip)
     endif
!     write(*,*)ip,radius(ip),dradiu(ip),densit(ip)
!         if(Alfven(ip-1).ge.2000.*sounds(ip-1) .or. iddns.eq.2)then
!            iddns=2
!            densit(ip) = densit(ip-1) *(radius(ip-1)/radius(ip))**(2)
!*            dxdt(ip) = dfloat(mod(ip,100))*.01
!*            densit(ip) = densit(ip) * (1.+dxdt(ip) )
!            prsgas(ip) = press / dinit * densit(ip) !*radius(ip)**(-200)
!c            prsgas(ip) = press * exp(grav*(1./radius(ip) - 1.) )
!         else
!             densit(ip) = dinit * exp(grav*(1./radius(ip) - 1.) )
!             prsgas(ip) = press / dinit * densit(ip)
!c             write(*,*)ip,densit(ip),prsgas(ip)
!          endif
     volume(ip) = 1.d0/densit(ip)
!         speene(ip) = prsgas(ip)*volume(ip)/gammi1
     sounds(ip) = sqrt(gamma*prsgas(ip)*volume(ip))
     clocal(ip) = dtmin
     dxlagr(ip) = dxplac(ip) * densit(ip)
     xlagra(ip) = xlagra(ip-1) + .5*(dxlagr(ip-1) + dxlagr(ip))
!     -------------------
     Bparal(ip) = Bplint/Aflow
     Alfven(ip) = fc1*Bparal(ip)/sqrt(densit(ip))
!         ran1 = ran1*102.+25.
     ettrb(ip) = 0.025d0*rdkm / dslmbd0 / sqrt(Aflow)
     ivtr = mod(ivtr * iran1 + iran2, iran3)
!         write(*,*)ip,ran1
!         vtrans(ip) = 2.2*(2.*dfloat(ivtr)/dfloat(iran3) -1.) 
!     &        / (densit(ip)**.2) 
!     vtrans(ip) =0.d0* sin(omgtr/Alfven(ip)*(1.d0-radius(ip))) &
!            &   / (densit(ip)**.1d0)
     vtrans(ip) = 0.d0
!     &        *(dfloat(mod(ip,2)) - .5)*2.
!         vtrans(ip) = 0.
     if(ip >= np)vtrans(ip) = 0.d0
!         Btrans(ip) = 0.
     Btrans(ip) = -vtrans(ip)*densit(ip)
!     -------------------
!     vtrant(ip) = 0.d0*(2.d0*dble(ivtr)/dble(iran3) -1.d0) &
!          & / (densit(ip)**.1d0)
     vtrant(ip) = 0.
     if(ip >= np)vtrant(ip) = 0.d0
!         Btrant(ip) = 0.
     Btrans(ip) = -vtrant(ip)*densit(ip)
!----------------------------
     Bprs = fc3 * (Btrans(ip)*Btrans(ip) + Btrant(ip)*Btrant(ip))
     pressu(ip) = prsgas(ip) + Bprs
     energy(ip) = speene(ip)+.5d0 * (dxdt(ip)*dxdt(ip) &
          &  + vtrans(ip)*vtrans(ip) + vtrant(ip)*vtrant(ip) ) &
          &  + Bprs / densit(ip) - grav / radius(ip)
!         write(*,*)ip,energy(ip),dxdt(ip),vtrans(ip),vtrant(ip)
!     &        ,Bprs,Bparal(ip),densit(ip),speene(ip),grav/radius(ip)
!         if(ip.eq.2)write(*,*)ip,Bparal(ip),Bplint,Aflow,rnow,fltb(rnow)
     velmax = sqrt ( sounds(ip)*sounds(ip) + fc2/densit(ip)  &
          & * (Bparal(ip)*Bparal(ip) + Btrans(ip)*Btrans(ip) & 
           & + Btrant(ip)*Btrant(ip) ) ) + abs(dxdt(ip))
     dttmp = dradiu(ip)/velmax
     dtmin = min(dtmin,dttmp)
  end do
!  stop
  dt    = dtmin*cflnum
  
!     ------------------------------------------------------------------
  write(*,'(/,4x,a,i8)')  'Total Number of Grid points =',np
!      write(11,'(1pe15.7,1p5e11.3)')
!     &                      ( xplace(i),dxdt(i),densit(i)
!     &                       ,energy(i),pressu(i),speene(i),i=1,np )
  write(11,'(1pe15.7,1p13e15.7)') &
       &  ( radius(i),xlagra(i),dxdt(i),densit(i),energy(i),prsgas(i) &
       & ,speene(i),Btrans(i),vtrans(i),Btrant(i),vtrant(i),dradiu(i) &
       & ,dxdtbd(i),Alfven(i),i=1,np )
!   ------------------------------------------------------------------
      
end subroutine initia
!======================================================================*
!            Prepare an Initial State by reading file                  *
!======================================================================*
subroutine initiaread
  use nrtype
  USE vars
  implicit none
  integer(I4B) :: ip,ivtr,iddns
  real(DP) :: dr,volum,press,uinit,sound
  real(DP) :: rnow,Aflow,Bprs,velmax,dtmin,drad,rad,fltb,dttmp &
       & ,r1,r2
  real(DP) :: dum2,dum3,dum4,dum10,dum11
  integer(I4B) :: i
!  real(SP) :: t_real
!  real(SP),allocatable :: varread(:,:)
  real(DP),allocatable :: varread(:,:)

  open(unit=14,file='Initial_index.dat')
  read(14,*) irstr
  close(14)
  
  read(12,rec=2)tprnt0

  allocate(varread(nd,nmvar))
!  read(12,rec=irstr)t_real,varread
!  t = dble(t_real)
!  radius(1:nd) = dble(varread(1:nd,1))
!  dradiu(1:nd) = dble(varread(1:nd,2))
!  densit(1:nd) = dble(varread(1:nd,3))
!  dxdt(1:nd) = dble(varread(1:nd,4))
!  vtrans(1:nd) = dble(varread(1:nd,5)) 
!  vtrant(1:nd) = dble(varread(1:nd,6))
!  Btrans(1:nd) = dble(varread(1:nd,7)) 
!  Btrant(1:nd) = dble(varread(1:nd,8))
!  speene(1:nd) = dble(varread(1:nd,9))
  read(12,rec=irstr)t,varread
  radius(1:nd) = varread(1:nd,1)
  dradiu(1:nd) = varread(1:nd,2)
  densit(1:nd) = varread(1:nd,3)
  dxdt(1:nd) = varread(1:nd,4)
  vtrans(1:nd) = varread(1:nd,5)
  vtrant(1:nd) = varread(1:nd,6)
  Btrans(1:nd) = varread(1:nd,7) 
  Btrant(1:nd) = varread(1:nd,8)
  speene(1:nd) = varread(1:nd,9)
  dxdtbd(1:nd) = varread(1:nd,10)

  irstr = irstr - 1
  
!  write(*,*)dradiu(1:4)
! ------------------------------------------------------------------
  dr = dradiu(1)
! -------------------------------------
  volum = 1.0d0/dinit
!      press = pinit + Bpre 
  press = pinit
  uinit = pinit*volum/gammi1
  sound = sqrt( gamma*pinit*volum )

!  dvelo = ddens*sound/(dinit + 2.d0*ddens) ! dvelo is input
  ddens = dvelo * dinit / max((sound - 2.d0*dvelo),0.1d0) 
  dpres = gamma*ddens*press/dinit ! p perturbation
      
  dvtrit = fc1*dBtrit/sqrt(dinit)
  
  dvttit = fc1*dBttit/sqrt(dinit)
! ------------ ghost mesh at ip = 0----------------------
  ip = 1
  drad = dr
  rad = 1.d0 - drad
  Aflow = rad * rad *fltb(rad)
  r1 = radius(ip) - 1.5d0 * drad
  r2 = radius(ip) - 0.5d0 * drad
!  dxplac0 = drad * Aflow 
  dxplac0 = 0.5d0 * dradiu(ip) *(r1*r1*fltb(r1) + r2*r2*fltb(r2))
!  dxplac0 = onovth * (r2*r2*r2*fltb(r2) - r1*r1*r1*fltb(r1))  
  densit0 =  dinit * exp(grav*(1.d0/rad - 1.d0) )
  pressu0 = press / dinit * densit0
  dxdt0 = 0.d0
  speene0 = pinit / (densit0*gammi1)
  Btrans0 = 0.d0
  vtrans0 = 0.d0
  Btrant0 = 0.d0
  vtrant0 = 0.d0
  energy0 = speene0 + .5d0*( dxdt0*dxdt0 + vtrans0*vtrans0 &
       & + vtrant0*vtrant0 ) + fc3*( Btrans0*Btrans0 &
       & + Btrant0*Btrant0 ) / densit0   - grav / rad
!-----------------------------------------------------------
  dtmin = 100.d0

  do ip=1,np+1
!     read(16,*) dum10,xlagra(ip),dxdt(ip),densit(ip) &
!          &  ,dum2,prsgas(ip),speene(ip),Btrans(ip),vtrans(ip) &
!          &  ,Btrant(ip),vtrant(ip),dum11,dum3,dum4

     prsgas(ip) = densit(ip) * speene(ip) * gammi1
     rnow = radius(ip)
     Aflow = rnow*rnow*fltb(rnow)
     r1 = radius(ip) - 0.5d0*dradiu(ip)
     r2 = radius(ip) + 0.5d0*dradiu(ip)
!     dxplac(ip) = dradiu(ip) * Aflow
     dxplac(ip) = 0.5d0 * dradiu(ip) *(r1*r1*fltb(r1) + r2*r2*fltb(r2))
!     dxplac(ip) = onovth*(r2*r2*r2*fltb(r2) - r1*r1*r1*fltb(r1))
     if(ip==1)then
        xplace(ip) = 1.d0
     else
        xplace(ip) = xplace(ip-1) + .5d0*(dxplac(ip-1)+dxplac(ip))
     endif
     radori(ip) = radius(ip)
     drador(ip) = dradiu(ip)
     xplori(ip) = xplace(ip)
     dxplor(ip) = dxplac(ip) 
     radbnd(ip) = radius(ip) + .5d0*dradiu(ip)
     rsqfbnd(ip) = radbnd(ip)*radbnd(ip) * fltb(radbnd(ip))
     volume(ip) = 1.d0/densit(ip)
     ddxdtt(ip) = 0.0d0
     sounds(ip) = sqrt(gamma*prsgas(ip)*volume(ip))
     dxlagr(ip) = dxplac(ip) * densit(ip) 
!     -------------------
     Bparal(ip) = Bplint/Aflow
     Alfven(ip) = fc1*Bparal(ip)/sqrt(densit(ip))
     ettrb(ip) = 0.025d0*rdkm / dslmbd0 / sqrt(Aflow)
     Bprs = fc3 * (Btrans(ip)*Btrans(ip) + Btrant(ip)*Btrant(ip))
     pressu(ip) = prsgas(ip) + Bprs
     energy(ip) = speene(ip)+.5d0 * (dxdt(ip)*dxdt(ip) &
          & +vtrans(ip)*vtrans(ip) + vtrant(ip)*vtrant(ip) ) &
          & + Bprs / densit(ip) - grav / radius(ip)
!         write(*,*)ip,energy(ip)
     velmax = sqrt ( sounds(ip)*sounds(ip) + fc2/densit(ip)  &
          & * (Bparal(ip)*Bparal(ip) + Btrans(ip)*Btrans(ip) &
          & + Btrant(ip)*Btrant(ip) ) ) + abs(dxdt(ip)) 
     dttmp = dradiu(ip)/velmax
     dtmin = min(dtmin,dttmp)
     clocal(ip) = dttmp
!      if(ip.eq.4)write(*,*)sounds(ip),densit(ip),Bparal(ip),Btrans(ip) &
!           & ,velmax,Aflow,rnow
  enddo
!      speene(np+1) = speene(np) ! for conduction
!      densit(np+1) = densit(np)

  dt    = dtmin*cflnum
!  if(dt>2.d-7)dt = 2.d-7
!      write(*,*)dt
! ------------------------------------------------------------------
  write(*,'(/,4x,a,i8)')  'Total Number of Grid points =',np
  write(11,'(1pe15.7,1p13e15.7)') &
       & ( radius(i),xlagra(i),dxdt(i),densit(i),energy(i),prsgas(i) &
       & ,speene(i),Btrans(i),vtrans(i),Btrant(i),vtrant(i),dradiu(i) &
       & ,dxdtbd(i),Alfven(i),i=1,np )

end subroutine initiaread

!======================================================================*
subroutine initiaread_rdexpnd
  use nrtype
  USE vars
  implicit none
  integer(I4B) :: ip,ivtr,iddns
  real(DP) :: dr,volum,press,uinit,sound
  real(DP) :: rnow,Aflow,Bprs,velmax,dtmin,drad,rad,fltb,dttmp &
       & ,r1,r2
  real(DP) :: dum2,dum3,dum4,dum10,dum11
  integer(I4B) :: i,irec
!  real(SP) :: t_real
!  real(SP),allocatable :: varread(:,:)
  real(DP),allocatable :: varread(:,:)

  open(unit=14,file='Initial_index.dat')
  read(14,*) irstr
  close(14)
  
  allocate(varread(nd_ori,nmvar))
!  read(12,rec=irstr)t_real,varread
!  t = dble(t_real)
!  radius(1:nd) = dble(varread(1:nd,1))
!  dradiu(1:nd) = dble(varread(1:nd,2))
!  densit(1:nd) = dble(varread(1:nd,3))
!  dxdt(1:nd) = dble(varread(1:nd,4))
!  vtrans(1:nd) = dble(varread(1:nd,5)) 
!  vtrant(1:nd) = dble(varread(1:nd,6))
!  Btrans(1:nd) = dble(varread(1:nd,7)) 
!  Btrant(1:nd) = dble(varread(1:nd,8))
!  speene(1:nd) = dble(varread(1:nd,9))
  read(12,rec=irstr)t,varread

  close(12)

  radius(1:nd_ori) = varread(1:nd_ori,1)
  dradiu(1:nd_ori) = varread(1:nd_ori,2)
  densit(1:nd_ori) = varread(1:nd_ori,3)
  dxdt(1:nd_ori) = varread(1:nd_ori,4)
  vtrans(1:nd_ori) = varread(1:nd_ori,5)
  vtrant(1:nd_ori) = varread(1:nd_ori,6)
  Btrans(1:nd_ori) = varread(1:nd_ori,7) 
  Btrant(1:nd_ori) = varread(1:nd_ori,8)
  speene(1:nd_ori) = varread(1:nd_ori,9)
  dxdtbd(1:nd_ori) = varread(1:nd_ori,10)

  deallocate(varread)

  inquire(iolength=irec)t,radius,dradiu,densit,dxdt,vtrans,vtrant &
       & ,Btrans,Btrant,speene(1:nd),Tpprt,real(Alfven),etaprt,eadprt,muprt &
       & ,real(densav),real(dxdtav),real(vtraav),real(vtrtav),real(Btraav) &
       & ,real(Btrtav),real(spenav),real(dvbaav),real(dvbtav)
  open(unit=12,file='printvars_rdexpnd.dat',access='direct',form='unformatted' &
       & ,recl = irec)
!       & ,recl=icmpl*(2+(2*nmvar+nmave)*nd))

!  irstr = irstr - 1
  irstr = 0

! ------------------------------------------------------------------
  dr = dradiu(1)
! -------------------------------------
  volum = 1.0d0/dinit
!      press = pinit + Bpre 
  press = pinit
  uinit = pinit*volum/gammi1
  sound = sqrt( gamma*pinit*volum )

!  dvelo = ddens*sound/(dinit + 2.d0*ddens) ! dvelo is input
  ddens = dvelo * dinit / max((sound - 2.d0*dvelo),0.1d0) 
  dpres = gamma*ddens*press/dinit ! p perturbation
      
  dvtrit = fc1*dBtrit/sqrt(dinit)
  
  dvttit = fc1*dBttit/sqrt(dinit)
! ------------ ghost mesh at ip = 0----------------------
  ip = 1
  drad = dr
  rad = 1.d0 - drad
  Aflow = rad * rad *fltb(rad)
  r1 = radius(ip) - 1.5d0 * drad
  r2 = radius(ip) - 0.5d0 * drad
  write(*,*)t
!  dxplac0 = drad * Aflow 
  dxplac0 = 0.5d0 * dradiu(ip) *(r1*r1*fltb(r1) + r2*r2*fltb(r2))
!  dxplac0 = onovth * (r2*r2*r2*fltb(r2) - r1*r1*r1*fltb(r1))  
  write(*,*)r1,r2
  densit0 =  dinit * exp(grav*(1.d0/rad - 1.d0) )
  pressu0 = press / dinit * densit0
  dxdt0 = 0.d0
  speene0 = pinit / (densit0*gammi1)
  Btrans0 = 0.d0
  vtrans0 = 0.d0
  Btrant0 = 0.d0
  vtrant0 = 0.d0
  energy0 = speene0 + .5d0*( dxdt0*dxdt0 + vtrans0*vtrans0 &
       & + vtrant0*vtrant0 ) + fc3*( Btrans0*Btrans0 &
       & + Btrant0*Btrant0 ) / densit0   - grav / rad
!-----------------------------------------------------------
  write(*,*)t
  dtmin = 100.d0

  do ip=1,np+1
!     read(16,*) dum10,xlagra(ip),dxdt(ip),densit(ip) &
!          &  ,dum2,prsgas(ip),speene(ip),Btrans(ip),vtrans(ip) &
!          &  ,Btrant(ip),vtrant(ip),dum11,dum3,dum4

     if(ip>np_ori-100)then
!        if(ip<=np-1000)then
        if(ip<=np_ori)then
           dradiu(ip) = dradiu(ip-1)
        else
!           dradiu(ip) = 1.01d0*dradiu(ip-1)
!           dradiu(ip) = dradiu(np-1000) * (1.0d0 + 999.0d0*(tanh(0.01d0 &
!                & *dble(ip-(np-200))) + 1.d0))
!           dradiu(ip) = dradiu(np-1000) * (1.0d0 + 1.d-6 &
!                & *(dble(ip-np+1000)**3.d0 ))
           if(ip>np-500)then
              dradiu(ip) = 1.01d0*dradiu(ip-1)
           else
              dradiu(ip) = dradiu(np_ori) * (1.0d0 + 0.001d0 &
                   & * (dble(ip-np_ori))**2.d0 / dble(np-np_ori))
           end if
        end if
        radius(ip) = radius(ip-1) + 0.5d0*(dradiu(ip-1) + dradiu(ip))
     end if
     rnow = radius(ip)
     Aflow = rnow*rnow*fltb(rnow)

     r1 = radius(ip) - 0.5d0*dradiu(ip)
     r2 = radius(ip) + 0.5d0*dradiu(ip)
!     dxplac(ip) = dradiu(ip) * Aflow
     dxplac(ip) = 0.5d0 * dradiu(ip) *(r1*r1*fltb(r1) + r2*r2*fltb(r2))
!     dxplac(ip) = onovth*(r2*r2*r2*fltb(r2) - r1*r1*r1*fltb(r1))
     if(ip==1)then
        xplace(ip) = 1.d0
     else
        xplace(ip) = xplace(ip-1) + .5d0*(dxplac(ip-1)+dxplac(ip))
     endif
     radori(ip) = radius(ip)
     drador(ip) = dradiu(ip)
     xplori(ip) = xplace(ip)
     dxplor(ip) = dxplac(ip) 
     radbnd(ip) = radius(ip) + .5d0*dradiu(ip)
     rsqfbnd(ip) = radbnd(ip)*radbnd(ip) * fltb(radbnd(ip))

     if(ip>np_ori-100)then
        densit(ip) = densit(np_ori-100) * radius(np_ori-100)*radius(np_ori-100)&
             & /(rnow*rnow)
        dxdt(ip) = dxdt(np_ori-100)
        dxdtbd(ip) = dxdt(np_ori-100)
        speene(ip) = speene(np_ori-100)
        btrans(ip) = 0.0d0
        btrant(ip) = 0.0d0
        vtrans(ip) = 0.0d0
        vtrant(ip) = 0.0d0
     end if

     prsgas(ip) = densit(ip) * speene(ip) * gammi1
     volume(ip) = 1.d0/densit(ip)
     ddxdtt(ip) = 0.0d0
     sounds(ip) = sqrt(gamma*prsgas(ip)*volume(ip))
     dxlagr(ip) = dxplac(ip) * densit(ip) 
!     -------------------
     Bparal(ip) = Bplint/Aflow
     Alfven(ip) = fc1*Bparal(ip)/sqrt(densit(ip))
     Bprs = fc3 * (Btrans(ip)*Btrans(ip) + Btrant(ip)*Btrant(ip))
     pressu(ip) = prsgas(ip) + Bprs
     energy(ip) = speene(ip)+.5d0 * (dxdt(ip)*dxdt(ip) &
          & +vtrans(ip)*vtrans(ip) + vtrant(ip)*vtrant(ip) ) &
          & + Bprs / densit(ip) - grav / radius(ip)
!         write(*,*)ip,energy(ip)
     velmax = sqrt ( sounds(ip)*sounds(ip) + fc2/densit(ip)  &
          & * (Bparal(ip)*Bparal(ip) + Btrans(ip)*Btrans(ip) &
          & + Btrant(ip)*Btrant(ip) ) ) + abs(dxdt(ip)) 
     dttmp = dradiu(ip)/velmax
     dtmin = min(dtmin,dttmp)
     clocal(ip) = dttmp
!      if(ip.eq.4)write(*,*)sounds(ip),densit(ip),Bparal(ip),Btrans(ip) &
!           & ,velmax,Aflow,rnow
  enddo
!      speene(np+1) = speene(np) ! for conduction
!      densit(np+1) = densit(np)

  dt    = dtmin*cflnum
!  if(dt>2.d-7)dt = 2.d-7
!      write(*,*)dt
! ------------------------------------------------------------------
  write(*,'(/,4x,a,i8)')  'Total Number of Grid points =',np
  write(11,'(1pe15.7,1p13e15.7)') &
       & ( radius(i),xlagra(i),dxdt(i),densit(i),energy(i),prsgas(i) &
       & ,speene(i),Btrans(i),vtrans(i),Btrant(i),vtrant(i),dradiu(i) &
       & ,dxdtbd(i),Alfven(i),i=1,np )

end subroutine initiaread_rdexpnd


!======================================================================*
!                 Right Hand Sides of Euler's Equations                *
!======================================================================*
subroutine rhs
  use nrtype
  USE vars
  implicit none

  integer(I4B) :: ip
  real(DP) :: speedm,speedm1,depend,rad1,rad2,d1,d2,p1,p2,v1,v2,plf,prt &
       & ,btr1,btr2,btt1,btt2,b1,b2,pa,va,da,db,dxdt1,dxdt2,press1,press2 &
       & ,r1,r2,f1,f2,A1,A2,Bpl1,Bpl2,prgas1,prgas2,btrav,bttav &
       & ,rmid,fmid,Amid,f,fltb,df
  real(DP) :: domain,radeg1,rsqf1,veg1,vlf,radeg2,rsqf2,veg2,vrt
! ------------------------------------------------------------------
  ip = 1
  ipcalc = ip
  speedm = sqrt(sounds(ip)*sounds(ip) + fc2/densit(ip) &
       & *(Btrans(ip)*Btrans(ip) + Btrant(ip)*Btrant(ip) &
       & + Bparal(ip)*Bparal(ip) ) )
  domain = speedm * dt
  depend = ( dradiu(ip) - domain )*0.5d0 !*facdep 
  rad1 = radius(ip) + depend
  radeg1 = rad1 - 0.5d0*depend
  depend = depend * facdep * 0.5d0 * (radeg1*radeg1*fltb(radeg1) &
       & + rsqfbnd(ip)) !* rad1 * rad1 * fltb(rad1)
!   --------------------------------------
  d1     = densit(ip) + depend*gradd(ip)
  p1     = pressu(ip) + depend*gradp(ip)
  v1     =   dxdt(ip) + depend*gradv(ip)
  btr1   = Btrans(ip) + depend*dBtran(ip)
  btt1   = Btrant(ip) + depend*dBtrat(ip)
  b1     = sqrt(btr1*btr1 + btt1*btt1)
!  --------------------------------------
  plf = p1 - .5d0*d1*speedm*dt*grav/(rad1*rad1) 
                                ! gravity
  if(plf<=0.d0 .or. plf<= .5d0*p1)plf = .5d0*p1
  rsqf1 = rad1 * rad1 * fltb(rad1)
  veg1 = dxdt(ip) + 0.5d0 * dxplac(ip) * gradv(ip)
  vlf = (2.d0 * rsqf1 * v1 - (rsqfbnd(ip) - rsqf1) * veg1) &
       & / (rsqf1 + rsqfbnd(ip)) 

!  write(*,*)ip,depend,speedm,dt,radius(ip:ip+1)
!-----------------------------------------------------------
  speedm = sqrt(sounds(ip+1)*sounds(ip+1)+fc2/densit(ip+1) &
       &  * (Btrans(ip+1)*Btrans(ip+1) + Btrant(ip+1)*Btrant(ip+1) &
       &  + Bparal(ip+1)*Bparal(ip+1) ) )
  domain = speedm*dt
  depend = ( dradiu(ip+1) - domain )*0.5d0 !*facdep 
  rad2 = radius(ip+1) - depend
  radeg2 = rad2 + 0.5d0*depend
  depend = depend * facdep * 0.5d0 * (radeg2*radeg2*fltb(radeg2) &
       & + rsqfbnd(ip) ) !* rad2 * rad2 * fltb(rad2)
!  ------------------------------------------
  d2     = densit(ip+1) - depend*gradd(ip+1)
  p2     = pressu(ip+1) - depend*gradp(ip+1)
  v2     =   dxdt(ip+1) - depend*gradv(ip+1)
  btr2   = Btrans(ip+1) - depend*dBtran(ip+1)
  btt2   = Btrant(ip+1) - depend*dBtrat(ip+1)
  b2     = sqrt(btr2*btr2 + btt2*btt2)
!   ------------------------------------------
  prt = p2 + .5d0*d2*speedm*dt*grav / (rad2*rad2) ! gravity
  if(prt>=1.5d0 * p2)prt = 1.5d0*p2
  rsqf2 = rad2 * rad2 * fltb(rad2)
  veg2 = dxdt(ip+1) - 0.5d0 * dxplac(ip+1) * gradv(ip+1)
  vrt = (2.d0 * rsqf2 * v2 + (rsqf2 - rsqfbnd(ip)) * veg2) &
       & / (rsqf2 + rsqfbnd(ip))
!---------------------------------------------------------
  call rieman(d1,d2,p1,plf,p2,prt,v1,v2,b1,b2,va,pa,da,db)   !<=======

!  write(*,*)ip,densit(ip)

  press1 = pa
  dxdt1 = va

!          press1 = p1
!          dxdt1 = v1

  r1 = radbnd(ip) + .5d0*dt*dxdt1 ! for spherical coordinate
  f1 = fltb(r1)
  A1 = r1*r1*f1
  Bpl1 = Bplint / A1

!          write(*,*)ip,dxdt1,d1,d2,p1,p2,v1,v2,b1,b2
!          write(*,*)ip,dxdt1,d1,d2,densit(ip),densit(ip+1)
          
  dxdtbd(ip)=dxdt1      !for derive xplace 

!-------- updata B_perp (necessary for longitudinal fast wave) ----          
  Btrrbd(ip)=btr1*da/d1
  Bttrbd(ip)=btt1*da/d1
  Btrlbd(ip+1)=btr2*db/d2
  Bttlbd(ip+1)=btt2*db/d2
!-------------------------------------------------------------------
  prgas1 = press1 - fc3 * ( Btrlbd(ip+1)*Btrlbd(ip+1) &
       & + Bttlbd(ip+1)*Bttlbd(ip+1) )

  speedm1 = speedm

  do ip = 2,np-1

     ipcalc = ip
     speedm = speedm1
     domain = speedm * dt
     depend = ( dradiu(ip) - domain )*0.5d0 !*facdep 
     rad1 = radius(ip) + depend
     radeg1 = rad1 - 0.5d0*depend
     depend = depend * facdep * 0.5d0 * (radeg1*radeg1*fltb(radeg1) &
          & + rsqfbnd(ip)) !* rad1 * rad1 * fltb(rad1)
!     depend = depend * facdep * rad1 * rad1 * fltb(rad1)
!     --------------------------------------
     d1     = densit(ip) + depend*gradd(ip)
     p1     = pressu(ip) + depend*gradp(ip)
     v1     =   dxdt(ip) + depend*gradv(ip)
     btr1   = Btrans(ip) + depend*dBtran(ip)
     btt1   = Btrant(ip) + depend*dBtrat(ip)
     b1     = sqrt(btr1*btr1 + btt1*btt1)
!     --------------------------------------
     plf = p1 - .5d0*d1*speedm*dt*grav/(rad1*rad1) 
           ! gravity
     if(plf<=0.d0 .or. plf<= .5d0*p1)plf = .5d0*p1
     rsqf1 = rad1 * rad1 * fltb(rad1)
     veg1 = dxdt(ip) + 0.5d0 * dxplac(ip) * gradv(ip)
     vlf = (2.d0 * rsqf1 * v1 - (rsqfbnd(ip) - rsqf1) * veg1) &
          & / (rsqf1 + rsqfbnd(ip)) 

!----------------------------------------------------------
     speedm = sqrt(sounds(ip+1)*sounds(ip+1)+fc2/densit(ip+1) &
          & * (Btrans(ip+1)*Btrans(ip+1) + Btrant(ip+1)*Btrant(ip+1) &
          & + Bparal(ip+1)*Bparal(ip+1) ))
     domain = speedm*dt
     depend = ( dradiu(ip+1) - domain )*0.5d0 !*facdep 
     rad2 = radius(ip+1) - depend
     radeg2 = rad2 + 0.5d0*depend
     depend = depend * facdep * 0.5d0 * (radeg2*radeg2*fltb(radeg2) &
          & + rsqfbnd(ip) ) !* rad2 * rad2 * fltb(rad2)
!     depend = depend * facdep * rad2 * rad2 * fltb(rad2)
!     ------------------------------------------
     d2     = densit(ip+1) - depend*gradd(ip+1)
     p2     = pressu(ip+1) - depend*gradp(ip+1)
     v2     =   dxdt(ip+1) - depend*gradv(ip+1)
     btr2   = Btrans(ip+1) - depend*dBtran(ip+1)
     btt2   = Btrant(ip+1) - depend*dBtrat(ip+1)
     b2     = sqrt(btr2*btr2 + btt2*btt2)
!     ------------------------------------------
     prt = p2 + .5d0*d2*speedm*dt*grav / (rad2*rad2) ! gravity
     if(prt>=1.5d0 * p2)prt = 1.5d0*p2
     rsqf2 = rad2 * rad2 * fltb(rad2)
     veg2 = dxdt(ip+1) - 0.5d0 * dxplac(ip+1) * gradv(ip+1)
     vrt = (2.d0 * rsqf2 * v2 + (rsqf2 - rsqfbnd(ip)) * veg2) &
          & / (rsqf2 + rsqfbnd(ip))
!----------------------------------------------------------
     if(ip==(np-1))then ! for outgoing boundary at the right
        if(dxdt1>0.d0)then
           dnout1=da
           pgout1=prgas2
        else
           dnout1=db
           pgout1=prgas1
        endif
        vxout1=dxdt1
        ptout1=press1
     endif

     call rieman(d1,d2,p1,plf,p2,prt,v1,v2,b1,b2,va,pa,da,db)    !<=======

!     if(ip==4)write(*,*)ip,real(depend),real(v1),real(v2),real(vlf),real(vrt)&
!          & ,real(p1),real(p2),real(plf),real(prt),real(va),real(pa) &
!          & ,real(densit(ip)+0.5d0*dxplac(ip)*gradd(ip)),real(densit(ip+1) &
!          & -0.5d0*dxplac(ip+1)*gradd(ip+1))
!          & ,real(speedm),real(sounds(ip+1))
!         if(ip.ge.(np-4))write(*,*)ip,p1,p2,plf,prt,pressu(ip)
!     &        ,pressu(ip+1),gradp(ip),gradp(ip+1)

     press2 = pa
     dxdt2 = va

!-------- updata B_perp (necessary for longitudinal fast wave) ----          
     Btrrbd(ip)=btr1*da/d1
     Bttrbd(ip)=btt1*da/d1
     Btrlbd(ip+1)=btr2*db/d2
     Bttlbd(ip+1)=btt2*db/d2
 ! 'Monotonicity' of B_pp for MOC
     f=0.d0 ! 1 : strong monotonicity; 0 : No Monotonicity
     if(btr1>btr2 .and. Btrrbd(ip)<Btrlbd(ip+1))then
        btrav = .5d0*( Btrrbd(ip) + Btrlbd(ip+1) )
!             Btrrbd(ip) = f*btrav + (1.-f)*Btrrbd(ip)
!             Btrlbd(ip+1) = f*btrav + (1.-f)*Btrlbd(ip+1)
        Btrrbd(ip) = f*btr1 + (1.d0-f)*Btrrbd(ip)
        Btrlbd(ip+1) = f*btr2 + (1.d0-f)*Btrlbd(ip+1)
     else if(btr1<btr2 .and. Btrrbd(ip)>Btrlbd(ip+1))then
        btrav = .5d0*( Btrrbd(ip) + Btrlbd(ip+1) )
!             Btrrbd(ip) = f*btrav + (1.-f)*Btrrbd(ip)
!             Btrlbd(ip+1) = f*btrav + (1.-f)*Btrlbd(ip+1)
        Btrrbd(ip) = f*btr1 + (1.d0-f)*Btrrbd(ip)
        Btrlbd(ip+1) = f*btr2 + (1.d0-f)*Btrlbd(ip+1)
     endif
     if(btt1>btt2 .and. Bttrbd(ip)<Bttlbd(ip+1))then
        bttav = .5d0*( Bttrbd(ip) + Bttlbd(ip+1) )
!             Bttrbd(ip) = f*bttav + (1.-f)*Bttrbd(ip)
!             Bttlbd(ip+1) = f*bttav + (1.-f)*Bttlbd(ip+1)
        Bttrbd(ip) = f*btt1 + (1.d0-f)*Bttrbd(ip)
        Bttlbd(ip+1) = f*btt2 + (1.d0-f)*Bttlbd(ip+1)
     else if(btt1<btt2 .and. Bttrbd(ip)>Bttlbd(ip+1))then
        bttav = .5d0*( Bttrbd(ip) + Bttlbd(ip+1) )
!             Bttrbd(ip) = f*bttav + (1.-f)*Bttrbd(ip)
!             Bttlbd(ip+1) = f*bttav + (1.-f)*Bttlbd(ip+1)
        Bttrbd(ip) = f*btt1 + (1.d0-f)*Bttrbd(ip)
        Bttlbd(ip+1) = f*btt2 + (1.d0-f)*Bttlbd(ip+1)
     endif
!          if(ip.le.35 .and. ip.ge.30)write(*,*)ip,v1,va,v2,btr1
!     &         ,Btrrbd(ip),Btrlbd(ip+1),btr2,d1,da,db,d2
!-------------------------------------------------------------------
     prgas2 = press2 - fc3 * ( Btrrbd(ip)*Btrrbd(ip) + Bttrbd(ip)*Bttrbd(ip) )
!          if(ip.eq.2)write(*,*)ip,depend,speedm,sounds(ip+1)
!     &         ,Btrans(ip+1),Bparal(ip+1)

     r2 = radbnd(ip) + .5d0*dt*dxdt2 ! for spherical coordinate
     f2 = fltb(r2)
     A2 = r2*r2*f2
     Bpl2 = Bplint / A2

     rmid = .5d0 * (r1 + r2)
     fmid = fltb(rmid)
     Amid = rmid*rmid*fmid

!    *****************************************************
!         dvoldt(ip) =   dxlagi*(  dxdt2 -  dxdt1 )
!         ddxdtt(ip) = - dxlagi*( press2 - press1 )
!         denedt(ip) = - dxlagi*( press2*dxdt2 - press1*dxdt1 )
     dvoldt(ip) =   ( A2*dxdt2 - A1*dxdt1 )/dxlagr(ip)
!     ddxdtt(ip) = - Amid*(press2 - press1)/dxlagr(ip)
     ddxdtt(ip) = - 0.5d0*(A2+A1)*(press2 - press1)/dxlagr(ip)
     denedt(ip) = (- ( A2*press2*dxdt2 - A1*press1*dxdt1 ) ) /dxlagr(ip)
!     &        + fc3 * Bplint * (Bpl2*dxdt2 - Bpl1*dxdt1) ) / dxlagr(ip)
!         if(ip.eq.(np-1))then
!            denedt(ip) = denedt(ip) !+ fc3 * Bplint * dxdt1
!     &           *(Bpl2 - Bpl1 ) / dxlagr(ip)            
!         else
!            denedt(ip) = denedt(ip) + fc3 * Bplint * 
!     &           (Bpl2*dxdt2 - Bpl1*dxdt1)  / dxlagr(ip)
!         endif
!         denedt(ip) = - .5*(prgas1+prgas2) * ( A2*dxdt2 - A1*dxdt1 )
!     &               /dxlagr(ip)
!     *****************************************************
!     if(ip==100)write(*,*)ip,dxplac(ip),dxplac(ip)+dt*(A2*dxdt2 - A1*dxdt1) &
!          & ,dxlagr(ip)*volume(ip),dxlagr(ip)*(volume(ip)+dvoldt(ip)*dt)
!         if(ip==4)write(*,*)ip,r2,radbnd(ip),dt,dxdt2 
!              & ,dxlagr(ip)
!         if(ip >= (np-1))write(*,*)ip,dxdt1,dxdt2,press1,press2
!     &        ,p1,p2,plf,prt,d1,d2,speedm,sounds(ip),sounds(ip+1)

     press1 = pa
     dxdt1 = va

     prgas1 = press1 - fc3 * ( Btrlbd(ip+1)*Btrlbd(ip+1) &
          & + Bttlbd(ip+1)*Bttlbd(ip+1) )
!          if(ip.eq.100)write(*,*)ip,A1,A2,Amid,dxlagr(ip)

     if(ip==(np-1))then ! for outgoing boundary at the right
        if(dxdt2>0.d0)then
           dnout2=da
           pgout2=prgas2
        else
           dnout2=db
           pgout2=prgas1
        endif
        vxout2=dxdt2
        ptout2=press2
     endif

     r1 = r2
     f1 = f2
     A1 = A2
     Bpl1 = Bpl2

     dxdtbd(ip)=dxdt1      !for derive xplace 

!        i = ip
!        write(13,'(1p8e11.3)') xplace(i),dxdt1,dxdt2,press1,press2
!    &                         ,dvoldt(i),ddxdtt(i),denedt(i)

     speedm1 = speedm

!      if(ip.le.10)write(*,*)ip,p1,pressu(ip),gradp(ip),speedm,dxplac(ip)
!     &     ,depend,dt
!      write(*,*)ip,dxdt1,press1,d1,d2,p1,p2,v1,v2,densit(ip)

  end do
  
end subroutine rhs
!======================================================================*
!                      Constraint for Monotonicity                     *
!======================================================================*
subroutine monoto(index)
  use nrtype
  USE vars
  implicit none

  integer(I4B),intent(in) :: index
  integer(I4B) :: ip
  real(DP) :: dx1inv,dx2inv,fcrd,agradd,xgradd,agradp,xgradp,agradv,xgradv &
       & ,agrade,xgrade,adBtran,xdBtran,advtran,xdvtran,adBtrat,xdBtrat &
       & ,advtrat,xdvtrat,fac

! index=0 : Usual Monotonicity
! index=1 : Weak Monotonicity for Euler Remapping
!  ------------------------------------------------------------------

  ip = 1
!         dx2inv = 2.0/dxplac(ip)
  dx1inv = 1.0d0/(dxplac(ip) + .5d0*(dxplac0 + dxplac(ip+1)) )
!         dx2inv = 4.0/(dxplac(ip)+dxplac(ip+1)) ! by stakeru
!         dx2inv = 2.0/dxplac(ip) ! by nakano
  dx2inv = min(2.0d0/dxplac(ip),2.0d0/dxplac(ip+1) &
       & ,4.0d0/(dxplac(ip)+dxplac(ip+1)) )
                                ! safer monotonicity
  fcrd = .5d0

  agradd       = ( densit(ip+1) - densit0    )*dx1inv
  xgradd       = ( densit(ip+1) - densit(ip) )*dx2inv
  gradd(ip) = min( abs(xgradd),abs(agradd) )
  gradd(ip) = gradd(ip)*sign(unity,xgradd) * fcrd
!          gradd(ip) = xgradd*fcrd
  gradd(ip+1) = xgradd 
  
  agradp       = ( pressu(ip+1) - pressu0    )*dx1inv
  xgradp       = ( pressu(ip+1) - pressu(ip) )*dx2inv
  gradp(ip) = min( abs(xgradp),abs(agradp) )
  gradp(ip) = gradp(ip)*sign(unity,xgradp) * fcrd
!          gradp(ip) = xgradp*fcrd
  gradp(ip+1) = xgradp

  agradv       = (   dxdt(ip+1) -   dxdt0    )*dx1inv
  xgradv       = (   dxdt(ip+1) -   dxdt(ip) )*dx2inv
  gradv(ip) = min( abs(xgradv),abs(agradv) )
  gradv(ip) = gradv(ip)*sign(unity,xgradv) * fcrd
!          gradv(ip) = xgradv*fcrd
  gradv(ip+1) = xgradv

  agrade       = ( energy(ip+1) - energy0    )*dx1inv
  xgrade       = ( energy(ip+1) - energy(ip) )*dx2inv
  grade(ip) = min( abs(xgrade),abs(agrade) )
  grade(ip) = grade(ip)*sign(unity,xgrade) * fcrd
!          grade(ip) = xgrade*fcrd
  grade(ip+1) = xgrade

  adBtran      = ( Btrans(ip+1) - Btrans0 )*dx1inv 
  xdBtran      = ( Btrans(ip+1) - Btrans(ip) )*dx2inv 
  dBtran(ip) = min( abs(xdBtran),abs(adBtran) )
  dBtran(ip) = dBtran(ip)*sign(unity,xdBtran) * fcrd
  dBtran(ip+1) = xdBtran

  advtran      = ( vtrans(ip+1) - vtrans0 )*dx1inv 
  xdvtran      = ( vtrans(ip+1) - vtrans(ip) )*dx2inv 
  dvtran(ip) = min( abs(xdvtran),abs(advtran) )
  dvtran(ip) = dvtran(ip)*sign(unity,xdvtran) * fcrd
  dvtran(ip+1) = xdvtran
! --- 3rd components -----------
  adBtrat      = ( Btrant(ip+1) - Btrant0 )*dx1inv 
  xdBtrat      = ( Btrant(ip+1) - Btrant(ip) )*dx2inv 
  dBtrat(ip) = min( abs(xdBtrat),abs(adBtrat) )
  dBtrat(ip) = dBtrat(ip)*sign(unity,xdBtrat) * fcrd
  dBtrat(ip+1) = xdBtrat

  advtrat      = ( vtrant(ip+1) - vtrant0 )*dx1inv 
  xdvtrat      = ( vtrant(ip+1) - vtrant(ip) )*dx2inv 
  dvtrat(ip) = min( abs(xdvtrat),abs(advtrat) )
  dvtrat(ip) = dvtrat(ip)*sign(unity,xdvtrat) * fcrd
  dvtrat(ip+1) = xdvtrat
! --- end of 3rd components----

  do ip = 2,np-1
!      do 1000 ip = 2,np
                  ! to derive grad at ip = np ; np+1 is ghost, unchanged 
                                !from the beginning 
!         dx1inv = 1.0/( dxplac(ip) + dxplac(ip-1) ) ! by stakeru
     dx1inv = 1.0d0/( dxplac(ip) + .5d0*(dxplac(ip-1)+dxplac(ip+1)) )
                                ! by stakeru 
!         dx2inv = 2.0/  dxplac(ip)
!         dx2inv = 4.0/(dxplac(ip)+dxplac(ip+1)) ! by stakeru
!         dx2inv = 2.0/dxplac(ip) ! by nakano
!         if(ip.eq.(np-1))then
!            dx2inv = dmin1(1.0/dxplac(ip),1.0/dxplac(ip+1)
!     &           ,2.0/(dxplac(ip)+dxplac(ip+1)) )            
!         else
     dx2inv = min(2.0d0/dxplac(ip),2.0d0/dxplac(ip+1) &
          & ,4.0d0/(dxplac(ip)+dxplac(ip+1)) )
                                ! safer monotonicity
     if(ip==(np-1))dx2inv = 0.3d0*dx2inv ! for 1st order outflow boud.
!         endif
!        xxgrad = ( 1.0 - ( dt/clocal(ip) )**3 )
!        xxgrad = 0.5
!        if( clocal(ip) .gt. 8.0*dt )
!    &       xxgrad = 1.0
         
     if(index==1)then          ! for Euler Remap     
        fac=abs(dxlg(ip)/dxlagr(ip))
        if(fac>=1.d0)then
           write(*,*)'fac > 1',fac, ip,dxlg(ip),dxlagr(ip)
           stop
        endif
     endif

     agradd       = ( densit(ip+1) - densit(ip-1) )*dx1inv
     xgradd       = ( densit(ip+1) - densit(ip  ) )*dx2inv
     if( sign(unity,xgradd) .eq. sign(unity,gradd(ip)) ) then
        if(index==1)then          ! for Euler Remap     
!               gradd(ip) = amin1( abs(gradd(ip)/fac)
!     &              ,abs(xgradd/(1.-fac)),abs(agradd) )             
           gradd(ip) = min( abs(gradd(ip)/fac),abs(xgradd/(1.d0-fac)) &
                & ,abs(agradd) )  
        else                ! Unsual Monotonicity
!               gradd(ip) = amin1( abs(gradd(ip)),abs(xgradd)
!     &              ,abs(agradd) )
           gradd(ip) = min( abs(gradd(ip)),abs(xgradd),abs(agradd) )
        endif
        gradd(ip) = gradd(ip)*sign(unity,xgradd)
!    &                   *xxgrad
     else
        gradd(ip) = 0.0d0
     endif
     gradd(ip+1) = xgradd

     if(index==2)then     ! only to derive dxlg(ip) in EULER
        if(ip==np)then
!            return
           exit
        else
!            goto 1001
           cycle
        endif
     endif

     agradp       = ( pressu(ip+1) - pressu(ip-1) )*dx1inv
     xgradp       = ( pressu(ip+1) - pressu(ip  ) )*dx2inv
     if( sign(unity,xgradp) == sign(unity,gradp(ip)) ) then
        if(index==1)then          ! for Euler Remap     
!               gradp(ip) = amin1( abs(gradp(ip)/fac)
!     &              ,abs(xgradp/(1.-fac)),abs(agradp) )
           gradp(ip) = min( abs(gradp(ip)/fac),abs(xgradp/(1.d0-fac)) &
                & ,abs(agradp) )
        else                   ! Unsual Monotonicity
!               gradp(ip) = amin1( abs(gradp(ip)),abs(xgradp)
!     &              ,abs(agradp) )
           gradp(ip) = min( abs(gradp(ip)),abs(xgradp),abs(agradp) )
        endif
        gradp(ip) = gradp(ip)*sign(unity,xgradp)
!    &                   *xxgrad
     else
        gradp(ip) = 0.0d0
     endif
     gradp(ip+1) = xgradp

!         if(ip.eq.610)write(*,*)ip,gradp(ip),dx1inv,dx2inv,pressu(ip-1)
!     &        ,pressu(ip),pressu(ip+1)

     agradv       = (   dxdt(ip+1) -   dxdt(ip-1) )*dx1inv
     xgradv       = (   dxdt(ip+1) -   dxdt(ip  ) )*dx2inv
     if( sign(unity,xgradv) == sign(unity,gradv(ip)) ) then
        if(index==1)then   ! for Euler Remap
!               gradv(ip) = amin1( abs(gradv(ip)/fac)
!     &              ,abs(xgradv/(1.-fac)),abs(agradv) )            
           gradv(ip) = min( abs(gradv(ip)/fac),abs(xgradv/(1.d0-fac)) &
                & ,abs(agradv) ) 
        else                ! Unsual Monotonicity
!               gradv(ip) = amin1( abs(gradv(ip)),abs(xgradv)
!     &              ,abs(agradv) )
           gradv(ip) = min( abs(gradv(ip)),abs(xgradv),abs(agradv) )
        endif
        gradv(ip) = gradv(ip)*sign(unity,xgradv)
!    &                   *xxgrad
     else
        gradv(ip) = 0.0d0
     endif
     gradv(ip+1) = xgradv
!         if(ip.eq.615)write(*,*)ip,gradv(ip),dxdt(ip-1),dxdt(ip)
!     &        ,dxdt(ip+1)

     agrade       = ( energy(ip+1) - energy(ip-1) )*dx1inv
     xgrade       = ( energy(ip+1) - energy(ip  ) )*dx2inv
     if( sign(unity,xgrade) == sign(unity,grade(ip)) ) then
        if(index==1)then   ! for Euler Remap
           grade(ip) = min( abs(grade(ip)/fac),abs(xgrade/(1.d0-fac)) &
                & ,abs(agrade) ) 
        else                ! Unsual Monotonicity
           grade(ip) = min( abs(grade(ip)),abs(xgrade),abs(agrade) )
        endif
        grade(ip) = grade(ip)*sign(unity,xgrade)
!    &                   *xxgrad
     else
        grade(ip) = 0.0d0
     endif
     grade(ip+1) = xgrade

     adBtran     = ( Btrans(ip+1) -  Btrans(ip-1) )*dx1inv
     xdBtran     = ( Btrans(ip+1) -  Btrans(ip  ) )*dx2inv
     if( sign(unity,xdBtran) == sign(unity,dBtran(ip)) ) then
        if(index==1)then   ! for Euler Remap
!               dBtran(ip) = amin1( abs(dBtran(ip)/fac)
!     &              ,abs(xdBtran/(1.-fac)),abs(adBtran) )            
           dBtran(ip) = min( abs(dBtran(ip)/fac),abs(xdBtran/(1.d0-fac)) &
                & ,abs(adBtran) )            
        else                ! Unsual Monotonicity
!               dBtran(ip) = amin1( abs(dBtran(ip)),abs(xdBtran)
!     &              ,abs(adBtran) )
           dBtran(ip) = min( abs(dBtran(ip)),abs(xdBtran),abs(adBtran) )
        endif
        dBtran(ip) = dBtran(ip)*sign(unity,xdBtran)
!    &                   *xxgrad
     else
        dBtran(ip) = 0.0d0
     endif
!         if(ip.eq.100)write(*,*)ip,dx1inv,dx2inv,dBtran(ip),adBtran
!     &        ,xdBtran
     dBtran(ip+1) = xdBtran

     advtran     = ( vtrans(ip+1) -  vtrans(ip-1) )*dx1inv
     xdvtran     = ( vtrans(ip+1) -  vtrans(ip  ) )*dx2inv
     if( sign(unity,xdvtran) == sign(unity,dvtran(ip)) ) then
        if(index==1)then   ! for Euler Remap
!               dvtran(ip) = amin1( abs(dvtran(ip)/fac)
!     &              ,abs(xdvtran/(1.-fac)),abs(advtran) )            
           dvtran(ip) = min( abs(dvtran(ip)/fac),abs(xdvtran/(1.d0-fac)) &
                & ,abs(advtran) )            
        else                ! Unsual Monotonicity
!               dvtran(ip) = amin1( abs(dvtran(ip)),abs(xdvtran)
!     &              ,abs(advtran) )
           dvtran(ip) = min( abs(dvtran(ip)),abs(xdvtran),abs(advtran) )
        endif
        dvtran(ip) = dvtran(ip)*sign(unity,xdvtran)
!    &                   *xxgrad
     else
        dvtran(ip) = 0.0d0
     endif
     dvtran(ip+1) = xdvtran

!---------- 3rd components ------------------------------------
     adBtrat     = ( Btrant(ip+1) -  Btrant(ip-1) )*dx1inv
     xdBtrat     = ( Btrant(ip+1) -  Btrant(ip  ) )*dx2inv
     if( sign(unity,xdBtrat) == sign(unity,dBtrat(ip)) ) then
        if(index==1)then   ! for Euler Remap
!               dBtrat(ip) = amin1( abs(dBtrat(ip)/fac)
!     &              ,abs(xdBtrat/(1.-fac)),abs(adBtrat) )  
           dBtrat(ip) = min( abs(dBtrat(ip)/fac),abs(xdBtrat/(1.d0-fac)) &
                & ,abs(adBtrat) )             
        else                ! Unsual Monotonicity
!               dBtrat(ip) = amin1( abs(dBtrat(ip)),abs(xdBtrat)
!     &              ,abs(adBtrat) )
           dBtrat(ip) = min( abs(dBtrat(ip)),abs(xdBtrat),abs(adBtrat) )
        endif
        dBtrat(ip) = dBtrat(ip)*sign(unity,xdBtrat)
!    &                   *xxgrad
     else
        dBtrat(ip) = 0.0d0
     endif
     dBtrat(ip+1) = xdBtrat

     advtrat     = ( vtrant(ip+1) -  vtrant(ip-1) )*dx1inv
     xdvtrat     = ( vtrant(ip+1) -  vtrant(ip  ) )*dx2inv
     if( sign(unity,xdvtrat) == sign(unity,dvtrat(ip)) ) then
        if(index==1)then   ! for Euler Remap
!               dvtrat(ip) = amin1( abs(dvtrat(ip)/fac)
!     &              ,abs(xdvtrat/(1.-fac)),abs(advtrat) )            
           dvtrat(ip) = min( abs(dvtrat(ip)/fac),abs(xdvtrat/(1.d0-fac)) &
                & ,abs(advtrat) )            
        else                ! Unsual Monotonicity
!               dvtrat(ip) = amin1( abs(dvtrat(ip)),abs(xdvtrat)
!     &              ,abs(advtrat) )
           dvtrat(ip) = min( abs(dvtrat(ip)),abs(xdvtrat),abs(advtrat) )
        endif
        dvtrat(ip) = dvtrat(ip)*sign(unity,xdvtrat)
!    &                   *xxgrad
     else
        dvtrat(ip) = 0.0d0
     endif
     dvtrat(ip+1) = xdvtrat
!--------- end of 3rd components ----------------------------
  end do

end subroutine monoto
!======================================================================*
!                     Calculate Shock Tube Parameter                   *
!======================================================================*
subroutine rieman(d1,d2,p1,plf,p2,prt,u1,u2,b1,b2,ua,pa,da,db)
  use nrtype
  USE vars
  implicit none

  integer(I4B) :: icount
  real(DP),intent(in) :: d1,d2,p1,plf,p2,prt,u1,u2,b1,b2
  real(DP),intent(out) :: ua,pa,da,db
  real(DP) :: bpre,gprt1,gprt2,qmplus,fac,qmplui,qmminu,qmmini,pasold
!  ------------------------------------------------------------------

  pa = p1
!  do icount = 1,nshock
  do icount = 1,10
!         pasovp = pa/p2
!         if( pasovp .lt. critif ) then
!             phipop = gammx2*( 1.0 - pasovp )/( 1.0 - pasovp**gammx1 )
!           else
!             phipop = sqrt( gammp1*pasovp + gammi1 )
!         endif
!         qmplus  = sqrt( p2*d2 )*phipop
     bpre  = fc2*gammm2*b2*b2 
     gprt1 = gammp3*pa + gammm3*prt
     gprt2 = gammi1*pa + gammp1*prt         
     qmplus = .5d0*sqrt(d2*(gprt1 - bpre &
          & + sqrt(gprt2*gprt2-bpre*(2.*gprt1-bpre))))
     if(prt/=p2)then
        fac = (gamma * p2 - bpre)/(gamma * prt - bpre)
        qmplus = qmplus*fac
     endif
     qmplui = 1.0d0/qmplus

!         pasovp = pa/p1
!         if( pasovp .lt. critif ) then
!             phipop = gammx2*( 1.0 - pasovp )/( 1.0 - pasovp**gammx1 )
!           else
!             phipop = sqrt( gammp1*pasovp + gammi1 )
!         endif
!         qmminu = sqrt( p1*d1 )*phipop
     bpre  = fc2*gammm2*b1*b1 
     gprt1 = gammp3*pa + gammm3*plf
     gprt2 = gammi1*pa + gammp1*plf         
     qmminu = .5d0*sqrt(d1*(gprt1 - bpre &
          & + sqrt(gprt2*gprt2-bpre*(2.*gprt1-bpre))))
     if(plf/=p1)then
        fac = (gamma * p1 - bpre)/(gamma * plf - bpre)
        qmminu = qmminu*fac
     endif
     qmmini = 1.0d0/qmminu

     pasold = pa
     pa = qmmini*plf + qmplui*prt + u1  - u2
     pa = pa/( qmplui + qmmini )
     if( pa <= 0.0d0 ) pa = pasold*0.7d0
!         if( abs( pa-pasold ) .lt. critit ) go to 3299
!      if( abs( (pa-pasold)/pa ) <= critit ) go to 3299 ! by sano 
     if( abs( (pa-pasold)/pa ) <= critit )exit
  end do

! 3299 ua = qmminu*u1 + qmplus*u2 + plf - prt
  ua = qmminu*u1 + qmplus*u2 + plf - prt
  ua = ua/( qmplus + qmminu )

!--- determine density from mass contimuity for B_perp -------
  da = 1.d0 / ( (ua - u1)*qmmini + 1.d0/d1 )
  db = 1.d0 / ( (u2 - ua)*qmplui + 1.d0/d2 )
!-------------------------------------------------------------
!      write(*,*)qmplus,qmminu
!      write(*,*)gprt1,gprt2,bpre


end subroutine rieman

!=====================================================================*
!           Updating B_perp (for longitudinal fast mode)              *
!=====================================================================*
! This subroutine derives time-centered B_Perp from Riemann solver. 
! Please do not use this subroutine for Alfven wave propagation. 
!(2nd order cannot be achieved for Alfven wave.) 
subroutine updBpp
  use nrtype
  USE vars
  implicit none

  integer(I4B) :: ip
  real(DP) :: dx1inv,dx2inv,adBtran,xdBtran,pdBtran,adBtrat,xdBtrat,pdBtrat
  
!      write(*,*)'100',Btrans(100),dBtran(100)

  do ip=2,np-1
!         if(ip.le.50.and.ip.ge.45)write(*,*)ip,xplace(ip),Btrans(ip)
!     &        ,Btrlbd(ip),Btrrbd(ip)
     Btrans(ip) = .5d0*(Btrlbd(ip)+Btrrbd(ip))
     Btrant(ip) = .5d0*(Bttlbd(ip)+Bttrbd(ip))
!         if(ip.le.50.and.ip.ge.45)write(*,*)ip,xplace(ip),Btrans(ip)
!     &        ,Btrlbd(ip),Btrrbd(ip)
!         if(ip.le.5)write(*,*)Btrans(ip)
  end do

!---- derive gradient for 2nd order scheme -----------------
  ip = 1
  dx2inv = 4.0d0/(dxplac(ip)+dxplac(ip+1)) ! by stakeru     
  xdBtran      = ( Btrans(ip+1) - Btrans(ip) )*dx2inv 
  dBtran(ip)   = xdBtran 
  dBtran(ip+1) = xdBtran
  xdBtrat      = ( Btrant(ip+1) - Btrant(ip) )*dx2inv 
  dBtrat(ip)   = xdBtrat 
  dBtrat(ip+1) = xdBtrat

  do ip=2,np-1
!         if(ip.eq.100)write(*,*)ip,dBtran(ip)
     dx1inv = 1.0d0/( dxplac(ip) + .5d0*(dxplac(ip-1)+dxplac(ip+1)) )
     dx2inv = 4.0d0/(dxplac(ip)+dxplac(ip+1)) 

     adBtran     = ( Btrans(ip+1) -  Btrans(ip-1) )*dx1inv
     xdBtran     = ( Btrans(ip+1) -  Btrans(ip  ) )*dx2inv
     pdBtran     = ( Btrrbd(ip) - Btrlbd(ip) ) / dxplac(ip)
     if( sign(unity,xdBtran) == sign(unity,dBtran(ip)) ) then
        dBtran(ip) = min( abs(dBtran(ip)),abs(xdBtran),abs(adBtran) &
             & ,abs(pdBtran) )
        dBtran(ip) = dBtran(ip)*sign(unity,xdBtran)
     else
        dBtran(ip) = 0.0d0
     endif
!         if(ip.eq.100)write(*,*)dBtran(ip),adBtran,xdBtran
     dBtran(ip+1) = xdBtran

     adBtrat     = ( Btrant(ip+1) -  Btrant(ip-1) )*dx1inv
     xdBtrat     = ( Btrant(ip+1) -  Btrant(ip  ) )*dx2inv
     pdBtrat     = ( Bttrbd(ip) - Bttlbd(ip) ) / dxplac(ip)
     if( sign(unity,xdBtrat) == sign(unity,dBtrat(ip)) ) then
        dBtrat(ip) = min( abs(dBtrat(ip)),abs(xdBtrat),abs(adBtrat) &
             & ,abs(pdBtrat) )
        dBtrat(ip) = dBtrat(ip)*sign(unity,xdBtrat)
     else
        dBtrat(ip) = 0.0d0
     endif
     dBtrat(ip+1) = xdBtrat
  end do

end subroutine updBpp
      
!=====================================================================*
!           MOC for transverse components                             *
!=====================================================================*
subroutine moc
  use nrtype
  USE vars
  implicit none

  real(DP) :: sq4d1,wid1,wid2,dalfv1,dalfv2,rp,rm,fp,fm,sqfp,sqfm,fmov &
       & ,r1,r2,f1,f2,Aflow1,Aflow2,dalhf1,dalhf2,dens1,dens2,dxdt1,dxdt2 &
       & ,dnssqi1,dnssqi2,r1p,r2p,r1m,r2m,rst1,rst2,Ap,Am,Bp,Bm &
       & ,vfrtra1,vfrtra2,frBtra1,frBtra2,vfrtrt1,vfrtrt2,frBtrt1,frBtrt2 &
       & ,fac1,fac2,rad,fsq,geo
  real(DP) :: fltb
  integer(I4B) :: ip

!-------- Derive B_perp & v_perp at cell boundary & time center ------      
      
  sq4d1 = fc1/sqrt(densit(1)) 
  do ip=1,np-1
! --------- for +characteristics-------------------------------------
     wid1 = dt * Alfven(ip)  
     dalfv1 = .5d0*(dradiu(ip) - wid1) ! at time n 
     rp = radbnd(ip) - .5d0*wid1 ! at time n for +chracteristics
     fp = fltb(rp)
     sqfp = sqrt(fp)
     dalfv1 = dalfv1*rp*rp*fp
     
     fmov = .5d0*dxdtbd(ip)*dt
     r1 = radbnd(ip) + .5d0*(fmov - .5d0*wid1) ! at time n+1/4
     f1 = fltb(r1)
!         fsq1 = sqrt(f1)
     Aflow1 = r1*r1*f1
!         Bprl1 = Bplint / (r1*r1*f1)

     dalhf1 = .5d0*(dradiu(ip) + fmov - .5d0*wid1)*Aflow1 ! at time n + 1/4  
     dens1 = densit(ip) + gradd(ip)*dalhf1 ! place is for n+1/4 
     if(dens1 < 0.d0)dens1=.9d0*abs(densit(ip) ) ! for safety
     dxdt1 = dxdt(ip) + gradv(ip)*dalhf1  ! but at time n
     dnssqi1 = fc1/sqrt(dens1)
     
     r1p = 1.00001d0*r1
     r1m =  .99999d0*r1
     rst1 = .5d0 * dxdt1 / (r1*r1*f1)*(r1p*r1p*fltb(r1p) -r1m*r1m*fltb(r1m) ) &
          & / (r1p - r1m)

!         write(*,*)ip,rst1,r1p,fltb(r1p) 

     Ap = 2.d0/dt + rst1
     Bp = 2.d0/dt - rst1
!         vAlfv1 = Bprl1 * dnssqi1  
!         vatr1 = .5 * vAlfv1 * dt / r1 
!         dnsr1 = dnssqi1 / r1

!         if(ip.le.50)write(*,*)ip,r1,Bprl1,Alfven(ip),vAlfv1


     vfrtra1 = (vtrans(ip) + dvtran(ip)*dalfv1)/(rp*sqfp)
     frBtra1 = sqfp*rp*(Btrans(ip) + dBtran(ip)*dalfv1)
!         Btran1 = Btrans(ip) + dBtran(ip)*dalfv1
!------3rd components ------------------------
     vfrtrt1 = (vtrant(ip) + dvtrat(ip)*dalfv1)/(rp*sqfp)
     frBtrt1 = sqfp*rp*(Btrant(ip) + dBtrat(ip)*dalfv1)
!         Btrat1 = Btrant(ip) + dBtrat(ip)*dalfv1
!------end of 3rd components------------------

!         if(ip.eq.1 .or. ip.eq.(np-1)) then
!            rBtrbd(ip) = rBtran1
!            vftrbd(ip) = vfrtra1
!            rBttbd(ip) = rBtrat1
!            vfttbd(ip) = vfrtrt1
!            goto 51
!         endif

!--------- for -characteristics--------------------------------------
     wid2 = dt * Alfven(ip+1)
     dalfv2 = .5d0*(dradiu(ip+1) - wid2)
     rm = radbnd(ip) + .5d0*wid2 ! at time n for -chracteristics
     fm = fltb(rm)
     sqfm = sqrt(fm)
     dalfv2 = dalfv2*rm*rm*fm

     r2 = radbnd(ip) + .5d0*(fmov + .5*wid2)
     f2 = fltb(r2)
!         fsq2 = sqrt(f2)
     Aflow2 = r2*r2*f2
!         Bprl2 = Bplint / (r2*r2*f2)

     dalhf2 = .5d0*(dradiu(ip+1) - fmov - .5d0*wid2)*Aflow2 ! at time n + 1/4  
     dens2 = densit(ip+1) - gradd(ip+1)*dalhf2 ! place is for n+1/4 
     if(dens2 < 0.d0)dens2=.9*abs(densit(ip+1) ) ! for safety
     dxdt2 = dxdt(ip+1) - gradv(ip+1)*dalhf2  ! but at time n
     dnssqi2 = fc1/sqrt(dens2)

     r2p = 1.00001d0*r2
     r2m =  .99999d0*r2
     rst2 = .5d0 * dxdt2 / (r2*r2*f2)*(r2p*r2p*fltb(r2p) -r2m*r2m*fltb(r2m) )&
          & / (r2p - r2m)

     Am = 2.d0/dt + rst2
     Bm = 2.d0/dt - rst2

!         write(*,*)ip,rst1,rst2

     vfrtra2 = (vtrans(ip+1) - dvtran(ip+1)*dalfv2)/(rm*sqfm)
     frBtra2 = sqfm*rm*(Btrans(ip+1) - dBtran(ip+1)*dalfv2)
!------3rd components ------------------------
     vfrtrt2 = (vtrant(ip+1) - dvtrat(ip+1)*dalfv2)/(rm*sqfm)
     frBtrt2 = sqfm*rm*(Btrant(ip+1) - dBtrat(ip+1)*dalfv2)
!------end of 3rd components------------------

     if(ip==(np-1)) then ! outflow boundary
        frBtra2 = 0.d0
        vfrtra2 = 0.d0
        frBtrt2 = 0.d0
        vfrtrt2 = 0.d0
     endif

!         if(ip.le.5)write(*,*)ip,Alfven(ip),dt,dradiu(ip)
!         if(ip.eq.1)write(*,*)Btrans(1),Btrans(2),Btran1,Btran2

     fac1 = dnssqi1/(r1*r1*f1)
     fac2 = dnssqi2/(r2*r2*f2)
     fBtrbd(ip) = (- vfrtra1*Bp*Am + vfrtra2*Bm*Ap +Ap*Am & 
          & * (frBtra1*fac1 + frBtra2*fac2)) / (Am*Bp*fac1 + Ap*Bm*fac2)
     vftrbd(ip) = (vfrtra1*Bp + fac1*(fBtrbd(ip)*Bp - frBtra1*Ap)) / Ap

     fBttbd(ip) = (- vfrtrt1*Bp*Am + vfrtrt2*Bm*Ap +Ap*Am &
          & * (frBtrt1*fac1 + frBtrt2*fac2)) / (Am*Bp*fac1 + Ap*Bm*fac2)
     vfttbd(ip) = (vfrtrt1*Bp + fac1*(fBttbd(ip)*Bp - frBtrt1*Ap)) / Ap

     if(ip==(np-2))then ! for outflow boundary
        rad = radbnd(ip) + .5d0*dt*dxdtbd(ip)
        fsq = sqrt(fltb(rad))
        geo = rad*fsq
        Btrbd1 = fBtrbd(ip) / geo
        vtrbd1 = vftrbd(ip) * geo
        Bttbd1 = fBttbd(ip) / geo
        vttbd1 = vfttbd(ip) * geo
     else if(ip==(np-1))then
        rad = radbnd(ip) + .5d0*dt*dxdtbd(ip)
        fsq = sqrt(fltb(rad))
        geo = rad*fsq
        Btrbd2 = fBtrbd(ip) / geo
        vtrbd2 = vftrbd(ip) * geo
        Bttbd2 = fBttbd(ip) / geo
        vttbd2 = vfttbd(ip) * geo            
     endif

!         if(ip.eq.70)write(*,*)ip,vfrtra1,vfrtra2,rBtran1,rBtran2
!     &        ,vftrbd(ip),rBtrbd(ip),fm,rm
!     &        ,rBtrbd(ip),Ap,Am,Bp,Bm,fac1,fac2,densit(ip),dens1
!         if(ip.ge.896 .and. ip.le.900)write(*,*)ip,dnssqi1,dens1
!     &        ,densit(ip),gradd(ip),dalhf1
!         tmp1 = .5*(Am*Bp*fac1 + Ap*Bm*fac2)
!         tmp2 = - vfrtra1*Bp*Am /tmp1
!         tmp3 = vfrtra2*Bm*Ap / tmp1
!         tmp4 = rBtran1*fac1*Ap*Am / tmp1
!         tmp5 = rBtran2*fac2*Ap*Am / tmp1
!         if(ip.le.6)write(*,*) ip,-vtrans(ip)/radius(ip)
!     &        ,vtrans(ip+1)/radius(ip+1),Btrans(ip)*radius(ip)
!     &        ,Btrans(ip+1)*radius(ip+1)
!     &        ,-vfrtra1,vfrtra2,rBtran1,rBtran2,rBtrbd(ip)
!         if(ip.le.4)write(*,*) ip,-vfrtra1*Bp*Am,vfrtra2*Bm*Ap
!     &        ,Ap*Am*rBtran1*fac1,Ap*Am*rBtran2*fac2
  end do

end subroutine moc

!///////////////// Latter half of MOC ///////////////////////////////
subroutine trans
  use nrtype
  USE vars
  implicit none

  integer(I4B) :: ip
  real(DP) :: delem1,rav,r1,r2,dr,fav,f1,f2,spexp,centrf,Bcent,fltb
  real(DP) :: Aflow1,Aflow2,etam,etap,eadm,eadp,Bsqttm,Bsqttp
!----------Derive dv_tr/dt, dB_tr/dt*dxplac, & (1/2)dB^2/dr -----------

  if(irss/=0.and.t>t_e2T_rss)call resistivity
!! dxplac is output from Riemann solver.
  do ip=2,np-1
!-- derive d(r v_perp)/dt --------------------------------------
     delem1 = fc2*Bplint/dxlagr(ip)
     dtvrtr(ip) = delem1 * ( fBtrbd(ip) - fBtrbd(ip-1) )
     dtvrtt(ip) = delem1 * ( fBttbd(ip) - fBttbd(ip-1) )
!------------------------------------------------------------
     rav = .5d0*(radius(ip)+radori(ip))
     r1 = rav - .25d0*(dradiu(ip)+drador(ip))
     r2 = rav + .25d0*(dradiu(ip)+drador(ip))
     dr=.5d0*(dradiu(ip) + drador(ip))
     fav = fltb(rav)
     f1 = fltb(r1)
     f2 = fltb(r2)

     spexp = (r2*sqrt(f2) - r1*sqrt(f1)) / (dr * sqrt(fav) )

!         write(*,*)ip,f1
         
     centrf = spexp * .5d0*rav*fav*(vftrbd(ip-1)*vftrbd(ip-1) &
          & + vftrbd(ip)*vftrbd(ip) + vfttbd(ip-1)*vfttbd(ip-1) &
          & + vfttbd(ip)*vfttbd(ip)) 
     Bcent = - spexp * fc2 / (rav*rav*rav*fav*densit(ip)) * .5d0 & 
          & *(fBtrbd(ip-1)*fBtrbd(ip-1) + fBtrbd(ip)*fBtrbd(ip) &
          & + fBttbd(ip-1)*fBttbd(ip-1) + fBttbd(ip)*fBttbd(ip) )

         
     Bcenrf(ip) = centrf + Bcent
!         if(ip.eq.617)write(*,*)ip,Bcenrf(ip),centrf,Bcent,densit(ip)
!     &        ,vftrbd(ip-1),vftrbd(ip),rBtrbd(ip-1),rBtrbd(ip)

     if(irss==0.or.ip==2.or.t<t_e2T_rss)then
        Btense(ip) = fc2 * Bplint * (  ( fBtrbd(ip)*vftrbd(ip) &
             & + fBttbd(ip)*vfttbd(ip) ) -  (fBtrbd(ip-1)*vftrbd(ip-1) &
             & + fBttbd(ip-1)*vfttbd(ip-1) ) ) / dxlagr(ip)
     else
        etam = min(0.5d0 * (eta(ip-1) + eta(ip)),etamax)       
        etap = min(0.5d0 * (eta(ip) + eta(ip+1)),etamax)
        if(irss==1)then
           eadm = 0.d0
           eadp = 0.d0
        else
           Aflow1 = r1 * r1 * f1
           Aflow2 = r2 * r2 * f2
           if(irss==2)then
              Bsqttm = (fBtrbd(ip-1)*fBtrbd(ip-1) + fBttbd(ip-1)*fBttbd(ip-1))&
                   & / Aflow1
              Bsqttp = (fBtrbd(ip)*fBtrbd(ip) + fBttbd(ip)*fBttbd(ip) ) &
                   & / Aflow2
           else
              Bsqttm = (fBtrbd(ip-1)*fBtrbd(ip-1) + fBttbd(ip-1)*fBttbd(ip-1) &
                   & + Bplint*Bplint/ Aflow1) / Aflow1
              Bsqttp = (fBtrbd(ip)*fBtrbd(ip) + fBttbd(ip)*fBttbd(ip) &
                   & + Bplint*Bplint/ Aflow2) / Aflow2
           end if
           eadm = min(0.5d0 * (ead(ip-1) + ead(ip)) * Bsqttm,etamax)
           eadp = min(0.5d0 * (ead(ip) + ead(ip+1)) * Bsqttp,etamax)
        end if
!        write(*,*)ip,eadm,eadp,etam,etap
        Btense(ip) = fc2 * ( Bplint * (  ( fBtrbd(ip)*vftrbd(ip) &
             & + fBttbd(ip)*vfttbd(ip) ) -  (fBtrbd(ip-1)*vftrbd(ip-1) &
             & + fBttbd(ip-1)*vfttbd(ip-1) ) ) + ((etap+eadp)*((fBtrbd(ip) &
             & *(fBtrbd(ip+1) - fBtrbd(ip-1)) + fBttbd(ip) &
             & * (fBttbd(ip+1) - fBttbd(ip-1))) )/ (dradiu(ip) + dradiu(ip+1)) &
             & - (etam+eadm) * ((fBtrbd(ip-1) * (fBtrbd(ip) - fBtrbd(ip-2)) &
             & + fBttbd(ip-1) * (fBttbd(ip) - fBttbd(ip-2))) ) &
             & / (dradiu(ip-1) + dradiu(ip)) )) / dxlagr(ip)
     end if

!         elm1 = fc2*Bplint/dxlagr(ip) 

!         Btense(ip) = centrf * .5*(dxdtbd(ip-1) + dxdtbd(ip))
!     &        + elm1 * .5* ( (rBtrbd(ip-1) + rBtrbd(ip)) 
!     &        * (vftrbd(ip) * sqrt(f2) - vftrbd(ip-1) * sqrt(f1) ) 
!     &        + (rBttbd(ip-1) + rBttbd(ip) ) 
!     &        * ( vfttbd(ip)* sqrt(f2)  - vfttbd(ip-1) * sqrt(f1) ) )  

!         if(ip.le.3)write(*,*)ip,centrf * .5*(dxdtbd(ip-1) + dxdtbd(ip))
!     &        ,elm1 * .5*  (rBtrbd(ip-1) + rBtrbd(ip)) 
!     &        , vftrbd(ip) , sqrt(f2), - vftrbd(ip-1) , sqrt(f1) 
!     &        ,vtrans(ip),vtrans(ip-1)

  end do

end subroutine trans


!=====================================================================*
!           Euler Remapping                                           *
!=====================================================================*
subroutine eulerrmp
  use nrtype
  USE vars
  implicit none

  real(DP) :: dxpl(nd),xlgnew(nd),dxlnew(nd)
  real(DP) :: dxdtmp(nd),engtmp(nd),dnstmp(nd)
  real(DP) :: vtrtmp(nd),vtttmp(nd)
!      real dxpldm(nd)
  real(DP) :: dxpl0,dxlg0,rnow,fcrd,dsav
  real(DP) :: fltb,df,sum1,sum2,sum3
  integer(I4B) :: ip

! Left and Right (ip=1,np) boundaries are fixed.

! derive dX^i-dX_i(=dxpl) & dxi_i-dxi^i(=dxlg) --------------------------

!      ip=1
!      dxpl(ip)=xplace(ip)-xplori(ip)+.5*(dxplac(ip)-dxplor(ip))
! ---------left boundary for wave generation -----------------
  dxpl0 = xplace(1)-xplori(1)-.5*(dxplac(1)-dxplor(1))
! -----------------------------------------------------------
  do ip=1,np-1
     dxpl(ip) = xplace(ip) - xplori(ip) + .5d0*(dxplac(ip)-dxplor(ip))
!     dxpl(ip) = dxdtbd(ip)*dt 
!         write(*,*)ip,dxpl(ip),xplace(ip),xplori(ip),dxplac(ip)
!         if(ip.eq.(np-1))write(*,*)ip,dxpl(ip),xplace(ip),dxplac(ip)
!     &        ,xplori(ip),dxplor(ip)
  end do
!  write(*,*)dxplac(100),radius(100)*radius(100)*dradiu(100)*fltb(radius(100))&
!       &, dxpl(100),radbnd(100)*radbnd(100)*fltb(radbnd(100)) &
!       & *(radori(100)+0.5d0*dradiu(100) - radbnd(100))

  call monoto(2)

! --------- left boundary for wave generation---------------
  if (abs(dxpl0/dxplac(1)) < 1.d-12)then ! No movement 
     dxlg0 = 0.d0
  else
     dxlg0 = dxpl0*( densit(1) -.5d0*gradd(1)*(dxplac(1)+dxpl0) )
  endif
! ----------------------------------------------------------
  do ip=1,np-1
     if (abs(dxpl(ip)/dxplac(ip)) < 1.d-12)then ! No movement 
        dxlg(ip)=0.d0
     else
        rnow = radius(ip)
!            sphf = 1./(rnow*rnow*fltb(rnow))
        if (dxpl(ip) > 0.d0)then  ! boundary =>
           dxlg(ip)=dxpl(ip)*( densit(ip) +.5d0*gradd(ip) &
                & * (dxplac(ip)-dxpl(ip)) )
        else                 ! boundary <=
           dxlg(ip)=dxpl(ip)*( densit(ip+1) +.5d0*gradd(ip+1) &
                & *(-dxplac(ip+1)-dxpl(ip)) )   
        endif
!            if(ip.eq.(np-1))write(*,*)ip,dxlg(ip),dxpl(ip),dxlagr(ip)
!     &           ,dxplac(ip),gradd(ip+1),densit(ip+1)
        if(abs(dxlg(ip)/dxlagr(ip))>=1.d0)then
           write(*,*)'overmove',ip,dxplac(ip),dxplac(ip+1),dxpl(ip) &
                & ,dxlagr(ip),dxlg(ip),densit(ip),densit(ip+1) &
                & ,gradd(ip),gradd(ip+1)
           dxlg(ip)=.99d0*dxlagr(ip)*sign(unity,dxlg(ip))
        endif
     endif
  end do

  ip=1
  xlgnew(ip) = xlagra(ip) - .5d0 *(dxlg0+dxlg(ip))
  dxlnew(ip) = dxlagr(ip) + dxlg0 - dxlg(ip)
  do ip = 2,np-1
     xlgnew(ip) = xlagra(ip) - .5d0*(dxlg(ip-1)+dxlg(ip))
     dxlnew(ip) = dxlagr(ip) + dxlg(ip-1) - dxlg(ip)
!         write(*,*)ip,xlagra(ip),xlgnew(ip),dxlagr(ip),dxlnew(ip)
!     &        ,dxpl(ip),dxlg(ip),xplace(ip),xplori(ip)
!         if(ip.eq.1)write(*,*)ip,dxlnew
  end do

!      write(*,*)densit(506),dxlagr(506),dxlnew(506),dxplor(506)

! Weak Monotonicity for Euler Remap
      
  call monoto(0)
  
! Remapping-----------------------Remapping-----------------------------
! ----- Remapping for left boundary  ---------------------------
  ip = 1
  rnow = radius(ip) 
!         sphf = 1./(rnow * rnow *fltb(rnow)) 
                                ! additional factor for gradient
  if( (abs(dxpl0/dxplor(ip))< 1.d-5) .and. &
       & (abs(dxpl(ip)/dxplor(ip)) < 1.d-5) )then 
     dxdtmp(ip)=dxdt(ip)
     engtmp(ip)=energy(ip)
     dnstmp(ip)=densit(ip)
     vtrtmp(ip)=vtrans(ip)
     vtttmp(ip)=vtrant(ip)
     goto 98
  endif

  fcrd = .5d0 ! reduction factor of slope at ip = 1
  if(dxpl(ip) > 0.d0) then
     dxdtmp(ip) = dxdt(ip) - .5d0 * fcrd * gradv(ip)*Volume(ip) &
          & * (dxlg0 + dxlg(ip))
     engtmp(ip) = energy(ip) - .5d0 * fcrd * grade(ip)*Volume(ip) &
          & * (dxlg0 + dxlg(ip))
     dnstmp(ip) = densit(ip) - .5d0 * fcrd * gradd(ip) &
          & * (dxpl0 + dxpl(ip))
     vtrtmp(ip) = vtrans(ip) - .5d0 * fcrd * dvtran(ip)*Volume(ip) &
          & * (dxlg0 + dxlg(ip))
     vtttmp(ip) = vtrant(ip) - .5d0 * fcrd * dvtrat(ip)*Volume(ip) &
          & * (dxlg0 + dxlg(ip))
  else
     dxdtmp(ip) = ( (dxdt(ip) - .5d0 * fcrd * gradv(ip) * Volume(ip)*dxlg0 ) &
          & *(dxlagr(ip) + dxlg0) + (dxdt(ip+1) - .5d0* gradv(ip+1) &
          & * Volume(ip+1) *(dxlagr(ip+1)+dxlg(ip)) )*(-dxlg(ip)) ) &
          & / dxlnew(ip)

     engtmp(ip) = ( (energy(ip) - .5d0 * fcrd * grade(ip) * Volume(ip)*dxlg0 )&
          & *(dxlagr(ip) + dxlg0) + (energy(ip+1) - .5d0 * grade(ip+1) &
          & * Volume(ip+1) * (dxlagr(ip+1)+dxlg(ip)) )*(-dxlg(ip)) ) &
          & / dxlnew(ip)

     dnstmp(ip) = ( (densit(ip) - .5d0 * fcrd * gradd(ip)*dxpl0 ) &
          & * (dxplac(ip) + dxpl0)  + (densit(ip+1) - .5d0* gradd(ip+1) &
          & *(dxplac(ip+1)+dxpl(ip)) )*(-dxpl(ip)) ) / dxplor(ip)
!         write(*,*)ip,densit(ip),dnstmp(ip),dxplac(ip),dxpl0,dxplor(ip)
!     &        ,densit(ip+1),dxpl(ip)

     vtrtmp(ip) = ( (vtrans(ip) - .5d0 * fcrd * dvtran(ip) * Volume(ip)*dxlg0)&
          & *(dxlagr(ip) + dxlg0) + (vtrans(ip+1) - .5* dvtran(ip+1) &
          & * Volume(ip+1) * (dxlagr(ip+1)+dxlg(ip)) )*(-dxlg(ip)) ) &
          & / dxlnew(ip)

     vtttmp(ip) = ( (vtrant(ip) - .5d0 * fcrd * dvtrat(ip) * Volume(ip)*dxlg0)&
          & * (dxlagr(ip) + dxlg0) + (vtrant(ip+1) - .5d0* dvtrat(ip+1) &
          & * Volume(ip+1) * (dxlagr(ip+1)+dxlg(ip)) )*(-dxlg(ip)) ) &
          & / dxlnew(ip)
  endif

!      write(*,*)ip,densit(ip),densit(ip+1),dnstmp(ip),dxpl0,dxpl(ip)

98 continue
!-----End of remapping of left boundary -----------------------------

!  do ip =2,np-1
  do ip = 2,np

! Remapping is not necessary
!         if (abs((xlnwrt-xlolrt)/dxlagr(ip)).lt.1.e-8 .and. 
!     &        abs((xlnwlf-xlollf)/dxlagr(ip)).lt.1.e-8) goto 99
     if( abs(dxpl(ip-1)/dxplor(ip-1)) < 1.d-12 )then !.and. &
!          &        (abs(dxpl(ip)/dxplor(ip)) < 1.d-12) )then 
        dxdtmp(ip) = dxdt(ip) * dxlagr(ip) 
        engtmp(ip) = energy(ip) * dxlagr(ip)
        dnstmp(ip) = densit(ip) * dxplac(ip)
        vtrtmp(ip) = vtrans(ip) * dxlagr(ip)
        vtttmp(ip) = vtrant(ip) * dxlagr(ip)
 
     else if(dxpl(ip-1) > 0.d0) then ! .and. (dxpl(ip) >= 0.d0))then
        df = ( dxdt(ip-1) + gradv(ip-1)*Volume(ip-1) *.5d0 &
             & *(dxlagr(ip-1)-dxlg(ip-1)) )*dxlg(ip-1)
        if(ip>2)dxdtmp(ip-1) = dxdtmp(ip-1) - df
        if(ip<np)dxdtmp(ip) = dxdt(ip) * dxlagr(ip) + df

        df = ( energy(ip-1) + grade(ip-1)*Volume(ip-1) *.5d0 &
             & *(dxlagr(ip-1)-dxlg(ip-1)) )*dxlg(ip-1)
        if(ip>2)engtmp(ip-1) = engtmp(ip-1) - df
        if(ip<np)engtmp(ip) = energy(ip) * dxlagr(ip) + df

        df = ( vtrans(ip-1) + dvtran(ip-1)*Volume(ip-1) *.5d0 &
             & *(dxlagr(ip-1)-dxlg(ip-1)) )*dxlg(ip-1)
        if(ip>2)vtrtmp(ip-1) = vtrtmp(ip-1) - df
        if(ip<np)vtrtmp(ip) = vtrans(ip) * dxlagr(ip) + df

        df = ( vtrant(ip-1) + dvtrat(ip-1)*Volume(ip-1) *.5d0 &
             & *(dxlagr(ip-1)-dxlg(ip-1)) )*dxlg(ip-1)
        if(ip>2)vtttmp(ip-1) = vtttmp(ip-1) - df
        if(ip<np)vtttmp(ip) = vtrant(ip) * dxlagr(ip) + df

        df = ( densit(ip-1) + gradd(ip-1)*.5d0 &
             & *(dxplac(ip-1)-dxpl(ip-1)) )*dxpl(ip-1)
        if(ip>2)dnstmp(ip-1) = dnstmp(ip-1) - df
        if(ip<np)dnstmp(ip) = densit(ip) * dxplac(ip) + df
        
!        if(ip==5992.or.ip==5993)write(*,*)'+',ip,dxlg(ip-1),df,dnstmp(ip)

     else if(dxpl(ip-1)< 0.d0) then ! .and. (dxpl(ip)<=0.d0))then

        df = (dxdt(ip) - gradv(ip)*Volume(ip) * .5d0*(dxlagr(ip)+dxlg(ip-1)) ) &
             & *(-dxlg(ip-1)) 
        if(ip>2)dxdtmp(ip-1) = dxdtmp(ip-1) + df
        if(ip<np)dxdtmp(ip) = dxdt(ip) * dxlagr(ip) - df

        df = (energy(ip) - grade(ip)*Volume(ip) &
             & * .5d0*(dxlagr(ip)+dxlg(ip-1)) ) *(-dxlg(ip-1)) 
        if(ip>2)engtmp(ip-1) = engtmp(ip-1) + df
        if(ip<np)engtmp(ip) = energy(ip) * dxlagr(ip) - df

        df = (vtrans(ip) - dvtran(ip)*Volume(ip) &
             & * .5d0*(dxlagr(ip)+dxlg(ip-1)) ) *(-dxlg(ip-1)) 
        if(ip>2)vtrtmp(ip-1) = vtrtmp(ip-1) + df
        if(ip<np)vtrtmp(ip) = vtrans(ip) * dxlagr(ip) - df

        df = (vtrant(ip) - dvtrat(ip)*Volume(ip) &
             & * .5d0*(dxlagr(ip)+dxlg(ip-1)) ) *(-dxlg(ip-1)) 
        if(ip>2)vtttmp(ip-1) = vtttmp(ip-1) + df
        if(ip<np)vtttmp(ip) = vtrant(ip) * dxlagr(ip) - df

        df = ( densit(ip) - gradd(ip) &
             & *.5d0*(dxplac(ip)+dxpl(ip-1)) )*(-dxpl(ip-1)) 
        if(ip>2) dnstmp(ip-1) = dnstmp(ip-1) + df
        if(ip<np)dnstmp(ip) = densit(ip) * dxplac(ip) - df

!        if(ip==5995)write(*,*)'-',dxlg(ip-1),df

     else 
        write(*,*)'Something wrong in Remapping'
        write(*,*)ip,dxpl(ip-1),dxpl(ip),xplace(ip-1),xplori(ip-1) &
             & ,xplace(ip),xplori(ip),dxplac(ip-1),dxplac(ip) &
             & ,radbnd(ip-1),radbnd(ip),dxdtbd(ip-1),dxdtbd(ip)
        stop
     endif

     if(ip>2)then
        dxdtmp(ip-1) = dxdtmp(ip-1) / dxlnew(ip-1)
        engtmp(ip-1) = engtmp(ip-1) / dxlnew(ip-1)
        vtrtmp(ip-1) = vtrtmp(ip-1) / dxlnew(ip-1)
        vtttmp(ip-1) = vtttmp(ip-1) / dxlnew(ip-1)
        dnstmp(ip-1) = dnstmp(ip-1) / dxplor(ip-1)
     end if

  end do
!---Remapping is finished-------------------------------------------------

!  sum1 = 0.d0; sum2 = 0.d0; sum3 = 0.d0
!  do ip = 2,np-1
!     sum1 = sum1 + dnstmp(ip)*dxplor(ip)
!     sum2 = sum2 + dxlnew(ip)
!     sum3 = sum3 + densit(ip)*dxplac(ip)
!  end do
!  write(*,*)'sum1,sum2,sum3',sum1,sum2,sum3
!  write(*,*)dxdt(3004:3005),dxdtmp(3004:3005),densit(3004:3005),dnstmp(3004:3005)
!---Substitution of the Remapped values-----------------------------------
!      write(*,*)dxpl0,dxpl(1),dnstmp(1),densit(1),dnstmp(2),densit(2)
!     &     ,dxpl(2),dxplac(2),dxplor(2)
!      write(*,*)dnstmp(1),densit(1)
  ip=1
  radbnd(ip) = radori(ip+1) - .5d0*drador(ip)
!///This part is necessary only for the wave generation at ip=1////////
  radius(ip) = radori(ip)
  dradiu(ip) = drador(ip)
  xplace(ip) = xplori(ip)
  dxplac(ip) = dxplor(ip)
  xlagra(ip) = xlgnew(ip)
!  dxlagr(ip) = dxlnew(ip)
  dxdt(ip)   = dxdtmp(ip)
  energy(ip) = engtmp(ip)
  densit(ip) = dnstmp(ip)
  volume(ip) = 1.d0/densit(ip)
  dxlagr(ip) = dxplac(ip)*densit(ip)
!  -----------------------
  vtrans(ip) = vtrtmp(ip)
  vtrant(ip) = vtttmp(ip)
!//////////////////////////////////////////////////////////////////////
  do ip=2,np-1
!         if(ip.eq.100)write(*,*)ip,densit(ip),dnstmp(ip)
!         if(ip.gt.400 .and. ip.lt.450)write(*,*)ip,dxpl(ip-1), dxpl(ip)
!     &        ,Btrans(ip),Btrtmp(ip)

     radius(ip) = radori(ip)
     dradiu(ip) = drador(ip)
     xplace(ip) = xplori(ip)
     dxplac(ip) = dxplor(ip)
     radbnd(ip) = radori(ip) + .5d0*drador(ip)
     dxdt(ip)   = dxdtmp(ip)
     energy(ip) = engtmp(ip)
     if(dnstmp(ip)<0.d0)then
!            write(*,*)ip,densit(ip),densit(ip-1),densit(ip),dnstmp(ip)
        dsav = .5d0*(densit(ip-1) + dnstmp(ip+1))
        if(dsav > 0.d0)then
           densit(ip) = dsav
        else
           densit(ip) = .9d0 * max(densit(ip-1),dnstmp(ip+1))
        endif
        write(*,*)'density < 0 in EULER', ip,densit(ip),dsav &
             & ,dxpl(ip-1),dxpl(ip),dnstmp(ip)
        stop
     else
        densit(ip) = dnstmp(ip)
     endif
     
     volume(ip) = 1.d0/densit(ip)

!         xlagra(ip) = xlgnew(ip)
!         dxlagr(ip) = dxlnew(ip)
     dxlagr(ip) = dxplac(ip)*densit(ip)
     xlagra(ip) = xlagra(ip-1) + .5d0*( dxlagr(ip-1)+dxlagr(ip) )
!    -----------------------
     vtrans(ip) = vtrtmp(ip)
     vtrant(ip) = vtttmp(ip)

  end do
  ip=np
  xplace(ip) = xplori(ip)
  dxplac(ip) = dxplor(ip)

!      write(*,*)dxlnew(506),densit(506)*dxplac(506)

end subroutine eulerrmp

!---------------------------------------------------------------------

!=====================================================================*
!           derive dBpp/dt by direct Eulerian                         *
!=====================================================================*
subroutine dBppcl
  use nrtype
  USE vars
  implicit none

  real(DP) :: x1,x2,dx,dx1,dx2,rmid,r1,r2,sqfmid,fltb
  integer(I4B) :: ip
  real(DP) :: f1,f2,Aflow1,Aflow2,etam,etap,eadm,eadp,Bsqttm,Bsqttp
! This subroutine derives dB_pp/dt from direct Euler method.   

  ip=1
  if (dxdtbd(ip)==0.d0)then
     fBtreb(ip) = fBtrbd(ip)
     fBtteb(ip) = fBttbd(ip)
     vftreb(ip) = vftrbd(ip)
     vftteb(ip) = vfttbd(ip)
     dxdteu(ip) = 0.d0
  else
     x1 = radbnd(ip) + .5d0*dt*dxdtbd(ip)
     x2 = radbnd(ip+1) + .5d0*dt*dxdtbd(ip+1)           
     dx = x2-x1
     dx1 = radbnd(ip) - x1
     dx2 = x2 - radbnd(ip)
     fBtreb(ip) = (fBtrbd(ip)*dx2 + fBtrbd(ip+1)*dx1)/dx
     fBtteb(ip) = (fBttbd(ip)*dx2 + fBttbd(ip+1)*dx1)/dx
     vftreb(ip) = (vftrbd(ip)*dx2 + vftrbd(ip+1)*dx1)/dx
     vftteb(ip) = (vfttbd(ip)*dx2 + vfttbd(ip+1)*dx1)/dx
     dxdteu(ip) = (dxdtbd(ip)*dx2 + dxdtbd(ip+1)*dx1)/dx
  endif
!      write(*,*)ip,rBtrbd(ip),rBtreb(ip),dxdtbd(ip)
  do ip=2,np-1
     if (dxdtbd(ip)==0.d0)then
        fBtreb(ip) = fBtrbd(ip)
        fBtteb(ip) = fBttbd(ip)
        vftreb(ip) = vftrbd(ip)
        vftteb(ip) = vfttbd(ip)
        dxdteu(ip) = 0.d0
     else if (dxdtbd(ip)>0.d0)then
        x1 = radbnd(ip-1) + .5d0*dt*dxdtbd(ip-1)
        x2 = radbnd(ip) + .5d0*dt*dxdtbd(ip)           
        dx = x2-x1
        dx1 = radbnd(ip) - x1
        dx2 = x2 - radbnd(ip)
        fBtreb(ip) = (fBtrbd (ip-1)*dx2 + fBtrbd(ip)*dx1)/dx
        fBtteb(ip) = (fBttbd (ip-1)*dx2 + fBttbd(ip)*dx1)/dx
        vftreb(ip) = (vftrbd (ip-1)*dx2 + vftrbd(ip)*dx1)/dx
        vftteb(ip) = (vfttbd (ip-1)*dx2 + vfttbd(ip)*dx1)/dx
        dxdteu(ip) = (dxdtbd (ip-1)*dx2 + dxdtbd(ip)*dx1)/dx
     else
        x1 = radbnd(ip) + .5d0*dt*dxdtbd(ip)
        x2 = radbnd(ip+1) + .5d0*dt*dxdtbd(ip+1)           
        dx = x2-x1
        dx1 = radbnd(ip) - x1
        dx2 = x2 - radbnd(ip)
        fBtreb(ip) = (fBtrbd(ip)*dx2 + fBtrbd(ip+1)*dx1)/dx
        fBtteb(ip) = (fBttbd(ip)*dx2 + fBttbd(ip+1)*dx1)/dx
        vftreb(ip) = (vftrbd(ip)*dx2 + vftrbd(ip+1)*dx1)/dx
        vftteb(ip) = (vfttbd(ip)*dx2 + vfttbd(ip+1)*dx1)/dx
        dxdteu(ip) = (dxdtbd(ip)*dx2 + dxdtbd(ip+1)*dx1)/dx
!            if(ip.eq.100)write(*,*)dx,dx1,dx2,Btrbnd(ip),Btrbnd(ip+1)
     endif

  end do
!         write(*,*)'100',Btrrbd(100),Btrlbd(101),Btrbnd(100)
!     &        ,Btreub(100),dxdtbd(100)

  if(irss/=0.and.t>t_e2T_rss)call resistivity
!---- derive dBpp/dt--------------------------------------------
  do ip=2,np-1
     rmid = radori(ip)
     r1 = rmid - .5*drador(ip)
     r2 = rmid + .5*drador(ip)
     sqfmid = sqrt(fltb(rmid))
!         sqf1 = sqrt(fltb(r1))
!         sqf2 = sqrt(fltb(r2))
     if(irss==0.or.ip==2.or.t<t_e2T_rss)then
        dtdBtr(ip) = ( Bplint * (vftreb(ip) -vftreb(ip-1) ) &
             & - (dxdteu(ip)*fBtreb(ip) - dxdteu(ip-1)*fBtreb(ip-1) ) ) &
             & / (drador(ip)*rmid*sqfmid)
        dtdBtt(ip) = ( Bplint * (vftteb(ip) -vftteb(ip-1) ) &
             & - (dxdteu(ip)*fBtteb(ip) - dxdteu(ip-1)*fBtteb(ip-1) ) ) &
             & / (drador(ip)*rmid*sqfmid) 
     else 
        etam = min(0.5d0 * (eta(ip-1) + eta(ip)), etamax)
        etap = min(0.5d0 * (eta(ip) + eta(ip+1)), etamax)
        if(irss==1)then
           eadm = 0.d0
           eadp = 0.d0
        else
           f1 = fltb(r1)
           f2 = fltb(r2)
           Aflow1 = r1 * r1 * f1
           Aflow2 = r2 * r2 * f2
           if(irss==2)then
              Bsqttm = (fBtreb(ip-1)*fBtreb(ip-1) + fBtteb(ip-1)*fBtteb(ip-1))&
                   &  / Aflow1
              Bsqttp = (fBtreb(ip)*fBtreb(ip) + fBtteb(ip)*fBtteb(ip) ) &
                   & / Aflow2
           else
              Bsqttm = (fBtreb(ip-1)*fBtreb(ip-1) + fBtteb(ip-1)*fBtteb(ip-1) &
                   & + Bplint*Bplint/ Aflow1) / Aflow1
              Bsqttp = (fBtreb(ip)*fBtreb(ip) + fBtteb(ip)*fBtteb(ip) &
                   & + Bplint*Bplint/ Aflow2) / Aflow2
           end if
           eadm = min(0.5d0 * (ead(ip-1) + ead(ip)) * Bsqttm,etamax)
           eadp = min(0.5d0 * (ead(ip) + ead(ip+1)) * Bsqttp,etamax)
        end if
        dtdBtr(ip) = ( Bplint * (vftreb(ip) -vftreb(ip-1) ) &
             & - (dxdteu(ip)*fBtreb(ip) - dxdteu(ip-1)*fBtreb(ip-1) ) &
             & + ( (etap+eadp) * (fBtreb(ip+1) - fBtreb(ip-1)) / (drador(ip+1) &
             & + drador(ip)) - (etam+eadm) * (fBtreb(ip) + fBtreb(ip-2)) &
             & / (drador(ip-1) + drador(ip)) ) ) / (drador(ip)*rmid*sqfmid)
!        if(ip==21)write(*,*)eta(ip-1:ip+1),eadp,eadm
        dtdBtt(ip) = ( Bplint * (vftteb(ip) -vftteb(ip-1) ) &
             & - (dxdteu(ip)*fBtteb(ip) - dxdteu(ip-1)*fBtteb(ip-1) ) &
             & + ( (etap+eadp) * (fBtteb(ip+1) - fBtteb(ip-1)) / (drador(ip+1) &
             & + drador(ip)) - (etam+eadm) * (fBtteb(ip) + fBtteb(ip-2)) &
             & / (drador(ip-1) + drador(ip)) ) ) / (drador(ip)*rmid*sqfmid)
     end if

!         if(ip.eq.611)write(*,*)ip,vftreb(ip-1),vftreb(ip),dtdBtr(ip)
!     &        ,dxdteu(ip),dxdteu(ip-1),rBtreb(ip),rBtreb(ip-1)
!         if(ip.eq.2)write(*,*)ip,Bplint * (vftreb(ip) -vftreb(ip-1) ) 
!     &        ,dxdteu(ip)*rBtreb(ip)*sqf1,dxdteu(ip-1)*rBtreb(ip-1)*sqf2
!     dtrBtr(ip),vftreb(ip),dxdteu(ip)
!     &        ,rBtreb(ip),dxplor(ip),vftreb(ip-1),dxdteu(ip-1)
!     &        ,rBtreb(ip-1),sqf1,sqf2
  end do

end subroutine dBppcl

!------------------------------------------------------------------------------------------
subroutine turbdisp
  use nrtype
  USE vars
  implicit none

  integer(i4b) :: ip
  real(dp) :: ztrans_p, ztrans_m, ztrant_p, ztrant_m
  
  do ip = 2, np-1
     ztrans_p = abs( vtrans(ip) - fc1*Btrans(ip)/sqrt(densit(ip)) )
     ztrans_m = abs( vtrans(ip) + fc1*Btrans(ip)/sqrt(densit(ip)) )
     ztrant_p = abs( vtrant(ip) - fc1*Btrant(ip)/sqrt(densit(ip)) )
     ztrant_m = abs( vtrant(ip) + fc1*Btrant(ip)/sqrt(densit(ip)) )

     vtrans(ip) = vtrans(ip) - dt*ettrb(ip) * ((ztrans_p + ztrans_m)*vtrans(ip) &
          & + (ztrans_p - ztrans_m)*fc1*Btrans(ip)/sqrt(densit(ip)) )
     Btrans(ip) = Btrans(ip) - dt*ettrb(ip) * ((ztrans_p + ztrans_m)*Btrans(ip) &
          & + (ztrans_p - ztrans_m)/fc1*vtrans(ip)*sqrt(densit(ip)) )
     vtrant(ip) = vtrant(ip) - dt*ettrb(ip) * ((ztrant_p + ztrant_m)*vtrant(ip) &
          & + (ztrant_p - ztrant_m)*fc1*Btrant(ip)/sqrt(densit(ip)) )
     Btrant(ip) = Btrant(ip) - dt*ettrb(ip) * ((ztrant_p + ztrant_m)*Btrant(ip) &
          & + (ztrant_p - ztrant_m)/fc1*vtrant(ip)*sqrt(densit(ip)) )
  end do
  
end subroutine turbdisp

!=====================================================================*
!           Integration of Evolution by Runge-Kutta Method            *
!=====================================================================*
subroutine evolve
  use nrtype
  USE vars
  implicit none

  real(DP) :: dtmin,wvlgsw,wvtrsw,Bprs,rnow,rnow1,rnow2&
       & ,riol,ricu,dxdtol,fsqold,fsqnew,spav &
       & ,velmax,dtminx,Aflow,trflux,pflux,grv,qkinpl,qkintr,entha,btr,ten &
       & ,fltb,r1,r2 ! sqdav,vij,xij,absxij,,thraw,th
  real(DP) :: perturblg, perturb1t,perturb2t
  real(DP) :: sum
  integer(I4B) :: in10,ip,k,i,ivtr,iprint,itrest1,irand!,input,mth
  real(DP) :: tsumav
  real(DP) :: dvcrs,dnscorr,specorr
  INTEGER(I4B),parameter :: ntime = 3900
  INTEGER(I4B) :: itime,iefave
  REAL(DP),parameter :: dtprnt = 1.d-2, dtmin_th=1.d-9
  REAL(DP) :: tprnt(0:ntime)
! ------------------------------------------------------------------
  write(*,'(/,a)')          ' 2. Subroutine Evolve'
! -----------------------------------------------------------------
!     xdmin = - ql1/float(np1)
!     dt    = xdmin/sounds(1)*cflnum 
!  t=0.d0
!  t=3.0d0
  call avmnmx(0)

  do itime = 0,ntime
     if(icont==2)then
        tprnt(itime) = dtprnt*dble(itime) + dtprnt*dble(int(t/dtprnt))
     else if(icont==1)then
        tprnt(itime) = dtprnt*dble(int(tprnt0/dtprnt)-1) + dtprnt*dble(itime)
     else
        tprnt(itime) = dtprnt*dble(itime)
     end if
  end do
  call printvars(0,0) ! initializing average values
  write(*,*)icont,t,tprnt0,tprnt(irstr),irstr
  if(t>tprnt(irstr))irstr = irstr+1

  do ip=1,np  ! initial Bpp for starting point of each time step 
     Btrprv(ip)=Btrans(ip)
     Bttprv(ip)=Btrant(ip)
  enddo
!     -----------
!      call monoto(0)
!      call rhs
!---- please comment out for Alfven wave propagation ---------------
!c      call updBpp ! for updating B_perp
!-------------------------------------------------------------------
!      call moc
!      call dBppcl               ! derive dBpp/dt by direct Eulerian
!     -----------

!      call pgssk(.FALSE.)

  ivtr = 1 ! randam perturbation

!  iprint = 0
!  irstr = 1 ! for Final_Data
!  irstr = 0
!  itrest1 = 10000 ! for average

  if(iwv==1)tsumav = 0.d0 ! for time average of wave energy flux
  iefave = 0

  do in10=1,100000000
!     if(t > (tfinal-tefave).and.iefave==0)then
!        call engflx(0)
!        iefave = 1
!     end if
! ////////////Output by pgplot//////////////////////// 
!         if(mod(in10,20).eq.1)then
!            call movie (iprint)
!         endif
!//////End of Output ///////////////////////////
!     write(*,*)dradiu(1:3)

!        xdmin = ql2 - ql1
!        xdmin = ql
!          dtmin = xdmin/sounds(1)*cflnum
     dtmin=tfinal ! initilization by sano
!       ----------
     t = t + dt
          
!     write(*,*)radbnd(1),dxdtbd(1)
!       ----------
     ip = 1
     radbnd(ip) = radbnd(ip) + dxdtbd(ip)*dt 
     drador(ip) = dradiu(ip)
     radori(ip) = radius(ip)
     xplori(ip) = xplace(ip)
     dxplor(ip) = dxplac(ip)
!//////////Wave Geration at ip=1////////////////////////////////////
     wvlgsw = 1.d0
     wvtrsw = 1.d0
!---- for soliton ----------------------------------------
!          if(omega*t .gt. 9.)wvlgsw=0. 
!          if(omgtr*t .gt. 9.)wvtrsw=0.
!-------------------------------------------------------
! smoothed random perterbation ; smoothing scale is trd(i)-trd(i-1)
     do irand=1,nrand
!             write(*,*)irand
        if(t<trd(irand))then
           if(irand>1)then
              perturblg = randlg(irand-1) + (t-trd(irand-1)) &
                   & /(trd(irand)-trd(irand-1)) &
                   & *(randlg(irand)-randlg(irand-1))
!                   write(*,*)t,irand,perturb,trd(irand),random(irand)
              perturb1t = random(irand-1) + (t-trd(irand-1)) &
                   & /(trd(irand)-trd(irand-1)) &
                   & *(random(irand)-random(irand-1))
              perturb2t = rand3c(irand-1) + (t-trd(irand-1)) &
                   & /(trd(irand)-trd(irand-1)) &
                   & *(rand3c(irand)-rand3c(irand-1))
           else 
              perturblg = randlg(irand)
              perturb1t = random(irand)
              perturb2t = rand3c(irand)
           endif
           exit
        endif
     enddo

!     densit(ip) = dinit + wvlgsw * ddens*sin(omega*t)
     densit(ip) = dinit
!          densit(ip) = densit(ip) + dt * omega * ddens 
!     &         * cos(omega*(t-.5*dt))
!          densit(ip) = densit(ip+1) * exp (grav *(1./radius(np) 
!     &         -1./radius(np+1)) )

!     dxdt(ip)   = vinit + wvlgsw * dvelo*sin(omega*t)
     dxdt(ip) = vinit + wvlgsw * dvelo*perturblg
!          dxdt(ip) = dxdt(ip) + dt * omega * dvelo 
!     &         * cos(omega*(t-.5*dt))
!          dxdt(ip) = dxdt(ip+1) ! for Alfven wave
!          write(*,*)dxdt(ip),vinit,dvelo

!     prsgas(ip) = pinit + wvlgsw * dpres*sin(omega*t)
     prsgas(ip) = pinit 
!          prsgas(ip) = prsgas(ip) + dt * omega * dpres 
!     &         * cos(omega*(t-.5*dt)) 
!          prsgas(ip) = prsgas(ip+1) * exp (grav *(1./radius(np) 
!     &         -1./radius(np+1)) )         
!          speene(ip) = speene(ip+1)
!          prsgas(ip) = speene(ip)*gammi1*densit(ip) 

     volume(ip) = 1.d0/densit(ip)    

! smoothed random perterbation ; smoothing scale is trd(i)-trd(i-1)
!     do irand=1,nrand
!!             write(*,*)irand
!        if(t<trd(irand))then
!           if(irand > 1)then
!              perturb = random(irand-1) + (t-trd(irand-1)) &
!                   & /(trd(irand)-trd(irand-1)) &
!                   & *(random(irand)-random(irand-1))
!!                   write(*,*)t,irand,perturb,trd(irand),random(irand)
!           else 
!              perturb = random(irand)
!           endif
!           exit
!        endif
!     enddo
!     write(*,*)irand,t,trd(irand),perturb

!     input = int(.5*omgtr*t / pi)

     vtrans(ip) = dvtrit * perturb1t

     Btrans(ip) = 0.d0

!     thraw = omgtr*t/pi
!     mth = int(thraw)
!     th = thraw - dble(mth)
!     if(th<= .5d0)then
!        Btrans(ip) = 1.0d0*Btrans(ip)
!        vtrans(ip) = 1.0d0*vtrans(ip)
!     endif

     Btrant(ip) = 0.d0
     vtrant(ip) = dvttit * perturb2t

!     thraw = (omgtr*t)/pi -.5d0
!     mth = int(thraw)
!     th = thraw - dble(mth)
!     if(th<= .5d0)then
!        Btrant(ip) = 1.0d0*Btrant(ip)
!        vtrant(ip) = 1.0d0*vtrant(ip)
!     endif
        
     Bprs = fc3 * (Btrans(ip)*Btrans(ip)+Btrant(ip)*Btrant(ip))
     pressu(ip) = prsgas(ip) + Bprs 

     dxplac(ip) = volume(ip)*dxlagr(ip)
     dradiu(ip) = dxplac(ip)         ! r=1
     radius(ip) = radbnd(ip) - .5d0*dradiu(ip)
     xplace(ip) = radius(ip)  !! might be improved

     speene(ip) = prsgas(ip)/(densit(ip)*gammi1)
     energy(ip) = speene(ip) + 0.5*( dxdt(ip)*dxdt(ip) &
          & + vtrans(ip)*vtrans(ip) + vtrant(ip)*vtrant(ip) ) &
          & + Bprs / densit(ip) - grav / radius(ip)
     sounds(ip) = sqrt( gamma*gammi1*speene(ip) )
     Alfven(ip) = fc1*Bparal(ip)/sqrt(densit(ip))
     engtot(ip) = energy(ip) + fc3*Bparal(ip)*Bparal(ip)/densit(ip) 

!     write(*,*)ip,radius(ip),radbnd(ip),dradiu(ip)
!///////////////////////////////////////////////////////////////////   


! acoustic waves in the corona
!          dxdt(400) = dxdt(400) + 10.*wvlgsw * dvelo*sin(omega*t)

!     -----------
     call monoto(0)
     call rhs
!---- please comment out for Alfven wave propagation ---------------
!      call updBpp ! for updating B_perp
!-------------------------------------------------------------------
     call moc
     call dBppcl               ! derive dBpp/dt by direct Eulerian
!     -----------

! ---- rest of ip=1 ------------------------------------------

!--------------------------------------------------------------
!     ip = 4
!     r1 = radius(ip) - 0.5d0*dradiu(ip)
!     r2 = radius(ip) + 0.5d0*dradiu(ip)
!     write(*,*)'before',real(densit(ip)),real(densit(ip) + 0.25d0*dt &
!          &/ dxplac(ip) *((densit(ip-1)+densit(ip))*(dxdt(ip-1)+dxdt(ip))&
!          & *r1*r1*fltb(r1) - (densit(ip)+densit(ip+1)) &
!          & *(dxdt(ip)+dxdt(ip+1))*r2*r2*fltb(r2) ) )

     do ip = 2,np-1
        drador(ip) = dradiu(ip)
        radori(ip) = radius(ip)
        dxplor(ip) = dxplac(ip)
        xplori(ip) = xplace(ip)
        radbnd(ip) = radbnd(ip)+dxdtbd(ip)*dt ! by stakeru
        radius(ip) = .5d0*(radbnd(ip-1)+radbnd(ip)) ! by stakeru
        dradiu(ip) = radbnd(ip)-radbnd(ip-1) ! by stakeru
        rnow = radius(ip)
!        dxplac(ip) = rnow*rnow*fltb(rnow)*dradiu(ip)
!        xplace(ip) = xplace(ip-1) + .5d0*(dxplac(ip-1) + dxplac(ip))
             ! dxplac is calclated below       
!             if(ip.ne.(np-1))write(*,*)ip,radbnd(ip),dxdtbd(ip)
!     &            ,radius(ip),xplace(ip)
     end do

     radori(np) = radius(np)
     drador(np) = dradiu(np)
     xplori(np) = xplace(np)
     dxplor(np) = dxplac(np)    

!     sum = 0.d0
!     do ip = 2,np-1
!        sum = sum + densit(ip)*dxplac(ip)
!     end do
!     write(*,*)'sum1=',sum
!-----------------------------------------------------------------------
! complete MOC by using improved xplace etc.
     call trans
!-----------------------------------------------------------------------
     do ip = 2,np-1
        riol = 1.d0/(radori(ip)*radori(ip))
        ricu = 1.d0/(radius(ip)*radius(ip))
        dxdtol = dxdt(ip)
!             drh = .5*dradiu(ip)
!             dns = densit(ip)
!             ddndr = gradd(ip)*radius(ip)*radius(ip)*fltb(radius(ip))
!             r1 = - drh
!             r2 = .5*drh
!             racc = 1.e-6
!             rgr = radius(ip)+rtnewt(centmass,r1,r2,racc)
!             write(*,*)radius(ip),rgr
!     ------------------------------------
!             if(ip.eq.2)write(*,*)ip,volume(ip),dvoldt(ip),ddxdtt(ip)
!     &            ,denedt(ip)
        volume(ip) = volume(ip) + dt*dvoldt(ip)
        dxplac(ip) = volume(ip) * dxlagr(ip)
        xplace(ip) = xplace(ip-1) + .5d0*(dxplac(ip-1) + dxplac(ip))
        dxdt(ip) =   dxdt(ip) + dt*ddxdtt(ip) + dt * (Bcenrf(ip) &
             & - grav *.5d0 * (riol + ricu) )
                                !/(rgr*rgr)
        energy(ip) = energy(ip) + dt*(denedt(ip) + Btense(ip))
!*     &            - grav*.5d0 * (dxdtol*riol + dxdt(ip)*ricu) )
!*             speene(ip) = speene(ip) + dt*denedt(ip) 
!        if(ip==10)write(*,*)ip,real(dxdt(ip)),real(ddxdtt(ip)),real(Bcenrf(ip))&
!             & ,real(-grav *.5d0 * (riol + ricu) ),real(riol),real(ricu)
!     &            ,Btense(ip),rBtrbd(ip),vftrbd(ip)
!     &            ,Btrans(ip+1),vtrans(ip+1)
!             if(ip==70)write(*,*)ip,denedt(ip),Btense(ip)
!     &            ,dradiu(ip),radius(ip),dxdtbd(ip),dxdtbd(ip-1)

        densit(ip) = 1.0d0/volume(ip)

!        if(ip==80)write(*,*)ip,densit(ip),speene(ip)
!     -------------------------
        fsqold = sqrt(fltb(radori(ip)) )
        fsqnew = sqrt(fltb(radius(ip)) )
        vtrans(ip) = (vtrans(ip)*radori(ip)*fsqold + dt*dtvrtr(ip)) &
             & / ( radius(ip)*fsqnew )
        vtrant(ip) = (vtrant(ip)*radori(ip)*fsqold + dt*dtvrtt(ip)) &
             & / ( radius(ip)*fsqnew )
!            if(ip==2)write(*,*)ip,vtrans(ip)
!             if(ip<=6)write(*,*)ip,dtvrtr(ip)/radius(ip)
     end do

!     sum = 0.d0
!     do ip = 2,np-1
!        sum = sum + densit(ip)*dxplac(ip)
!     end do
!     write(*,*)'sum2=',sum
!     write(*,*)'after',dxlagr(100),densit(100)*dxplac(100),densit(100) &
!          & ,dxplac(100)

! ------------ outflow boundary condition -------------------
!              densit(np) = densit(np-1) * exp (grav * speene0 
!     &         / speene(np-1) * (1./radius(np) -1./radius(np-1)) )
!                    !Evaluate at the local temperature 
!                                !! speene(np-1) is at previous step 
!              dxdt(np)   = dxdt(np-1)
!              energy(np) = energy(np-1)
!              vtrans(np) = vtrans(np-1)
!              vtrant(np) = vtrant(np-1)

!     write(*,*)'before remap',real(densit(4)),real(dxdt(4)),real(dxdtbd(3:4))
!-----------------------------------------------------------------------
! Euler Remapping for update values
     call eulerrmp
!-----------------------------------------------------------------------
!     write(*,*)'after remap',real(densit(4)),real(dxdt(4))
!          write(*,*)'third',energy(32)
!          write(*,*)densit(614),densit(615),speene(614),speene(615)
!          write(*,*)densit(617),dxdt(617)

     do ip=2,np-1
        Btrans(ip) = Btrprv(ip) + dt*dtdBtr(ip)
        Btrant(ip) = Bttprv(ip) + dt*dtdBtt(ip)
!                 if(ip.le.6)write(*,*)ip,dtdBtr(ip)
     end do
     if(itrb==1)call turbdisp ! turbulent dissipation
!     write(*,*)Btrans(21),dtdBtr(21)
     
     do ip=2,np-1
        Btrprv(ip) = Btrans(ip)
        Bttprv(ip) = Btrant(ip)
!                 if(ip.eq.2)write(*,*)ip,Btrprv(ip),dtrBtr(ip)
!                 if(ip.le.50.and.ip.ge.45)write(*,*)ip,Btrans(ip)

        Bprs = fc3 * (Btrans(ip)*Btrans(ip) + Btrant(ip)*Btrant(ip))
!                 speene(ip) = energy(ip) - 0.5*dxdt(ip)*dxdt(ip)
!     &                - Bprs / densit(ip)
!                 if(speene(ip).le.0.)write(*,*)'sp',ip,energy(ip)
!     &                ,dxdt(ip),Bprs
!                 th=0.
!                 speene(ip) = dmax1(speene(ip),th)

        speene(ip) = energy(ip) - 0.5 * (dxdt(ip)*dxdt(ip) &
             & + vtrans(ip)*vtrans(ip) + vtrant(ip)*vtrant(ip) ) &
             & - Bprs / densit(ip) + grav / radius(ip)

!                 if(ip.eq.32)write(*,*)ip,'fourth',energy(ip),speene(ip)
!             if(mod(ip,400).eq.32)write(*,*)ip,energy(ip),speene(ip),dt
!     &           ,denedt(ip),Btense(ip),dxdt(ip),vtrans(ip),Bprs
!     &           ,Bparal(ip),densit(ip), grav / radius(ip)    
     end do

     do ip = 2, np-1
        if(speene(ip)<0.d0)then
           spav = .5d0*(speene(ip-1) + speene(ip+1))
!                    write(*,*)'e < 0', ip,speene(ip),spav,energy(ip)
!     &                   ,dxdt(ip),vtrans(ip),vtrant(ip),Bprs
!     &                   ,densit(ip)
           if(spav>0.d0)then
              speene(ip) = .5d0*spav  ! reduce by 1/2
           else
              speene(ip) = .5d0*max(speene(ip-1),speene(ip+1))
           endif
        endif        
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!        if(radius(ip)>1.1)speene(ip) = max(speene(ip),1.d-1)
        if(radius(ip)>1.1)then
           speene(ip) = max(speene(ip),spemin)
           speene(ip) = min(speene(ip),spemax)
        end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        prsgas(ip) = densit(ip)*gammi1*speene(ip)
        pressu(ip) = prsgas(ip) + Bprs
     end do

!     write(*,*)energy(2),speene(2)

!-----------outflow boundary --------------------------------------
!              call outgobdsc
!              Btrans(np)=0.
!              vtrans(np)=0.
!              Btrant(np)=0.
!              vtrant(np)=0.
!----------- conduction & radiation cooling for speene (Euler)-------
     if(t>=0.d0)then
!                 if(t.lt.1.e-1)then
!                    radfac = radfac0 * .5
!                 else
!        radfac = radfac0 * .8
!                 endif
!        if(t>1.194d-2)write(*,*)'before cond'
        call conduction
!        if(t>1.194d-2)write(*,*)'after cond'
     endif
!---------------------------------------------------------------------
     
     do ip = 2,np
        Bprs = fc3 * (Btrans(ip)*Btrans(ip) + Btrant(ip)*Btrant(ip))
        energy(ip) = speene(ip) + 0.5d0 * (dxdt(ip)*dxdt(ip) &
             & + vtrans(ip)*vtrans(ip) + vtrant(ip)*vtrant(ip) ) &
             & + Bprs / densit(ip) - grav / radius(ip)
                                ! total energy is recalculate from  
                                ! corrected speene
!                 if(ip.eq.32)write(*,*)ip,' fifth',energy(ip),speene(ip)

        prsgas(ip) = densit(ip)*gammi1*speene(ip)
        pressu(ip) = prsgas(ip) + Bprs
        sounds(ip) = sqrt( gamma*gammi1*speene(ip) )
        Alfven(ip) = fc1*Bparal(ip)/sqrt(densit(ip))

        engtot(ip) = energy(ip) + fc3*Bparal(ip)*Bparal(ip)/densit(ip) 
!                 if(ip.eq.2)write(*,*)ip,t,speene(ip)
!                 if(ip.eq.70)write(*,*)ip,energy(ip),speene(ip)
!     &                ,prsgas(ip),densit(ip),denedt(ip),Btense(ip)
!     &                ,- grav / radius(ip),(Bprs + fc3 
!     &                * Bparal(ip)*Bparal(ip) ) / densit(ip)
!                --------------------------------------------------
!        xij        = dradiu(ip)
                !xij        = xplace(ip) - xplace(ip-1)
!        vij        =   dxdt(ip) -   dxdt(ip-1)
!        absxij     = abs( xij )
!                 xdmin      = dmin1( xdmin , absxij )
!                 velmax     = sounds(ip) + abs(dxdt(ip)) ! by sano
        velmax     = sqrt(sounds(ip)*sounds(ip) + fc2/densit(ip) &
             & *(Btrans(ip)*Btrans(ip) + Btrant(ip)*Btrant(ip) &
             & + Bparal(ip)*Bparal(ip) )) + abs(dxdt(ip))  
        dtminx = dradiu(ip) / velmax
!        dtminx     = absxij/velmax
!                 if(Alfven(ip) .ge. 2.*sounds(ip) .and. 
!     &               .1*Alfven(ip)*t .le. radius(ip) )
!     &              dtminx = .7*(abs(Alfven(ip)/sounds(ip)))**.5*dtminx  
        if(dtminx<dtmin_th)then ! density floor 
           dvcrs = dradiu(ip)/dtmin_th
           if(sounds(ip)>0.5d0*dvcrs)then
              sounds(ip) = 0.5d0*dvcrs
              speene(ip) = sounds(ip)*sounds(ip)/(gamma*gammi1)
           end if
           dnscorr = (densit(ip)*Alfven(ip)*Alfven(ip) + Btrans(ip)*Btrans(ip) &
                & + Btrant(ip)*Btrant(ip)) / ((dvcrs - abs(dxdt(ip)))**2.d0&
                & - sounds(ip)*sounds(ip))
!           specorr = min(10.d0*speene(ip) * densit(ip) / dnscorr,speene(ip))
           dnscorr = max(dnscorr,0.5d0*densit(ip-1),densit(ip))
           specorr = max(speene(ip) * densit(ip) / dnscorr,spemin)
           specorr = min(specorr,spemax)
           write(*,*)'density floor',ip,densit(ip),dnscorr,speene(ip),specorr,velmax,dtminx &
                & ,sounds(ip),Btrans(ip),Btrant(ip),Bparal(ip),dxdt(ip)           
           densit(ip) = dnscorr
           speene(ip) = specorr
           prsgas(ip) = densit(ip)*gammi1*speene(ip)
           pressu(ip) = prsgas(ip) + Bprs
           sounds(ip) = sqrt( gamma*gammi1*speene(ip) )
           Alfven(ip) = fc1*Bparal(ip)/sqrt(densit(ip))
           energy(ip) = speene(ip) + 0.5d0 * (dxdt(ip)*dxdt(ip) &
                & + vtrans(ip)*vtrans(ip) + vtrant(ip)*vtrant(ip) ) &
                & + Bprs / densit(ip) - grav / radius(ip)
           engtot(ip) = energy(ip) + fc3*Bparal(ip)*Bparal(ip)/densit(ip) 
           volume(ip) = 1.d0/densit(ip) 
           dxlagr(ip) = dxplac(ip) / volume(ip) 
           velmax     = sqrt(sounds(ip)*sounds(ip) + fc2/densit(ip) &
                & *(Btrans(ip)*Btrans(ip) + Btrant(ip)*Btrant(ip) &
                & + Bparal(ip)*Bparal(ip) )) + abs(dxdt(ip))  
           dtminx = dradiu(ip) / velmax
        end if
        dtmin      = min( dtmin , dtminx )
!                 write(*,*)ip,dtmin,velmax,sounds(ip),dxdt(ip)
!     &                ,densit(ip)
        clocal(ip) = dtminx

     end do
           
! ------------ outflow boundary condition -------------------
     call outgobdsc ! sholud be put above 
!     engtot(np) = energy(np) + fc3*Bparal(np)*Bparal(np) / densit(np) 

!        if( t .gt. tfinal ) go to 9000
     
     call tavevars
!     call engflx !(ef)
!     if(iefave==1)call engflx(1)
     call avmnmx(1)
     
!          if(t.gt.1.145e0 .and. t.lt.1.148e0)then
!             if(dt.gt.3.e-8)dt = 3.e-8 ! for conduction 
!          else if(t.gt.1.148e0)then
!             if(dt.gt.1.e-8)dt = 1.e-8
!          endif
!         -----------
!          call monoto(0)
!          call rhs
!---- please comment out for Alfven wave propagation ---------------
!c          call updBpp ! for updating B_perp
!-------------------------------------------------------------------
!          call moc
!          call dBppcl       ! derive dBpp/dt by direct Eulerian
!         -----------
     if(mod(in10,50)==1)then
        write(6,'(5x,a,i8,a,1p2e11.3)') 'Step',in10,' ended.  t =',t,dt
        call flush(6)
     endif
     if(t>=tprnt(irstr))then
        call printvars(irstr+1,1)
        write(43,'(f12.6,1p5e12.4)')t,(radems(k),k=1,5)
        irstr = irstr + 1
     endif
!     write(*,*)dxplac(1)*densit(1),dxlagr(1),dxplac(2)*densit(2),dxlagr(2)
!10 continue
     if(iwv==1)then
        call wveng
        do ip = 1,np
           Faloav(ip) = Faloav(ip) + Falout(ip) * dt
           Faliav(ip) = Faliav(ip) + Falin(ip) * dt
           Faloav3(ip) = Faloav3(ip) + Falout3(ip) * dt
           Faliav3(ip) = Faliav3(ip) + Falin3(ip) * dt
           Fcsoav(ip) = Fcsoav(ip) + Fcsout(ip) * dt
        enddo
        tsumav = tsumav + dt
     end if

     if( t > tfinal ) exit
     dt = dtmin*cflnum
  end do
!     ------------------------------------------------------------------
9000 write(6,'(5x,a)')         'The End of Integration'
  call printvars(irstr+1,1)
  write(43,'(f12.6,1p5e12.4)')t,(radems(k),k=1,5)
  call avmnmx(2)
  
  open(unit=14,file='Initial_index.dat')
  write(14,*) irstr+1
  close(14)

!--- print out energy flux ----------------------------
  do ip=1,np
!     Aflow = radius(ip) * radius(ip) * fltb(radius(ip))
!     trflux = -Aflow * fc2*Bparal(ip)* (Btrans(ip)*vtrans(ip) &
!          & + Btrant(ip)*vtrant(ip) )
!     pflux = Aflow * densit(ip) * dxdt(ip)
!     grv = - pflux * grav / radius(ip) 
!     qkinpl = pflux * .5*dxdt(ip)*dxdt(ip)
!     qkintr = pflux * .5 * (vtrans(ip)*vtrans(ip) + vtrant(ip)*vtrant(ip) )
!     entha = pflux * (speene(ip) + prsgas(ip) / densit(ip) )
!     btr = fc2*(Btrans(ip)*Btrans(ip)+Btrant(ip)*Btrant(ip)) * dxdt(ip)*Aflow
!     ten= qkinpl + qkintr + entha + btr + trflux + grv
!     write(13,'(1p8e13.5)')radius(ip),trflux,qkinpl,qkintr,entha,grv,btr,ten
     write(13,'(1p11e13.5)')radius(ip),ef(1:10,ip)/tefsum
  enddo

  if(iwv==1)then
     do ip = 1,np
        Faloav(ip) = Faloav(ip)/tsumav
        Faliav(ip) = Faliav(ip)/tsumav
        Faloav3(ip) = Faloav3(ip)/tsumav
        Faliav3(ip) = Faliav3(ip)/tsumav
        Fcsoav(ip) = Fcsoav(ip)/tsumav
        write(61,'(1p6e13.5)')radius(ip),Faloav(ip),Faliav(ip),Fcsoav(ip) &
             & ,Faloav3(ip),Faliav3(ip)
     enddo
  end if

end subroutine evolve
