module vars 
  implicit none
  integer,parameter :: nmvar=10,nmave=13,ntime=3901
  real,parameter :: Ysun = 0.2485, Zsun = 0.0134, dYdZ = 0.125
  real,parameter :: Thigh = 1.d4, Tlow = 2.d3
  real,parameter :: Rcon=8.26e7,pi=3.14159265,Rsun=6.96e5,prtnms=1.67e-24 &
       & ,Teffsun=5780.d0 ! Rcon=kB/mH
  real(8),parameter :: fmax=1.0d0 ,re1=1.2d0,sgm=0.15d0, sgmcsun = 0.015d0
  real :: Mstar,Teff,akm,rhocgs,gamma,gammi1,spetkl,Rstar,rdkm,taucsc,dnnH1,clmfct &
       & ,zmetal,YHe,YoX,XHyd,mufull,muzero,mufrac,coefT1,coefT0,akmsun
  integer :: np,nd

  real(8),allocatable :: var(:,:)
  real,allocatable :: varave(:,:)
  real,allocatable :: rd(:),zd(:),dn(:),vl(:),Tprt(:),Bt(:),eta(:),ead(:),mul(:)&
       & ,vt(:),vt3(:),Bt3(:),Br(:),BtBr_tav(:),vrot(:),vtave(:),vtave3(:) &
       & ,Btave(:),Btave3(:),va(:),dr(:),dvsqtm2(:),dvsqtm3(:) &
       & ,dvbcrt2(:),dvbcrt3(:),beta(:),mmw(:)
  real(8) :: fmaxc, re1c ,sgmc 
end module vars

program main
  use vars
  implicit none

  real :: spvar,Temperature,mu,xe,Mdot,fcMd,rsq,Mdot_ave
  real :: vrot0,Tmfrac
  real :: pfsum = 0.0
  real :: dum1,dum2,dum3
  integer :: ip,itime!,icmpl
  real(8) :: t,fltb
  integer :: pgopen
  integer, parameter :: nMdav = 100
  integer :: datanum,irec
  character*3 :: cnum 
  character*100 ::inputfile,datafile

  inputfile = 'InputParameters_MHDstwind'
  datafile = 'printvars.dat'
  
  open(unit= 8,file=inputfile,status='old' )
  read(8,*)  
  read(8,*)  np             ! Number of Grid points
  read(8,*)  
  read(8,*)  
  read(8,*)  dum1,dum2,gamma 
  read(8,*)  
  read(8,*)  Mstar,Teff,rhocgs,Rstar,dum2,zmetal 
  read(8,*)  fmaxc
  read(8,*)  vrot0

  rdkm = Rstar * Rsun
  nd = np + 1
  gammi1 = gamma - 1.0
  
  YHe = Ysun + dYdZ*(zmetal - 1.d0)*Zsun
  XHyd = 1.d0 - YHe - zmetal*Zsun
  YoX = YHe / XHyd
  mufull = 1.d0 / (2.d0*XHyd + 0.75d0*YHe + 0.45d0*zmetal*Zsun)
  muzero = 1.d0 / (XHyd + 0.25d0*YHe + 0.05d0*zmetal*Zsun)
  mufrac = mufull / muzero
  Tmfrac = Thigh - Tlow*mufrac
  coefT1 = (Thigh - Tlow) / Tmfrac
  coefT0 = Thigh*Tlow*(1.d0 - mufrac) / Tmfrac
  dnnH1 = rhocgs *XHyd / prtnms
  clmfct = rdkm * 1.d5 * dnnH1 !rhocgs / (1.4d0*prtnms) ! mass per H = 1.4  
!  write(*,*)clmfct,rdkm,dnnH1  

  call Teff2akm2
  call Teff2akm_sun
  
!  taucsc = Rstar/Mstar * akm*akm / (akmsun*akmsun) ! normalization for wave period
  taucsc = Rstar/Mstar * Teff / Teffsun ! normalization for wave period
  sgmc = sgmcsun * taucsc 
  re1c = 1.d0 + 2.d0*sgmc ! H/R + R

  write(*,*)"YHe=",YHe,"nufull=",mufull,"muzero=",muzero,"taucsc=",taucsc&
       & ,"re1c=",re1c,"sgmc=",sgmc
  
  spetkl = 1.e10* akm*akm * mufull/Rcon * gammi1
  fcMd = 4.e15/2.e33*3.15576e7 * pi * rhocgs * akm *rdkm*rdkm / real(nMdav + 1) 

  allocate(var(nd,nmvar),varave(nd,nmave)) 
  allocate(rd(np),zd(np),dn(np),vl(np),Tprt(np),Bt(np),eta(np),ead(nd),mul(nd)&
       & ,vt(nd),vt3(np),Bt3(np),Br(np),BtBr_tav(np),vrot(np)&
       & ,vtave(np),vtave3(np),Btave(np),Btave3(np),va(np),dr(nd) &
       & ,dvsqtm2(np),dvsqtm3(np),dvbcrt2(np),dvbcrt3(np),beta(np),mmw(np))

  inquire(iolength=irec)t,var,varave
  open(unit=11,file=datafile,access='direct',form='unformatted',recl=irec)

  IF (pgopen('?') <= 0) STOP
  do itime = 1,ntime
     read(11,rec=itime)t,var,varave
     write(*,*)itime,t 
     do ip = 1,np
        rd(ip) = log10(real(var(ip,1)))
        zd(ip) = log10(real(var(ip,1))-1.e0)
        dr(ip) = real(var(ip,2))
        dn(ip) = log10(real(var(ip,3))*rhocgs)
        vl(ip) = real(var(ip,4))*akm
        vt3(ip) = real(var(ip,6))*akm
        Bt3(ip) = real(var(ip,8))
        vt(ip) = real(var(ip,5))*akm
        Bt(ip) = real(var(ip,7)/var(ip,3)**0.25)
        Br(ip) = varave(ip,1) * real(sqrt(var(ip,3)))
        va(ip) = varave(ip,1) !* akm
        beta(ip) = 2.0*gammi1*real(var(ip,9))/va(ip)/va(ip)
        Tprt(ip) = log10(real(var(ip,10)))

        mmw(ip) = varave(ip,4)
        vtave(ip) = varave(ip,7)*akm
        vtave3(ip) = varave(ip,8)*akm
        Btave(ip) = varave(ip,9)/sqrt(varave(ip,5))*akm
        Btave3(ip) = varave(ip,10)/sqrt(varave(ip,5))*akm
        dvsqtm2(ip) = vtave(ip)*vtave(ip) + Btave(ip)*Btave(ip)
        dvsqtm3(ip) = vtave3(ip)*vtave3(ip) + Btave3(ip)*Btave3(ip)
        dvbcrt2(ip) = -2.0*varave(ip,12)/sqrt(varave(ip,5))*akm*akm
        dvbcrt3(ip) = -2.0*varave(ip,13)/sqrt(varave(ip,5))*akm*akm
        vrot(ip) = real(var(ip,1)) * vrot0 * akm
        BtBr_tav(ip) = (vtave3(ip) - vrot(ip) ) / (sign(1.0,vl(ip)) &
             & *max(abs(vl(ip)),1.e-12))
     end do
     eta(1:np) = log10(varave(1:np,2))
     ead(1:np) = log10(varave(1:np,3))
     mul(1:np) = log10(varave(1:np,4))

     Mdot = 0.0
     Mdot_ave = 0.0
     do ip = np-1000-nMdav,np-1000
!     do ip = np*2/3-100-nMdav,np*2/3-100
        rsq = real(var(ip,1)*var(ip,1))
!        Mdot = Mdot + real(var(ip,3) * var(ip,4) * var(ip,1)*var(ip,1))
        Mdot = Mdot + real(var(ip,3) * var(ip,4)) * rsq
        Mdot_ave = Mdot_ave + varave(ip,5)*varave(ip,6) * rsq 
     end do
     Mdot = Mdot * fcMd
     Mdot_ave = Mdot_ave * fcMd

     call movie(real(t),Mdot,Mdot_ave)
  end do
  call pgclos

end program main

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
  tbalmer = 3000.d0!max(3000.d0, 0.85d0*Teff) 
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
  write(*,*)mu,akm
  
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
  tbalmer = 3000.d0!max(3000.d0, 0.85d0*Teffsun) 
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
  
!---------------------------------------------------------------------
function fltb(r)
  USE vars
  implicit none   

  real(8) :: fltb
  real(8),intent(in) :: r
  real(8) :: ps,f1,ex,fltbc

! This subroutine calculate area expansion factor, f, on r.

! In the corona
  ps = (r-re1)/sgm
  if(ps>10.d0)then
     fltb=fmax
  else
     f1 = 1.d0 - (fmax-1.d0)*exp((1.d0-re1)/sgm)
     ex = exp(ps)
     fltb=(fmax*ex+f1)/(ex+1.)
  end if

! In the chromosphere
  ps=(r-re1c)/sgmc
  if(ps>10.d0)then
     fltbc=fmaxc
  else
     f1 = 1.d0 - (fmaxc-1.d0)*exp((1.d0-re1c)/sgmc)
     ex = dexp(ps)
     fltbc = (fmaxc*ex+f1) / (ex+1.d0)
  end if

  fltb = fltb * fltbc

end function fltb

!--------------------------------------------------------------------
subroutine movie(t,Mdot,Mdot_ave)
  use vars
  implicit none

  real,intent(in) :: t,Mdot,Mdot_ave
  character*8 ::time
  character*12 :: massloss,masslossave

  write(time,'(f8.4)')t
  write(massloss,'(1pe12.4)')Mdot
  write(masslossave,'(1pe12.4)')Mdot_ave

  call pgpage
  call pgslw(2)
  call pgsch(1.3)
  call pgsvp(0.07,0.32,0.55,0.95)       
  call pgmtxt('T',0.5,0.6,0.,time)
!  call pgswin(rd(1),rd(np-70),3.,6.3)
  call pgswin(zd(3),zd(np-30),2.7,6.3)
  call pgbox('BLNST',0.,1,'BLNTS',0.,1)
  call pgmtxt('L',2.,0.49,0.,'T')
!  call pgline(np,rd,Tp)
  call pgline(np,zd,Tprt)

  call pgsvp(0.07,0.32,0.55,0.95)       
  call pgswin(zd(3),zd(np-30),-3.,1.5)
  call pgbox('CLST',0.,1,'CLMTS',0.,1)
  call pgmtxt('R',1.5,0.56,0.,'\gb')
  call pgsls(2)
  call pgline(np,zd,log10(beta))
  
  call pgsls(1)
  call pgslw(2)
  call pgsch(1.3)
  call pgsvp(0.07,0.32,0.1,0.5)       
!  call pgswin(rd(1),rd(np-70),log10(rhocgs)-15.,log10(rhocgs))
  call pgswin(zd(3),zd(np-30),log10(rhocgs)-15.,log10(rhocgs))
  call pgbox('BLNST',0.,1,'BNTS',0.,0)
  call pgmtxt('L',2.,0.49,0.,'log(\gr)')
  call pgmtxt('B',2.,0.3,0.,'r/R\d*\u-1')
!  call pgline(np,rd,dn)
  call pgline(np,zd,dn)

  call pgsvp(0.07,0.32,0.1,0.5)       
  call pgswin(zd(3),zd(np-30),0.5,1.35)
  call pgbox('CLNST',0.,1,'CMTS',0.2,2)
  call pgmtxt('R',1.5,0.48,0.,'\gm')
  call pgsls(2)
  call pgline(np,zd,mmw)

  call pgsls(1)  
  call pgslw(2)
  call pgsch(1.3)
  call pgsvp(0.42,0.67,0.1,0.5)       
  call pgswin(rd(10),rd(np-30),-1.e2,9.e2)
!  call pgswin(zd(10),zd(np-30),-1.e2,1.5e3)
  call pgbox('BCLNST',0.,1,'BCNTS',5.e2,5)
  call pgmtxt('L',1.5,0.2,0.,'v\dr\u & v\dA')
  call pgline(np,rd(:),vl(:))
!  call pgline(np,zd,vl)
  call pgsls(4)
  call pgline(np,rd(:),va(:)*akm)
  call pgmtxt('B',2.,0.4,0.,'r/R\d*')
!  call pgline(np,zd,va)
  call pgsls(1)

  call pgslw(2)
  call pgsch(1.3)
  call pgsvp(0.42,0.67,0.55,0.95) 
  call pgmtxt('T',0.5,0.01,0.,"snap="//massloss)
!  call pgswin(rd(1),rd(np-70),-1.e2,1.e2)
!  call pgswin(rd(1),rd(np-70),0.,4.e2)
  call pgswin(zd(10),zd(np-30),0.,2.e2)
  call pgbox('BCLNST',0.,1,'BCNTS',1.e2,2)
  call pgmtxt('L',1.5,0.38,0.,'v\dt')
!  call pgline(np,rd,vt)
!  call pgline(np,rd,vtave)
  call pgline(np,zd,vtave)
  call pgsls(2)
!  call pgline(nd,rd,Btave)
  call pgline(np,zd,Btave)
  call pgsls(3)
  call pgline(np,zd,vtave3)
  call pgsls(4)
  call pgline(np,zd,Btave3)
  call pgsls(1)

  call pgslw(2)
  call pgsch(1.3)
  call pgsvp(0.74,0.99,0.55,0.95)
  call pgmtxt('T',0.5,0.01,0.,"ave="//masslossave)
  call pgswin(zd(3),zd(np-70),0.,7.)
!  call pgswin(zd(10),zd(np-30),-1.e2,2.e5)
!  call pgswin(1.,1.01,-1.e1,1.e2)
!  call pgbox('BCLNST',0.,1,'BCNTS',2.e4,2)
  call pgbox('BCLNST',0.,1,'BCLNTS',0.,5)
  call pgmtxt('L',2.,0.1,0.,'dB\u2\d+dv\u2\d,-2dv*dB,')
  call pgsls(3)
!  call pgline(np,zd,dvsqtm2)
  call pgsls(4)
!  call pgline(np,zd,dvsqtm3)
  call pgsls(1)
!  call pgline(np,zd,dvbcrt2)
  call pgsls(2)
!  call pgline(np,zd,dvbcrt3)
  call pgsls(1)
  call pgline(np,zd,log10(dvsqtm2+dvsqtm3+dvbcrt2+dvbcrt3))
  call pgsls(2)
  call pgline(np,zd,log10(dvsqtm2+dvsqtm3-dvbcrt2-dvbcrt3))

  call pgsls(1)

  call pgslw(2)
  call pgsch(1.3)
  call pgsvp(0.74,0.99,0.1,0.5) 
!  call pgswin(rd(1),rd(np-70),-3.e0,1.e0)
!  call pgswin(zd(2),zd(np-30),-1.e0,3.5e0)
  call pgswin(zd(3),zd(np-30),-1.e0,3.5e0) 
  call pgbox('BCLNST',0.,1,'BCLNTS',2.e0,2)
  call pgmtxt('L',2.,0.49,0.,'v\dA\u*\gt/dr')
  call pgmtxt('B',2.,0.3,0.,'r/R\d*\u-1')
!  call pgline(np,rd,Bt3/Br)
!  call pgline(np,zd,Bt3/Br)
  call pgline(np,zd(:),log10(va(:)/dr(:)*2.1e-4*taucsc))
  call pgsls(2)
  call pgline(np,zd(:),log10(max(va(:)+vl(:),va(:))/dr(:)*2.1e-4*taucsc))
  call pgsls(4)
  call pgsch(0.0)
  call pgarro(-3.,0.,3.,0.)
  call pgsls(1)

end subroutine movie
