!     **************************************************************
!     * MUSCL for MHD wave propagation in Spherical Coordinate      *
!     **************************************************************
!
! MHD waves are injected from the left boundary (ip=1) 
! in subroutine 'evolve'
! MHD Riemann solver (Only WITH B_perp)
! Enerby eq. : Internal Energy
module vars
  use nrtype
  implicit none
!  integer,parameter :: nd=30000,ndb=4000 
!  ------------------------------------------------------------------
  real(DP) :: fc1,fc2,fc3
  real(DP),parameter :: Gslm=1.327d26,Rsun= 6.96d5, cdcon=1.d-6 &
       & ,Rcon=8.26d7,prtnms=1.67d-24,Teffsun=5780.d0      !pi=3.141529d0,Rcon=kB/mH 
  real(DP),parameter :: fmax=1.d0 ,re1=1.2d0,sgm=0.15d0 
  real(DP) :: fmaxc, re1c ,sgmc
  real(DP) :: Mstar,Rstar,zmetal
!  real(DP) :: cofrad(1:4) = (/-127.3d0,56.57d0,-9.84d0,0.5548d0/) 
!                         !solar-metal
!  real(DP) :: cofrad(1:4) = (/-119.2d0,53.44d0,-9.628d0,0.5633d0/) 
                         !zero-metal 
  real(DP),parameter :: unity = 1.d0, onovth = 1.d0/3.d0
  real(DP),parameter :: Ysun = 0.2485d0, Zsun = 0.0134d0, dYdZ = 0.125d0
  real(DP),parameter :: Thigh = 1.0d4, Tlow = 4.d3, sgmcsun = 0.015d0
  real(DP) :: drin, drout, einit, rhoturn, drfac, drcnst, dtcd0, Tclso
  real(DP) ::  gamma,gammi1,dxlagc,dxlagi,facdep,critit !ql
  real(DP) ::  critif,rtgamm,gammp1,gammx1,gammx2,gammp3,gammm3,gammm2
  real(DP) ::  YoX,YHe,XHyd,mufull,muzero,mufrac,coefT1,coefT0
  integer(I4B) ::  np,nd,ipcalc,iwv,icont,irss,ibeam,irstr,itrb 
  Integer(I4B) :: np_ori,nd_ori
  real(DP) :: tprnt0
  integer(I4B),parameter :: nmvar=10,nmave=13
  real(DP) :: cflnum,t,dt,tfinal
  real(DP),parameter:: t_e2T_rss = 0.05d0
  real(DP),allocatable :: radius(:),dradiu(:),xplace(:),dxplac(:),ddxdtt(:) &
       & ,densit(:),volume(:),pressu(:),speene(:),energy(:) &
       & ,sounds(:),dxdt(:),dvoldt(:),denedt(:) &
       & ,gradd(:), gradp(:), gradv(:), grade(:) &
       & ,xlagra(:),dxlagr(:), prsgas(:) &
       & ,Bparal(:),Alfven(:),Btense(:),Bcenrf(:),Btrans(:),vtrans(:)&
       & ,dBtran(:),dvtran(:),dtdBtr(:),dtvrtr(:),fBtrbd(:),vftrbd(:)&
       & ,Btrant(:),vtrant(:),dBtrat(:),dvtrat(:),dtdBtt(:),dtvrtt(:) &
       & ,fBttbd(:),vfttbd(:),clocal(:)
  real(DP),allocatable :: radbnd(:),dxdtbd(:),rsqfbnd(:)
  real(DP),allocatable :: Btrrbd(:),Btrlbd(:),Bttrbd(:),Bttlbd(:)
  real(DP),allocatable :: xplori(:),dxplor(:),radori(:),drador(:)
  real(DP),allocatable :: dxlg(:)
  real(DP),allocatable :: fBtreb(:),fBtteb(:),vftreb(:),vftteb(:),dxdteu(:)
  real(DP),allocatable :: eta(:),ead(:),ettrb(:),Tpprt(:)
  real(SP),allocatable :: etaprt(:),eadprt(:),muprt(:) !,xeprt(:)
  real(DP) :: dinit,vinit,pinit,ddens,dvelo,dpres !,omega
  real(DP) :: Bplint,Btrint,dBtrit,dvtrit,dBttit,dvttit,wvnm !,omgtr
  real(DP) :: grav,cond,akm,rdkm,Teff,dnnH1,clmfct,etafct,eadfct,akmsun
  real(DP) :: etamax = 1.d-2
  real(DP) :: radfac,spetkl,rhocgs,rhocrt,Tcrt,lambdacrt !,radfac0
  INTEGER(I4B),parameter:: nradcol = 91
  real(DP) :: Ttab(nradcol),radtab(nradcol)
  real(DP) :: Ttab0,dTtab
  real(DP) :: dxplac0,densit0,pressu0,dxdt0,speene0,energy0,Btrans0,vtrans0 &
       & ,Btrant0,vtrant0
  INTEGER(I4B),parameter :: nrand = 262144
  REAL(DP),parameter :: trandsun = 2.1d-4 ! Pmin=18sec for akm_Rsun=8.5 
  real(DP) :: trd(nrand),random(nrand),randlg(nrand),rand3c(nrand) !&
!       & ,trd3c(nrand),trdlg(nrand),
  real(DP),allocatable :: dxdtav(:),densav(:),enrgav(:),prgsav(:),spenav(:) &
       & ,Btraav(:),vtraav(:),Btrtav(:),vtrtav(:),zplus(:),zmins(:) &
       & ,dvbaav(:),dvbtav(:) 
  real(DP) :: dtavtt
  real(DP) :: vxout1,vxout2,ptout1,ptout2,dnout1,dnout2,pgout1,pgout2 &
       &,Btrbd1,Btrbd2,vtrbd1,vtrbd2,Bttbd1,Bttbd2,vttbd1,vttbd2
  real(DP),allocatable :: radcol(:)
  REAL(DP) :: radems(5)

  real(DP) :: vsum,vmin,vmax,rvsum,rvmin,rvmax,Tsum,Tmin,Tmax,tmsum
  real(DP),allocatable :: Falout(:),Falin(:),Fcsout(:),Faloav(:),Faliav(:) &
       & ,Fcsoav(:),Falout3(:),Falin3(:),Faloav3(:),Faliav3(:)
  real(DP),allocatable :: ef(:,:) !,efav(:,:)
  real(DP) :: tefsum !,tefave
  real(DP),allocatable :: Btrprv(:),Bttprv(:),engtot(:)
  real(DP) :: vrot0
  real(DP) :: dslmbd0
  real(DP),parameter :: spemin = 1.d-1,spemax=1.d4

end module vars
