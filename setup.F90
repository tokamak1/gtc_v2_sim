!========================================================================

module particle_decomp

!========================================================================
  integer  :: ntoroidal,npartdom
  integer  :: toroidal_comm,partd_comm
  integer  :: nproc_partd,myrank_partd
  integer  :: nproc_toroidal,myrank_toroidal
  integer  :: left_pe,right_pe
  integer  :: toroidal_domain_location,particle_domain_location
end module particle_decomp


!========================================================================

    Subroutine setup

!========================================================================

  use global_parameters
  use particle_decomp
  use particle_array
  use field_array
  use diagnosis_array
  use particle_tracking
  implicit none

  integer i,j,k,ierror,ij,mid_theta,ip,jt,indp,indt,mtest,micell
  integer mi_local
  real(wp) r0,b0,temperature,tdum,r,q,sint,dtheta_dx,rhoi,b,zdum,&
       edensity0,delr,delt,rmax,rmin,wt
  CHARACTER(LEN=10) date, time
  namelist /run_parameters/ numberpe,mi,mgrid,mid_theta,mtdiag,delr,delt,&
       ulength,utime,gyroradius

#define FLUSH flush

! total # of PE and rank of PE
  call mpi_comm_size(mpi_comm_world,numberpe,ierror)
  call mpi_comm_rank(mpi_comm_world,mype,ierror)

! Read the input file that contains the run parameters
  call read_input_params(micell,r0,b0,temperature,edensity0)

  irest=1
  FileExit=0

! numerical constant
  pi=4.0_wp*atan(1.0_wp)
  mstep=max(2,mstep)
  msnap=min(msnap,mstep/ndiag)
  isnap=mstep/msnap
  idiag1=mpsi/2
  idiag2=mpsi/2
  if(nonlinear < 0.5)then
     paranl=0.0_wp
     mode00=0
     idiag1=1
     idiag2=mpsi
  endif
  rc=rc*(a0+a1)
  rw=1.0_wp/(rw*(a1-a0))

! Set up the particle decomposition within each toroidal domain
  call set_particle_decomp

! equilibrium unit: length (unit=cm) and time (unit=second) unit
  ulength=r0
  utime=1.0_wp/(9580._wp*b0) ! time unit = inverse gyrofrequency of proton 
! primary ion thermal gyroradius in equilibrium unit, vthermal=sqrt(T/m)
  gyroradius=102.0_wp*sqrt(aion*temperature)/(abs(qion)*b0)/ulength
  tstep=tstep*aion/(abs(qion)*gyroradius*kappati)
  
! allocate memory
  allocate (qtinv(0:mpsi),itran(0:mpsi),mtheta(0:mpsi),&
     deltat(0:mpsi),rtemi(0:mpsi),pfluxpsi(0:mpsi),rdtemi(0:mpsi),&
     rden(0:mpsi),igrid(0:mpsi),pmarki(0:mpsi),phi00(0:mpsi),phip00(0:mpsi),&
     hfluxpsi(0:mpsi),zonali(0:mpsi),gradt(mpsi),&
     eigenmode(m_poloidal,num_mode,mpsi),STAT=mtest)
  if (mtest /= 0) then
     write(0,*)mype,'*** Cannot allocate qtinv: mtest=',mtest
     call MPI_ABORT(MPI_COMM_WORLD,1,ierror)
  endif

! --- Define poloidal grid ---
! grid spacing
  deltar=(a1-a0)/real(mpsi)

! grid shift associated with fieldline following coordinates
  tdum=pi*a/real(mthetamax)
  do i=0,mpsi
     r=a0+deltar*real(i)
     mtheta(i)=2*max(1,int(pi*r/tdum+0.5_wp))
     deltat(i)=2.0_wp*pi/real(mtheta(i))
     q=q0+q1*r/a+q2*r*r/(a*a)
     itran(i)=int(real(mtheta(i))/q+0.5_wp)
     qtinv(i)=real(mtheta(i))/real(itran(i)) !q value for coordinate transformation
     qtinv(i)=1.0/qtinv(i) !inverse q to avoid divide operation
     itran(i)=itran(i)-mtheta(i)*(itran(i)/mtheta(i))
  enddo
! un-comment the next two lines to use magnetic coordinate
!  qtinv=0.0
!  itran=0

! When doing mode diagnostics, we need to switch from the field-line following
! coordinates alpha-zeta to a normal geometric grid in theta-zeta. This
! translates to a greater number of grid points in the zeta direction, which
! is mtdiag. Precisely, mtdiag should be mtheta/q but since mtheta changes
! from one flux surface to another, we use a formula that gives us enough
! grid points for all the flux surfaces considered.
  mtdiag=(mthetamax/mzetamax)*mzetamax
  mthetamax=mtheta(mpsi)

! starting point for a poloidal grid
  igrid(0)=1
  do i=1,mpsi
     igrid(i)=igrid(i-1)+mtheta(i-1)+1
  enddo

! number of grids on a poloidal plane
  mgrid=sum(mtheta+1)
  mi_local=micell*(mgrid-mpsi)*mzeta          !# of ions in toroidal domain
  mi=micell*(mgrid-mpsi)*mzeta/npartdom       !# of ions per processor
  if(mi<mod(mi_local,npartdom))mi=mi+1
  mimax=mi+100*ceiling(sqrt(real(mi))) !ions array upper bound

  delr=deltar/gyroradius
  delt=deltat(mpsi/2)*(a0+deltar*real(mpsi/2))/gyroradius
  mid_theta=mtheta(mpsi/2)
  if(mype == 0) then
	write(stdout,run_parameters)
	if(stdout /= 6 .and. stdout /= 0)close(stdout)
  end if	

! allocate memory
  allocate(pgyro(4,mgrid),tgyro(4,mgrid),markeri(mzeta,mgrid),&
     densityi(0:mzeta,mgrid),phi(0:mzeta,mgrid),evector(3,0:mzeta,mgrid),&
     jtp1(2,mgrid,mzeta),jtp2(2,mgrid,mzeta),wtp1(2,mgrid,mzeta),&
     wtp2(2,mgrid,mzeta),dtemper(mgrid,mzeta),heatflux(mgrid,mzeta),&
     STAT=mtest)
  pgyro=0.0
  tgyro=0.0
  markeri=0.0
  densityi=0.0
  phi=0.0
  evector=0.0
  jtp1=0.0
  jtp2=0.0
  wtp1=0.0
  wtp2=0.0
  dtemper=0.0
  heatflux=0.0
  ptracer=0.0
  etracer=0.0

  if (mtest /= 0) then
     write(0,*)mype,'*** setup: Cannot allocate pgyro: mtest=',mtest
     call MPI_ABORT(MPI_COMM_WORLD,1,ierror)
  endif

! temperature and density on the grid, T_i=n_0=1 at mid-radius
  rtemi=1.0
  rden=1.0
  phi=0.0
  phip00=0.0
  pfluxpsi=0.0
  rdtemi=0.0
  zonali=0.0
 
! # of marker per grid, Jacobian=(1.0+r*cos(theta+r*sin(theta)))*(1.0+r*cos(theta))
  pmarki=0.0
!$omp parallel do private(i,j,k,r,ij,zdum,tdum,rmax,rmin)
  do i=0,mpsi
     r=a0+deltar*real(i)
     do j=1,mtheta(i)
        ij=igrid(i)+j
        do k=1,mzeta
           zdum=zetamin+real(k)*deltaz
           tdum=real(j)*deltat(i)+zdum*qtinv(i)
           markeri(k,ij)=(1.0+r*cos(tdum))**2
           pmarki(i)=pmarki(i)+markeri(k,ij)
        enddo
     enddo
     rmax=min(a1,r+0.5*deltar)
     rmin=max(a0,r-0.5*deltar)
     tdum=real(mi*npartdom)*(rmax*rmax-rmin*rmin)/(a1*a1-a0*a0)
     do j=1,mtheta(i)
        ij=igrid(i)+j
        do k=1,mzeta
           markeri(k,ij)=tdum*markeri(k,ij)/pmarki(i)
           markeri(k,ij)=1.0/markeri(k,ij) !to avoid divide operation
        enddo
     enddo
     pmarki(i)=1.0/(real(ntoroidal)*tdum)
     markeri(:,igrid(i))=markeri(:,igrid(i)+mtheta(i))
  enddo

  if(track_particles == 1)then
     nparam=8
     nspec=1
     allocate(ptrackedi(nparam,max(mimax,1)))
     ptrackedi=0.0
     ntrackp=0
  else
     nparam=6

  endif

! allocate memory
  allocate(zion(nparam,mimax),zion0(nparam,mimax),jtion0(4,mimax),&
     jtion1(4,mimax),kzion(mimax),wzion(mimax),wpion(4,mimax),&
     wtion0(4,mimax),wtion1(4,mimax),STAT=mtest)
!XY intitialize
  zion=0.0
  zion0=0.0
  jtion0=0.0
  jtion1=0.0
  kzion=0.0
  wzion=0.0
  wpion=0.0
  wtion0=0.0
  wtion1=0.0
  ddeni=0.0
  if (mtest /= 0) then
     write(0,*)mype,'*** Cannot allocate zion: mtest=',mtest
     call MPI_ABORT(MPI_COMM_WORLD,1,ierror)
  endif

! 4-point gyro-averaging for sqrt(mu)=gyroradius on grid of magnetic coordinates
! rho=gyroradius*sqrt(2/(b/b_0))*sqrt(mu/mu_0), mu_0*b_0=m*v_th^2
! dtheta/delta_x=1/(r*(1+r*cos(theta))), delta_x=poloidal length increase
!$omp parallel do private(i,j,r,ij,tdum,q,b,dtheta_dx,rhoi)
  do i=0,mpsi
     r=a0+deltar*real(i)
     do j=0,mtheta(i)
        ij=igrid(i)+j
        tdum=deltat(i)*real(j)
        q=q0+q1*r/a+q2*r*r/(a*a)
        b=1.0/(1.0+r*cos(tdum))
        dtheta_dx=1.0/r
! first two points perpendicular to field line on poloidal surface            
        rhoi=sqrt(2.0/b)*gyroradius
        pgyro(1,ij)=-rhoi
        pgyro(2,ij)=rhoi
! non-orthorgonality between psi and theta: tgyro=-rhoi*dtheta_dx*r*sin(tdum)
        tgyro(1,ij)=0.0
        tgyro(2,ij)=0.0

! the other two points tangential to field line
        tgyro(3,ij)=-rhoi*dtheta_dx
        tgyro(4,ij)=rhoi*dtheta_dx
        pgyro(3:4,ij)=rhoi*0.5*rhoi/r
     enddo
  enddo

! initiate radial interpolation for grid
  do k=1,mzeta
     zdum=zetamin+deltaz*real(k)
!$omp parallel do private(i,ip,j,indp,indt,ij,tdum,jt,wt)
     do i=1,mpsi-1
        do ip=1,2
           indp=min(mpsi,i+ip)
           indt=max(0,i-ip)
           do j=1,mtheta(i)
              ij=igrid(i)+j
! upward
              tdum=(real(j)*deltat(i)+zdum*(qtinv(i)-qtinv(indp)))/deltat(indp)
              jt=floor(tdum)
              wt=tdum-real(jt)
              jt=mod(jt+mtheta(indp),mtheta(indp))
              if(ip==1)then
                 wtp1(1,ij,k)=wt
                 jtp1(1,ij,k)=igrid(indp)+jt
              else
                 wtp2(1,ij,k)=wt
                 jtp2(1,ij,k)=igrid(indp)+jt
              endif
! downward
               
              tdum=(real(j)*deltat(i)+zdum*(qtinv(i)-qtinv(indt)))/deltat(indt)
              jt=floor(tdum)
              wt=tdum-real(jt)
              jt=mod(jt+mtheta(indt),mtheta(indt))
              if(ip==1)then
                 wtp1(2,ij,k)=wt
                 jtp1(2,ij,k)=igrid(indt)+jt
              else
                 wtp2(2,ij,k)=wt
                 jtp2(2,ij,k)=igrid(indt)+jt
              endif
           enddo
        enddo
     enddo
  enddo

end subroutine setup


!=============================================================================

  Subroutine read_input_params(micell,r0,b0,temperature,edensity0)

!=============================================================================

  use global_parameters
  use particle_decomp
  use particle_tracking
  use diagnosis_array
  implicit none

  logical file_exist
  integer ierror,micell
  real(wp),intent(INOUT) :: r0,b0,temperature,edensity0
  CHARACTER(LEN=10) date, time

#ifdef _OPENMP
  integer nthreads,omp_get_num_threads
#endif

  namelist /input_parameters/ irun,mstep,msnap,ndiag,nonlinear,paranl,&
       mode00,tstep,micell,mpsi,mthetamax,mzetamax,npartdom,&
       a,a0,a1,q0,q1,q2,rc,rw,&
       aion,qion,kappati,kappate,kappan,tite,flow0,&
       flow1,flow2,r0,b0,temperature,edensity0,stdout,nbound,umax,iload,&
       track_particles,nptrack,rng_control,nmode,mmode
!
! Since it is preferable to have only one MPI process reading the input file,
! we choose the master process to set the default run parameters and to read
! the input file. The parameters will then be broadcast to the other processes.
!

  if(mype==0) then
! Default control parameters
    irun=0                 ! 0 for initial run, any non-zero value for restart
    mstep=1500             ! # of time steps
    msnap=1                ! # of snapshots
    ndiag=4                ! do diag when mod(istep,ndiag)=0
    nonlinear=1.0          ! 1.0 nonlinear run, 0.0 linear run
    paranl=0.0             ! 1: keep parallel nonlinearity
    mode00=1               ! 1 include (0,0) mode, 0 exclude (0,0) mode

! run size (both mtheta and mzetamax should be multiples of # of PEs)
    tstep=0.2              ! time step (unit=L_T/v_th), tstep*\omega_transit<0.1 
    micell=2               ! # of ions per grid cell
    mpsi=90                ! total # of radial grid points
    mthetamax=640          ! poloidal grid, even and factors of 2,3,5 for FFT
    mzetamax=64            ! total # of toroidal grid points, domain decomp.
    npartdom=1             ! number of particle domain partitions per tor dom.
     
! run geometry
    a=0.358                ! minor radius, unit=R_0
    a0=0.1                 ! inner boundary, unit=a
    a1=0.9                 ! outer boundary, unit=a
    q0=0.854               ! q_profile, q=q0 + q1*r/a + q2 (r/a)^2
    q1=0.0
    q2=2.184
    rc=0.5                 ! kappa=exp{-[(r-rc)/rw]**6}
    rw=0.35                ! rc in unit of (a1+a0) and rw in unit of (a1-a0)

! species information
    aion=1.0               ! species isotope #
    qion=1.0               ! charge state

! equilibrium unit: R_0=1, Omega_c=1, B_0=1, m=1, e=1
    kappati=6.9            ! grad_T/T
    kappate=6.9
    kappan=kappati*0.319   ! inverse of eta_i, grad_n/grad_T
    tite=1.0               ! T_i/T_e
    flow0=0.0              ! d phi/dpsi=gyroradius*[flow0+flow1*r/a+
    flow1=0.0              !                              flow2*(r/a)**2]
    flow2=0.0

! physical unit
    r0=93.4                ! major radius (unit=cm)
    b0=19100.0             ! on-axis vacuum field (unit=gauss)
    temperature=2500.0     ! electron temperature (unit=ev)
    edensity0=0.46e14      ! electron number density (1/cm^3)

! standard output: use 0 or 6 to terminal and 11 to file 'stdout.out'
    stdout=0  
    nbound=4               ! 0 for periodic, >0 for zero boundary 
    umax=4.0               ! unit=v_th, maximum velocity in each direction
    iload=0                ! 0: uniform, 1: non-uniform
    track_particles=0      ! 1: keep track of some particles
    nptrack=0              ! track nptrack particles every time step
    rng_control=1          ! controls seed and algorithm for random num. gen.
                           ! rng_control>0 uses the portable random num. gen.

! mode diagnostic: 8 modes.
    nmode=(/5,7,9,11,13,15,18,20,25,30,35,40,45/) ! n: toroidal mode number
    mmode=(/7,10,13,15,18,21,25,28,35,42,49,56,63/) ! m: poloidal mode number

! Test if the input file gtc.input exists
    inquire(file='gtc.input',exist=file_exist)
    if (file_exist) then
       open(55,file='gtc.input',status='old')
       read(55,nml=input_parameters)
       close(55)
    else
       write(0,*)'******************************************'
       write(0,*)'*** NOTE!!! Cannot find file gtc.input !!!'
       write(0,*)'*** Using default run parameters...'
       write(0,*)'******************************************'
    endif

! Changing the units of a0 and a1 from units of "a" to units of "R_0"
    a0=a0*a
    a1=a1*a

! open file for standard output, record program starting time
    if(stdout /= 6 .and. stdout /= 0)open(stdout,file='stdout.out',status='replace')
    call date_and_time(date,time)
    write(stdout,*) 'Program starts at DATE=', date, 'TIME=', time
    write(stdout,input_parameters)

#ifdef _OPENMP
!$omp parallel private(nthreads)
    nthreads=omp_get_num_threads()  !Get the number of threads if using OMP
!$omp single
    write(stdout,'(/,"===================================")')
    write(stdout,*)' Number of OpenMP threads = ',nthreads
    write(stdout,'("===================================",/)')
!$omp end single nowait
!$omp end parallel
#else
    write(stdout,'(/,"===================================")')
    write(stdout,*)' Run without OpenMP threads'
    write(stdout,'("===================================",/)')
#endif
  endif

! Now send the parameter values to all the other MPI processes
  call broadcast_input_params(micell,r0,b0,temperature,edensity0)

end Subroutine read_input_params


!=============================================================================

  Subroutine broadcast_input_params(micell,r0,b0,temperature,edensity0)

!=============================================================================

  use global_parameters
  use particle_tracking
  use particle_decomp
  use diagnosis_array

  integer,parameter :: n_integers=16+2*num_mode,n_reals=25
  integer  :: integer_params(n_integers)
  real(wp) :: real_params(n_reals)
  integer ierror,micell
  real(wp),intent(INOUT) :: r0,b0,temperature,edensity0

! The master process, mype=0, holds all the input parameters. We need
! to broadcast their values to the other processes. Instead of issuing
! an expensive MPI_BCAST() for each parameter, it is better to pack
! everything in a single vector, broadcast it, and unpack it.

  if(mype==0)then
!   Pack all the integer parameters in integer_params() array
    integer_params(1)=irun
    integer_params(2)=mstep
    integer_params(3)=msnap
    integer_params(4)=ndiag
    integer_params(5)=mode00
    integer_params(6)=micell
    integer_params(7)=mpsi
    integer_params(8)=mthetamax
    integer_params(9)=mzetamax
    integer_params(10)=npartdom
    integer_params(11)=stdout
    integer_params(12)=nbound
    integer_params(13)=iload
    integer_params(14)=track_particles
    integer_params(15)=nptrack
    integer_params(16)=rng_control
    integer_params(17:17+num_mode-1)=nmode(1:num_mode)
    integer_params(17+num_mode:17+2*num_mode-1)=mmode(1:num_mode)

!   Pack all the real parameters in real_params() array
    real_params(1)=nonlinear
    real_params(2)=paranl
    real_params(3)=tstep
    real_params(4)=a
    real_params(5)=a0
    real_params(6)=a1
    real_params(7)=q0
    real_params(8)=q1
    real_params(9)=q2
    real_params(10)=rc
    real_params(11)=rw
    real_params(12)=aion
    real_params(13)=qion
    real_params(14)=kappati
    real_params(15)=kappate
    real_params(16)=kappan
    real_params(17)=tite
    real_params(18)=flow0
    real_params(19)=flow1
    real_params(20)=flow2
    real_params(21)=r0
    real_params(22)=b0
    real_params(23)=temperature
    real_params(24)=edensity0
    real_params(25)=umax
  endif

! Send input parameters to all processes
  call MPI_BCAST(integer_params,n_integers,MPI_INTEGER,0,MPI_COMM_WORLD,ierror)
  call MPI_BCAST(real_params,n_reals,mpi_Rsize,0,MPI_COMM_WORLD,ierror)

  if(mype/=0)then
!   Unpack integer parameters
    irun=integer_params(1)
    mstep=integer_params(2)
    msnap=integer_params(3)
    ndiag=integer_params(4)
    mode00=integer_params(5)
    micell=integer_params(6)
    mpsi=integer_params(7)
    mthetamax=integer_params(8)
    mzetamax=integer_params(9)
    npartdom=integer_params(10)
    stdout=integer_params(11)
    nbound=integer_params(12)
    iload=integer_params(13)
    track_particles=integer_params(14)
    nptrack=integer_params(15)
    rng_control=integer_params(16)
    nmode(1:num_mode)=integer_params(17:17+num_mode-1)
    mmode(1:num_mode)=integer_params(17+num_mode:17+2*num_mode-1)

!   Unpack real parameters
    nonlinear=real_params(1)
    paranl=real_params(2)
    tstep=real_params(3)
    a=real_params(4)
    a0=real_params(5)
    a1=real_params(6)
    q0=real_params(7)
    q1=real_params(8)
    q2=real_params(9)
    rc=real_params(10)
    rw=real_params(11)
    aion=real_params(12)
    qion=real_params(13)
    kappati=real_params(14)
    kappate=real_params(15)
    kappan=real_params(16)
    tite=real_params(17)
    flow0=real_params(18)
    flow1=real_params(19)
    flow2=real_params(20)
    r0=real_params(21)
    b0=real_params(22)
    temperature=real_params(23)
    edensity0=real_params(24)
    umax=real_params(25)
  endif


end subroutine broadcast_input_params


!=============================================================================

    Subroutine set_particle_decomp

!=============================================================================

  use global_parameters
  use particle_decomp
  implicit none

  integer :: ierror

! ----- First we verify the consistency of ntoroidal and npartdom -------
! The number of toroidal domains (ntoroidal) times the number of particle
! "domains" (npartdom) needs to be equal to the number of processor "numberpe".
! numberpe cannot be changed since it is given on the command line.

! numberpe must be a multiple of npartdom so change npartdom accordingly
  do while (mod(numberpe,npartdom) /= 0)
     npartdom=npartdom-1
     if(npartdom==1)exit
  enddo
  ntoroidal=numberpe/npartdom
  if(mype==0)then
    write(stdout,*)'*******************************************************'
    write(stdout,*)'  Using npartdom=',npartdom,' and ntoroidal=',ntoroidal
    write(stdout,*)'*******************************************************'
    write(stdout,*)
  endif

! make sure that mzetamax is a multiple of ntoroidal
  mzetamax=ntoroidal*max(1,int(real(mzetamax)/real(ntoroidal)+0.5))

! Make sure that "mpsi", the total number of flux surfaces, is an even
! number since this quantity will be used in Fast Fourier Transforms
  mpsi=2*(mpsi/2)

! We now give each PE (task) a unique domain identified by 2 numbers: the
! particle and toroidal domain numbers.
!    particle_domain_location = rank of the particle domain holding mype
!    toroidal_domain_location = rank of the toroidal domain holding mype
! 
! Assign successive ranks to particle domains first. This keeps the particle
! decomposition local within each toroidal domain when possible.

  particle_domain_location=mod(mype,npartdom)
  toroidal_domain_location=mype/npartdom

! Domain decomposition in toroidal direction.
  mzeta=mzetamax/ntoroidal
  zetamin=2.0_wp*pi*real(toroidal_domain_location)/real(ntoroidal)
  zetamax=2.0_wp*pi*real(toroidal_domain_location+1)/real(ntoroidal)

! grid spacing in the toroidal direction
  deltaz=(zetamax-zetamin)/real(mzeta)

! ---- Create particle domain communicator and toroidal communicator -----
! We now need to create a new communicator which will include only the
! processes located in the same toroidal domain. The particles inside
! each toroidal domain are divided equally between "npartdom" processes.
! Each one of these processes will do a charge deposition on a copy of
! the same grid, requiring a toroidal-domain-wide reduction after the
! deposition. The new communicator will allow the reduction to be done
! only between those processes located in the same toroidal domain.
!
! We also need to create a purely toroidal communicator so that the
! particles with the same particle domain id can exchange with their
! toroidal neighbors.
!
! Form 2 subcommunicators: one that includes all the processes located in
! the same toroidal domain (partd_comm), and one that includes all the
! processes part of the same particle domain (toroidal_comm).
! Here is how to create a new communicator from an old one by using
! the MPI call "MPI_COMM_SPLIT()".
! All the processes passing the same value of "color" will be placed in
! the same communicator. The "rank_in_new_comm" value will be used to
! set the rank of that process on the communicator.

! particle domain communicator (for communications between the particle
! domains WITHIN the same toroidal domain)
  call MPI_COMM_SPLIT(MPI_COMM_WORLD,toroidal_domain_location,&
                      particle_domain_location,partd_comm,ierror)

! toroidal communicator (for communications BETWEEN toroidal domains of same
! particle domain number)
  call MPI_COMM_SPLIT(MPI_COMM_WORLD,particle_domain_location,&
                      toroidal_domain_location,toroidal_comm,ierror)

  call mpi_comm_size(partd_comm,nproc_partd,ierror)
  call mpi_comm_rank(partd_comm,myrank_partd,ierror)

  call mpi_comm_size(toroidal_comm,nproc_toroidal,ierror)
  call mpi_comm_rank(toroidal_comm,myrank_toroidal,ierror)

  if(nproc_partd/=npartdom)then
    write(0,*)'*** nproc_partd=',nproc_partd,' NOT EQUAL to npartdom=',npartdom
    call MPI_ABORT(MPI_COMM_WORLD,1,ierror)
  endif

  if(nproc_toroidal/=ntoroidal)then
    write(*,*)'*** nproc_toroidal=',nproc_toroidal,' NOT EQUAL to ntoroidal=',&
              ntoroidal
    call MPI_ABORT(MPI_COMM_WORLD,1,ierror)
  endif

! We now find the toroidal neighbors of the current toroidal domain and
! store that information in 2 easily accessible variables. This information
! is needed several times inside the code, such as when particles cross
! the domain boundaries. We will use the toroidal communicator to do these
! transfers so we don't need to worry about the value of myrank_partd.
! We have periodic boundary conditions in the toroidal direction so the
! neighbor to the left of myrank_toroidal=0 is (ntoroidal-1).

  left_pe=mod(myrank_toroidal-1+ntoroidal,ntoroidal)
  right_pe=mod(myrank_toroidal+1,ntoroidal)

end subroutine set_particle_decomp
