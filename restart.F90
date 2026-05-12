subroutine restart_io(iop)
  use global_parameters
  use particle_array
  use field_array
  use diagnosis_array
  use particle_decomp
  use particle_tracking
  implicit none

  integer merror
  integer i,j
  integer restart_mflux,restart_mpsi,restart_mzeta,restart_mi,restart_mimax
  integer mquantity,mflx,n_mode,mstepfinal,noutputs
  real(wp) dum
  character(*),intent(in)::iop
  character(len=50) :: restart_fname
  real(wp),dimension(:),allocatable::zion0_read
  logical sheareb_open

  if(iop/="read" .and. iop/="write")then
     write(*,*)'*** subroutine restart_io (iop <> "read" or "write")',iop
     call MPI_ABORT(MPI_COMM_WORLD,1,merror)
     return
  endif

  if(mod(irest,2)==0)then
     write(restart_fname,'(a,i5.5,".bp")')"restart_dir1/restart_",myrank_toroidal
  else
     write(restart_fname,'(a,i5.5,".bp")')"restart_dir2/restart_",myrank_toroidal
  endif

  if(iop=="write")then
     open(901,file=trim(restart_fname),form='unformatted',access='stream',status='replace')
     write(901)mflux,mpsi,mzeta,mi,mimax,mgrid,nparam,ntracer,etracer
     write(901)rdtemi,ptracer,pfluxpsi,phi00,phip00,zonali
     write(901)zion0(6,:)
     write(901)phi
     write(901)zion
     close(901)

     if(mype==0 .and. istep<=mstep)then
        write(*,*)'written ',trim(restart_fname)

        if(mod(irest,2)==0)then
           restart_fname="restart_dir1/history_restart.out"
        else
           restart_fname="restart_dir2/history_restart.out"
        endif

        flush(ihistory)
        open(776,file='history.out',status='old',action='read')
        open(777,file=trim(restart_fname),status='replace')
        read(776,101)j
        write(777,101)j
        read(776,101)mquantity
        write(777,101)mquantity
        read(776,101)mflx
        write(777,101)mflx
        read(776,101)n_mode
        write(777,101)n_mode
        read(776,101)mstepfinal
        noutputs=mstepfinal-mstep/ndiag+istep/ndiag
        write(777,101)noutputs
        do i=0,(mquantity+4*n_mode)*noutputs
           read(776,102)dum
           write(777,102)dum
        enddo
        close(776)
        close(777)
        write(*,*)'written ',trim(restart_fname)

        if(mod(irest,2)==0)then
           restart_fname="restart_dir1/sheareb_restart.out"
        else
           restart_fname="restart_dir2/sheareb_restart.out"
        endif

        inquire(unit=444,opened=sheareb_open)
        if(sheareb_open)flush(444)
        open(778,file='sheareb.out',status='old',action='read')
        open(777,file=trim(restart_fname),status='replace')
        read(778,101)j
        write(777,101)j
        read(778,101)mflx
        write(777,101)mflx
        do
           read(778,102,iostat=merror)dum
           if(merror/=0)exit
           write(777,102)dum
        enddo
        close(778)
        close(777)

        open(345,file="FileExit.dat",status="replace")
        write(345,"(A9,i1)")"FileExit=",FileExit
        write(345,"(A9,i5)")"irest   =",irest+1
        if(mod(irest,2)==0)then
           write(345,"(A12)")"restart_dir1"
        else
           write(345,"(A12)")"restart_dir2"
        endif
        close(345)
     endif

     irest=irest+1
  else
     if(mype==0)write(*,*)'read in',restart_fname
     open(901,file=trim(restart_fname),form='unformatted',access='stream',status='old')
     read(901)restart_mflux,restart_mpsi,restart_mzeta,restart_mi,restart_mimax,&
          mgrid,nparam,ntracer,etracer
     if(restart_mflux/=mflux .or. restart_mpsi/=mpsi .or. restart_mzeta/=mzeta .or. &
          restart_mi/=mi .or. restart_mimax/=mimax)then
        write(*,*)'*** restart file dimensions do not match current run'
        write(*,*)'restart:',restart_mflux,restart_mpsi,restart_mzeta,restart_mi,restart_mimax
        write(*,*)'current :',mflux,mpsi,mzeta,mi,mimax
        call MPI_ABORT(MPI_COMM_WORLD,1,merror)
     endif
     read(901)rdtemi,ptracer,pfluxpsi,phi00,phip00,zonali
     allocate(zion0_read(mimax))
     read(901)zion0_read
     read(901)phi
     read(901)zion
     close(901)
     zion0(6,:)=zion0_read
     deallocate(zion0_read)
     irest=irest+1
     if(mype==0)write(*,*)trim(restart_fname),' read over'
  endif

101  format(i6)
102  format(e12.6)
end subroutine restart_io
