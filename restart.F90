subroutine restart_io(iop)
  use global_parameters
  use particle_array
  use field_array
  use diagnosis_array
  use particle_decomp
  use data_type
  use particle_tracking
  implicit none
  integer merror, mrequest, mfmode
  integer i,j,k,subsize,startidx,endidx
  integer restart_mflux,restart_mpsi,restart_mzeta,restart_mi,restart_mimax
  integer mquantity,mflx,n_mode,mstepfinal,noutputs
  real(wp) dum
  character(*),intent(in)::iop
  character(len=50) :: restart_fname
  real(wp),dimension(:),allocatable::zion0_read
  logical sheareb_open
  INTEGER(KIND=MPI_OFFSET_KIND) mype_filesize, sum_filesize

#if ADIOS
#define ADIOS_WRITE(a,b) call adios_write(a,'b'//char(0),b)
#define ADIOS_WRITE_PATH(a,b,c) call adios_write_path(a,'b'//char(0),b,c//char(0))
#define ADIOS_READ(a,b) call adios_read(a,'b'//char(0),b)
  integer*8 buf_id, group_id
  character(len=50) :: dirstr
#endif

  if(iop/="read" .and. iop/="write")then
     write(*,*)'*** subroutine restart_io (iop <> "read" or "write")',iop
     call MPI_ABORT(MPI_COMM_WORLD,1,merror)
     return
  endif

 ! if(mype==0)write(*,*)"MFLUX,MPSI,mstepall:",mflux,mpsi,mstepall

!!xy rolling restart
  if(mod(irest,2)==0)then
        write(restart_fname,'(a,i5.5,".bp")')"restart_dir1/restart_",myrank_toroidal 
  else
        write(restart_fname,'(a,i5.5,".bp")')"restart_dir2/restart_",myrank_toroidal 
  endif

!!!  write(restart_fname,'(a,i5.5,".bp")')"restart_dir/restart_",myrank_toroidal !!!,mstepall+istep
!  #ifdef __TIMER
!  call MPI_BARRIER(MPI_COMM_WORLD,merror)
!  start_time=MPI_WTIME()
!  #endif

#if !ADIOS
  if(iop=="write")then
     open(901,file=trim(restart_fname),form='unformatted',access='stream',status='replace')
     write(901)mflux,mpsi,mzeta,mi,mimax,mgrid,nparam,ntracer,etracer
     write(901)rdtemi,ptracer,pfluxpsi,phi00,phip00,zonali
     write(901)zion0(6,:)
     write(901)phi
     write(901)zion
     close(901)

! S.Ethier 01/30/04 Save a copy of history.out and sheareb.out for restart
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

101  format(i6)
102  format(e12.6)
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
  return
#endif

#if ADIOS
     ! setup the element path for this node
     write(dirstr,'("/node",i5.5,"/param")')mype
     dirstr=trim(dirstr)//char(0)
!     call MPI_BARRIER(MPI_COMM_WORLD,merror)
!     start_time = MPI_WTIME()
     call adios_get_group (group_id, "restart"//char(0))
     ! set the path for all vars in the type for proper sizing
     call adios_set_path (group_id,dirstr//char(0));
     restart_fname=trim(restart_fname)//char(0)
     if (iop=="read") then
        call adios_open_read (buf_id, group_id, restart_fname)
     else
        call adios_open (buf_id, group_id, restart_fname)
     endif
     ! write the sizing paramters for both reading and writing
     ADIOS_WRITE(buf_id,partd_comm)
     ADIOS_WRITE(buf_id,mflux)
     ADIOS_WRITE(buf_id,mpsi+1)
     ADIOS_WRITE(buf_id,mzeta+1)
     ADIOS_WRITE(buf_id,nparam)
     ADIOS_WRITE(buf_id,mimax)
     ADIOS_WRITE(buf_id,mgrid)
  if(iop=="write")then
     ADIOS_WRITE(buf_id,mzeta)
     ADIOS_WRITE(buf_id,mi)
     ADIOS_WRITE(buf_id,ntracer)
     ADIOS_WRITE(buf_id,etracer)
     ADIOS_WRITE(buf_id,rdtemi)
     ADIOS_WRITE(buf_id,ptracer)
     ADIOS_WRITE(buf_id,pfluxpsi)
     ADIOS_WRITE(buf_id,phi00)
     ADIOS_WRITE(buf_id,phip00)
     ADIOS_WRITE(buf_id,zonali)
     call adios_write(buf_id,"zion0"//char(0),zion0(6,:))
     ADIOS_WRITE(buf_id,phi)
     ADIOS_WRITE(buf_id,zion)
!     call adios_get_data_size (buf_id, mype_filesize)
     call adios_close (buf_id)
!   #ifdef __TIMER
!     call MPI_BARRIER(MPI_COMM_WORLD,merror)
!     end_time=MPI_WTIME()
!     call MPI_REDUCE(mype_filesize,sum_filesize,1,MPI_INTEGER4,MPI_SUM,0,MPI_COMM_WORLD,merror)
!     if(mype==0)then
!        write(stdout,'("NtoM Time(s), MBytes, MB/s ",a)')iop
!        write(stdout,*)sum_filesize/1024/1024,mype_filesize/1024/1024
!        write(stdout,*)end_time-start_time,mype_filesize*numberpe/1024/1024,mype_filesize*numberpe/((end_time-start_time)*1024*1024)
!     endif
!   #endif

! S.Ethier 01/30/04 Save a copy of history.out and sheareb.out for restart
     if(mype==0 .and. istep<=mstep)then
!!diag XY
        write(*,*)'written ',trim(restart_fname)
   
        !!xy rolling restart
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
        !!diag 
      !!!  print *,'mstepfinal=',mstepfinal
      !!!  print *,'mstep=     ', mstep
      !!!  print *,'istep=     ', istep
      !!!  print *,'noutputs=  ',noutputs

        write(777,101)noutputs
        do i=0,(mquantity+4*n_mode)*noutputs
           read(776,102)dum
           write(777,102)dum
        enddo
        close(776)
        close(777)
        write(*,*)'written ',trim(restart_fname) !!diag XY

        ! Now do sheareb.out

        !!xy rolling restart
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
!        write(*,*)'written ',trim(restart_fname) !!diag XY

        open(345,file="FileExit.dat",status="replace")
        write(345,"(A9,i1)")"FileExit=",FileExit
        write(345,"(A9,i5)")"irest   =",irest+1
        if(mod(irest,2)==0)then
           write(345,"A12")"restart_dir1"
        else
           write(345,"A12")"restart_dir2"
        endif
        close(345)
     endif

     irest=irest+1 !!XY rolling rstart

101  format(i6)
102  format(e12.6)


  else
!!     if(mype==0)write(*,*)"param::",mpsi,mzeta
     if(mype==0)write(*,*)'read in',restart_fname
     allocate(zion0_read(mimax))
     ADIOS_READ(buf_id,mzeta)
     ADIOS_READ(buf_id,mi)
     ADIOS_READ(buf_id,ntracer)
     ADIOS_READ(buf_id,etracer)
     ADIOS_READ(buf_id,rdtemi)
     ADIOS_READ(buf_id,ptracer)
     ADIOS_READ(buf_id,pfluxpsi)
     ADIOS_READ(buf_id,phi00)
     ADIOS_READ(buf_id,phip00)
     ADIOS_READ(buf_id,zonali)
     call adios_read(buf_id,"zion0"//char(0),zion0_read)
     ADIOS_READ(buf_id,phi)
     ADIOS_READ(buf_id,zion)

!     call adios_get_data_size (buf_id, mype_filesize)
     call adios_close (buf_id)
     zion0(6,:)=zion0_read
     deallocate(zion0_read)
     irest=irest+1
!     if(myrank_partd<nproc_partd-1)call MPI_Wait(mrequest,mstatus,merror)
     if(mype==0)write(*,*)restart_fname,'read over'
  endif  ! end of read


#endif
end subroutine restart_io
