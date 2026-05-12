!========================================
!! originally written by s. either, modified by Y. Xiao April 2008
subroutine tag_particles

!========================================

  use global_parameters
  use particle_array
  use particle_tracking
  implicit none
  integer :: m,np
  real(wp) :: ainside,aoutside,thickratio

  do m=1,mi
     zion(7,m)=-real(m,wp)
     zion(8,m)=real(mype+1,wp)
     zion0(7,m)=-real(m,wp)
     zion0(8,m)=real(mype+1,wp)
  enddo

  thickratio=0.60
  ainside=0.5*(1.0-thickratio)*a1+0.5*(1.0+thickratio)*a0
  ainside=ainside*ainside*0.5

  aoutside=0.5*(1.0+thickratio)*a1+0.5*(1.0-thickratio)*a0
  aoutside=aoutside*aoutside*0.5

  np=0
  if(mype==0)write(*,*)'np=',np,'ainside=',ainside,'aoutside=',aoutside

end subroutine tag_particles

!========================================

subroutine locate_tracked_particles

!========================================

  use global_parameters
  use particle_array
  use particle_tracking
  implicit none
  integer :: npp

  ntrackp=0
  ptrackedi=0.0

  npp=0
  ntrackp(1)=npp

end subroutine locate_tracked_particles

!========================================

subroutine write_tracked_particles

!========================================

  use global_parameters
  use particle_tracking
  implicit none
 
  integer :: j,ntpart(0:numberpe-1)
  character(len=10) :: cdum
#if ADIOS
    character(len=50),SAVE:: fname
    character(len=50)::dirstr
    integer*8 :: group_id, buf_id,comm
    #define ADIOS_WRITE(a,b) call adios_write(a,'b'//char(0),b)
#endif


#if ADIOS
    comm=MPI_COMM_WORLD
    write(fname,'("trackp_dir/TRACKP_",i5.5,".bp")')mstepall+istep
    fname=trim(fname)//char(0)
    write(dirstr,'("node_",i5.5)')mype
    dirstr=trim(dirstr)//char(0)
    call adios_get_group (group_id, "particles"//char(0))
    call adios_set_path (group_id,dirstr//char(0))
    call adios_open (buf_id, group_id, fname)
    ADIOS_WRITE(buf_id,comm)
    ADIOS_WRITE(buf_id,mype)
    ADIOS_WRITE(buf_id,nparam)
    ADIOS_WRITE(buf_id,nspec)
    call adios_close (buf_id)
    if(mype==0)write(*,*)"tracking filename ",fname
#else
  
     write(cdum,'("TRACKP.",i5.5)')mype
     open(57,file=cdum,status='unknown',position='append')
     write(57,*)istep
     write(57,*)ntrackp(1:nspec)

     do j=1,ntrackp(1)
        write(57,*)ptrackedi(1:nparam,j)
     enddo

     close(57)
#endif

end subroutine write_tracked_particles

!========================================
