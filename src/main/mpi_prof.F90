!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module mpiprof
!
! Temporary per-rank profiling of the MPI cell exchange. Not for upstream.
!
! Every counter is per rank, summed over the whole run and over all OpenMP
! threads, so times are thread-seconds. Each rank writes mpiprof_NNNN.txt
! at the end of the run.
!
! :References: None
!
! :Owner: None
!
! :Runtime parameters: None
!
! :Dependencies: io, omputils
!
 implicit none

 integer, parameter, public :: iprof_dens = 1, iprof_force = 2

 ! counts
 integer(kind=8), public :: ncell_local(2)   = 0  ! active cells computed in the local loop
 integer(kind=8), public :: npart_local(2)   = 0  ! particles in those cells
 integer(kind=8), public :: ncell_export(2)  = 0  ! cells exported
 integer(kind=8), public :: nmsg_export(2)   = 0  ! messages sent for exports (cells x targets)
 integer(kind=8), public :: ncell_remote(2)  = 0  ! cells computed on behalf of other ranks
 integer(kind=8), public :: ncalls(2)        = 0  ! calls to densityiterate / force
 integer(kind=8), public :: nremote_its      = 0  ! density remote iterations
 integer(kind=8), public :: ntarget_hist(2,0:64) = 0  ! histogram of targets per exported cell
 integer(kind=8), public :: nopen_geom(2) = 0   ! global leaves opened by r < xsize_i+xsize_j+rcut_i
 integer(kind=8), public :: nopen_hj(2)   = 0   ! ... only because the node hmax enlarged rcut
 integer(kind=8), public :: nopen_grav(2) = 0   ! ... only by the gravity opening criterion
 integer,         public :: refine_min = 1000, refine_max = -1
 integer(kind=8), public :: npart_sum = 0, nbuild = 0
 ! per-call record for force: local active particles and local compute time
 integer, parameter, public :: maxcallrec = 20000
 integer(kind=8), public :: callrec_np(maxcallrec) = 0
 real(kind=8),    public :: callrec_t(maxcallrec)  = 0.

 ! thread-seconds
 real(kind=8), public :: t_compute(2)  = 0.  ! compute_cell on local cells
 real(kind=8), public :: t_sendwait(2) = 0.  ! spinning until previous sends complete
 real(kind=8), public :: t_rww(2)      = 0.  ! inside recv_while_wait (waiting for everyone)
 real(kind=8), public :: t_remote(2)   = 0.  ! computing cells for other ranks
 real(kind=8), public :: t_finish(2)   = 0.  ! finish_cell_exchange
 real(kind=8), public :: t_tglobal = 0., t_tlocal = 0., t_trefine = 0., t_tbal = 0.  ! tree build parts

 public :: write_mpiprof, wtime

 private

contains

real(kind=8) function wtime()
!$ use omp_lib, only:omp_get_wtime
 integer(kind=8) :: c, r
 wtime = 0.
!$ wtime = omp_get_wtime()
!$ return
 call system_clock(c,r)
 wtime = real(c,8)/real(r,8)
end function wtime

subroutine write_mpiprof()
 use io, only:id,nprocs
 integer :: lu, k, i
 character(len=32) :: fname
 character(len=5), parameter :: lab(2) = (/'dens ','force'/)

 write(fname,"(a,i4.4,a)") 'mpiprof_',id,'.txt'
 open(newunit=lu,file=fname,status='replace')
 write(lu,"(a,i6,a,i6)") '# rank ',id,' of ',nprocs
 do k=1,2
    write(lu,"(a,1x,a,i10)")       lab(k),'ncalls      ',ncalls(k)
    write(lu,"(a,1x,a,i14)")       lab(k),'ncell_local ',ncell_local(k)
    write(lu,"(a,1x,a,i14)")       lab(k),'npart_local ',npart_local(k)
    write(lu,"(a,1x,a,i14)")       lab(k),'ncell_export',ncell_export(k)
    write(lu,"(a,1x,a,i14)")       lab(k),'nmsg_export ',nmsg_export(k)
    write(lu,"(a,1x,a,i14)")       lab(k),'ncell_remote',ncell_remote(k)
    write(lu,"(a,1x,a,f14.3)")     lab(k),'t_compute   ',t_compute(k)
    write(lu,"(a,1x,a,f14.3)")     lab(k),'t_sendwait  ',t_sendwait(k)
    write(lu,"(a,1x,a,f14.3)")     lab(k),'t_rww       ',t_rww(k)
    write(lu,"(a,1x,a,f14.3)")     lab(k),'t_remote    ',t_remote(k)
    write(lu,"(a,1x,a,f14.3)")     lab(k),'t_finish    ',t_finish(k)
    write(lu,"(a,1x,a)",advance='no') lab(k),'targets     '
    do i=0,min(nprocs,64)
       write(lu,"(1x,i10)",advance='no') ntarget_hist(k,i)
    enddo
    write(lu,*)
 enddo
 write(lu,"(a,i14)") 'dens  nremote_its ',nremote_its
 do k=1,2
    write(lu,"(a,1x,a,3i14)") lab(k),'open g/hj/gr',nopen_geom(k),nopen_hj(k),nopen_grav(k)
 enddo
 write(lu,"(a,2i6)") 'tree  refinelevels min/max',refine_min,refine_max
 write(lu,"(a,4f10.3)") 'tree  t_global(excl bal) t_bal t_local t_refine',t_tglobal-t_tbal,t_tbal,t_tlocal,t_trefine
 write(lu,"(a,f14.1)") 'tree  mean npart per build',real(npart_sum)/max(nbuild,1_8)
 write(lu,"(a)") '# percall  icall  npart_active_local  t_compute_local'
 do i=1,int(min(ncalls(2),int(maxcallrec,8)))
    write(lu,"(a,i8,i12,es14.5)") 'percall',i,callrec_np(i),callrec_t(i)
 enddo
 close(lu)

end subroutine write_mpiprof

end module mpiprof
