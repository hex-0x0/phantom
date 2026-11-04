!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module moddump
!
! Convert a dump from the barotropic eos (ieos=8, no thermal energy) to an ideal gas (ieos=2):
! u = (P/rho)/(gamma-1) with P/rho taken from the barotropic law, so the temperature is unchanged.
! Use with an ISOTHERMAL=no build.
!
! :References: None
!
! :Owner: Daniel Price
!
! :Runtime parameters: None
!
! :Dependencies: eos, eos_barotropic, part
!
 implicit none
 character(len=*), parameter, public :: moddump_flags = ''

contains

subroutine modify_dump(npart,npartoftype,massoftype,xyzh,vxyzu)
 use eos,           only:ieos,gamma,polyk,polyk2,init_eos
 use eos_barotropic, only:get_eos_barotropic
 use part,          only:rhoh,igas,isdead_or_accreted
 integer, intent(inout) :: npart
 integer, intent(inout) :: npartoftype(:)
 real,    intent(inout) :: massoftype(:)
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:)
 integer :: i,ierr
 real    :: rhoi,ponrhoi,spsoundi,gammai,tmax

 ieos = 8
 call init_eos(ieos,ierr)
 tmax = 0.
 do i = 1,npart
    if (isdead_or_accreted(xyzh(4,i))) cycle
    rhoi = rhoh(xyzh(4,i),massoftype(igas))
    call get_eos_barotropic(rhoi,polyk,polyk2,ponrhoi,spsoundi,gammai)
    vxyzu(4,i) = ponrhoi/(5./3.-1.)
    tmax = max(tmax,ponrhoi)
 enddo
 ieos  = 2
 gamma = 5./3.
 print*,' set u = (P/rho)/(gamma-1) from the barotropic eos, max P/rho = ',tmax,'; output dump has ieos = 2, gamma = 5/3'

end subroutine modify_dump

end module moddump
