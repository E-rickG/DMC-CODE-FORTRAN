module dmc_module
use nr_module
use nrutil_module
implicit none

logical, parameter :: TRUE = .true., FALSE = .false.
real(8), parameter :: PI = 3.141592653589793d0
real(8), parameter :: hartree = 2.195d5
real(8), parameter :: harttoev = 27.21d0

real(8) :: vref, delta, rhh, dt = 0.1d0
real(8) :: crhh = 1.398d0, rhh2 = 2.0d0
integer(8) :: nact

integer(8) :: Npsips = 500, Nmax = 4000
integer(8) :: Nchec = 2000, Ncol = 2000
integer(8) :: Nwf = 100, seed, Ndim
real(8) :: mrf, ri, stwf

interface

subroutine initpsips(psips, m, n, potnum, Npsips, Nmax, rhh)
integer(8), intent(in) :: m, n, potnum, Npsips, Nmax
real(8), allocatable, intent(inout) :: psips(:,:)
real(8), intent(in) :: rhh
end subroutine initpsips

subroutine count(wz, psips, m, n, potnum, ri, Nwf, stwf)
integer(8), intent(in) :: m, n, potnum, Nwf
real(8), intent(in) :: ri, stwf
real(8), allocatable, intent(inout) :: psips(:,:), wz(:)
end subroutine count

subroutine branch(psips, m, n, potnum, Npsips, seed, dt, Ndim, rhh, vref)
integer(8), intent(in) :: m, n, potnum, Npsips, Ndim
integer(8), intent(inout) :: seed
real(8), intent(in) :: dt, rhh
real(8), intent(inout) :: vref
real(8), allocatable, intent(inout) :: psips(:,:)
end subroutine branch

subroutine walk(psips, m, n, potnum, dt, seed)
integer(8), intent(in) :: m, n, potnum
integer(8), intent(inout) :: seed
real(8), intent(in) :: dt
real(8), allocatable, intent(inout) :: psips(:,:)
end subroutine walk

end interface

end module dmc_module


