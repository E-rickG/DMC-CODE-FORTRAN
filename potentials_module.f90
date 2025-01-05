module potentials_module

	!***********************************************************************!
	!**  Potentials                                                       **!
	!**                                                                   **!
	!** In order to calculate these we must simply use the potitions of   **!
	!** our replicas and place them into the corresponding potential      **!
	!** function.                                                         **!
	!**                                                                   **!
	!** Example:  Hydrogen atom                                           **!
	!**                                                                   **!
	!** Ndim = 3                                                          **!
	!** psips[i][1] = (indicator of existance)                            **!
	!** psips[i][2] = x coordinate of electron                            **!
	!** psips[i][3] = y coordinate of electron                            **!
	!** psips[i][4] = z coordinate of electron                            **!
	!**                                                                   **!
	!** V = Sqrt(psips[i][2]**2 + psips[i][3]**2 + psips[i][4]**2         **!
	!**                                                                   **!
	!** Other Potentials follow in a similar form                         **!
	!***********************************************************************!

  use nrutil_module
  implicit none
  private
  public :: vho, vmp, vh, vh2p, vh2

contains

	!* Harmonic Oscillator *!
double precision function vho(rint2, nr)
double precision, dimension(:), intent(in) :: rint2
integer(8), intent(in) :: nr                           
double precision :: vtmp

vtmp = 0.5d0 * DSQR(rint2(nr))

vho = vtmp
end function vho

	!* Morse Potential *!
real(8) function vmp(rint2, nr)
implicit none
real(8), intent(in) :: rint2(:)
integer(8), intent(in) :: nr

real(8) :: vtmp

vtmp = 0.5d0 * (exp(-2.0d0 * rint2(nr)) - 2.0d0 * exp(-1.0d0 * rint2(nr)))

vmp = vtmp

end function vmp

	!* H *!
real(8) function vh(rint2, nr, Ndim)
  implicit none
  real(8), intent(in) :: rint2(:)
  integer(8), intent(in) :: nr, Ndim

  real(8) :: vtmp, r1p1
  integer(8) :: i
  real(8) :: dsqrarg

  r1p1 = 0.0d0

  do i = 1, Ndim
    dsqrarg = rint2(i)
    r1p1 = r1p1 + (dsqrarg * dsqrarg)
    	!* r from proton to electron *!
  end do

  r1p1 = sqrt(r1p1)

  vtmp = -1.0d0 / r1p1

  vh = vtmp

end function vh

	!* H2+ *!
double precision function vh2p(rint2, nr, Ndim, rhh)
implicit none
double precision, intent(in) :: rint2(:)
integer(8), intent(in) :: nr, Ndim
real(8), intent(in) :: rhh

double precision :: vtmp, r12, r1p1, r1p2
real(8), allocatable :: rp1(:), rp2(:)
integer :: i

rp1 = dvector(1_8, 3_8)
rp2 = dvector(1_8, 3_8)

        !* Function Body *!
r1p1 = 0.0d0
r1p2 = 0.0d0
rp1(1) = 0.0d0
rp1(2) = 0.0d0
rp1(3) = rhh / 2.0d0
rp2(1) = 0.0d0
rp2(2) = 0.0d0
rp2(3) = -rhh / 2.0d0

do i = 1, Ndim
  	r1p1 = r1p1 + dsqr(rint2(i) - rp1(i))
  	!*  proton 1 *!
  	r1p2 = r1p2 + dsqr(rint2(i) - rp2(i))
  	!*  proton 2 *!
end do

r1p1 = sqrt(r1p1)
r1p2 = sqrt(r1p2)

vtmp = -1.0d0 / r1p1 - 1.0d0 / r1p2 + 1.0d0 / rhh

call free_dvector(rp1, 1_8, 3_8)
call free_dvector(rp2, 1_8, 3_8)

vh2p = vtmp
end function vh2p

	!* H2 *!
double precision function vh2(rint2, nr, rhh)
	!* find the potential that the replicas are in *!
implicit none
double precision, intent(in) :: rint2(:)
integer(8), intent(in) :: nr
real(8), intent(in) :: rhh                

double precision :: vtmp, r12, r1p1, r1p2, r2p1, r2p2
double precision, allocatable :: rp1(:), rp2(:)
integer :: i 

rp1 = dvector(1_8, 3_8)
rp2 = dvector(1_8, 3_8)

	!* Function Body *!
r1p1 = 0.0d0
	!* radius from electron to proton 1 *!
r1p2 = 0.0d0
	!* radius from electron to proton 2 *!
r2p1 = 0.0d0
	!* radius from electron to proton 1 *!
r2p2 = 0.0d0
r12  = 0.0d0
rp1(1) = 0.0d0
rp1(2) = 0.0d0
rp1(3) = rhh / 2.0d0
rp2(1) = 0.0d0
rp2(2) = 0.0d0
rp2(3) = -rhh / 2.0d0

	!* now we step through our array declared in branch();.  It holds 		*!
       	!* the positions of our particles.  the following calc's the V for these	*!
       	!* electrons 									*!

do i = 1, 3
	r1p1 = r1p1 + dsqr(rint2(i) - rp1(i))
	r1p2 = r1p2 + dsqr(rint2(i) - rp2(i))
	r2p1 = r2p1 + dsqr(rint2(3 + i) - rp1(i))
	r2p2 = r2p2 + dsqr(rint2(3 + i) - rp2(i))
	r12  = r12  + dsqr(rint2(3 + i) - rint2(i))
end do

r1p1 = sqrt(r1p1)
	!* proton 1, electron 1 distance *!
r1p2 = sqrt(r1p2)
	!* proton 1, electron 2 distance*!
r2p1 = sqrt(r2p1)
	!* proton 2, electron 1 distance *!
r2p2 = sqrt(r2p2)
	!* proton 2, electron 2 distance *!
r12  = sqrt(r12)
	!* electron 1, electron 2 distance *!

vtmp = -1.0d0 / r1p1 - 1.0d0 / r1p2 - 1.0d0 / r2p1 - 1.0d0 / r2p2 + 1.0d0 / r12 + 1.0d0 / rhh

call free_dvector(rp1, 1_8, 3_8)
call free_dvector(rp2, 1_8, 3_8)

vh2 = vtmp

end function vh2

end module potentials_module

