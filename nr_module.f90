module nr_module
implicit none

type :: fcomplex
real(8) :: r, i
end type fcomplex

type :: arithcode
integer(8), pointer :: ilob(:), iupb(:), ncumfq(:)
integer(8) :: jdif, nc, minint, nch, ncum, nrad
end type arithcode

type :: huffcode
integer(8), pointer :: icode(:), ncod(:), left(:), right(:)
integer(8) :: nch, nodemax
end type huffcode

contains

subroutine asolve(n, b, x, itrnsp)
integer(8), intent(in) :: n
real(8), intent(in) :: b(n)
real(8), intent(out) :: x(n)
integer(8), intent(in) :: itrnsp

integer :: i

if (itrnsp == 0) then
	x = b
else
	x = b
end if
end subroutine asolve

subroutine atimes(n, x, r, itrnsp)
integer(8), intent(in) :: n
real(8), intent(in) :: x(n)
real(8), intent(out) :: r(n)
integer(8), intent(in) :: itrnsp

if (itrnsp == 0) then
	r = x
else
	r = x
end if
end subroutine atimes

function bessi(n, x) result(res)
integer, intent(in) :: n
real(8), intent(in) :: x
real(8) :: res
res = 0.0 
end function bessi

end module nr_module

