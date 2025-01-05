module nrutil_module
implicit none

real :: sqrarg, maxarg1, maxarg2, minarg1, minarg2
double precision :: dsqrarg, dmaxarg1, dmaxarg2, dminarg1, dminarg2
real, parameter :: zero = 0.0

contains

function vector(nl, nh) result(v)
integer(8), intent(in) :: nl, nh
real :: v(nl:nh)
end function vector

function ivector(nl, nh) result(v)
integer(8), intent(in) :: nl, nh
integer :: v(nl:nh)
end function ivector

function dvector(nl, nh) result(v)
integer(8), intent(in) :: nl, nh
double precision, allocatable :: v(:)

if (allocated(v)) deallocate(v)
allocate(v(nl:nh))
v = 0.0d0
end function dvector

function matrix(nrl, nrh, ncl, nch) result(m)
integer(8), intent(in) :: nrl, nrh, ncl, nch
real :: m(nrl:nrh, ncl:nch)
end function matrix

function dmatrix(nl, nh, ml, mh) result(v)
integer(8), intent(in) :: nl, nh, ml, mh
real(8), allocatable :: v(:,:)

allocate(v(nl:nh, ml:mh))
end function dmatrix

function f3tensor(nrl, nrh, ncl, nch, ndl, ndh) result(t)
integer(8), intent(in) :: nrl, nrh, ncl, nch, ndl, ndh
real(8) :: t(nrl:nrh, ncl:nch, ndl:ndh)
end function f3tensor

subroutine free_vector(v, nl, nh)
real(8), intent(inout) :: v(:)
integer(8), intent(in) :: nl, nh
end subroutine free_vector

subroutine free_ivector(v, nl, nh)
integer(8), intent(inout) :: v(:)
integer(8), intent(in) :: nl, nh
end subroutine free_ivector

subroutine free_dvector(v, nl, nh)
double precision, intent(inout) :: v(:)
integer(8), intent(in) :: nl, nh
end subroutine free_dvector

subroutine free_matrix(m, nrl, nrh, ncl, nch)
real(8), intent(inout) :: m(:,:)
integer(8), intent(in) :: nrl, nrh, ncl, nch
end subroutine free_matrix

subroutine free_dmatrix(m, nrl, nrh, ncl, nch)
double precision, intent(inout) :: m(:,:)
integer(8), intent(in) :: nrl, nrh, ncl, nch
end subroutine free_dmatrix

subroutine free_f3tensor(t, nrl, nrh, ncl, nch, ndl, ndh)
real(8), intent(inout) :: t(:,:,:)
integer(8), intent(in) :: nrl, nrh, ncl, nch, ndl, ndh
end subroutine free_f3tensor

function SQR(a)
real(8), intent(in) :: a
real(8) :: SQR
sqrarg = a
if (sqrarg == zero) then
	SQR = zero
else
      	SQR = sqrarg * sqrarg
end if
end function SQR

function DSQR(a)
double precision, intent(in) :: a
double precision :: DSQR
dsqrarg = a
if (dsqrarg == zero) then
	DSQR = zero
else
      	DSQR = dsqrarg * dsqrarg
end if
end function DSQR

function DMAX(a, b)
double precision, intent(in) :: a, b
double precision :: DMAX
dmaxarg1 = a
dmaxarg2 = b
if (dmaxarg1 > dmaxarg2) then
      	DMAX = dmaxarg1
else
      	DMAX = dmaxarg2
end if
end function DMAX

function DMIN(a, b)
double precision, intent(in) :: a, b
double precision :: DMIN
dminarg1 = a
dminarg2 = b
if (dminarg1 < dminarg2) then
      	DMIN = dminarg1
else
      	DMIN = dminarg2
end if
end function DMIN

function FMAX(a, b)
real(8), intent(in) :: a, b
real(8) :: FMAX
maxarg1 = a
maxarg2 = b
if (maxarg1 > maxarg2) then
      	FMAX = maxarg1
else
      	FMAX = maxarg2
end if
end function FMAX

function FMIN(a, b)
real(8), intent(in) :: a, b
real(8) :: FMIN
minarg1 = a
minarg2 = b
if (minarg1 < minarg2) then
      	FMIN = minarg1
else
      	FMIN = minarg2
end if
end function FMIN

function SIGN(a, b)
real(8), intent(in) :: a, b
real(8) :: SIGN
if (b >= zero) then
      	SIGN = abs(a)
else
SIGN = -abs(a)
end if
end function SIGN

end module nrutil_module

