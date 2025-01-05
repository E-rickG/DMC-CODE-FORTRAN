module gasdev_module
use ran3_module, only: ran3
implicit none
private
public :: gasdev

contains

function gasdev(idum) result(gasdev_value)
implicit none
integer(8), intent(inout) :: idum
real(8) :: gasdev_value
real(8), save :: gset
integer, save :: iset = 0
real(8) :: fac, rsq, v1, v2

if (iset == 0) then
	do
        	v1 = 2.0d0 * ran3(idum) - 1.0d0
        	v2 = 2.0d0 * ran3(idum) - 1.0d0
        	rsq = v1 * v1 + v2 * v2
        	if (rsq < 1.0d0 .and. rsq /= 0.0d0) exit
      	end do
	fac = sqrt(-2.0d0 * log(rsq) / rsq)
      	gset = v1 * fac
      	iset = 1
      	gasdev_value = v2 * fac
else
      	iset = 0
      	gasdev_value = gset
end if

end function gasdev

end module gasdev_module

