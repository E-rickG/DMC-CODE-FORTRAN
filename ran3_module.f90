module ran3_module
  implicit none
  integer, parameter :: MBIG = 1000000000
  integer, parameter :: MSEED = 161803398
  integer, parameter :: MZ = 0
  real(8), parameter :: FAC = 1.0d0 / MBIG
contains

  function ran3(idum) result(rand_num)
    implicit none
    integer(8), intent(inout) :: idum
    real(8) :: rand_num
    integer :: i, ii, k
    integer, save :: ma(56) = 0, inext = 0, inextp = 0, iff = 0
    integer(8) :: mj, mk

    ! Inicialização
    if (idum < 0 .or. iff == 0) then
      iff = 1
      mj = MSEED - abs(idum)
      mj = mod(mj, MBIG)
      ma(55) = mj
      mk = 1
      do i = 1, 54
        ii = mod(21 * i, 55)
        ma(ii) = mk
        mk = mj - mk
        if (mk < MZ) mk = mk + MBIG
        mj = ma(ii)
      end do

      do k = 1, 4
        do i = 1, 55
          ma(i) = ma(i) - ma(1 + mod(i + 30, 55))
          if (ma(i) < MZ) ma(i) = ma(i) + MBIG
        end do
      end do
      inext = 0
      inextp = 31
      idum = 1
    end if

    ! Geração do número aleatório
    inext = inext + 1
    if (inext == 56) inext = 1
    inextp = inextp + 1
    if (inextp == 56) inextp = 1
    mj = ma(inext) - ma(inextp)
    if (mj < MZ) mj = mj + MBIG
    ma(inext) = mj
    rand_num = mj * FAC

  end function ran3

end module ran3_module

