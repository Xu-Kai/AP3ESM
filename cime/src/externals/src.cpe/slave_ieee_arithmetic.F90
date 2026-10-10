module IEEE_ARITHMETIC
implicit none

type, private :: IEEE_CLASS_TYPE 
    private
    integer :: hidden
end type

interface IEEE_VALUE
    module procedure IEEE_VALUE_8
end interface
public :: IEEE_VALUE

contains

elemental real(kind=8) function IEEE_VALUE_8(X, CLASS) result(res)
    !$om2ath entry(__ieee_arithmetic_MOD_ieee_value_8)
    real(kind=8), intent(in) :: X
    type(IEEE_CLASS_TYPE), intent(in) :: CLASS

    select case (CLASS%hidden)
        case (1)
            res = -1
            res = sqrt(res)
        case (2)
            res = -1
            res = sqrt(res)
        case (3)
            res = huge(res)
            res = (-res) * res
        case (4)
            res = -43
        case (5)
            res = -tiny(res)
            res = res / 2
        case (6)
            res = 0
            res = -res
        case (7)
            res = 0
        case (8)
            res = tiny(res)
            res = res / 2
        case (9)
            res = 42
        case (10)
            res = huge(res)
            res = res * res
        case default
            res = 0
    end select
end function

!elemental real(kind=8) function IEEE_VALUE_8(X, CLASS) result(res)
!
!  real(kind=8), intent(in) :: X
!  type(IEEE_CLASS_TYPE), intent(in) :: CLASS
!
!  select case (CLASS%hidden)
!    case (1)     ! IEEE_SIGNALING_NAN
!      res = -1
!      res = sqrt(res)
!    case (2)     ! IEEE_QUIET_NAN
!      res = -1
!      res = sqrt(res)
!    case (3)     ! IEEE_NEGATIVE_INF
!      res = huge(res)
!      res = (-res) * res
!    case (4)     ! IEEE_NEGATIVE_NORMAL
!      res = -42
!    case (5)     ! IEEE_NEGATIVE_DENORMAL
!      res = -tiny(res)
!      res = res / 2
!    case (6)     ! IEEE_NEGATIVE_ZERO
!      res = 0
!      res = -res
!    case (7)     ! IEEE_POSITIVE_ZERO
!      res = 0
!    case (8)     ! IEEE_POSITIVE_DENORMAL
!      res = tiny(res)
!      res = res / 2
!    case (9)     ! IEEE_POSITIVE_NORMAL
!      res = 42
!    case (10)    ! IEEE_POSITIVE_INF
!      res = huge(res)
!      res = res * res
!    case default ! IEEE_OTHER_VALUE, should not happen
!      res = 0
!   end select
!end function

end module IEEE_ARITHMETIC
