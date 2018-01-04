module globvar
  intrinsic :: selected_real_kind
  integer,  parameter, public :: dp = selected_real_kind(15)
  real(dp), parameter, public :: zero = 0.0_dp, one = 1.0_dp
end module globvar
