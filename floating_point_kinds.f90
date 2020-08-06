module floating_point_kinds
  ! intrinsic :: selected_real_kind
  ! integer,  parameter, public :: dp = selected_real_kind(15)
  use,intrinsic :: iso_fortran_env
  real(real64), parameter :: zero = 0.0d0, one = 1.0d0
end module floating_point_kinds