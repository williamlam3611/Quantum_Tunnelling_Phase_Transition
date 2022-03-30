module constants
use iso_fortran_env, only: real32, real64
implicit none

! Precision
    integer, parameter :: sp = real32
    integer, parameter :: dp = real64

! Mathematical
    real(dp), parameter :: pi = 4.0_dp * datan(1.0_dp)

! Physical
    real(dp), parameter :: e = 1.60217662e-19_dp
    real(dp), parameter :: eps = 8.85418782e-12_dp
    real(dp), parameter :: c = 	299792458.0_dp

end module constants
