module hqt_constants
use iso_fortran_env, only: real32, real64
implicit none

! Precision
    integer, parameter :: hqt_sp = real32
    integer, parameter :: hqt_dp = real64

! Mathematical
    real(hqt_dp), parameter :: hqt_pi = 4.0_hqt_dp * datan(1.0_hqt_dp)

! Physical
    real(hqt_dp), parameter :: hqt_e = 1.60217662e-19_hqt_dp
    real(hqt_dp), parameter :: hqt_eps = 8.85418782e-12_hqt_dp
    real(hqt_dp), parameter :: hqt_c = 	299792458.0_hqt_dp

end module hqt_constants
