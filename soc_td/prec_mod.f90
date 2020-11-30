!===================================================================================================================================
! Date (last update): 15/07/2014                                                                Status:   Finished + Tested        |
!-----------------------------------------------------------------------------------------------------------------------------------
! This module explicitly sets up symbolic kind types for most used variables, i.e. integers, reals, complex, logical and           |
! user-defined...                                                                                                                  |
!--------------------------------------------------------- C 132 -------------------------------------------------------------------
module prec_mod
  implicit none
  public
! SYMBOLIC NAMES FOR KIND TYPES OF INTEGERS:
!    integer, parameter :: i1b = selected_int_kind(2)                      ! 1-byte
!    integer, parameter :: i2b = selected_int_kind(4)                      ! 2-bytes
!    integer, parameter :: i4b = selected_int_kind(9)                      ! 4-bytes
! SYMBOLIC NAMES FOR KIND TYPES OF REALS:
    integer, parameter :: spr = kind(1.0)                                 ! single precision
    integer, parameter :: dpr = kind(1.d0)                                ! double precision
    integer, parameter :: qpr = selected_real_kind(2*precision(1.0_dpr))  ! quadruple precision
! SYMBOLIC NAMES FOR KIND TYPES OF COMPLEX:
    integer, parameter :: spc = kind((1.0_spr, 1.0_spr))                  ! single precision complex
    integer, parameter :: dpc = kind((1.0_dpr, 1.0_dpr))                  ! double precision complex
    integer, parameter :: qpc = kind((1.0_qpr, 1.0_qpr))                  ! quadruple precision complex
! SYMBOLIC NAME FOR KIND TYPE OF DEFAULT LOGICAL:
!     integer, parameter :: lgc = kind(.true.) 
! USER-DEFINED TYPE OF VARIABLES:
    type spin_orb
      integer          :: numb
      character(len=1) :: spin
    end type spin_orb

    real(dpr), parameter  :: fine_stru = 7.297352568D-3
    !au2debye 1a.u.=h_bar**2/(m_e*e)     
    real(dpr), parameter  :: au2debye = 2.54174780119995
    !hartree2wavenumbers = 219474.6
    real(dpr), parameter  :: au2wavnum = 219474.6   
    real(dpr), parameter  :: au2ev = 27.21138505
    real(dpr), parameter  :: pi = 3.14159265359
    real(dpr), parameter  :: r_pi = 0.28209479177 !0.5 * 1.0/sqrt(pi)

end module prec_mod
