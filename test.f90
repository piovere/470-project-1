PROGRAM TEST

use iso_fortran_env
implicit none

real :: t, g, pe
real :: m(3)
integer :: n=3
real, allocatable :: k(:)
integer :: l=3

allocate (k(n))

t = 75
g = 63
pe = 0

call error(t, g, pe)

print *,pe

m(1) = 1.0
m(2) = 2.0
m(3) = 3.0

print*,m

k = m

print*,k

m = 3*m

print*,m

print*,k

END PROGRAM TEST

! User functions and subroutines
subroutine Error(truth, guess, e)
  implicit none
  real :: truth, guess, e

  e = abs(truth-guess)/truth
END subroutine Error
