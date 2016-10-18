program desperatev2
! Neutron Flux desperation Group 3
! DGETRS - solve a system of linear equations A * X = B or 
! See: http://www.netlib.org/lapack/explore-html/d6/d49/dgetrs_8f.html

use iso_fortran_env
implicit none

real(real64), allocatable :: A(:,:)                 ! matrix to LU decompose
real(real64), allocatable :: b(:)                   ! b vector
real(real64), allocatable :: ipiv(:)                ! pivoting array for LAPACK
integer                   :: n, i, G                ! problem size
integer                   :: info                   ! success flag for LAPACK
real(real64)              :: sigma_a, D, dx, w, L   !Declared below
!real(real64), parameter   :: pi=4.0*atan(1.0)

!------------------------Used for Iteration loop--------------------------!
integer :: j, MAX_ITERATIONS = 1000                  
real(real64), allocatable :: b_old(:)
real(real64) :: k_error, b_error
real(real64) :: m, m_old                            ! container for all our magnitude functions
real(real64) :: k_old
real(real64) :: min_error=0.01
real(real64) :: k=0                                 ! guess initial k

! External procedures defined in LAPACK
external DGETRF, DGETRS
open( unit = 55, file= "Analytic.dat")
open( unit = 77, file= "Numeric.dat")
!-------------------------------------------------------------------------!

print *,"Number of nodes? "
read *,n
allocate (A(n,n),b(n),ipiv(n),b_old(n))

print *,"What is you geometry? cartesian=1, cylindrical=2, spherical=3"
read *,G

D=3.62E-2                    ! Sigma Transport
sigma_a = 0.15               ! Sigma 
D = 1.0/(3.0*D)              ! Diffusion Coef
w=100.0                      ! Slab size
dx = w/n                     ! Step Size
A=0
b=0

L= sqrt(D/sigma_a)

!--------------------------Matrix Declaration-----------------------------!

   A(1,1) = (0.5*sigma_a)+(D/(dx**2.0))
   A(1,2) = -D/(dx**2.0)
   A(n,n) = 0.5*sigma_a+D/(dx**2.0)
   A(n,n-1) =-D/(dx**2.0)

   do i=2 , n-1

      A(i,i) = sigma_a+2.0*D/(dx**2.0)
      A(i,i+1) = -D/(dx**2.0)
      A(i,i-1) = -D/(dx**2.0)

   enddo

!LU factorization of a general M-by-N matrix A
! See: http://www.netlib.org/lapack/explore-html/d3/d6a/dgetrf_8f.html
call DGETRF(n, n, A, n, ipiv, info)
if (info /= 0) stop 'Matrix is numerically singular!'



!-------------------------Business loop-----------------------------------!
!-------------------------------------------------------------------------!
! guess initial b vector
b = 1.0

k_old = k
b_old = b

j=0

k_error = 1.0
b_error = 1.0

! while (k_error > min_error OR phi_error > min_error) and i < MAX_ITERATIONS:
do while (((k_error .gt. min_error) .or. (b_error .gt. min_error)) .AND. (j .lt. MAX_ITERATIONS))

  b_old = b      ! don't throw away old b
  k_old = k      ! don't throw away old k

  call DGETRS('N', n, 1, A, n, ipiv, k*b, n, info)  ! b_prime = A_inverse * k * b
  if (info /= 0) stop 'Solution of the linear system failed!'

  call magnitude(b, n, m)
  call magnitude(b_old, n, m_old)
  k = m / m_old                  ! k_prime = magnitude(b) / magnitude(b_old)

  call error(k_old, k, k_error)  ! k_error = error_function(k, k_prime)
  call error(b_old, b, b_error)  ! b_error = error_function(b, b_prime)

  j = j+1        ! i++
enddo

! VOMIT RESULTS
print *,'Number of iterations until convergence: Max=1000 ', j
if (j .eq. MAX_ITERATIONS) print *,'WHAT HAVE YOU DONE'
write(77, "(1e10.4)" ) b

end program desperatev2



subroutine magnitude(vector, n, m)
  implicit none
  real :: m ! contains magnitude variable
  integer :: n, i
  real :: vector(n) ! the vector
  m=0.0
  do i=1,n
    m = m+vector(i)*vector(i)
  enddo

  m=sqrt(m)
end subroutine magnitude

subroutine Error(truth, guess, e)
  implicit none
  real :: truth, guess, e

  e = abs(truth-guess)/truth
END subroutine Error
