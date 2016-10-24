program desperatev2
! Neutron Flux desperation Group 3
! DGETRS - solve a system of linear equations A * X = B or 
! See: http://www.netlib.org/lapack/explore-html/d6/d49/dgetrs_8f.html

use iso_fortran_env
implicit none

real(real64), allocatable :: A(:,:)                 ! matrix to LU decompose
real(real64), allocatable :: b(:)                   ! b vector
real(real64), allocatable :: ipiv(:)                ! pivoting array for LAPACK
real(real64), allocatable :: S(:)                   ! Constants array
integer                   :: n, i               ! problem size
integer                   :: info                   ! success flag for LAPACK
real(real64)              :: sigma_a, D, dx, w, nusig  !Declared below
!real(real64), parameter   :: pi=4.0*atan(1.0)

!------------------------Used for Iteration loop--------------------------!
integer :: j, MAX_ITERATIONS = 1000                
real(real64), allocatable :: b_old(:)
real(real64) :: k_error, b_error
real(real64) :: m, m_old                            ! container for all our magnitude functions
real(real64) :: k_old
real(real64) :: min_error=0.001
real(real64) :: k=1                                 ! guess initial k

real(real64) :: mag 
! External procedures defined in LAPACK
external DGETRF, DGETRS
!open( unit = 55, file= "Analytic.dat")
open( unit = 77, file= "Numeric.dat")
!-------------------------------------------------------------------------!

print *,"Number of nodes? "
read *,n
n=n-1
allocate (A(n,n),b(n),ipiv(n),b_old(n),S(n))



sigma_a = 0.15               ! Sigma absorb
nusig = 0.157                ! Sigma fission*nu
D=3.62E-2                    ! Sigma Transport
D = 1.0/(3.0*D)              ! Diffusion Coef

w=10                      ! Slab size
dx = w/n                     ! Step Size
A=0
b=0

S=nusig
S(1)=S(1)/2
!--------------------------Matrix Declaration-----------------------------!

   A(1,1) = (0.5*sigma_a)+(D/(dx**2.0))
   A(1,2) = -D/(dx**2.0)
   A(n,n) = sigma_a+2.0*D/(dx**2.0)
   A(n,n-1) = -D/(dx**2.0)

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
b = 10.0

do i=1, n
    b(i)=b(i)*S(i)
enddo


j=0

k_error = 1.0
b_error = 1.0

! while (k_error > min_error OR phi_error > min_error) and i < MAX_ITERATIONS:
do while ((k_error .gt. min_error) .or. (b_error .gt. min_error)) 

  b=S*b/k

  b_old = b      ! don't throw away old b
  k_old = k      ! don't throw away old k

  call DGETRS('N', n, 1, A, n, ipiv, b, n, info)  ! b_prime = A_inverse * k * b
  if (info /= 0) stop 'Solution of the linear system failed!'

!  call magnitude(b, n, m)
  mag=0.0
  do i=1,n
    mag = mag+(b(i)**2)
  enddo
  m=sqrt(mag)

! call magnitude(b_old, n, m_old)

    mag=0.0
  do i=1,n
    mag = mag+(b_old(i)**2)
  enddo
  m_old=sqrt(mag)

  k =  m / m_old                  ! k_prime = magnitude(b) / magnitude(b_old)
      
!  call error(k_old, k, k_error)  ! k_error = error_function(k, k_prime)
  k_error = (k-k_old)/k

!  call error(b_old, b, b_error)  ! b_error = error_function(b, b_prime)
  mag=0
  do i=1,n 
    mag=mag+(b(i)-b_old(i))
  enddo
  b_error=mag/n

  j = j+1        ! i++
enddo



! VOMIT RESULTS
print *,'Number of iterations until convergence:  ', j
if (j .eq. MAX_ITERATIONS) print *,'WHAT HAVE YOU DONE'
write(77, "(1e10.4)" ) b
write(77, "(1e10.4)" ) 0.0

end program desperatev2



subroutine magnitude(vector, n, mag)
  implicit none
  real :: mag ! contains magnitude variable
  integer :: n, i
  real, dimension (1:n) :: vector ! the vector
  mag=0.0
  do i=1,n
    mag = mag+(vector(i)*vector(i))
  enddo

  mag=sqrt(mag)
end subroutine magnitude

subroutine Error(truth, guess, e)
  implicit none
  real :: truth, guess, e

  e = (truth-guess)/truth
END subroutine Error
