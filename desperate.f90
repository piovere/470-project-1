program desperate
! Neutron Flux desperation

use iso_fortran_env
implicit none

real(real64), allocatable :: A(:,:)  ! matrix to LU decompose
real(real64), allocatable :: b(:)    ! b vector
real(real64), allocatable :: ipiv(:) ! pivoting array for LAPACK
integer                   :: n, i       ! problem size
integer                   :: info    ! success flag for LAPACK
character(len=200)        :: fmt     ! formatting string for pretty printing
real(real64)              :: sigma_a, D, dx, w
! External procedures defined in LAPACK
external DGETRF, DGETRS

! Input data
print *,"Size of the problem n? "
read *,n
n=n-1
allocate (A(n,n),b(n),ipiv(n))
A=0
b=0

sigma_a = 0.0924
D = 1000.0 
w=300
dx = w/n
print *,dx
A(1,1) = (0.5*sigma_a)+(D/(dx**2.0))
A(1,2) = -D/(dx**2.0)
A(n,n) = sigma_a+2.0*D/(dx**2.0)
A(n,n-1) =-D/(dx**2.0)

b(1)=1.0E8

do i=2 , n-1
   A(i,i) = sigma_a+2.0*D/(dx**2.0)
   A(i,i+1) = -D/(dx**2.0)
   A(i,i-1) = -D/(dx**2.0)
enddo

!print *,"Elements of the b vector ?"
!read *,b

!LU factorization of a general M-by-N matrix A
! See: http://www.netlib.org/lapack/explore-html/d3/d6a/dgetrf_8f.html
call DGETRF(n, n, A, n, ipiv, info)
if (info /= 0) stop 'Matrix is numerically singular!'

! DGETRS - solve a system of linear equations A * X = B or 
! A' * X = B with a general N-by-N matrix A using the LU 
! See: http://www.netlib.org/lapack/explore-html/d6/d49/dgetrs_8f.html
call DGETRS('N', n, 1, A, n, ipiv, b, n, info)  
if (info /= 0) stop 'Solution of the linear system failed!'

! Print solution
print *, "Solved vector x is:"
write (fmt,"(a,i3,a)") "(",n,"(f10.4,1x))"
print fmt, b
open( unit = 77, file= "desperate.dat")
write(77, "(e10.4)" ) b

end program desperate

