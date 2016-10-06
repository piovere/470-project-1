program desperate
! Neutron Flux desperation

use iso_fortran_env
implicit none

real(real64), allocatable :: A(:,:)  ! matrix to LU decompose
real(real64), allocatable :: b(:)    ! b vector
real(real64), allocatable :: ipiv(:) ! pivoting array for LAPACK
integer                   :: n, i       ! problem size
integer                   :: info    ! success flag for LAPACK
!character(len=200)        :: fmt     ! formatting string for pretty printing
real(real64)              :: sigma_a, D, dx, w, L, x, psi, C, f 
! External procedures defined in LAPACK
external DGETRF, DGETRS
open( unit =55, file="realsol.dat")

! Input data
print *,"Size of the problem n? "
read *,n
n=n
allocate (A(n,n),b(n),ipiv(n))
A=0
b=0
print *,"Relative Strength of other side of slab (ie f*s) f=?"
read *,f
D=3.62E-2                    ! Sigma Transport
sigma_a = 0.15               ! Sigma 
D = 1.0/(3.0*D)              ! Diffusion Coef
w=100.0                       ! Slab size
dx = w/(n)                     ! Step Size


L= sqrt(D/sigma_a)

A(1,1) = (0.5*sigma_a)+(D/(dx**2.0))
A(1,2) = -D/(dx**2.0)
A(n,n) = 0.5*sigma_a+D/(dx**2.0)
A(n,n-1) =-D/(dx**2.0)

b(1)=1.0E8                   ! Source Strength
b(n)=f*b(1)

   C =  b(1)*L/(2.0*D*(1.0+exp(-2.0*w/L)))
   print *,"C is", C

   psi = C*(1-exp((-2.0*w)/L))
   write(55, "(1e10.4)") psi
!   print *, "Psi(1) is ",psi
 
   psi = C*(exp(-dx/L)-exp((dx-2.0*w)/L))
   write(55, "(1e10.4)") psi


do i=2 , n-1
   x = dx*i
   A(i,i) = sigma_a+2.0*D/(dx**2.0)
   A(i,i+1) = -D/(dx**2.0)
   A(i,i-1) = -D/(dx**2.0)

   psi = C*(exp(-x/L)-exp((x-2.0*w)/L))
   write(55, "(1e10.4)") psi
enddo

b=b/(2.0*dx)
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
!print *, "Solved vector x is:"
!write (fmt,"(a,i3,a)") "(",n,"(f10.4,1x))"
!print fmt, b
open( unit = 77, file= "desperate.dat")
write(77, "(1e10.4)" ) b
print *,b(1)

end program desperate

