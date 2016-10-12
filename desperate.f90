program desperate
! Neutron Flux desperation Group 3

use iso_fortran_env
implicit none

real(real64), allocatable :: A(:,:)  ! matrix to LU decompose
real(real64), allocatable :: b(:)    ! b vector
real(real64), allocatable :: ipiv(:) ! pivoting array for LAPACK
integer                   :: n, i, G       ! problem size
integer                   :: info    ! success flag for LAPACK
!character(len=200)        :: fmt     ! formatting string for pretty printing
real(real64)              :: sigma_a, D, dx, w, L, x, psi, C, f 
real(real64), parameter   :: pi=4.0*atan(1.0)
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
print *,"What is you geometry? cartesian=1, cylindrical=2, spherical=3"
read *,G
print *,"Relative Strength of other side of slab (ie f*s) f=?"
read *,f

D=3.62E-2                    ! Sigma Transport
sigma_a = 0.15               ! Sigma 
D = 1.0/(3.0*D)              ! Diffusion Coef
w=100.0                       ! Slab size
dx = w/(n)                     ! Step Size

b(1)=1.0E8                  !Source Strength
b(n)=f*b(1)

L= sqrt(D/sigma_a)

if (G .eq. 1) then              !----------CARTESIAN-------------!

   A(1,1) = (0.5*sigma_a)+(D/(dx**2.0))
   A(1,2) = -D/(dx**2.0)
   A(n,n) = 0.5*sigma_a+D/(dx**2.0)
   A(n,n-1) =-D/(dx**2.0)

   !C =  b(1)*L/(2.0*D*(1.0+exp(-2.0*w/L)))

   !psi = C*(1-exp((-2.0*w)/L))
   !write(55, "(1e10.4)") psi
 
   !psi = C*(exp(-dx/L)-exp((dx-2.0*w)/L))
   !write(55, "(1e10.4)") psi


   do i=2 , n-1

      A(i,i) = sigma_a+2.0*D/(dx**2.0)
      A(i,i+1) = -D/(dx**2.0)
      A(i,i-1) = -D/(dx**2.0)

      !x = dx*i
      !psi = C*(exp(-x/L)-exp((x-2.0*w)/L))
      !write(55, "(1e10.4)") psi
   enddo

elseif (G .eq. 3) then           !------------SPHERE--------------------!

   A(1,1) = (0.5*sigma_a)+(D/(dx**2.0))
   A(1,2) = -(D/(dx**2.0))*(1.0-(2.0/(2.0*1.0-1.0)))
   A(n,n) = 0.5*sigma_a+D/(dx**2.0)
   A(n,n-1) =-(D/(dx**2.0))*(1.0-(2.0/(2.0*n-1.0)))

   C =  b(1)/(4*pi*D)

   psi = 0
   write(55, "(1e10.4)") psi
 
   psi = C*(exp(-dx/L)/(dx))
   write(55, "(1e10.4)") psi


    do i=2 ,n-1

       A(i,i) = sigma_a+2.0*D/(dx**2.0)
       A(i,i+1) = -(D/(dx**2.0))*(1.0+(2.0/(2.0*i+1.0)))
       A(i,i-1) = -(D/(dx**2.0))*(1.0-(2.0/(2.0*i-1.0)))


      x = dx*i
      psi = C*(exp(-x/L)/(x))
      write(55, "(1e10.4)") psi
    enddo


elseif (G .eq. 2) then          !------------------Cylinder---------------------!
    
   A(1,1) = (0.5*sigma_a)+(D/(dx**2.0))
   A(1,2) = -(D/(dx**2.0))*(1.0-(1.0/(2.0*1.0-1.0)))
   A(n,n) = 0.5*sigma_a+D/(dx**2.0)
   A(n,n-1) =-(D/(dx**2.0))*(1.0-(1.0/(2.0*n-1.0)))

!   C =  b(1)*L/(2.0*D*(1.0+exp(-2.0*w/L)))

!   psi = C*(1-exp((-2.0*w)/L))
!   write(55, "(1e10.4)") psi
 
!   psi = C*(exp(-dx/L)-exp((dx-2.0*w)/L))
!   write(55, "(1e10.4)") psi

    do i=2 ,n-1

        A(i,i) = sigma_a+2.0*D/(dx**2.0)
        A(i,i+1) = -(D/(dx**2.0))*(1.0+(1.0/(2.0*i+1.0)))
        A(i,i-1) = -(D/(dx**2.0))*(1.0-(1.0/(2.0*i-1.0)))
!        x = dx*i
!        psi = C*(exp(-x/L)-exp((x-2.0*w)/L))
!        write(55, "(1e10.4)") psi
    enddo

else
   print *,"Thats not a valid coordinate system"
endif

b=b/(2.0*dx)
!print *,"Elements of the b vector ?"
!read *,b

!LU factorization of a general M-by-N matrix A
! See: http://www.netlib.org/lapack/explore-html/d3/d6a/dgetrf_8f.html
call DGETRF(n, n, A, n, ipiv, info)
if (info /= 0) stop 'Matrix is numerically singular!'

   ! guess initial b vector
   ! guess initial k
   ! i = 0
   ! while (k_error > min_error OR phi_error > min_error) and i < MAX_ITERATIONS:
   !     b_prime = A_inverse * CONSTANTS * 1/k * b
   !     k_prime = magnitude(b_prime) / magnitude(b)
   !     k_error = error_function(k, k_prime)
   !     b_error = error_function(b, b_prime)
   !     k = k_prime
   !     b = b_prime
   !     i++
   ! VOMIT RESULTS

! DGETRS - solve a system of linear equations A * X = B or 
! A' * X = B with a general N-by-N matrix A using the LU 
! See: http://www.netlib.org/lapack/explore-html/d6/d49/dgetrs_8f.html
call DGETRS('N', n, 1, A, n, ipiv, b, n, info)  
if (info /= 0) stop 'Solution of the linear system failed!'

! Print solution
open( unit = 77, file= "desperate.dat")
write(77, "(1e10.4)" ) b

end program desperate
