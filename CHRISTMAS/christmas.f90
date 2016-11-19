program christmas

! CHRISTMAS: CRIticality Search for neuTron Multi-group reActor Slab

! Neutron Flux desperation Group 3
! DGETRS - solve a system of linear equations A * X = B or 
! See: http://www.netlib.org/lapack/explore-html/d6/d49/dgetrs_8f.html

use iso_fortran_env
implicit none

real(real64), allocatable :: A(:,:,:)               ! matrix to LU decompose
real(real64), allocatable :: b(:,:)                 ! b vector
real(real64), allocatable :: ipiv(:)                ! pivoting array for LAPACK
real(real64), allocatable :: S(:)                   ! Constants array
integer                   :: n, g                   ! n=number of spatial nodes g=energy groups
integer                   :: info                   ! success flag for LAPACK
real(real64), allocatable :: sigma_a(:), sigma_s(:,:), sigma_tr(:)       !specific crosssections for g groups
real(real64)              :: D, dx, w, sigma_fnu  !Declared below

!------------------------Used for Iteration loop--------------------------!
integer :: j, j_width ,j_tot, i, MAX_ITERATIONS = 1000                
real(real64), allocatable :: b_old(:)
real(real64) :: k_error, b_error
real(real64) :: m, m_old                            ! container for all our magnitude functions
real(real64) :: k_old
real(real64) :: min_error=0.0001
real(real64) :: k                                 ! guess initial k

! External procedures defined in LAPACK
external DGETRF, DGETRS
open( unit = 77, file= "Flux.dat")
!-------------------------------------------------------------------------!

print *,"Number of nodes? "
read *,n
print *,'Number of Energy groups?'
read *,g
n=n-1
allocate (A(n,n,g),b(n,g),ipiv(n),b_old(n,g),S(n))   ! Used in diff eqs

!---------------------Matrix Constants Declaration------------------------!

allocate (sigma_a(g), sigma_s(g,g), sigma_tr(g))
sigma_a = 0.1532               
sigma_fnu = 0.157                
sigma_tr=3.62E-2                           ! Sigma Transport
sigma_tr = 1.0/(3.0*sigma_tr)              ! Diffusion Coef
j_tot = 0
j_width=0

100 dx = w/n                               ! Step Size
A=0
b=0

!--------------------------Matrix Declaration-----------------------------!

   A(1,1) = (0.5*sigma_a)+(D/(dx**2))
   A(1,2) = -D/(dx**2.0)
   A(n,n) = sigma_a+2.0*D/(dx**2)
   A(n,n-1) = -D/(dx**2)

   do i=2 , n-1

      A(i,i) = sigma_a+2.0*D/(dx**2)
      A(i,i+1) = -D/(dx**2)
      A(i,i-1) = -D/(dx**2)

   enddo

! Pretty matrix printing
! From http://jblevins.org/log/array-write
! do i=1,n
!     write(*,"(100g15.5)") ( A(i,j), j=1,n )
! enddo

    !LU factorization of a general M-by-N matrix A
    ! See: http://www.netlib.org/lapack/explore-html/d3/d6a/dgetrf_8f.html
    call DGETRF(n, n, A, n, ipiv, info)
    if (info /= 0) stop 'Matrix is numerically singular!'

!-------------------------Business loop-----------------------------------!
!-------------------------------------------------------------------------!
!guess initial
b = 1.0
k=1
mag = 0 
j=0

k_error = 1.0
b_error = 1.0


do while ( ((b_error .gt. min_error) .or. (k_error .gt. min_error)) .and. (j .lt. MAX_ITERATIONS))
 
  b_old = b      ! don't throw away old b
  k_old = k      ! don't throw away old k
  b=b*S/k

  call DGETRS('N', n, 1, A, n, ipiv, b, n, info)  ! b_prime = A_inverse * b
  if (info /= 0) stop 'Solution of the linear system failed!'

  m=norm2(b)
  m_old=norm2(b_old)

  k = k_old* m / m_old                 
  k_error = abs((k-k_old)/k)


  do i=1 ,n
       if ( abs(b(i)-b_old(i)) .gt. mag) you=i
  enddo

  b_error=(b(you)-b_old(you))/b_old(you)
  b = b / m

  j = j+1        ! i++
enddo


!NOTE: I moved the k calculations outside the spatial loop since we only want those to happen 
!during the fast flux calculation (Not quite sure where to put them

!-------------------------width adjustments for k=1-----------------------!
  if(k .lt. 1.0-min_error)then
     w=1.1*w
    j_width=j_width+1
     goto 100
  else if(k .gt. 1.0+min_error)then  
     w=0.9*w
    j_width=j_width+1
     goto 100
  else
     print *,'Number of width adjustments:',j_width
  endif


!-------------------------VOMIT RESULTS-----------------------------------!
print *,'Number of iterations until convergence:  ', j_tot
if (j .eq. MAX_ITERATIONS) print *,'WHAT HAVE YOU DONE'
write(77, "(1e10.4)" ) b
write(77, "(1e10.4)" ) 0.0
print *,'Critical width:  ', 2.0*w
call system('gnuplot -p flux.plt')
call system('xdg-open Fluxdis.png')
 close(77)
end program christmas
