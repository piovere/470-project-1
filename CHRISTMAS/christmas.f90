program christmas

! CHRISTMAS: CRIticality Search for neuTron Multi-group reActor Slab

! Neutron Flux desperation Group 3
! DGETRS - solve a system of linear equations A * X = B or 
! See: http://www.netlib.org/lapack/explore-html/d6/d49/dgetrs_8f.html

use iso_fortran_env
implicit none

real(real64), allocatable :: A(:,:,:)               ! matrix to LU decompose
real(real64), allocatable :: b(:,:)                 ! b vector (flux array)
real(real64), allocatable :: ipiv(:)                ! pivoting array for LAPACK
real(real64), allocatable :: chi(:)                 ! fission population birthrate (only fast)
integer                   :: n, g                   ! n=number of spatial nodes, g=energy groups
integer                   :: info            
real(real64), allocatable :: sigma_a(:), sigma_s(:,:), sigma_tr(:),sigma_fnu(:)      
real(real64)              :: D, dx, w   
real(real64), allocatable :: b_old(:,:)
!------------------------Used for Iteration loop--------------------------!
integer                   :: j, j_width ,j_tot, i, h , y ,p ,x               
real(real64)              :: k_error, b_error
real(real64)              :: m, m_old                            
real(real64)              :: k , k_old    
real(real64)              :: min_error=0.0001
integer                   :: MAX_ITERATIONS = 1000 
real(real64)              :: mag
! External procedures defined in LAPACK
external DGETRF, DGETRS
open( unit = 77, file= "Flux.dat")
!-------------------------------------------------------------------------!

print *,"Number of nodes? "
read *,n
print *,'Number of Energy groups?'
read *,g
n=n-1
allocate (A(n,n,g),b(n,g),ipiv(n),b_old(n,g))   ! Used in diff eqs
w=10
!---------------------Matrix Constants Declaration------------------------!

allocate (sigma_a(g), sigma_s(g,g), sigma_tr(g),chi(g),sigma_fnu(g))
sigma_a = 0.1532               
sigma_fnu = 0.157                
sigma_tr=3.62E-2                           ! Sigma Transport
D = 1.0/(3.0*sigma_tr(1))                     ! Diffusion Coef
j_tot = 0
j_width=0

chi=0.0
chi(1)=1.0

100 dx = w/n                               ! Step Size
A=0

!--------------------------Matrix Declaration-----------------------------!
do p=1,g
   A(1,1,p) = (0.5*sigma_a(1))+(D/(dx**2))
   A(1,2,p) = -D/(dx**2.0)
   A(n,n,p) = sigma_a(1)+2.0*D/(dx**2)
   A(n,n-1,p) = -D/(dx**2)

   do i=2 , n-1

      A(i,i,p) = sigma_a(1)+2.0*D/(dx**2)
      A(i,i+1,p) = -D/(dx**2)
      A(i,i-1,p) = -D/(dx**2)

   enddo
enddo
! Pretty matrix printing
! From http://jblevins.org/log/array-write
! do i=1,n
!     write(*,"(100g15.5)") ( A(i,j), j=1,n )
! enddo

    !LU factorization of a general M-by-N matrix A
    ! See: http://www.netlib.org/lapack/explore-html/d3/d6a/dgetrf_8f.html
    do p=1,g    
        call DGETRF(n, n, A(:,:,p), n, ipiv, info)
        if (info /= 0) stop 'Matrix is numerically singular!'
    enddo
!-------------------------Business loop-----------------------------------!
!-------------------------------------------------------------------------!
!guess initial
b = 1.0
k=1.0
mag = 0 
j=0

k_error = 1.0
b_error = 1.0


do while ( ((b_error .gt. min_error) .or. (k_error .gt. min_error)) .and. (j .lt. MAX_ITERATIONS))

    b_old = b     
    k_old = k     
    m_old=0

!------------------Eigenvalue Search--------------------------------------!
    do p=1,g
        do h=1,g
            b(:,p) = (chi(h)*sigma_fnu(h) *b(:,h) + b(:,p) ) /k
        enddo

        call DGETRS('N', n, 1, A(:,:,p), n, ipiv, b(:,p), n, info) 
        if (info /= 0) stop 'Solution of the linear system failed!'

    enddo
    j = j+1        
!-------------------------------------------------------------------------!
!----------------------Error tests----------------------------------------!

    
    do p=1,g
        m=norm2(b(:,p))/norm2(b_old(:,p))
        m_old = m+m_old                             !add individual flux increase?
    enddo
    

    do i=1 ,n
        do p=1,g
            if ( abs(b(i,p)-b_old(i,p)) .gt. mag)then
            mag= abs(b(i,p)-b_old(i,p))
            x=i
            y=p
            endif
        enddo
    enddo

    b_error=(b(x,y)-b_old(x,y))/b_old(x,y)
    b = b / norm2(b)

    print *,'M_OLD=',m_old
    print *,'k=     ', k
    k = k_old * m_old                 
    k_error = abs((k-k_old)/k)

enddo

j_tot=j_tot+j


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
print *,'Critical width:  ', 2.0*w

write(77, "(1e10.4)" ) b
write(77, "(1e10.4)" ) 0.0
call system('gnuplot -p flux.plt')
call system('xdg-open Fluxdis.png')
 close(77)

end program christmas
