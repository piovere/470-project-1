program christmas

! CHRISTMAS: CRIticality Search for neuTron Multi-group reActor Slab

! Neutron Flux desperation Group 3
! DGETRS - solve a system of linear equations A * X = B or 
! See: http://www.netlib.org/lapack/explore-html/d6/d49/dgetrs_8f.html

use iso_fortran_env
implicit none

real(real64), allocatable :: A(:,:,:)               ! matrix to LU decompose
real(real64), allocatable :: b(:,:)                 ! b vector (flux array)
real(real64), allocatable :: ipiv(:,:)                ! pivoting array for LAPACK
real(real64), allocatable :: chi(:)                 ! fission population birthrate (only fast)
integer                   :: n, g,n_flector, n_tot         ! n=number of spatial nodes, g=energy groups
integer                   :: info            
real(real64), allocatable :: sigma_a(:,:), sigma_s(:,:), sigma_tr(:,:),sigma_fnu(:,:), sigma_r(:,:)      
real(real64)              :: w, ratio, dr1, dr2   
real(real64), allocatable :: b_old(:,:), S(:)
!------------------------Used for Iteration loop--------------------------!
integer                   :: j, j_width ,j_tot, i,p                
real(real64)              :: k_error, b_error
real(real64)              :: m, m_old                           
real(real64)              :: k , k_old    
real(real64)              :: min_error=0.001
integer                   :: MAX_ITERATIONS = 500 

! External procedures defined in LAPACK
external DGETRF, DGETRS

!-------------------------------------------------------------------------!

print *,"Number of core nodes? "
read *,n

print *,'Number of reflector nodes?'
read *,n_flector

print *,'What is the desired ratio between core and reflector width?( dw/W) real number please'
read *, ratio
print *,'Number of Energy groups?'
read *,g
n=n-1
n_tot=(n+n_flector)
allocate (A(n_tot,n_tot,g),b(n_tot,g),ipiv(n_tot,g),b_old(n_tot,g),S(n_tot),sigma_fnu(n_tot,g))   ! Used in diff eqs
w=10_real64
!---------------------Matrix Constants Declaration------------------------!

allocate (sigma_a(g,2), sigma_s(g,g), sigma_tr(g,2),chi(g),sigma_r(g,2))

j_tot = 0
j_width=0

sigma_s=0.0
sigma_fnu=0.0
!------------------------G=4 declarations---------------------------------!
if(g .eq. 4)then
    sigma_s(1,1)=0.37045_real64
    sigma_s(1,2)=0.04152_real64
    sigma_s(1,3)=0.00001_real64

    sigma_s(2,2)=0.98285_real64
    sigma_s(2,3)=0.07459_real64                !scattering cross sections(From, to)
    sigma_s(2,4)=0.01371_real64

    sigma_s(3,3)=0.76110_real64
    sigma_s(3,4)=0.31856_real64

    sigma_s(4,3)=0.00085_real64
    sigma_s(4,4)=1.96607_real64

    sigma_a(1,2)=0.00051_real64
    sigma_a(2,2)=0.00354_real64
    sigma_a(3,2)=0.01581_real64
    sigma_a(4,2)=0.04637_real64

                    ! D already
    sigma_tr(1,2)=1_real64/(3.0_real64*0.20608_real64)
    sigma_tr(2,2)=1_real64/(3.0_real64*0.60215_real64)
    sigma_tr(3,2)=1_real64/(3.0_real64*0.56830_real64)
    sigma_tr(4,2)=1_real64/(3.0_real64*1.21110_real64)


             !------------ Core is below-------------!              
    do i=1,n
        sigma_fnu(i,1)=0.009572_real64 
        sigma_fnu(i,2)=0.001193_real64
        sigma_fnu(i,3)=0.01768_real64
        sigma_fnu(i,4)=0.18514_real64
    enddo

    sigma_a(1,1)=0.004946_real64
    sigma_a(2,1)=0.002840_real64
    sigma_a(3,1)=0.03053_real64
    sigma_a(4,1)=0.1210_real64

                    ! D already
    sigma_tr(1,1)=2.1623_real64
    sigma_tr(2,1)=1.0867_real64
    sigma_tr(3,1)=0.6318_real64
    sigma_tr(4,1)=0.3543_real64

    sigma_r(1,1)=0.08795_real64
    sigma_r(2,1)=0.06124_real64
    sigma_r(3,1)=0.09506_real64
    sigma_r(4,1)=0.1210_real64
!------------------------G=2 declarations---------------------------------!
elseif(g .eq. 2)then
    sigma_s(1,2)=0.0494_real64

    sigma_a(1,2)=0.0004_real64
    sigma_a(2,2)=0.0197_real64

    sigma_tr(1,2)=1.13_real64
    sigma_tr(2,2)=0.16_real64


    sigma_r(1,2)=0.0494_real64
    sigma_r(2,2)=0.0197_real64


             !------------ Core is below-------------!              
    do i=1,n    
        sigma_fnu(i,1)=0.008476_real64 
        sigma_fnu(i,2)=0.18514_real64
    enddo

    sigma_a(1,1)=0.01207_real64
    sigma_a(2,1)=0.121_real64


                    ! D already
    sigma_tr(1,1)=1.2627_real64
    sigma_tr(2,1)=0.3543_real64


    sigma_r(1,1)=0.02619_real64
    sigma_r(2,1)=0.121_real64

else
    print *,'I cant do that Dave'
    goto 12
endif


!--------------------------Matrix Declaration-----------------------------!

100 continue
 dr1 = w/(n-1)                               ! Step Size
dr2 = w*ratio/(n_flector-1)
A=0
ipiv=0
do p=1,g
   A(1,1,p) = (0.5_real64*sigma_r(p,1)*dr1)+(sigma_tr(p,1)/(dr1))
   A(1,2,p) = -sigma_tr(p,1)/(dr1)

   A(n_tot,n_tot,p) = sigma_a(1,2)*dr2 + 2.0_real64   *sigma_tr(p,2)/(dr2)
   A(n_tot,n_tot-1,p) = -(1_real64-1_real64/(2_real64*n_tot-1_real64))*sigma_tr(p,2)/(dr2)

   A(n,n-1,p) = -sigma_tr(p,1)/(dr1)  
   A(n,n+1,p) =  -sigma_tr(p,2)/(dr2) 
   A(n,n,p)=(sigma_tr(p,1)/(dr1)+sigma_tr(p,2)/(dr2))+(sigma_a(p,2)*dr2+sigma_a(p,1)*dr1)/2.0_real64
!((2*n+1)/(n-1))*
   do i=2 , n-1
 
      A(i,i,p) = sigma_r(p,1) *dr1 + 2.0_real64   *sigma_tr(p,1)/(dr1)
      A(i,i+1,p) = -(1_real64+1_real64/(2.0_real64*i-1_real64))*sigma_tr(p,1)/(dr1)
      A(i,i-1,p) = -(1_real64-1_real64/(2.0_real64*i-1_real64))*sigma_tr(p,1)/(dr1)

   enddo
   do i=n+1 , n_tot-1

      A(i,i,p) = sigma_a(p,2)*dr2+ 2.0_real64   *sigma_tr(p,2)/(dr2)
      A(i,i+1,p) = -sigma_tr(p,2)/(dr2)*(1_real64+1_real64/(2.0_real64*i-1_real64))
      A(i,i-1,p) = -sigma_tr(p,2)/(dr2)*(1_real64-1_real64/(2.0_real64*i-1_real64))

   enddo
enddo
! Pretty matrix printing
! From http://jblevins.org/log/array-write


    !LU factorization of a general M-by-N matrix A
    ! See: http://www.netlib.org/lapack/explore-html/d3/d6a/dgetrf_8f.html
    do p=1,g    
        call DGETRF(n_tot, n_tot, A(:,:,p), n_tot, ipiv(:,p), info)
        if (info /= 0) stop 'Matrix is numerically singular!'

    enddo
!-------------------------Business loop-----------------------------------!
!-------------------------------------------------------------------------!
!guess initial

b = 1.0_real64
k=1.0_real64
S=0.0_real64
b_old = 1.0_real64
k_old = 1.0_real64
do p=1,g
    do i=1,n_tot
        S(i)= S(i)+dr1*sigma_fnu(i,p)*b(i,p)
    enddo
enddo

j=0

k_error = 1.0_real64
b_error = 1.0_real64
m=0
m_old=0

do while ( ((b_error .gt. min_error) .or. (k_error .gt. min_error)) .and. (j .lt. MAX_ITERATIONS))

    


!------------------Eigenvalue Search--------------------------------------!

    do p=1,g
        if(p .gt. 1)then
            
            do i=1, g
                b(:,p)= b(:,p) + sigma_s(i,p)*b_old(:,i)
            enddo
        endif
        if(p .eq. 1) b(:,1)=S/k

        call DGETRS('N', n_tot, 1, A(:,:,p), n_tot, ipiv(:,p), b(:,p), n_tot, info) 
        if (info /= 0) stop 'Solution of the linear system failed!'

    enddo
    j = j+1        
!-------------------------------------------------------------------------!
!----------------------Error tests----------------------------------------!
    m_old=norm2(S)


    S=0.0_real64
    do p=1,g
        do i=1,n_tot
            S(i)= S(i)+dr1*sigma_fnu(i,p)*b(i,p)
        enddo
    enddo
    m= norm2(S)

    k =k_old* m/m_old                 
    k_error = abs((k-k_old)/k_old)

    do p=1,g
        chi(p)=maxval(abs(b(:,p) - b_old(:,p))/b_old(:,p))
    enddo
    b_error=maxval(chi)

    do i=1,g
        b(:,i)  = b(:,i) / norm2(b(:,i) )
    enddo

    b_old = b     
    k_old = k 

    !print *,b_error
enddo

print *,'Iteration number:  ', j_width
print *,'Mid width:  ', 2.0_real64*w
print *,(b_error .gt. min_error) ,(k_error .gt. min_error)
print *,'final k value:  ', k
print *,''

j_tot=j_tot+j
!-------------------------width adjustments for k=1-----------------------!
!if (j_width .eq. 10)goto 12
 
  if(k .lt. 1.0_real64-min_error)then
     w=1.1_real64*w
    j_width=j_width+1
     goto 100
  else if(k .gt. 1.0_real64+min_error)then  
     w=0.9_real64*w
    j_width=j_width+1
     goto 100
  else
     print *,'Number of width adjustments:',j_width
  endif

12 continue
!-------------------------VOMIT RESULTS-----------------------------------!
print *,'Number of iterations until convergence:  ', j_tot
print *,'Critical width:  ', 2.0_real64*(1_real64+ratio)*w
print *,'K-value is',k
print *,'Final b error is',b_error
print *,'Final k error is',k_error

open( unit = 77, file= "Flux1Group.dat")
open( unit = 66, file= "Flux2Group.dat")

write(77, "(1e10.4)" ) b(:,1)
write(77, "(1e10.4)" ) 0.0
write(66, "(1e10.4)" ) b(:,2)
write(66, "(1e10.4)" ) 0.0


if(g .eq. 2)goto 14

!call system('xdg-open Fluxdis.png')
open( unit = 55, file= "Flux3Group.dat")
open( unit = 44, file= "Flux4Group.dat")

write(55, "(1e10.4)" ) b(:,3)
write(55, "(1e10.4)" ) 0.0
write(44, "(1e10.4)" ) b(:,4)
write(44, "(1e10.4)" ) 0.0
 close(55)
 close(44)
14 continue 
 close(77)
 close(66)

call system('gnuplot -p flux.plt')
end program christmas


