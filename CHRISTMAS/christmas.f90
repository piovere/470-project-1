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
integer                   :: n, g,n_flector, n_tot         ! n=number of spatial nodes, g=energy groups
integer                   :: info            
real(real64), allocatable :: sigma_a(:,:), sigma_s(:,:), sigma_tr(:,:),sigma_fnu(:), sigma_r(:,:)      
real(real64)              :: w, ratio, dr1, dr2   
real(real64), allocatable :: b_old(:,:), S(:)
!------------------------Used for Iteration loop--------------------------!
integer                   :: j, j_width ,j_tot, i,p                
real(real64)              :: k_error, b_error
real(real64)              :: m_old                           
real(real64)              :: k , k_old    
real(real64)              :: min_error=0.0001
integer                   :: MAX_ITERATIONS = 1000 

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
allocate (A(n_tot,n_tot,g),b(n_tot,g),ipiv(n_tot),b_old(n_tot,g),S(n_tot))   ! Used in diff eqs
w=10
!---------------------Matrix Constants Declaration------------------------!

allocate (sigma_a(g,2), sigma_s(g,g), sigma_tr(g,2),chi(g),sigma_fnu(g),sigma_r(g,2))

j_tot = 0
j_width=0


chi=1_real64   
chi(1)=0.5_real64   

sigma_s=0.0

!------------------------G=4 declarations---------------------------------!
if(g .eq. 4)then
    sigma_s(1,1)=0.37045
    sigma_s(1,2)=0.04152
    sigma_s(1,3)=0.00001

    sigma_s(2,2)=0.98285
    sigma_s(2,3)=0.07459                !scattering cross sections(From, to)
    sigma_s(2,4)=0.01371

    sigma_s(3,3)=0.76110
    sigma_s(3,4)=0.31856

    sigma_s(4,3)=0.00085
    sigma_s(4,4)=1.96607

    sigma_a(1,2)=0.00051
    sigma_a(2,2)=0.00354
    sigma_a(3,2)=0.01581
    sigma_a(4,2)=0.04637

                    ! D already
    sigma_tr(1,2)=1/(3.0*0.20608)
    sigma_tr(2,2)=1/(3.0*0.60215)
    sigma_tr(3,2)=1/(3.0*0.56830)
    sigma_tr(4,2)=1/(3.0*1.21110)


             !------------ Core is below-------------!              
    sigma_fnu(1)=0.009572 
    sigma_fnu(2)=0.001193
    sigma_fnu(3)=0.01768
    sigma_fnu(4)=0.18514

    sigma_a(1,1)=0.004946
    sigma_a(2,1)=0.002840
    sigma_a(3,1)=0.03053
    sigma_a(4,1)=0.1210

                    ! D already
    sigma_tr(1,1)=2.1623
    sigma_tr(2,1)=1.0867
    sigma_tr(3,1)=0.6318
    sigma_tr(4,1)=0.3543

    sigma_r(1,1)=0.08795
    sigma_r(2,1)=0.06124
    sigma_r(3,1)=0.09506
    sigma_r(4,1)=0.1210
!------------------------G=2 declarations---------------------------------!
elseif(g .eq. 2)then
    sigma_s(1,2)=0.0494

    sigma_a(1,2)=0.0004
    sigma_a(2,2)=0.0197

    sigma_tr(1,2)=1.13
    sigma_tr(2,2)=0.16


    sigma_r(1,2)=0.0494
    sigma_r(2,2)=0.0197


             !------------ Core is below-------------!              
    sigma_fnu(1)=0.008476 
    sigma_fnu(2)=0.18514


    sigma_a(1,1)=0.01207
    sigma_a(2,1)=0.121


                    ! D already
    sigma_tr(1,1)=1.2627
    sigma_tr(2,1)=0.3543


    sigma_r(1,1)=0.02619
    sigma_r(2,1)=0.121

else
    print *,'I cant do that John'
    goto 12
endif


!--------------------------Matrix Declaration-----------------------------!


100 dr1 = w/n                               ! Step Size
dr2 = w*ratio/n_flector
A=0
!print *,sigma_tr(1,2)

do p=1,g
   A(1,1,p) = (0.5_real64   *sigma_r(p,1))+(2.0*sigma_tr(p,1)/(dr1**2))
   A(1,2,p) = -sigma_tr(p,1)/(dr1**2.0)

   A(n_tot,n_tot,p) = sigma_r(1,2)+2.0_real64   *sigma_tr(p,2)/(dr2**2)
   A(n_tot,n_tot-1,p) = -sigma_tr(p,2)/(dr2**2)*(1-1/(2*n_tot-1))

   A(n,n-1,p) = -sigma_tr(p,1)/(dr1)  -sigma_tr(p,1)/(2.0*w)
   A(n,n+1,p) =  -sigma_tr(p,2)/(dr1)  -sigma_tr(p,2)/(2.0*w)
   A(n,n,p)=(-sigma_tr(p,2)/(dr2)+sigma_tr(p,2)/(2.0*w))+(-A(n,n-1,p))+(sigma_r(p,2)*dr2+sigma_r(p,1)*dr1)/2.0

   do i=2 , n-1
 
      A(i,i,p) = sigma_r(p,1) + 2.0_real64   *sigma_tr(p,1)/(dr1**2)
      A(i,i+1,p) = -(1+1/(2.0*i-1))*sigma_tr(p,1)/(dr1**2)
      A(i,i-1,p) = -(1-1/(2.0*i-1))*sigma_tr(p,1)/(dr1**2)

   enddo
   do i=n+1 , n_tot

      A(i,i,p) = sigma_r(p,2)+2.0_real64   *sigma_tr(p,2)/(dr2**2)
      A(i,i+1,p) = -sigma_tr(p,2)/(dr2**2)*(1+1/(2.0*i-1))
      A(i,i-1,p) = -sigma_tr(p,2)/(dr2**2)*(1-1/(2.0*i-1))

   enddo
enddo
! Pretty matrix printing
! From http://jblevins.org/log/array-write
! do i=1,n_tot
!     write(*,"(10g15.5)") ( A(i,j,1), j=1,n_tot )
! enddo

    !LU factorization of a general M-by-N matrix A
    ! See: http://www.netlib.org/lapack/explore-html/d3/d6a/dgetrf_8f.html
    do p=1,g    
        call DGETRF(n_tot, n_tot, A(:,:,p), n_tot, ipiv, info)
!        if (info /= 0) stop 'Matrix is numerically singular!'
!        print *,'HELP ME'
    enddo
!-------------------------Business loop-----------------------------------!
!-------------------------------------------------------------------------!
!guess initial

b = 1.0
k=1.0
S=0
    do p=1,g
        do i=1,n_tot
        if(i .lt. n)then
            S(i)= S(i)+sigma_fnu(p)*b(i,p)
        else
            S(i)=0.0   
        endif
        enddo
    enddo
S(1)=S(1)/2.0
j=0

k_error = 1.0
b_error = 1.0


do while ( ((b_error .gt. min_error) .or. (k_error .gt. min_error)) .and. (j .lt. MAX_ITERATIONS))

    b_old = b     
    k_old = k     


!------------------Eigenvalue Search--------------------------------------!
    do p=1,g
        if(p .gt. 1)then
            do i=p, g
                b(:,p)= b(:,p) + sigma_s(i,p)*b(:,i)
            enddo
        endif
        if(p .eq. 1) b(:,p)=b(:,p)+S*b(:,p)/k
        call DGETRS('N', n_tot, 1, A(:,:,p), n_tot, ipiv, b(:,p), n_tot, info) 
        if (info /= 0) stop 'Solution of the linear system failed!'

    enddo
    j = j+1        
!-------------------------------------------------------------------------!
!----------------------Error tests----------------------------------------!
    m_old=norm2(S)

 

    S=0
    do p=1,g
        do i=1,n_tot
        if(i .lt. n+1)then
            S(i)= S(i)+sigma_fnu(p)*b(i,p)

        else
            S(i)=0.0   
        endif
        enddo
    enddo
    S(1)=S(1)/2.0

    k = k_old *norm2(S)/ m_old                 
    k_error = abs((k-k_old)/k)
    
    do p=1,g
        chi(p)=maxval(abs(b(:,p) - b_old(:,p)))/maxval(b_old(:,p))
    enddo
    b_error=maxval(chi)

   do i=1,g
        b(:,i)  = b(:,i) / norm2(b(:,i) )
    enddo
enddo
!print *, (b_error .gt. min_error)  
!print *, (k_error .gt. min_error)
!print *, (j .lt. MAX_ITERATIONS)
j_tot=j_tot+j
if(j_width .eq. 500)goto 12
print *,j_width
!-------------------------width adjustments for k=1-----------------------!
print *,''
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

12 continue
!-------------------------VOMIT RESULTS-----------------------------------!
print *,'Number of iterations until convergence:  ', j_tot
print *,'Critical width:  ', 2.0*w

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


