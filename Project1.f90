program Mock
! One-speed diffusion neutron flux using 
! LAPACK linear system solver example, Ondrej Chvala <ochvala@utk.edu> 2014-10-31
!
! Uses LAPACK's DGETRF to compute an LU factorization of a matrix A
! using partial pivoting with row interchanges, and then DGETRS to 
! solve a system of linear equations A*x=b by the LU factorization 
! already computed by DGETRF.
! 
use iso_fortran_env
implicit none

real(real64), allocatable :: A(:,:)  ! matrix to LU decompose
real(real64), allocatable :: S(:)    ! b vector (SOURCE)
real(real64), allocatable :: ipiv(:) ! pivoting array for LAPACK
integer                   :: N       ! problem size
integer                   :: info    ! success flag for LAPACK
!character(len=200)        :: fmt     ! formatting string for pretty printing
integer                   :: i       ! loop constants
real                      :: m, D, sig, W, DeltaX  ! D is the diffusion coef and
                                                   ! sig is macro cross section

! External procedures defined in LAPACK
external DGETRF, DGETRS

! Input data
print *,'How many nodes would you like?' 
read *, N                           ! Number of discrete points
allocate (ipiv(1:N-1), S(1:N-1), A(1:N-1,1:N-1))
A=0.0
S=0.0
print *, 'What would you like the geometry to be?(1 for spherical, 2 for cylindrical, 3 for cartesian)'
!read *, G
print * ,'Too Bad its going to be a plane'
sig=0.15
m=N-1
W=10.0                                     !Length
DeltaX = W/(m)
D=9.0 / (DeltaX**2.0)                               !! This is the Diffusion coeffecient / x step size ^2
!!------------------------Source Initializing----------------!!
S(1)=1.0E8 /(2.0*DeltaX)
!S(N)=0.0
!print *,S(1)
!!------------------------Declaring A------------------------!!
A(1,1)=D+0.5*sig
A(m,m)=D+0.5*sig
A(1,2)=-D
A(m,m-1)=-D
do i=2, m-1

    A(i,i) = sig+2.0*D
    A(i-1,i) = -D
    A(i+1,i) = -D

enddo

! DGETRF computes an LU factorization of a general M-by-N matrix A
! using partial pivoting with row interchanges.
! See: http://www.netlib.org/lapack/explore-html/d3/d6a/dgetrf_8f.html
call DGETRF(m, m, A, m, ipiv, info)
if (info /= 0) stop 'Matrix is numerically singular!'

! DGETRS - solve a system of linear equations A * X = B or 
! A' * X = B with a general N-by-N matrix A using the LU 
! factorization computed by DGETRF. 
! See: http://www.netlib.org/lapack/explore-html/d6/d49/dgetrs_8f.html
call DGETRS('N', m, 1, A, m, ipiv, S, m, info)  
if (info /= 0) stop 'Solution of the linear system failed!'

! Print solution
open( unit = 77, file= "Project1.dat")
write(77,'(E10.4)') S

end program Mock

