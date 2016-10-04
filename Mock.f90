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
real                      :: D, sig  ! D is the diffusion coef and sig is macro cross section

! External procedures defined in LAPACK
external DGETRF, DGETRS

! Input data
print *,'How many nodes would you like?' 
read *, N                           ! Number of discrete points
allocate (ipiv(1:N), S(1:N), A(1:N,1:N))
A=0.0
S=0.0
print *, 'What would you like the geometry to be?(1 for spherical, 2 for cylindrical, 3 for cartesian)'
!read *, G
print * ,'Too Bad its going to be a plane'
sig=3.0E-2
D=1.0E-2                                !! This is the Diffusion coeffecient / x step size ^2
!!------------------------Source Initializing----------------!!
S(1)=1.0E8
S(N)=0.0
!!------------------------Declaring A------------------------!!
do i=1, N

    A(i,i) = sig+2.0*D
    if (i .ne. 1) then
       A(i,i-1) = -D
    endif
    if (i .ne. N) then
       A(i,i+1) = -D
    endif
	
enddo


! DGETRF computes an LU factorization of a general M-by-N matrix A
! using partial pivoting with row interchanges.
! See: http://www.netlib.org/lapack/explore-html/d3/d6a/dgetrf_8f.html
call DGETRF(N, N, A, N, ipiv, info)
if (info /= 0) stop 'Matrix is numerically singular!'

! DGETRS - solve a system of linear equations A * X = B or 
! A' * X = B with a general N-by-N matrix A using the LU 
! factorization computed by DGETRF. 
! See: http://www.netlib.org/lapack/explore-html/d6/d49/dgetrs_8f.html
call DGETRS('N', N, 1, A, N, ipiv, S, N, info)  
if (info /= 0) stop 'Solution of the linear system failed!'

! Print solution
open( unit = 77, file= "Project1.dat")
write(77,'(E10.4)') S

end program Mock

