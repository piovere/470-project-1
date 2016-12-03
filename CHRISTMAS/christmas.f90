program christmas

! CHRISTMAS: CRIticality Search for neuTron Multi-group reActor Slab

! Neutron Flux desperation Group 3
! DGETRS - solve a system of linear equations A * X = B or 
! See: http://www.netlib.org/lapack/explore-html/d6/d49/dgetrs_8f.html

use iso_fortran_env
implicit none

! Guess core geometry and composition

! Guess initial fission source S(0) and k(0)

! Inner iterations
    ! matrix * flux(n+1) = (1 / k) * S(n)
    ! S(n+1) = F * flux(n+1)
    ! k(n+1) = k * S(n+1) / S(n)
    ! Convergence test
        ! abs(k(n+1) - k(n)) / k(n) < max_error
        ! abs(mag(S(n+1) - S(n)) / S(n)) < max_error

! abs(k - 1.00000) < max_error

end program christmas