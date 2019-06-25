! *********************************************************
!   Matrix inverse: order 3
! *********************************************************
	subroutine Inverse3(jac, ijac, det_)
      
      implicit none

      real(8), intent(in) :: jac(3,3)
      real(8), intent(out) :: ijac(3,3)
      real(8), intent(out) :: det_
      real(8) :: idet

!Compute for determinant
      det_ = jac(1,1)*jac(2,2)*jac(3,3) + jac(2,1)*jac(3,2)*jac(1,3) + jac(3,1)*jac(1,2)*jac(2,3) &
           -(jac(3,1)*jac(2,2)*jac(1,3) + jac(2,1)*jac(1,2)*jac(3,3) + jac(1,1)*jac(3,2)*jac(2,3))

!Compute for inverse matrix
      ijac(1,1) = jac(2,2)*jac(3,3) - jac(2,3)*jac(3,2)
	  ijac(2,1) = jac(2,3)*jac(3,1) - jac(2,1)*jac(3,3)
	  ijac(3,1) = jac(2,1)*jac(3,2) - jac(3,1)*jac(2,2)
	
	  ijac(1,2) = jac(1,3)*jac(3,2) - jac(1,2)*jac(3,3)
	  ijac(2,2) = jac(1,1)*jac(3,3) - jac(1,3)*jac(3,1)
	  ijac(3,2) = jac(1,2)*jac(3,1) - jac(1,1)*jac(3,2)

	  ijac(1,3) = jac(1,2)*jac(2,3) - jac(2,2)*jac(1,3)
	  ijac(2,3) = jac(1,3)*jac(2,1) - jac(1,1)*jac(2,3)
	  ijac(3,3) = jac(1,1)*jac(2,2) - jac(1,2)*jac(2,1)

      idet = 1.0 / det_
	  ijac = ijac * idet

	  return

	end subroutine Inverse3
