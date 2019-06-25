
	subroutine FEForce(tf)
! ---------------------------------------------------------------------
! -                                                                   -
! -  Purpose                                                          -
! -    Update FE nodal force                                          -
! -                                                                   -
! -  Input                                                            -
! -    tf - value of time function at current time                    -
! -                                                                   -
! ---------------------------------------------------------------------
	use ParticleData
	use MaterialData
	use ElementData
	use Simulation
	implicit none

	real(8), intent(in) :: tf

	real(8):: vol, den			! element volume and density 
	real(8):: xyz(3,8), v(3,8)	! element nodal coordinates and velocities 
	real(8):: ijac(3,3)			! element inverse Jacobian matrix
	real(8):: de(6), vort(3)	! incremental strain and spin
	real(8):: pkxj(4,3)			! derivatives of shape function with respect to coordinates
	integer:: ie, i, j, p		! loop counter

	real(8) :: dte, at			! element time step and distortion
	real(8) :: c				! element wave speed and hourglass control constant

	type(Element),  POINTER :: el
	type(Material), POINTER :: mat
	type(History),  POINTER :: his, phis

	EleDistortion = 1.0
	
	Acc = Fp*tf

!	loop over all elements
	do ie = 1, nb_element
		el => element_list(ie)

		his => history_list(el%nHistory)	! element history variables

!		retrieve the nodal coordinates and velocities of the element
		do j = 1, 8
			p = el%nNode(j)
			xyz(:,j) = Pos(:,p)		! Element nodal coordinates
			v(:,j) = Vel(:,p)		! Element nodal velocities
		end do

!		evaluate the element inverse Jacobian matrix and volume
		call Jacobian(xyz, vol, ijac)

!		the element is severely distorted
		if (vol .le. 0) then
			write(*, 100) ie
100			format(1x, "*** Error *** Jacobian determinant of element ", i6, " is negative ! ")
			stop
		end if

!		evaluate the increament strain and spin
		call Strain(DTk, v, ijac, de, vort, pkxj)

		mat => mat_list(el%mat)				! element material data
		den = mat%density*his%VOL/vol		! Current density
		c   = mat%WaveSpeed					! wave speed

!		evaluate element time step size and distortion
		call ElementTimeStep(xyz, vol, c, DTk1, EleDistortion)

!		update stress by constitution law
		call Constitution(de, vort, den, his, mat, vol)

!		evaluate and accumulate the equivalent nodal forces due to element stress
		call NodalForce(el%nNode, pkxj, his, vol)

!		evaluate and accumulate the hourglass-resisting forces
		call HourGlassForce(el%nNode, v, xyz, pkxj, vol, den, c)

	end do

!	determine the new time step size
	if (abs(DTk1-1.0e6) .gt. 1.0e-5) then
		DTk1 = DTk1*DTScale
		DTk1 = min(DTk1, 1.05*DTk1Old)
		DTk1Old = DTk1
	end if

	return
	end subroutine FEForce


	subroutine NodalForce(nNode, pkxj, his, vol)
! -------------------------------------------------------------------------
! - Purpose                                                               -
! -    Calculate and accumulate the nodal forces due to element stress    -
! -                                                                       -
! - Input                                                                 -
! -    nNode - element nodes                                              -
! -    pkxj  - derivatives of shape function with respect to coordinates  -
! -    his   - element history variable                                   -
! -    vol   - element volume                                             -
! -                                                                       -
! -------------------------------------------------------------------------
	use ParticleData
	use MaterialData
	implicit none

	real(8), intent(in) :: pkxj(4,3), vol
	integer, intent(in) :: nNode(8)
	type(History), intent(in) :: his

	real(8) :: sm, sig(6)	! mean stress and stresses 
	real(8) :: f(3,8)		! nodal forces due to element stresses
	integer :: i, k, p		! loop counter

	sm = his%SM  ! Mean stress

	sig(1) = (his%SDS(1) + sm)*vol
	sig(2) = (his%SDS(2) + sm)*vol
	sig(3) = (his%SDS(3) + sm)*vol
	sig(4) = (his%SDS(4))*vol
	sig(5) = (his%SDS(5))*vol
	sig(6) = (his%SDS(6))*vol

!	evaluate nodal forces
	do k = 1, 4
		f(1,k) = -(sig(1)*pkxj(k,1) + sig(6)*pkxj(k,2) + sig(5)*pkxj(k,3))
		f(2,k) = -(sig(2)*pkxj(k,2) + sig(6)*pkxj(k,1) + sig(4)*pkxj(k,3))
		f(3,k) = -(sig(3)*pkxj(k,3) + sig(4)*pkxj(k,2) + sig(5)*pkxj(k,1))
	end do

	do i = 1, 3
		f(i,5) = -f(i,3)
		f(i,6) = -f(i,4)
		f(i,7) = -f(i,1)
		f(i,8) = -f(i,2)
	end do

!	accumulate nodal forces
	do k = 1, 8
		p = nNode(k)
		Acc(:,p) = Acc(:,p) + f(:,k)
	end do

	end subroutine NodalForce


	subroutine UpdateFEGeometry()
! -------------------------------------------------------------------------
! - Purpose                                                               -
! -    Compute accelerations, apply displacement b.c.'s, and update       -
! -       velocity and position for FE nodes                              -
! -                                                                       -
! -------------------------------------------------------------------------
	use Simulation
	use ElementData
	use ParticleData
	implicit none

	integer p, nDir, nUnit
	real(8) :: DT2

	DT2 = (Dtk+DTk1)*0.5

!	Compute accelerations of FE nodes

	do p = 1, nb_node
		Acc(:,p) = Acc(:,p) / Mp(p)
	end do

!	Apply displacement boundary conditions - rigid wall

	nDir = plane%nDir		! unit normal vector of the rigid wall

	if (nDir .ne. 0) then
		nUnit = sign(1, nDir)

		do p = 1, nb_node

!			the rigid wall is perpendicular to x axis
			if (abs(nDir).eq.1 .and. nUnit*(plane%coor-Pos(1,p)).ge.0) then
				if (nUnit*Vel(1,p).lt.0) Vel(1,p) = 0
				if (nUnit*Acc(1,p).lt.0) Acc(1,p)  = 0
			end if

!			the rigid wall is perpendicular to y axis
			if (abs(nDir).eq.2 .and. nUnit*(plane%coor-Pos(2,p)).ge.0) then
				if (nUnit*Vel(2,p).lt.0) Vel(2,p) = 0
				if (nUnit*Acc(2,p).lt.0) Acc(2,p)  = 0
			end if

!			the rigid wall is perpendicular to z axis
			if (abs(nDir).eq.3 .and. nUnit*(plane%coor-Pos(3,p)).ge.0) then
				if (nUnit*Vel(3,p).lt.0) Vel(3,p) = 0
				if (nUnit*Acc(3,p).lt.0) Acc(3,p)  = 0
			end if

		end do
	end if

!	Update velocity and position for FE nodes
	Vel = Vel + Acc*DT2
	Pos = Pos + Vel*DTk1

	end subroutine UpdateFEGeometry


	subroutine HourGlassForce(nNode, v, xyz, pkxj, vol, den, c)
! -------------------------------------------------------------------------
! - Purpose                                                               -
! -    evaluate and accumulate hourglass-resisting nodal forces           -
! -                                                                       -
! - Input                                                                 -
! -    nNode - element nodes                                              -
! -    v     - element nodeal velocities                                  -
! -    xyz   - same as before                                             -
! -    pkxj  - same as before                                             -
! -    vol   - element volume                                             -
! -    den   - element density                                            -
! -    c     - material sound speed                                       -
! -                                                                       -
! -------------------------------------------------------------------------
	use ParticleData
	use MaterialData
	implicit none

	integer, intent(in) :: nNode(8)
	real(8), intent(in) :: v(3,8), xyz(3,8), pkxj(4,3) 
	real(8), intent(in) :: vol, den, c

	real(8) :: v3478(3), v2358(3), v1467(3), v1256(3)
	real(8) :: hgr(3,4)		! projection of velocity filed on to hourglass modes 
	real(8) :: f(3,8)		! coefficient and hourglass-resisting forces

	real(8) :: hap(3), ham(3), hbp(3), hbm(3)

	real(8) :: a, v32, ah, al
	integer :: i, k, p		! loop counter

	if (HourGlass%method.eq.0 .or. HourGlass%Qhg.le.1.0e-6) return

	if (HourGlass%method .eq. 1) then
		call HGForceStandard(v, vol, den, c, f)
	else if (HourGlass%method .eq. 2) then
		call HGForceFlangan(v, xyz, pkxj, vol, den, c, f)
	else
		stop "***Error*** Invalid Hourglass control method !"
	end if

!	accumulate the hourglass-resisting forces
	do k = 1, 8
		p = nNode(k)
		Acc(:,p) = Acc(:,p) + f(:,k)
	end do

	end subroutine HourGlassForce


	subroutine HGForceStandard(v, vol, den, c, f)
! -------------------------------------------------------------------------
! - Purpose                                                               -
! -    evaluate hourglass-resisting nodal forces using standard method    -
! -                                                                       -
! - Input                                                                 -
! -    v     - element nodeal velocities                                  -
! -    vol   - element volume                                             -
! -    den   - element density                                            -
! -    c     - material sound speed                                       -
! -                                                                       -
! - Output                                                                -
! -    f     - hourglass-resisting nodal forces                           -
! -                                                                       -
! -------------------------------------------------------------------------
	use ParticleData
	use MaterialData
	implicit none

	real(8), intent(in) :: v(3,8), vol, den, c
	real(8), intent(out):: f(3,8)

	real(8) :: v3478(3), v2358(3), v1467(3), v1256(3)
	real(8) :: hgr(3,4)		! projection of velocity filed on to hourglass modes 

	real(8) :: hap(3), ham(3), hbp(3), hbm(3)

	real(8) :: qh, a, v32, ah, al
	integer :: i, k, p		! loop counter

	qh = HourGlass%Qhg

	a  = qh*den/4.0
	v32= vol**(2.0/3.0)
	ah = a*v32*c
	al = a*v32*(100*qh)

	do i = 1, 3
		v3478(i) = v(i,3) - v(i,4) - v(i,7) + v(i,8)
		v2358(i) = v(i,2) - v(i,3) - v(i,5) + v(i,8)
		v1467(i) = v(i,1) - v(i,4) - v(i,6) + v(i,7)
		v1256(i) = v(i,1) - v(i,2) - v(i,5) + v(i,6)
	end do

	do i = 1, 3
		hgr(i,1) = v1467(i) - v2358(i)
		hgr(i,2) = v1467(i) + v2358(i)
		hgr(i,3) = v1256(i) - v3478(i)
		hgr(i,4) = v1256(i) + v3478(i)
	end do

	do i = 1, 3
		do k = 1, 4
			hgr(i,k) = hgr(i,k)*(ah + abs(al*hgr(i,k)))
		end do
	end do

	do i = 1, 3
		hap(i) = hgr(i,1) + hgr(i,4)
		ham(i) = hgr(i,1) - hgr(i,4)
		hbp(i) = hgr(i,2) + hgr(i,3)
		hbm(i) = hgr(i,2) - hgr(i,3)
	end do

!	evaluate the hourglass-resisting forces
	do i = 1, 3
		f(i,1) = -hap(i) - hbp(i)
		f(i,2) =  hap(i) - hbm(i)
		f(i,3) = -hap(i) + hbp(i)
		f(i,4) =  hap(i) + hbm(i)
		f(i,5) = -ham(i) + hbp(i)
		f(i,6) =  ham(i) + hbm(i)
		f(i,7) = -ham(i) - hbp(i)
		f(i,8) =  ham(i) - hbm(i)
	end do

	end subroutine HGForceStandard


	subroutine HGForceFlangan(v, xyz, pkxj, vol, den, c, f)
! -------------------------------------------------------------------------
! - Purpose                                                               -
! -    evaluate hourglass-resisting nodal forces using Flangan-Belytschko -
! -       method                                                          -
! -                                                                       -
! - Input                                                                 -
! -    v     - element nodeal velocities                                  -
! -    vol   - element volume                                             -
! -    den   - element density                                            -
! -    c     - material sound speed                                       -
! -                                                                       -
! - Output                                                                -
! -    f     - hourglass-resisting nodal forces                           -
! -                                                                       -
! -------------------------------------------------------------------------
	use ParticleData
	use MaterialData
	implicit none

	real(8), intent(in) :: v(3,8), xyz(3,8), pkxj(4,3) 
	real(8), intent(in) :: vol, den, c
	real(8), intent(out):: f(3,8)

	integer, parameter :: h(4,4) = (/ 1,-1, 1,-1, & 
									  1, 1,-1,-1, &
									  1,-1,-1, 1, &
									  1,-1, 1,-1  /)

	real(8), parameter :: ss(4,4) = (/ 2.0,-2.0, 2.0,-2.0,  &
									  -2.0,-2.0, 2.0, 2.0,  &
									  -2.0, 2.0, 2.0,-2.0,  &
									   0.0, 0.0, 0.0, 0.0  /)

	real(8) :: x3478(3), x2358(3), x1467(3), x1256(3)
	real(8) :: hgr(3,4), beta(3,4), gama(4,8), g(3,4)

	real(8) :: ah, qh
	integer :: i, j, k, p		! loop counter

	qh = HourGlass%Qhg

	ah = -qh*c*den*vol**(2.0/3.0)/4.0

	do i = 1, 3
		x3478(i) = xyz(i,3) - xyz(i,4) - xyz(i,7) + xyz(i,8)
		x2358(i) = xyz(i,2) - xyz(i,3) - xyz(i,5) + xyz(i,8)
		x1467(i) = xyz(i,1) - xyz(i,4) - xyz(i,6) + xyz(i,7)
		x1256(i) = xyz(i,1) - xyz(i,2) - xyz(i,5) + xyz(i,6)
	end do

	do i = 1, 3
		beta(i,1) = x1467(i) - x2358(i)
		beta(i,2) = x1467(i) + x2358(i)
		beta(i,3) = x1256(i) - x3478(i)
		beta(i,4) = x1256(i) + x3478(i)
	end do

	do j = 1, 4
		do k = 1, 4
			gama(j,k) = h(k,j)
			do i = 1, 3
				gama(j,k) = gama(j,k) - beta(i,j)*pkxj(k,i) 
			end do
		end do 
	end do

	do j = 1, 4
		gama(j,5) = ss(1,j) - gama(j,3)
		gama(j,6) = ss(2,j) - gama(j,4)
		gama(j,7) = ss(3,j) - gama(j,1)
		gama(j,8) = ss(4,j) - gama(j,2)
	end do

	do i = 1, 3
		do j= 1, 4
			g(i,j) = 0.0
			do k = 1, 8
				g(i,j) = g(i,j) + v(i,k)*gama(j,k)
			end do
		end do
	end do

!	evaluate the hourglass-resisting forces
	do i = 1, 3
		do k = 1, 8
			f(i,k) = 0.0
			do j = 1, 4
				f(i,k) = f(i,k) + g(i,j)*gama(j,k)
			end do
			f(i,k) = ah*f(i,k)
		end do
	end do

	end subroutine HGForceFlangan


	subroutine Jacobian(xyz, vol, ijac)
! ---------------------------------------------------------------------
! -  Purpose                                                          -
! -    Evaluate Jacobian determinant and inverse Jacobian matrix at   -
! -       the center of the element                                   -
! -                                                                   -
! -  Input                                                            -
! -    xyz  - Coordinates of vertexes of the element                  -
! -                                                                   -
! -  Output                                                           -
! -    vol  - Volume of the element                                   -
! -    ijac - Inverse Jacobian matrix                                 -
! -                                                                   -
! ---------------------------------------------------------------------
!	use lin_sol_gen_int		! IMSL Fortran 90 MP LIBRARY 
							! need to link with SF90MP.LIB SMATHD.LIB
	use imsl
	implicit none

	real(8), intent(in)  :: xyz(3,8)  !  Coordinates of vertexes of the element
	real(8), intent(out) :: vol       !  Volume of the element
	real(8), intent(out) :: ijac(3,3) !  Inverse Jacobian

	real(8):: det_       !  Determinant of jacobian
	real(8):: jac(3,3)  !  Jacobian matrix

	real(8):: x17(3), x28(3), x35(3), x46(3)
	real(8):: b(3,0), x(3,0), dett(2)   !  Required for lin_sol_gen
	integer:: j

!   JACOBIAN MATRIX at the center of a element (one-point quadrature)
	do j = 1, 3
		x17(j) = xyz(j,7) - xyz(j,1)
		x28(j) = xyz(j,8) - xyz(j,2)
		x35(j) = xyz(j,5) - xyz(j,3)
		x46(j) = xyz(j,6) - xyz(j,4)
	end do

	do j = 1, 3
		jac(1,j) = (x17(j) - x28(j) - x35(j) + x46(j))/8.0
		jac(2,j) = (x17(j) + x28(j) - x35(j) - x46(j))/8.0
		jac(3,j) = (x17(j) + x28(j) + x35(j) + x46(j))/8.0
	enddo

!	Compute the matrix inverse and its determinant. 
 	call lin_sol_gen(jac, b, x, nrhs=0, ainv=ijac, det=dett) 

	det_ = abs(dett(1))**dett(2) * (dett(1))/abs(dett(1))

!	element volume
	vol = 8.0*det_

	return
	end subroutine Jacobian


	subroutine Strain(DTk, v, ijac, de, vort, pkxj)
! ---------------------------------------------------------------------
! - Purpose                                                           -
! -    Calculate incremental strain and vorts                         -
! -                                                                   -
! - Input                                                             -
! -    DTk   - Time step size                                         -
! -    v     - element nodal velocities                               -
! -    ijac  - Inverse Jacobian                                       -
! -                                                                   -
! - Output                                                            -
! -    de    - Velocity strain rates                                  -
! -    vort  - Non zero entries of incremental spin tensor            -
! -    pkxj  - derivatives of shape function d(Phi_k)/d(x_j)          -
! -                                                                   -
! ---------------------------------------------------------------------
	implicit none

	real(8), intent(in)  :: DTk        ! Current time step size
	real(8), intent(in)  :: v(3,8)     ! element nodal velocity
	real(8), intent(in)  :: ijac(3,3)  ! Inverse Jacobian
	real(8), intent(out) :: de(6)      ! Velocity strain rates
	real(8), intent(out) :: vort(3)    ! incremental spin
	real(8), intent(out) :: pkxj(4,3)  ! derivatives of shape function d(Phi_k)/d(x_j)

	real(8) :: v17(3), v28(3), v35(3), v46(3)
	real(8) :: dij(3,3)
	integer i, j

!   evaluate the derivatives of shape function with respect to spatial coordinates x_j
	do j = 1, 3
		pkxj(1,j) = (-ijac(j,1) - ijac(j,2) - ijac(j,3))/8.0
		pkxj(2,j) = ( ijac(j,1) - ijac(j,2) - ijac(j,3))/8.0
		pkxj(3,j) = ( ijac(j,1) + ijac(j,2) - ijac(j,3))/8.0
		pkxj(4,j) = (-ijac(j,1) + ijac(j,2) - ijac(j,3))/8.0
	end do
!
	do i = 1, 3
		v17(i) = v(i,1) - v(i,7)
		v28(i) = v(i,2) - v(i,8)
		v35(i) = v(i,3) - v(i,5)
		v46(i) = v(i,4) - v(i,6)
	end do

	do i = 1, 3
		do j = 1, 3
			dij(i,j) = v17(i)*pkxj(1,j)+v28(i)*pkxj(2,j)+v35(i)*pkxj(3,j)+v46(i)*pkxj(4,j)
		end do
		de(i) = dij(i,i)   ! Normal velocity strain rates
	end do

	de(4) = dij(2,3) + dij(3,2)   ! Engineering shear velocity strain rates
	de(5) = dij(1,3) + dij(3,1)
	de(6) = dij(1,2) + dij(2,1)

	vort(1) = (dij(2,3) - dij(3,2))*0.5   ! spin : wyz
	vort(2) = (dij(1,3) - dij(3,1))*0.5   ! spin : wxz
	vort(3) = (dij(1,2) - dij(2,1))*0.5   ! spin : wxy

	de = de * DTk		! Incremental strain
	vort = vort *DTk	! Incremental vort

	end subroutine Strain


	subroutine ElementTimeStep(xyz, vol, c, DTk1, EleDistortion)
! -------------------------------------------------------------------------
! - Purpose                                                               -
! -    evaluate element time step size and distortion                     -
! -                                                                       -
! - Input                                                                 -
! -    xyz  - element nodal coordinates                                   -
! -    vol  - element volume                                              -
! -    c    - wave speed                                                  -
! -    DTk1          - current time step size                             -
! -    EleDistortion - current element distortion                         -
! -                                                                       -
! - Output                                                                -
! -    DTk1          - current time step size                             -
! -    EleDistortion - current element distortion                         -
! -                                                                       -
! -------------------------------------------------------------------------
	implicit none

	real(8), intent(in) :: xyz(3,8), vol, c
	real(8), intent(inout) :: DTk1, EleDistortion

	integer :: fac(4,6) = (/ 1,2,3,4, 5,6,7,8, 1,2,6,5, &	! element surface definition
							 2,3,7,6, 3,4,8,7, 4,1,5,8 /)

	real(8) ::dt, at
	real(8) :: e, f, g, atest, areal, aream, x13(3), x24(3), fs(3), ft(3)
	integer :: i, j, k1, k2, k3, k4

	areal = 1.0e20
	aream = 0.0

!	Loop over all six surfaces of the element
	do j = 1, 6
		do i = 1, 3
			k1 = fac(1,j)
			k2 = fac(2,j)
			k3 = fac(3,j)
			k4 = fac(4,j)

			x13(i) = xyz(i,k3) - xyz(i,k1)
			x24(i) = xyz(i,k4) - xyz(i,k2)

			fs(i) = x13(i) - x24(i)
			ft(i) = x13(i) + x24(i)
		end do

		e = fs(1)*fs(1) + fs(2)*fs(2) + fs(3)*fs(3)
		f = fs(1)*ft(1) + fs(2)*ft(2) + fs(3)*ft(3)
		g = ft(1)*ft(1) + ft(2)*ft(2) + ft(3)*ft(3)

		atest = e*g - f*f	! area/4

		aream = max(atest, aream)
		areal = min(atest, areal)

	end do

	at = areal/aream
	dt = 4*vol/sqrt(aream)
	dt = dt/c

	DTk1 = min(DTk1, dt)
	EleDistortion = min(at, EleDistortion)

	end subroutine ElementTimeStep


	subroutine KineticE
! ------------------------------------------------------------------
! -                                                                -
! -  purpose:                                                      -
! -     calculate kinematic energy                                 -
! -                                                                -
!-------------------------------------------------------------------
	use ParticleData
	use Simulation
	implicit none

	integer :: p	! loop counter

	EngKinetic = 0.0

	do p = 1, nb_node
		EngKinetic = EngKinetic + dot_product(Vel(:,p), Vel(:,p))*Mp(p)*0.5d0
	end do
	end subroutine KineticE
