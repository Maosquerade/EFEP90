module ParticleData

	integer, parameter :: nDim = 3;		! Number of dimension
	integer :: nb_node = 0				! number of nodes
	integer :: nb_node_root = 0

	real(8), allocatable :: Acc(:,:)	! particle acceleration
	real(8), allocatable :: Pos(:,:), Pos_root(:,:)	! particle position
	real(8), allocatable :: Vel(:,:), Vel_root(:,:)	! particle velocity
	real(8), allocatable :: Fp(:,:), Fp_root(:,:)		! external load amplitude
	real(8), allocatable :: Mp(:), Mp_root(:)		! Mass of material point

end module ParticleData
