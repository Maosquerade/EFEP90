module ElementData

!	Type ELEMENT is used to describe a element
	type Element
		integer:: nNode(8)		!  Element connectivity data
		integer:: mat =0		! material set number (1 ~ nb_mat)
		integer:: nHistory =0	!  Index of element history in HistoryList
	end type Element

	integer:: nb_element	! number of elements
	integer:: nb_element_root
	type(Element), TARGET, allocatable:: element_list_root(:)  ! Element list
	type(Element), TARGET, allocatable:: element_list(:)

	real(8):: EleDistortion	! element distortion. 1 for no distortion

!	Type RIGIDPLANE is used to describe a rigid plane (wall)
	type RigidPlane
		integer:: nDir	! The unit normal vector of the rigid plane (¡À1,¡À2,¡À3)
						!     Must point towards the impacting body
		real(8):: coor  ! x(¡À1), y(¡À2) or z(¡À3) coordinate of a point on the plane
	end type RigidPlane

	type(RigidPlane) :: plane = RigidPlane(0,0)

	integer:: nb_comp = 0					! number of components
	integer, allocatable:: CompLen(:)		! nb_comp
	integer, allocatable:: CompMember(:,:)	! nb_comp * nb_node

end module ElementData
