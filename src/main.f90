
! ----------------------------------------------------------------------
! -                          E F E P                                   -
! - Explicit Finite Element Program for high velocity impact analysis  -
! -                                                                    -        
! ----------------------------------------------------------------------

	program main
	use DataIn
	use DataOut
	use PreProcessing
	use FFI
	use mpi

	implicit none
	real:: time_begin, time_end
	integer:: ierr
!
	call MPI_INIT(ierr)
	call MPI_COMM_SIZE(MPI_COMM_WORLD, procs, ierr)
	call MPI_COMM_RANK(MPI_COMM_WORLD, pid, ierr)

!   use root processor to read data
	if (pid .eq. 0) then
		call cpu_time( time_begin )
!	
		call InputPara()
	end if
!
	call PreData()

	call MPI_Barrier(MPI_COMM_WORLD, ierr)

	if (pid .eq. 0) then
		call cpu_time( time_end )
		print *, '** Time for preprocessing is', time_end - time_begin, ' seconds'
	end if

	if (pid .eq. 0) then 
		call cpu_time( time_begin )
		print *, 'solving...'
	end if

	call Integration()	
	
	if (pid .eq. 0) then
		call cpu_time( time_end )
		print *, '** Time for computation is  ', time_end - time_begin, ' seconds'
		
		close(ior)
		close(iow1)
		close(iow2)
	
		close(iores)
	end if

	call MPI_FINALIZE(ierr)
	end program main

	subroutine integration()
! -------------------------------------------------------------------
! -                                                                 -
!   Integrate momentum equations
! -                                                                 -
! -------------------------------------------------------------------
	use Simulation
	use ElementData
	use ParticleData
	use DataOut
	use TimeFunctionData
	use MPIVariable
	use mpi
	use PreProcessing
	implicit none

	real(8) :: prt = 0.0
	real(8) :: plt = 0.0

	real(8) :: tf
	real(8), allocatable:: buf_Pos(:, :)
	integer:: ierr, i, j, status, node

	allocate(buf_Pos(nDim, nb_node_root))
	do while (CurrentTime .le. EndTime)

		DTk  = DTk1
		DTk1 = 1.0d6
!	Calculate kinetic energy
		call KineticE()

!	Time function at the current time step
		tf = CurrentTMFunction(CurrentTime)
		
		call FEForce(tf)

		! call updateInformation(tf)

		call UpdateFEGeometry()

		if (CurrentTime.ge.prt) then
			prt = prt + OutTimeStep

			if (pid .eq. 0) then
				do i = 2, procs
					call MPI_Recv(buf_Pos, nodeCounter(i) * 3, MPI_REAL8, i - 1, 21, MPI_COMM_WORLD, status, ierr)
					do j = 1, nodeCounter(i)
						Pos_root(:, nodeGlobal(j, i)) = buf_Pos(:, j)
					end do
				end do
				do j = 1, nodeCounter(1)
					Pos_root(:, nodeGlobal(j, 1)) = Pos(:, j)
				end do
			else
				call MPI_Send(Pos, nb_node * 3, MPI_REAL8, 0, 21, MPI_COMM_WORLD, ierr)
			end if
			if (pid .eq. 0) then
				call OutCurve(CurrentTime)
				call OutAnim(CurrentTime)
			end if
		end if
		
		if (CurrentTime.ge.plt) then
			plt = plt + RepTimeStep
			if (pid .eq. 0) then
			write(*,100) istep, CurrentTime, DTk1, EleDistortion, EngKinetic
100			format(1x, "Step=", i6, 1p, "  T=", e9.3, "  DT=", e9.3,  &
					   "  EDist=", e9.3, "  K.E.=", e9.3)
			end if
		end if

		istep = istep+1
		CurrentTime = CurrentTime + DTk1
	end do

	end subroutine integration
