! ---------------------------------------------------------------------
! -                                                                   -
! -  Result output procedures                                         -
! -                                                                   -
! ---------------------------------------------------------------------

module DataOut

	integer, parameter :: MaxCurves = 10	! Limited 10 curves
	integer:: nCurves = 0					! Number of curves defined
	integer:: CurveOption(MaxCurves)		! Output variable for each curve
	integer:: CurvePoint(MaxCurves)			! Particle number for curve output

	integer, parameter :: MaxAnim = 10		! Limited 10 animation variables
	integer:: nAnimate = 0					! Number of animation variables defined
	integer:: AnimateOption(MaxAnim)		! Animation variables

	real(8):: OutTimeStep = 0.0		! time interval for plot data output
	real(8):: RepTimeStep = 0.0		! time interval for status report

	integer, parameter:: nbname = 5
	character(4), parameter:: OutputName(nbname) = &
	(/'seqv', 'epef','engk', 'engs', 'ener'/)

	integer:: nRes = 0
	integer:: nPts = 0

	integer, allocatable :: nElem(:)  ! Number of elements connected with a particle

contains


	subroutine OutCurve(Time)
! ---------------------------------------------------------------------
! -                                                                   -
! -  Output result for plotting time-history curve using TecPlot      -
! -                                                                   -
! ---------------------------------------------------------------------
	use FFI
	use ElementData, only: Element, element_list
	use MaterialData, only: History, history_list
	use Simulation, only: Title
	implicit none

	type(Element),  POINTER :: el
	type(History),  POINTER :: his

	real(8):: i, Time

	nPts = nPts + 1

	if (nPts .eq. 1) then
		write(iow2,10) Title
		write(iow2,20) ('"', OutputName(CurveOption(i)), CurvePoint(i), '"', i=1,nCurves)
		write(iow2,*)
10		format('TITLE = "', a60, '"')
20		format('VARIABLES= "Time"', 50(A2, A4, ' at', i6, A1))
	endif

	write(iow2,"(d15.6)", advance='no') Time

	do i = 1, nCurves

		el => element_list(CurvePoint(i))
		his => history_list(el%nHistory)

		call OutVariable(iow2, CurveOption(i), el%mat, his)

	end do
	write(iow2,*)

	end subroutine OutCurve


	subroutine OutAnim(Time)
! ---------------------------------------------------------------------
! -                                                                   -
! -  Output result for post processing using TecPlot                  -
! -                                                                   -
! ---------------------------------------------------------------------
	use ParticleData
	use MaterialData
	use ElementData
	use FFI, only: iow1
	use Simulation
	implicit none

	type(Element),  POINTER :: el
	type(History),  POINTER :: his

	real(8):: Time
	integer p, nh, i, j

	nRes = nRes + 1

	if (nRes .eq. 1) then
		write(iow1,10) Title
		write(iow1,20)
		write(iow1,*)

		write(iow1, 30) Time, nb_node, nb_element
	else
		write(iow1, 40) Time, nb_node, nb_element
	endif

10	format('TITLE = "', a60, '"')
20	format('VARIABLES= "X"   "Y"   "Z"  ')
30	format('ZONE T="Time =',1p, E12.4,'"', ' F=FEPOINT N=', I10,   &
		   ' E=', I10, ' ET=BRICK C=CYAN')
40	format('ZONE T="Time =',1p, E12.4,'"', ' F=FEPOINT N=', I10,   &
		   ' E=', I10, ' ET=BRICK C=CYAN D=(FECONNECT)')

	do p = 1, nb_node
		write(iow1,"(1p, 3e12.4)") Pos_root(:,p)
	end do

	if (nRes .eq. 1) then
		write(iow1,"(8i7)") ((element_list_root(i)%nNode(j), j=1,8), i=1,nb_element)
	endif

	end subroutine OutAnim


	subroutine OutVariable(iow, opt, mat, his)
! ---------------------------------------------------------------------
! -                                                                   -
! -  Output the value of variable name[opt]                           -
! -                                                                   -
! -  Input                                                            -
! -     iow  - File unit                                              -
! -     opt  - Variable index in name list                            -
! -     pt   - Particle                                               -
! -     his  - History varialbe                                       -
! -                                                                   -
! ---------------------------------------------------------------------
	use MaterialData, only: History
	use Simulation
	implicit none

	integer, intent(in) :: iow, opt, mat
	type(History),  intent(in) :: his

	select case(opt)

		case(1) ! seqv
			write(iow,"(d15.6)", advance='no') his%seqv
		case(2) ! epeff
			write(iow,"(d15.6)", advance='no') his%epeff
		case(3) ! engk
			write(iow,"(d15.6)", advance='no') EngKinetic
		case(4) ! engs
			write(iow,"(d15.6)", advance='no') EngStrain
		case(5) ! ener
			write(iow,"(d15.6)", advance='no') EngStrain+EngKinetic
		case default
			stop "***Error*** Invalid output variable !"

	end select

	end subroutine OutVariable


	integer function SetResOption()
! ---------------------------------------------------------------------
! -                                                                   -
! -  Purpose                                                          -
! -     Read and set result output option                             -
! -                                                                   -
! ---------------------------------------------------------------------
	use FFI
	implicit none

	SetResOption = keyword(OutputName,nbname)

	if(SetResOption.le.0 .or. SetResOption.gt.nbname) then
		call ErrorMsg()
		stop 'Result option error'
	end if

	end function SetResOption

end module DataOut