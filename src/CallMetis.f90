module Calling_METIS
    use iso_c_binding
    implicit none
    integer(c_int) :: ne, nn 

    integer(c_int), dimension(:), allocatable :: eptr, eind

    type(c_ptr) :: vwgt, vsize 

    integer(c_int) :: nparts

    type(c_ptr) :: tpwgt

    integer(c_int), dimension(0:100) :: opts 

    integer(c_int) :: objval

    integer(c_int), dimension(:), allocatable :: epart, npart 
    
    interface
        subroutine METIS_PartMeshNodal( ne, nn, eptr, eind, vwgt, vsize, nparts, tpwgt, opts, objval, epart, npart) bind(c)
        use iso_c_binding
        implicit none

            integer (c_int), intent(in):: ne, nn
            integer (c_int), dimension(*), intent(in) :: eptr, eind

            type(c_ptr), value :: vwgt, vsize 
            integer (c_int), intent(in) :: nparts

            type(c_ptr), value :: tpwgt
            integer (c_int), dimension(0:100) :: opts
            integer (c_int), intent(out) :: objval
            integer (c_int), dimension(*), intent(out) :: epart
            integer (c_int), dimension(*), intent(out) :: npart
        end subroutine METIS_PartMeshNodal 
    end interface
end module