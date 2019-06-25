module MPIVariable
    integer:: procs, pid,  elementSharedNumber, elementInteriorNumber, commonMax, commonNumber
    integer, allocatable:: commonNum(:), curveNum(:), curveSequence(:)
    integer, allocatable:: Neighbour(:, :), nodeType(:), nodeGlobal(:, :), requestArray(:), elementShared(:), elementInterior(:), statusArray(:)
    real(8), allocatable:: test(:, :), DtkTest(:), buf_Acc1(:, :, :), buf_Acc2(:, :, :)

    contains

    subroutine updateInformation(tf)
    use ParticleData
    use TimeFunctionData
    use Simulation
    use ElementData
    use mpi
    implicit none

        real(8), intent(in):: tf
        integer:: i, j, k, ierr, status(MPI_STATUS_SIZE), requestNumber
        ! real(8), allocatable:: buf_Acc1(:, :, :), buf_Acc2(:, :, :)

        ! allocate(buf_Acc1(nDim, commonMax, procs))
        ! allocate(buf_Acc2(nDim, commonMax, procs))
        requestNumber = 0

        buf_Acc1 = 0
        buf_Acc2 = 0
        do i = 1, procs
            if (commonNum(i) .gt. 0) then
                requestNumber = requestNumber + 1
                do j = 1, commonNum(i)
                    buf_Acc1(:, j, i) = Acc(:, Neighbour(j, i))
                end do
                call MPI_Isend(buf_Acc1(1, 1, i), commonNum(i) * 3, MPI_REAL8, i - 1, 10, MPI_COMM_WORLD, requestArray(requestNumber), ierr)
                ! call MPI_Waitall(1, request(requestNumber), status, ierr)
                ! call MPI_Send(buf_Acc1(1, 1, i), commonNum(i) * 3, MPI_REAL8, i - 1, 10, MPI_COMM_WORLD, ierr)
            end if
        end do

        do i = 1, procs
            if (commonNum(i) .gt. 0) then
                ! print *, pid, requestNumber
                requestNumber = requestNumber + 1
                call MPI_Irecv(buf_Acc2(1, 1, i), commonNum(i) * 3, MPI_REAL8, i - 1, 10, MPI_COMM_WORLD, requestArray(requestNumber), ierr)
                ! call MPI_Waitall(1, request(requestNumber), status, ierr)
                ! call MPI_Recv(buf_Acc2(1, 1, i), commonNum(i) * 3, MPI_REAL8, i - 1, 10, MPI_COMM_WORLD, status, ierr)
            end if
            ! do j = 1, commonNum(i)
            !     Acc(:, Neighbour(j, i)) = buf_Acc(:, j) + Acc(:, Neighbour(j, i))
            ! end do
        end do

        ! print *, pid, 'before'
        ! call MPI_Waitall(2 * commonNumber, requestArray, status, ierr)

        requestNumber = commonNumber
        do i = 1, procs
            if (commonNum(i) .gt. 0) then
                requestNumber = requestNumber + 1
                call MPI_Wait(requestArray(requestNumber), status, ierr)
                do j = 1, commonNum(i)
                    Acc(:, Neighbour(j, i)) = buf_Acc2(:, j, i) + Acc(:, Neighbour(j, i))
                end do
            end if
        end do

        ! deallocate(buf_Acc1)
        ! deallocate(buf_Acc2)
        do i = 1, nb_node
            Acc(:, i) = Acc(:, i) + Fp(:, i) * tf
        end do
    end subroutine updateInformation

    subroutine updateDTk()
    use ParticleData
    use TimeFunctionData
    use Simulation
    use ElementData
    use mpi
    implicit none

        integer:: status, i, j, ierr
        real(8):: buff1(4), buff2(4)
        real(8):: buf_DTK, buf_EleDist, buf_EngKine, buf_Strain

        buff1(1) = DTk1
        buff1(2) = EleDistortion
        buff1(3) = EngKinetic
        buff1(4) = EngStrain
        ! call MPI_Reduce(DTk1, buf_DTK, 1, MPI_REAL8, MPI_MIN, 0, MPI_COMM_WORLD, ierr)
        ! call MPI_Reduce(EleDistortion, buf_EleDist, 1, MPI_REAL8, MPI_MIN, 0, MPI_COMM_WORLD, ierr)
        call MPI_Reduce(buff1, buff2, 2, MPI_REAL8, MPI_MIN, 0, MPI_COMM_WORLD, ierr)
        call MPI_Reduce(buff1(3), buff2(3), 2, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    
        ! call MPI_Reduce(EngKinetic, buf_EngKine, 1, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
        ! call MPI_Reduce(EngStrain, buf_Strain, 1, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
        call MPI_Bcast(buff2, 4, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)

        DTk1 = buff2(1)
        EleDistortion = buff2(2)
        EngKinetic = buff2(3)
        EngStrain = buff2(4)
        
        ! do i = 1, procs
        !     if (i - 1 .ne. pid) then
        !         buff(1) = DTk1
        !         buff(2) = EleDistortion
        !         buff(3) = EngKinetic
        !         buff(4) = EngStrain
        !         call MPI_Send(buff, 4, MPI_REAL8, i - 1, 11, MPI_COMM_WORLD, ierr)
        !     end if
        ! end do

        ! do i = 1, procs
        !     if (i - 1 .ne. pid) then
        !         call MPI_Recv(buff, 4, MPI_REAL8, i - 1, 11, MPI_COMM_WORLD, status, ierr)
        !         EngKinetic = EngKinetic + buff(3)
        !         EngStrain = EngStrain + buff(4)
        !         if (DTk1 .ge. buff(1)) then
        !             DTk1 = buff(1)
        !         end if
        !         if (EleDistortion .ge. buff(2)) then
        !             EleDistortion = buff(2)
        !         end if
        !     end if
        ! end do
    end subroutine updateDTk

end module MPIVariable