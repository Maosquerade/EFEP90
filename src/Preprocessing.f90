module PreProcessing
use DataIn
use MPIVariable
use TimeFunctionData
use mpi
integer, allocatable:: partE(:, :), partN(:, :), partSize(:), nodeCounter(:), nodeLocal(:), nodeType_root(:), nodeShared(:)

contains 
subroutine PreData()
    use iso_c_binding
    ! use Calling_METIS

    implicit none
    
    integer:: i, j, k, m, node, numPerProc, ierr, status, count
    integer, allocatable:: commNum(:), bufNeighbor(:), buf_NodeType(:)
    type(Element), allocatable:: buf_Element(:)
    real(8), allocatable:: buf_Pos(:, :), buf_Vel(:, :), buf_Mp(:), buf_Fp(:, :)
    type(History), allocatable:: buf_History(:)
    integer:: ioMetis = 21, ie, sharedFlag
    integer, allocatable:: epart(:)
    real:: start, end

    type(Element):: el
    character*20:: str

    allocate(commonNum(procs))

    call MPI_Bcast(LenTMFunction, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(Jaum, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(EngStrain, 1, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(nb_comp, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(nb_mat, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(EndTime, 1, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(OutTimeStep, 1, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(RepTimeStep, 1, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(DTScale, 1, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(LenTMFunction, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    if (LenTMFunction .gt. 0) then
        if (pid .ne. 0) then
            allocate(TMFunction(LenTMFunction))
        end if
        call MPI_Bcast(TMFunction, 2 * LenTMFunction, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
    end if
    call MPI_Bcast(HourGlass%method, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(HourGlass%Qhg, 1, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(plane%nDir, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(plane%coor, 1, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)

    if (pid .eq. 0) then
        allocate(partE(procs, nb_element_root))
        allocate(partN(procs, nb_node_root))
        allocate(partSize(procs))
        allocate(nodeCounter(procs))
        allocate(nodeLocal(nb_node_root))
        allocate(nodeType_root(nb_node_root))
        allocate(nodeGlobal(nb_node_root, procs))

        allocate(buf_Pos(nDim, nb_node_root))
        allocate(buf_Vel(nDim, nb_node_root))
        allocate(buf_Mp(nb_node_root))
        allocate(buf_Fp(nDim, nb_node_root))
        allocate(buf_NodeType(nb_node_root))

        allocate(buf_Element(nb_element_root))
        allocate(buf_History(nb_element_root))

        allocate(bufNeighbor(nb_node_root))
        allocate(commNum(procs))

        allocate(nodeShared(nb_node_root))
        allocate(elementShared(nb_element_root))
        allocate(elementInterior(nb_element_root))

        partSize = 0

        open(ioMetis, file='../make/element.mesh.epart.8', status='old')
        allocate(epart(nb_element_root))
        do i = 1, nb_element_root
            read(ioMetis, *) epart(i)
            epart(i) = epart(i) + 1     ! The default index is from 0, add 1 to adapt fortran array
        end do
        close(ioMetis)

        numPerProc = 1 + nb_element_root / procs
        do i = 1, nb_element_root
            j = epart(i)
            partSize(j) = partSize(j) + 1
            partE(j, partSize(j)) = i
        end do

        nodeType_root = 0
        nodeShared = 0

        do i = 1, procs
            nodeCounter(i) = 0
            do j = 1, nb_node_root
                nodeLocal(j) = -1
                partN(i, j) = 0
            end do
            do j = 1, partSize(i)
                sharedFlag = 0
                do k = 1, 8
                    node = element_list_root(partE(i, j))%nNode(k)
                    if (nodeLocal(node) .eq. -1) then
                        nodeCounter(i) = nodeCounter(i) + 1
                        nodeLocal(node) = nodeCounter(i)
                    end if
                    if (partN(i, node) == 0) then
                        nodeType_root(node) = i
                        partN(i, node) = 1
                        nodeShared(node) = nodeShared(node) + 1
                    end if
                end do
            end do
        end do

        do i = 2, procs
            do j = 1, nb_node_root
                nodeLocal(j) = -1
            end do

            call MPI_Send(nodeCounter(i), 1, MPI_INTEGER, i - 1, i - 1, MPI_COMM_WORLD, ierr)
            call MPI_Send(partSize(i), 1, MPI_INTEGER, i - 1, i - 1 + procs, MPI_COMM_WORLD, ierr)

            nodeCounter(i) = 0
            buf_NodeType = 0

            do j = 1, partSize(i)
                do k = 1, 8
                    node = element_list_root(partE(i, j))%nNode(k)
                    if (nodeLocal(node) .eq. -1) then
                        nodeCounter(i) = nodeCounter(i) + 1
                        nodeLocal(node) = nodeCounter(i) 
                        buf_Pos(:, nodeCounter(i)) = Pos_root(:, node)
                        buf_Vel(:, nodeCounter(i)) = Vel_root(:, node)
                        buf_Mp(nodeCounter(i)) = Mp_root(node)
                        buf_Fp(:, nodeCounter(i)) = Fp_root(:, node)
                        if (nodeType_root(node) .eq. i) then
                            buf_NodeType(nodeCounter(i)) = 1
                        end if
                        nodeGlobal(nodeCounter(i), i) = node
                    end if
                end do
            end do

            call MPI_Send(buf_Pos, nodeCounter(i) * 3, MPI_REAL8, i - 1, i - 1 + procs * 2, MPI_COMM_WORLD, ierr)
            call MPI_Send(buf_Vel, nodeCounter(i) * 3, MPI_REAL8, i - 1, i - 1 + procs * 8, MPI_COMM_WORLD, ierr)
            call MPI_Send(buf_Mp, nodeCounter(i), MPI_REAL8, i - 1, i - 1 + procs*9, MPI_COMM_WORLD, ierr)
            call MPI_Send(buf_Fp, nodeCounter(i) * 3, MPI_REAL8, i - 1, i - 1 + procs * 10, MPI_COMM_WORLD, ierr)
            call MPI_Send(buf_NodeType, nodeCounter(i), MPI_INTEGER, i - 1, i - 1 + procs * 20, MPI_COMM_WORLD, ierr)

            count = 0
            do j = 1, partSize(i)
                do k = 1, 8
                    node = element_list_root(partE(i, j))%nNode(k)
                    buf_Element(j)%nNode(k) = nodeLocal(node)
                end do
                buf_Element(j)%mat = element_list_root(partE(i, j))%mat
                buf_Element(j)%nHistory = j
                buf_History(j) = history_list_root(partE(i, j))
                if (buf_History(j)%VOL .le. 1e-6) then
                    count = count + 1
                end if
            end do

            call MPI_Send(buf_Element, partSize(i) * 10, MPI_INTEGER, i - 1, i - 1 + procs * 3, MPI_COMM_WORLD, ierr)
            call MPI_Send(buf_History, partSize(i) * 11, MPI_REAL8, i - 1, i - 1 + procs*7, MPI_COMM_WORLD, ierr)

            do k = 1, procs
                commNum(k) = 0
                if (i .ne. k) then
                    do j = 1, nb_node_root
                        if (partN(i, j) > 0 .and. partN(k, j) > 0) then
                            commNum(k) = commNum(k) + 1
                        end if
                    end do
                end if
            end do

            call MPI_Send(commNum, procs, MPI_INTEGER, i - 1, i - 1 + procs*4, MPI_COMM_WORLD, ierr)

            do k = 1, procs
                if (commNum(k) > 0) then
                    commNum(k) = 0
                    do j = 1, nb_node_root
                        if (partN(i, j) > 0 .and. partN(k, j) > 0) then
                            commNum(k) = commNum(k) + 1
                            bufNeighbor(commNum(k)) = nodeLocal(j)
                        end if
                    end do
                end if
                call MPI_Send(bufNeighbor, commNum(k), MPI_INTEGER, i - 1, i - 1 + k + procs*5, MPI_COMM_WORLD, ierr)
            end do

            elementSharedNumber = 0
            elementInteriorNumber = 0
            
            do j = 1, partSize(i)
                sharedFlag = 0
                do k = 1, 8
                    node = element_list_root(partE(i, j))%nNode(k)
                    do m = 1, procs
                        if (partN(m, node) > 0 .and. m .ne. i) then
                            sharedFlag = 1
                        end if
                    end do
                end do
                if (sharedFlag .eq. 1) then
                    elementSharedNumber = elementSharedNumber + 1
                    elementShared(elementSharedNumber) = j
                else 
                    elementInteriorNumber = elementInteriorNumber + 1
                    elementInterior(elementInteriorNumber) = j
                end if
            end do

            call MPI_Send(elementSharedNumber, 1, MPI_INTEGER, i - 1, i - 1 + procs*30, MPI_COMM_WORLD, ierr)
            call MPI_Send(elementShared, elementSharedNumber, MPI_INTEGER, i - 1, i - 1 + procs*31, MPI_COMM_WORLD, ierr)
            call MPI_Send(elementInteriorNumber, 1, MPI_INTEGER, i - 1, i - 1 + procs*32, MPI_COMM_WORLD, ierr)
            call MPI_Send(elementInterior, elementInteriorNumber, MPI_INTEGER, i - 1, i - 1 + procs*33, MPI_COMM_WORLD, ierr)
        end do

!   allocate root processor data
        do j = 1, nb_node_root
            nodeLocal(j) = -1
        end do

        nb_node = nodeCounter(1)
        nb_element = partSize(1)
        allocate(Pos(nDim, nb_node))
        allocate(Vel(nDim, nb_node))
        allocate(Mp(nb_node))
        allocate(Fp(nDim, nb_node))
        allocate(nodeType(nb_node))
        allocate(element_list(nb_element))
        allocate(history_list(nb_element))

        nodeCounter(1) = 0
        nodeType = 0
        do j = 1, partSize(1)
            do k = 1, 8
                node = element_list_root(partE(1, j))%nNode(k)
                if (nodeLocal(node) .eq. -1) then
                    nodeCounter(1) = nodeCounter(1) + 1
                    nodeLocal(node) = nodeCounter(1) 
                    Pos(:, nodeCounter(1)) = Pos_root(:, node)
                    Vel(:, nodeCounter(1)) = Vel_root(:, node)
                    Mp(nodeCounter(1)) = Mp_root(node)
                    Fp(:, nodeCounter(1)) = Fp_root(:, node)
                    if (nodeType_root(node) .eq. 1) then
                        nodeType(nodeCounter(1)) = 1
                    end if
                    nodeGlobal(nodeCounter(1), 1) = node
                end if
            end do
        end do

        do j = 1, partSize(1)
            do k = 1, 8
                node = element_list_root(partE(1, j))%nNode(k)
                element_list(j)%nNode(k) = nodeLocal(node)
            end do
            element_list(j)%mat = element_list_root(partE(1, j))%mat
            element_list(j)%nHistory = j
            history_list(j) = history_list_root(partE(1, j))
        end do

        commonMax = 0
        do k = 1, procs
            commonNum(k) = 0
            if (k .ne. 1) then
                do j = 1, nb_node_root
                    if (partN(1, j) > 0 .and. partN(k, j) > 0) then
                        commonNum(k) = commonNum(k) + 1
                    end if
                end do
            end if
            if (commonMax < commonNum(k)) then
                commonMax = commonNum(k)
            end if
        end do

        allocate(Neighbour(commonMax, procs))
        
        do k = 1, procs
            if (commonNum(k) > 0) then
                commonNum(k) = 0
                do j = 1, nb_node_root
                    if (partN(1, j) > 0 .and. partN(k, j) > 0) then
                        commonNum(k) = commonNum(k) + 1
                        Neighbour(commonNum(k), k) = nodeLocal(j)
                    end if
                end do
            end if
        end do

        do i = 1, nb_node_root
            if (partN(1, i) > 0 .and. partN(2, i) > 0) then
                print *, 'The first common node between 1 and 2 is', i
                exit
            end if
        end do

        elementSharedNumber = 0
        elementInteriorNumber = 0
            
        do j = 1, partSize(1)
            sharedFlag = 0
            do k = 1, 8
                node = element_list_root(partE(1, j))%nNode(k)
                do m = 1, procs
                    if (partN(m, node) > 0 .and. m .ne. 1) then
                        sharedFlag = 1
                    end if
                end do
            end do
            if (sharedFlag .eq. 1) then
                elementSharedNumber = elementSharedNumber + 1
                elementShared(elementSharedNumber) = j
            else 
                elementInteriorNumber = elementInteriorNumber + 1
                elementInterior(elementInteriorNumber) = j
            end if
        end do
    else 
        call MPI_Recv(nb_node, 1, MPI_INTEGER, 0, pid, MPI_COMM_WORLD, status, ierr)
        call MPI_Recv(nb_element, 1, MPI_INTEGER, 0, pid + procs, MPI_COMM_WORLD, status, ierr)
        allocate(Pos(nDim, nb_node))
        allocate(Vel(nDim, nb_node))
        allocate(Mp(nb_node))
        allocate(Fp(nDim, nb_node))
        allocate(nodeType(nb_node))
        allocate(element_list(nb_element))
        allocate(history_list(nb_element))
        allocate(elementShared(nb_element))
        allocate(elementInterior(nb_element))
        call MPI_Recv(Pos, nb_node * 3, MPI_REAL8, 0, pid + procs * 2, MPI_COMM_WORLD, status, ierr)
        call MPI_Recv(Vel, nb_node * 3, MPI_REAL8, 0, pid + procs * 8, MPI_COMM_WORLD, status, ierr)
        call MPI_Recv(Mp, nb_node, MPI_REAL8, 0, pid + procs * 9, MPI_COMM_WORLD, status, ierr)
        call MPI_Recv(Fp, nb_node * 3, MPI_REAL8, 0, pid + procs * 10, MPI_COMM_WORLD, status, ierr)
        call MPI_Recv(nodeType, nb_node, MPI_INTEGER, 0, pid + procs * 20, MPI_COMM_WORLD, status, ierr)
        call MPI_Recv(element_list, nb_element * 10, MPI_INTEGER, 0, pid + procs*3, MPI_COMM_WORLD, status, ierr)
        call MPI_Recv(history_list, nb_element * 11, MPI_REAL8, 0, pid + procs*7, MPI_COMM_WORLD, status, ierr)
        call MPI_Recv(commonNum, procs, MPI_INTEGER, 0, pid + procs*4, MPI_COMM_WORLD, status, ierr)
        commonMax = 0
        do i = 1, procs
            if (commonMax < commonNum(i)) then
                commonMax = commonNum(i)
            end if
        end do
        allocate(Neighbour(commonMax, procs))
        do i = 1, procs
            call MPI_Recv(Neighbour(1, i), commonNum(i), MPI_INTEGER, 0, pid + i + procs*5, MPI_COMM_WORLD, status, ierr)
        end do

        call MPI_Recv(elementSharedNumber, 1, MPI_INTEGER, 0, pid + procs*30, MPI_COMM_WORLD, status, ierr)
        call MPI_Recv(elementShared, elementSharedNumber, MPI_INTEGER, 0, pid + procs*31, MPI_COMM_WORLD, status, ierr)
        call MPI_Recv(elementInteriorNumber, 1, MPI_INTEGER, 0, pid + procs*32, MPI_COMM_WORLD, status, ierr)
        call MPI_Recv(elementInterior, elementInteriorNumber, MPI_INTEGER, 0, pid + procs*33, MPI_COMM_WORLD, status, ierr)

        allocate(Acc(nDim, nb_node))
        allocate(mat_list(nb_mat))
    end if

    do i = 1, nb_mat
        call MPI_Bcast(mat_list(i)%MatType, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        call MPI_Bcast(mat_list(i)%EosType, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        call MPI_Bcast(mat_list(i)%Density, 1, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
        call MPI_Bcast(mat_list(i)%Young, 1, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
        call MPI_Bcast(mat_list(i)%Poisson, 1, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
        call MPI_Bcast(mat_list(i)%WaveSpeed, 1, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
        call MPI_Bcast(mat_list(i)%Yield0, 1, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
        call MPI_Bcast(mat_list(i)%TangMod, 1, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
        call MPI_Bcast(mat_list(i)%B_jnc, 1, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
        call MPI_Bcast(mat_list(i)%n_jnc, 1, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
        call MPI_Bcast(mat_list(i)%C_jnc, 1, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
    end do
    j = 0
    do i = 1, nb_element
        if (abs(history_list(i)%VOL) .le. 1e-6) then
            j = j + 1
        end if
    end do

    commonNumber = 0
    do i = 1, procs
        if (commonNum(i) > 0) then
            commonNumber = commonNumber + 1
        end if
    end do
    allocate(requestArray(2 * procs))
    allocate(statusArray(2 * procs))
    allocate(buf_Acc1(nDim, commonMax, procs))
    allocate(buf_Acc2(nDim, commonMax, procs))

    print *, pid, commonNumber, commonMax

end subroutine PreData
end module PreProcessing