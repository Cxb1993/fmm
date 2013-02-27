subroutine getNewIt(Zex, nodes, sN, i, j, k)
    use precision
    use utils
    implicit none

    ! Input
    logical, intent(in), dimension(:,:,:)   :: Zex
    integer, intent(in), dimension(:,:)     :: nodes
    integer, intent(in)                     :: sN
    
    ! Output
    integer, intent(out)                    :: i, j, k
        
    ! Local variables
    integer                                 :: nZ, mZ, pZ, ii, jj, kk
    integer                                 :: exs, total, v
    real(dp)                                :: rand
    logical                                 :: newIt
    
    ! Dimensions
    nZ = size(Zex,1); mZ = size(Zex,2); pZ = size(Zex,3)       
        
    ! Iterate till new iteration is found
    newIt = .false. 
    do while (.not.newIt)
    
        ! Choose random voxel
        call random_number(rand)
        i = floor(1.0_dp + (nZ * rand))
        call random_number(rand)
        j = floor(1.0_dp + (mZ * rand))
        call random_number(rand)
        k = floor(1.0_dp + (pZ * rand))

        ! Check feasible patterns in the neighborhood
        exs = 0
        if (Zex(i,j,k)) exs = exs + 1    
        total = sN + 2
        do v = 1, sN
            ii = i + nodes(v,1)
            jj = j + nodes(v,2)
            kk = k + nodes(v,3)

            ! If neighboring voxel exists
            if (ii >= 1 .and. ii <= nZ .and. &
                jj >= 1 .and. jj <= mZ .and. &
                kk >= 1 .and. kk <= pZ) then

                ! If pattern of neighboring voxel exists
                if (Zex(ii,jj,kk)) exs = exs + 1                
            else
                total = total - 1                
            end if
        end do
        
        ! Check if the iteration should be accepted
        call random_number(rand)
        if (rand > real(exs, kind = dp)/real(total, kind = dp)) newIt = .true.
    end do

end subroutine getNewIt