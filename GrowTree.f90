subroutine GrowTree(Z, zold, Zex, i, j, k, nodes, sV, sN, Ttrain, T)
    use TypeDef
    use utils
    use Interfaces, only : getNeighborhood, CenterCount, AddCount, SubtractCount, wrapUpdateTreeBoundary
    implicit none
    
    ! Input
    integer, intent(in), dimension(:,:,:)   :: Z
    integer, intent(in)                     :: zold, i, j, k, sV, sN
    integer, intent(in), dimension(:,:)     :: nodes
    type(streenode), intent(in)             :: Ttrain
    
    ! Input/Output
    logical, intent(inout), allocatable, dimension(:,:,:) :: Zex
    type(streenode), intent(inout)          :: T
    
    ! Local variables
    integer                                 :: nZ, mZ, pZ, znew, zcen
    integer                                 :: ijk, ii, jj, kk, v
    integer, dimension(sN)                  :: Zvec
    logical                                 :: innervoxel, exist
    real(dp)                                :: cpdf
     
    ! Dimensions
    nZ = size(Z,1); mZ = size(Z,2); pZ = size(Z,3)
    
    ! New value
    znew = Z(i,j,k)
        
    ! UPDATE THE CONTRIBUTION FROM THE CENTER VOXEL
    T%repl(zold+1) = T%repl(zold+1) - 1 
    T%repl(znew+1) = T%repl(znew+1) + 1     
    
    ! Neighborhood
    call getNeighborhood(Z, i, j, k, nodes, sN, Zvec, innervoxel)       
    
    ! Update the tree 
    if (innervoxel) then
        exist = .false.
        call CenterCount(T, sN, znew, zold, Zvec, exist)    
        Zex(i,j,k) = exist
    else                
        ! Add new counts
        cpdf = 1.0_dp 
        call wrapUpdateTreeBoundary(T, Ttrain, cpdf, Zvec, znew+1, sV, sN, exist)
        Zex(i,j,k) = exist

        ! Subtract the old count
        cpdf = -1.0_dp
        call wrapUpdateTreeBoundary(T, Ttrain, cpdf, Zvec, zold+1, sV, sN, exist)
    end if

    ! UPDATE THE CONTRIBUTIONS FROM IT BEING A NEIGHBORING
    
    ! For each potential neighbor position
    do v = 1, sN

        ! Coordinates of the center for this pattern
        ii = i - nodes(v,1)
        jj = j - nodes(v,2)
        kk = k - nodes(v,3)
        
        ! If center exists update the tree
        if (ii >= 1 .and. ii <= nZ .and. &
            jj >= 1 .and. jj <= mZ .and. &
            kk >= 1 .and. kk <= pZ) then
            
            ! Color of the current neighbor
            zcen = Z(ii,jj,kk)

            ! Find their neighborhood
            call getNeighborhood(Z, ii, jj, kk, nodes, sN, Zvec, innervoxel)
                        
            ! If neighbor was an inner voxel
            if (innervoxel) then
                                
                exist = .false.
                
                ! Add new counts
                call AddCount(T, sN, v, Zvec, zcen, exist)
                Zex(ii,jj,kk) = exist

                ! Subtract the old count
                Zvec(v) = zold       
                call SubtractCount(T, sN, v, Zvec, zcen)
            else
                
                ! Add new counts
                cpdf = 1.0_dp 
                call wrapUpdateTreeBoundary(T, Ttrain, cpdf, Zvec, zcen+1, sV, sN, exist)
                Zex(ii,jj,kk) = exist
                
                ! Subtract the old count
                cpdf = -1.0_dp
                Zvec(v) = zold
                call wrapUpdateTreeBoundary(T, Ttrain, cpdf, Zvec, zcen+1, sV, sN, exist)
            end if
        end if    
    end do
end subroutine