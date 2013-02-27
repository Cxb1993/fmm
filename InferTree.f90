subroutine InferTree(Z, nodes, Ttrain, sV, sN, T, Zex)
    use utils
    use precision
    use TypeDef
    use Interfaces, only : ShapeTree, getNeighborhood, wrapUpdateTree, wrapUpdateTreeBoundary
    implicit none
        
    ! Input
    integer, intent(in), dimension(:,:,:)   :: Z
    integer, intent(in), dimension(:,:)     :: nodes
    type(streenode), intent(in)             :: Ttrain
    integer, intent(in)                     :: sV, sN
    
    
    ! Output
    type(streenode), intent(out)            :: T
    logical, intent(out), allocatable, dimension(:,:,:) :: Zex
        
    ! Local variables
    integer                                 :: nZ, mZ, pZ, i, j, k, ii, jj, kk, cv
    real(dp)                                :: cpdf
    integer, dimension(sN)                  :: Zvec
    logical                                 :: innervoxel, exist
    
    ! Dimensions
    nZ = size(Z,1);  mZ = size(Z,2);  pZ = size(Z,3)
        
    ! Initialization
    allocate(Zex(nZ,mZ,pZ), source = .false.)
    T%depth = 0; allocate(T%repl(sV+1), source = 0.0_dp); nullify(T%next)  
    call ShapeTree(T, Ttrain, sV, sN)
    
    ! Loop over all voxels
    do k = 1, pZ
        do j = 1, mZ
            do i = 1, nZ    
                
                ! Current value of the voxel
                cv = Z(i,j,k) + 1
                                
                ! Add the count to the 0-p statistics
                T%repl(cv) = T%repl(cv) + 1
                
                ! Get neighborhood
                call getNeighborhood(Z, i, j, k, nodes, sN, Zvec, innervoxel)
                
                ! If an inner voxel
                if (innervoxel) then
                    
                    ! Update the search tree
                    call wrapUpdateTree(T, Ttrain, Zvec, cv, sV, sN, exist)
                    Zex(i,j,k) = exist
                    
                ! If boundary voxel    
                else
                    
                    ! Update the search tree
                    cpdf = 1.0_dp
                    call wrapUpdateTreeBoundary(T, Ttrain, cpdf, Zvec, cv, sV, sN, exist)
                    Zex(i,j,k) = exist
                end if
            end do
        end do
    end do
    
end subroutine InferTree
                
                

