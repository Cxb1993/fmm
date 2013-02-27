subroutine InferTrainTree(ZTrain, nodes, sV, sN, Ttrain)
    use TypeDef
    use Interfaces, only : getNeighborhood, UpdateTrainTree
    implicit none
    ! InferTrainTree: Infers the search tree for a training image
    
    ! Input
    integer, intent(in), dimension(:,:,:)   :: ZTrain
    integer, intent(in), dimension(:,:)     :: nodes
    integer, intent(in)                     :: sV, sN
    
    ! Output
    type(streenode), intent(inout)          :: Ttrain
    
    ! Local variables
    integer                                 :: nZt, mZt, pZt, i, j, k, cv
    integer, dimension(sN)                  :: Zvec
    logical                                 :: innervoxel
    
    ! Dimensions
    nZt = size(ZTrain,1);  mZt = size(ZTrain,2);  pZt = size(ZTrain,3)
            
    ! Loop over all voxels 
    do k = 1, pZt
        do j = 1, mZt
            do i = 1, nZt
                
                ! Get neighborhood
                call getNeighborhood(Ztrain, i, j, k, nodes, sN, Zvec, innervoxel)
                
                ! Only inner voxels contribute
                if (innervoxel) then
                    
                    ! Current value of the voxel
                    cv = ZTrain(i,j,k) + 1
                
                    ! Add the count to the 0-p statistics
                    Ttrain%repl(cv) = Ttrain%repl(cv) + 1
                
                    ! Update the tree
                    call UpdateTrainTree(Ttrain, Zvec, cv, sV, sN)
                end if
            end do
        end do
    end do
    
end subroutine InferTrainTree           