recursive subroutine UpdateTrainTree(Ttrain, Zvec, cv, sV, sN)
    use precision
    use TypeDef
    implicit none
    ! InferTrainTree: Infers the search tree for a training image
    
    ! Input
    integer, intent(in), dimension(:)       :: Zvec
    integer, intent(in)                     :: cv, sV, sN
    
    ! Input/Output
    type(streenode), intent(inout)          :: Ttrain
    
    ! Local variables
    integer                                 :: nZ, mZ, pZ, cvnext, v
    
    ! Current voxel is the vth neighboring voxel
    v = Ttrain%depth+1
     
    ! Value of current voxel
    cvnext = Zvec(v) + 1             

    ! If the search tree node corresponding to the current conditioning
    ! configuration does not exit, create it:       
    if (.not.associated(Ttrain%next)) call ExtendTree(sV, Ttrain)

    ! Update the node corresponding to the current central value: 
    Ttrain%next(cvnext)%repl(cv) = Ttrain%next(cvnext)%repl(cv) + 1.0_dp

    ! If locations need still be visited in the data template, visit them
    ! and update the search tree accordingly:  
    if (v < sN) call UpdateTrainTree(Ttrain%next(cvnext), Zvec, cv, sV, sN) 
                
end subroutine UpdateTrainTree
    
    
                
                

