recursive subroutine UpdateTree(T, Ttrain, Zvec, cv, sV, sN, exist)
    use precision
    use TypeDef
    implicit none
    ! InferTrainTree: Infers the search tree for a training image
    
    ! Input
    integer, intent(in), dimension(:)       :: Zvec
    integer, intent(in)                     :: cv, sV, sN
    
    ! Input/Output
    type(streenode), intent(inout)          :: T
    type(streenode), intent(inout)          :: Ttrain
    logical, intent(inout)                  :: exist
    
    ! Local variables
    integer                                 :: nZ, mZ, pZ, cvnext, v
   
    ! Current voxel is the vth neighboring voxel
    v = Ttrain%depth+1
       
    ! Value of current voxel
    cvnext = Zvec(v) + 1             

    ! If the current pattern is found in the training image
    if (associated(Ttrain%next)) then

        ! Update the node corresponding to the current central value:
        T%next(cvnext)%repl(cv) = T%next(cvnext)%repl(cv) + 1.0_dp

        ! If locations need still be visited in the data template, visit them
        ! and update the search tree accordingly:  
        if (v < sN) then
            call UpdateTree(T%next(cvnext), Ttrain%next(cvnext), Zvec, cv, sV, sN, exist) 
        else
            
            ! Otherwise we made it to the end and the voxel has an existing pattern
            exist = .true.
        end if        
    end if   
end subroutine UpdateTree             