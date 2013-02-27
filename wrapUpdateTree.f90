subroutine wrapUpdateTree(T, Ttrain, Zvec, cv, sV, sN, exist)
    use TypeDef
    use Interfaces, only : UpdateTree
    implicit none
    ! InferTrainTree: Infers the search tree for a training image
    
    ! Input
    type(streenode), intent(in)             :: Ttrain
    integer, intent(in), dimension(:)       :: Zvec
    integer, intent(in)                     :: cv, sV, sN
    
    ! Input/Output
    type(streenode), intent(inout)          :: T
    
    ! Output
    logical, intent(out)                    :: exist
    
    ! Local variables
    type(streenode)                         :: Ttmp
    Ttmp = Ttrain
    exist = .false.
        
    ! Update the tree
    call UpdateTree(T, Ttmp, Zvec, cv, sV, sN, exist) 

end subroutine wrapUpdateTree             

