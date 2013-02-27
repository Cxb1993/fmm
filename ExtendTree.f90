subroutine ExtendTree(sV, T) 
    use precision
    use TypeDef
    implicit none
    ! ExtendTree: Extends the search tree 'T' by adding a new node
    !
    ! Input
    integer, intent(in)                 :: sV
    
    ! Input/Output
    type(streenode), intent(inout)      :: T
    
    ! Local variables
    integer                             :: c
    
    ! Create new node
    allocate(T%next(sV+1))
       
    ! Initialize
    do c = 1, sV+1
        
        allocate(T%next(c)%repl(sV+1), source = 0.0_dp)
        T%next(c)%depth = T%depth + 1
        nullify(T%next(c)%next)
    end do
    
end subroutine ExtendTree
    
    
    
    

