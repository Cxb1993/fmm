recursive subroutine DeallocateTree(sV, T) 
    use TypeDef
    implicit none
    ! DeallocateTree: Deallocates recursively the search tree 'stree'
    ! starting from leaves to root.
    !
    ! Input
    integer, intent(in)                 :: sV
    
    ! Input/Output
    type(streenode), intent(inout)      :: T
    
    ! Local variables
    integer                             :: c
   
    ! If more nodes to be deallocated
    if (associated(T%next)) then
        
        do c = 1, sV+1
            
            if (associated(T%next(c)%next)) then
                call DeallocateTree(sV, T%next(c))
            else
                deallocate(T%next(c)%repl)
            end if
        end do
        
        deallocate(T%next)
    end if
    
    deallocate(T%repl)
           
end subroutine DeallocateTree
    
    
    
    



