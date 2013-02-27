recursive subroutine CopyTree(Told, sV, sN, Tnew)
    use TypeDef
    use Interfaces, only : ExtendTree
    use precision
    implicit none
    
    ! Input
    type(streenode), intent(in)     :: Told
    integer                         :: sV, sN    
    
    ! Output
    type(streenode), intent(inout)  :: Tnew
    
    ! Local variables
    integer                         :: d, c
        
    if (Told%depth == 0) Tnew%repl = Told%repl
       
    ! Go deeper
    if (associated(Told%next)) then
        
        ! Extend the new tree
        if (.not.associated(Tnew%next)) call ExtendTree(sV, Tnew)
        
        ! Fill in the values
        do c = 1, sV+1
            Tnew%next(c)%repl(1:sV+1) = Told%next(c)%repl(1:sV+1)
        end do
        
        ! Current depth
        d = Told%depth + 1
        
        if (d < sN) then
            do c = 1, sV+1
                call CopyTree(Told%next(c), sV, sN, Tnew%next(c))        
            end do
        end if
    end if
    
end subroutine CopyTree