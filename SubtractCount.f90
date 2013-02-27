recursive subroutine SubtractCount(T, sN, level, Zvec, zcen)
    use precision
    use TypeDef
    implicit none
    
    ! Input
    integer, intent(in)                 :: sN, level, zcen
    integer, intent(in), dimension(:)   :: Zvec
    
    ! Input/Output
    type(streenode), intent(inout)      :: T
    
    ! Local variables
    integer                             :: d, c
         
    ! If represented in the tree
    if (associated(T%next)) then
            
        ! Current depth
        d = T%depth+1
                
        ! Value of current neighbor of the neighbor
        c = Zvec(d) + 1
        
        ! Go to the level of change
        if (d < level) then       
            
            call SubtractCount(T%next(c), sN, level, Zvec, zcen)
        else

            ! Add the count
            T%next(c)%repl(zcen+1) = T%next(c)%repl(zcen+1) - 1.0_dp
            
            ! If more leaves needs to be updated
            if (level < sN)then     
                call SubtractCount(T%next(c), sN, level+1, Zvec, zcen)
            end if
        end if
    end if
    
end subroutine SubtractCount