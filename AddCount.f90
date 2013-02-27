recursive subroutine AddCount(T, sN, level, Zvec, zcen, exist)
    use precision
    use TypeDef
    implicit none
    
    ! Input
    integer, intent(in)                 :: sN, level, zcen
    integer, intent(in), dimension(:)   :: Zvec
    
    ! Input/Output
    type(streenode), intent(inout)      :: T
    logical, intent(inout)              :: exist
    
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
            
            call AddCount(T%next(c), sN, level, Zvec, zcen, exist)
        else
            
            ! Add the count
            T%next(c)%repl(zcen+1) = T%next(c)%repl(zcen+1) + 1.0_dp

            ! If more leaves needs to be updated
            if (level < sN) then     
                 call AddCount(T%next(c), sN, level+1, Zvec, zcen, exist)
            else
                exist = .true.
            end if
        end if
    end if    
end subroutine AddCount