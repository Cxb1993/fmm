recursive subroutine CenterCount(T, sN, znew, zold, Zvec, exist)
    use TypeDef
    implicit none
    ! CenterCount: Updates the tree wrt. the center voxel when it is flipped
    
    ! Input
    integer, intent(in)                 :: sN, znew, zold
    integer, intent(in), dimension(:)   :: Zvec
    
    ! Input/Output
    type(streenode), intent(inout)      :: T
    logical, intent(inout)              :: exist
    
    ! Local variables
    integer                             :: c, v
    
    ! Current voxel is the vth neighboring voxel
    v = T%depth + 1
       
    ! If represented in the tree
    if (associated(T%next)) then
    
        ! Color of current neighbor
        c = Zvec(v) + 1
         
        ! Add/subtract the count for new/old value
        T%next(c)%repl(znew+1) = T%next(c)%repl(znew+1) + 1
        T%next(c)%repl(zold+1) = T%next(c)%repl(zold+1) - 1
        
        ! If more leaves needs to be updated
        if (v < sN) then     
            call CenterCount(T%next(c), sN, znew, zold, Zvec, exist)
        else
            exist = .true.
        end if
    end if
    
end subroutine CenterCount