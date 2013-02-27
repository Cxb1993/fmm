recursive subroutine ShapeTree(Ttest, Ttrain, sV, sN)
    use precision
    use TypeDef
    use Interfaces, only : ExtendTree
    implicit none

    ! Input
    type(streenode), intent(in)     :: Ttrain
    integer, intent(in)             :: sV, sN
    
    ! Input/Output
    type(streenode), intent(inout)  :: Ttest
    
    ! Local variables
    integer                 :: c
        
    ! If Ttrain has a branch but Ttest not    
    if (associated(Ttrain%next)) then

        ! Extend Ttest if necessary
        if (.not.associated(Ttest%next)) call ExtendTree(sV, Ttest)

        ! Check further branches
        do c = 1, sV+1
            call ShapeTree(Ttest%next(c), Ttrain%next(c), sV, sN)
        end do                
    end if

end subroutine ShapeTree