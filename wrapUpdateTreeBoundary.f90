subroutine wrapUpdateTreeBoundary(T, Ttrain, cpdf, Zvec, cv, sV, sN, exist)
    use precision
    use TypeDef
    use Interfaces, only : UpdateTreeBoundary
    implicit none
    
    ! Input
    type(streenode), intent(in)             :: Ttrain
    real(dp), intent(in)                    :: cpdf
    integer, intent(in), dimension(:)       :: Zvec
    integer, intent(in)                     :: cv, sV, sN
    
    ! Input/Output
    type(streenode), intent(inout)          :: T
    
    ! Output
    logical, intent(out)                    :: exist 
    
    ! Local variables
    type(streenode)                         :: Ttmp
    real(dp)                                :: cpdftmp
    integer, dimension(:), allocatable      :: Zvectmp
    
    Ttmp = Ttrain  
    cpdftmp = cpdf
    allocate(Zvectmp(size(Zvec,1)), source = Zvec)
    exist = .false.
    
    call UpdateTreeBoundary(T, Ttmp, cpdftmp, Zvectmp, cv, sV, sN, exist)
    
end subroutine wrapUpdateTreeBoundary
