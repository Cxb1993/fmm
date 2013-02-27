recursive subroutine getCPDF(Ttrain, Zvec, sV, sN, cpdf)
    use utils
    use precision
    use TypeDef
    implicit none
    
    ! Input
    type(streenode), intent(in)                         :: Ttrain
    integer, intent(in), dimension(:)                   :: Zvec
    integer, intent(in)                                 :: sV, sN
    
    ! Output
    real(dp), intent(inout), dimension(sV+1)            :: cpdf
    
    ! Local variables
    integer                                             :: cv, c, v
    
    v = Ttrain%depth + 1
    
    ! If more known indexes
    if (v <= sN) then
        
        ! Value of current voxel in the pattern Zvec
        cv = Zvec(v) + 1
        
        ! If present in the training image
        if (associated(Ttrain%next)) then
            
            ! If z is known
            if (cv > 0) then        
                call getCPDF(Ttrain%next(cv), Zvec, sV, sN, cpdf)
            
            ! If z is unknown it can belong to all sV+1 categories     
            else
                do c = 1, sV+1
                    call getCPDF(Ttrain%next(c), Zvec, sV, sN, cpdf)
                end do
            end if
        end if
        
    ! If only unknown indexes are left    
    else
        cpdf = cpdf + Ttrain%repl
    end if
end subroutine getCPDF

