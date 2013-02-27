recursive subroutine UpdateTreeBoundary(T, Ttrain, cpdfold, Zvec, cv, sV, sN, exist)
    use utils
    use precision
    use TypeDef
    use Interfaces, only : getCPDF
    implicit none
    
    ! Input
    integer, intent(inout), dimension(:)    :: Zvec
    integer, intent(in)                     :: cv, sV, sN
    type(streenode), intent(in)             :: Ttrain
    
    ! Input/Output
    type(streenode), intent(inout)          :: T
    real(dp), intent(inout)                 :: cpdfold
    logical, intent(inout)                  :: exist
    
    ! Local variables
    integer                                 :: cvnext, c, v
    real(dp), dimension(sV+1)               :: cpdf, ratio
    real(dp)                                :: cpdfnew
    integer, dimension(sN)                  :: Ztmp
       
    ! Current voxel is the vth neighboring voxel
    v = Ttrain%depth+1
    
    ! Value of current voxel
    cvnext = Zvec(v) + 1
    
    ! If the current pattern is found in the training image
    if (associated(Ttrain%next)) then
                            
        ! If the current value is known
        if (cvnext > 0) then
            
            ! Update the node corresponding to the current central value:
            T%next(cvnext)%repl(cv) = T%next(cvnext)%repl(cv) + cpdfold
        
            ! If locations still need to be visited in the data template, visit them
            if (v < sN) then
                call UpdateTreeBoundary(T%next(cvnext), Ttrain%next(cvnext), cpdfold, Zvec, cv, sV, sN, exist) 
            else
                ! Otherwise we made it to the end and the voxel has an existing pattern
                exist = .true.
            end if    
        else
                                               
            ! For each possible category
            Ztmp = Zvec
            do c = 1, sV+1
            
                ! Current neighborhood
                Ztmp(v) = c-1
            
                ! Get ratio
                cpdf = 0.0_dp   
                call getCPDF(Ttrain, Ztmp, sV, sN, cpdf) 
                ratio(c) = cpdf(cv)
            end do
            
            ! For each possible category   
            Ztmp = Zvec
            do c = 1, sV+1
                Ztmp(v) = c-1      
                
                if (ratio(c) > 0) then
                    cpdfnew = cpdfold * ratio(c) / sum(ratio)
                    T%next(c)%repl(cv) = T%next(c)%repl(cv) + cpdfnew
                        
                    ! If locations still need to be visited in the data template, visit them
                    if (v < sN) then
                        call UpdateTreeBoundary(T%next(c), Ttrain%next(c), cpdfnew, Ztmp, cv, sV, sN, exist) 
                    else
                        ! Otherwise we made it to the end and the voxel has an existing pattern
                        exist = .true.
                    end if
                end if
            end do
        end if
    end if 
end subroutine UpdateTreeBoundary             