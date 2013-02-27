recursive subroutine Tree2Hist(sV, sN, T, H)
    use precision
    use TypeDef
    implicit none
    ! Input
    integer, intent(in)                                 :: sV, sN
    type(streenode), intent(in)                         :: T
    
    ! Input/Output
    real(dp), intent(inout), allocatable, dimension(:,:):: H
    
    ! Local variables
    integer                                             :: c, nh
    real(dp), allocatable, dimension(:,:)               :: Htmp
            
    ! If more nodes
    if (associated(T%next)) then
        
        ! Go to the last node
        do c = 1, sV+1
            call Tree2Hist(sV, sN, T%next(c), H)
        end do
        
    else
        
        ! If at an end node
        if (T%depth == sN) then 

            ! Update the histogram
            if (allocated(H)) then

                ! Number of patterns already found
                nh = size(H,2)        

                ! Temporarily variable
                allocate(Htmp(sV+1,nh+1))
                Htmp(:,1:nh) = H
                Htmp(:,nh+1) = T%repl
                
                ! Updated histogram
                deallocate(H)
                allocate(H(sV+1,nh+1), source = Htmp)
                deallocate(Htmp)
            else

                ! First pattern found
                allocate(H(sV+1,1))
                do c = 1, sV+1
                    H(c,1) = T%repl(c)
                end do

            end if
        end if
    end if
    
end subroutine Tree2Hist
    

