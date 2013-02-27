subroutine SimVoxel(Zvec, Ttrain, sV, sN, D, hc)
    use precision
    use utils
    use Interfaces, only : getCPDF
    use TypeDef
    implicit none
    
    ! Input
    integer, intent(in), dimension(:)                   :: Zvec, D
    type(streenode), intent(in)                         :: Ttrain
    integer, intent(in)                                 :: sV, sN
            
    ! Output
    real(dp), intent(out), allocatable, dimension(:)    :: hc

    ! Local variables
    integer                                             :: v, cmin
    integer, allocatable, dimension(:)                  :: nknown, Dknown, vknown 
    integer, dimension(sN)                              :: Ztmp
    logical                                             :: simed
    real(dp), dimension(sV+1)                           :: htmp, cpdf
    real(dp)                                            :: rand
            
    ! Initialization
    simed = .false.
    cmin = 0
    cpdf = 0
    Ztmp = Zvec
    
    ! While not yet simulated
    do while (.not. simed)
    
        ! Index of known neighbors
        call find(Ztmp > -1, nknown)
        
        ! If no known neighbors use 1p statistics
        if (size(nknown)==0) then
          
            htmp = Ttrain%repl                         
            simed = .true.
            
        else
            
            ! Get the conditional PDF
            cpdf = 0.0_dp
            
            call getCPDF(Ttrain, Ztmp, sV, sN, cpdf)
                        
            ! If any neighborhoods found in the TI
            if (sum(cpdf) > 0.0_dp) then                
                htmp = cpdf
                if (sum(htmp) > cmin) simed = .true.
            end if
            
            ! If not, remove voxel that is conditioned upon   
            if (.not. simed) then
                if (size(nknown) > 1) then
                                       
                    ! Distance to known neighbors
                    call flatten(D, nknown, Dknown)
                                        
                    ! Remove random voxel with highest distance
                    call find(Ztmp > -1 .and. D==maxval(Dknown), vknown)
                    call random_number(rand)
                    v = floor(1.0_dp + (size(vknown) * rand))
                    Ztmp(vknown(v)) = - 1   
                    
                else 
                    htmp = Ttrain%repl
                    simed = .true.
                end if
            end if
        end if   
    end do
    
    ! Return the cumulative cpdf
    allocate(hc(sV+1), source = cumsum(htmp))
        
end subroutine SimVoxel