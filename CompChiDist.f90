function CompChiDist(H,Htrain,N) result(X)
    use precision
    use utils
    implicit none
    real(dp), intent(in), dimension(:)             :: H, Htrain
    integer                                        :: N
    real(dp)                                       :: X
    ! Local Variables
    real(dp)                                       :: np, npi
    integer, allocatable, dimension(:)             :: ids
    real(dp), allocatable, dimension(:)            :: p, pi      
           
    call find(Htrain > 0, ids)
    p = H(ids)
    pi = Htrain(ids)
        
    np = real(N, kind = dp)
    npi = sum(pi)
    
    X = sum(( sqrt(np/npi)*pi - sqrt(npi/np)*p )**2 / (pi + p)) + npi - npi/np * sum(p)  
    
    if (X < 0) then
        print *, 'ERROR in CompChiDist: result is negative.'
        stop
    end if
    
end function CompChiDist