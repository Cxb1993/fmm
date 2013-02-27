function CompDataFit(Z,invProb) result(L)
    use TypeDef
    use precision
    use utils
    implicit none
    
    ! Input
    integer, intent(in), dimension(:,:,:)       :: Z
    type(inverseproblem), intent(in)            :: invProb
    real(dp)                                    :: L
    
    ! Local Variables    
    real(dp), allocatable, dimension(:,:)       :: G
    real(dp), allocatable, dimension(:)         :: dobs
    real(dp), allocatable, dimension(:,:)       :: invCov
    real(dp), allocatable, dimension(:)         :: cat
    integer                                     :: mZ, nZ, pZ, sV, c, nobs
    real(dp), allocatable, dimension(:,:,:)     :: m
    real(dp), allocatable, dimension(:)         :: Gm, res, invCovr, rinvCovr 
    real(dp), allocatable, dimension(:,:)       :: resT 
    
    ! Inverse problem parameters
    G       = invProb%G
    dobs    = invProb%dobs
    invCov  = invProb%invCov
    cat     = invProb%cat
    
    ! Dimensions
    nZ = size(Z,1); mZ = size(Z,2); pZ = size(Z,3)
    sV = size(cat)-1
    nobs = size(dobs)
    
    ! From categorical image to model paramters
    allocate(m(nZ,mZ,pZ))
    do c = 1, sV+1    
        where (Z == c-1)
            m = cat(c)
        end where
    end do
       
    ! Compute Gmat * m
    if (allocated(Gm)) deallocate(Gm)
    allocate(Gm(1:size(dobs)))
    call multMatVec(G, ravel(m), Gm)
    
    ! Residual
    if (allocated(res)) deallocate(res)
    allocate(res(1:nobs), source = dobs-Gm)
    if (allocated(resT)) deallocate(resT)
    allocate(resT(1:1,1:nobs), source = reshape(res, shape=[1,nobs]))    
    
    ! invCov * res
    if (allocated(invCovr)) deallocate(invCovr)
    allocate(invCovr(1:nobs))
    call multMatVec(invCov, res, invCovr)
    
    ! res' * (invCov*res)
    if (allocated(rinvCovr)) deallocate(rinvCovr)
    allocate(rinvCovr(1:nobs))
    call multMatVec(resT, invCovr, rinvCovr)
            
    L = 0.5_dp * rinvCovr(1)
    
    if (L < 0) then
        print *, 'ERROR in CompDataFit: result is negative.'
    end if
             
end function CompDataFit