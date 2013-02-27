!     
! File:   test_FMM.f90
! Author: Dr. Jan Frydendall
!
! Created on May 17, 2011, 9:12 AM
program test_FMM
    use utils
    use interfaces
    use TypeDef
    use precision
    implicit none
    
    ! Parameters
    integer, parameter                                   :: nZt = 100, mZt = 100, pZt = 1               ! training image
    integer, parameter                                   :: nZ = 50, mZ = 50, pZ = 1                    ! image 
    integer, parameter                                   :: sV = 2                                      ! categories
    integer, parameter                                   :: nN = 5, mN = 5, pN = 1                      ! neighborhood mask
    integer, parameter                                   :: nD = 15, mD = 15, pD = 1                    ! perturbation domain
    integer, parameter                                   :: nobs = 1575, npar = 2500                    ! data fit
    real(dp), parameter                                  :: iter = 0.10_dp                              ! iterations per voxel
            
    ! Local variables
    integer, dimension(nZt,mZt,pZt)                      :: Ztrain
    integer, dimension(nZ,mZ,pZ)                         :: Z0
    logical, dimension(nZ,mZ,pZ)                         :: Zcond
    integer, dimension(nN,mN,pN)                         :: N
    integer, allocatable, dimension(:)                   :: D
    type(inverseproblem)                                 :: invProb
    real(dp), dimension(nobs,npar)                       :: Gmat
    real(dp), dimension(nobs)                            :: dobs
    real(dp)                                             :: lambda, alpha
    real(dp), dimension(nobs,nobs)                       :: invCov
    real(dp), dimension(sV+1)                            :: cat
    integer, allocatable, dimension(:)                   :: tmp1D
    real(dp), allocatable, dimension(:)                  :: tmp1Dr
    integer                                              :: nrseed, clock, i, v, aval
    real(dp), allocatable, dimension(:)                  :: avec, t0
    integer, allocatable, dimension(:)                   :: nseed
    type(option)                                         :: options
    type(solution)                                       :: solutions
    real(dp), allocatable, dimension(:,:)                :: Ps
    type(NeighborMask)                                   :: Nmask
    type(DomainMask)                                     :: Dblock
   
    ! Initialize seed of random generator
    call random_seed(size = nrseed)
    if (.True.)then
        allocate(nseed(1:nrseed))
        nseed = 4
    else
        CALL SYSTEM_CLOCK(COUNT = clock)
        allocate(nseed(1:nrseed))
        nseed = clock
    end if
    call random_seed(PUT = nseed)

    ! --------------------- TRAINING IMAGE ----------------------
    
    ! Read from file
    open(010,file='Sara/TI_2.txt', status='old')
    if (allocated(tmp1D)) deallocate(tmp1D)
    allocate(tmp1D(nZt*mZt*pZt))
    do i = 1,nZt*mZt*pZt
        read(010,*) tmp1D(i)
    end do
    close(010)    
    Ztrain = reshape(tmp1D, shape = [nZt, mZt, pZt])
                 
    
    ! --------------------- STARTING IMAGE ----------------------

    ! Read from file
    open(011,file='Sara/Z0_2.txt', status='old')
    if (allocated(tmp1D)) deallocate(tmp1D)
    allocate(tmp1D(nZ*mZ*pZ))        
    do i = 1,nZ*mZ*pZ
        read(011,*) tmp1D(i)
    end do
    close(011)    
    Z0 = reshape(tmp1D, shape = [nZ, mZ, pZ])
        
    ! Condition on hard data
    Zcond = .False.
    
    ! ---------------------- SEISMIC PART ------------------------
    
    ! System matrix read from file
    open(020,file='Sara/Gmat_2.txt',status='old')
    if (allocated(tmp1Dr)) deallocate(tmp1Dr)
    allocate(tmp1Dr(nobs*npar))
    do i = 1, nobs*npar
        read(020,*) tmp1Dr(i)
    end do
    close(020)     
    Gmat = reshape(tmp1Dr, shape = [nobs, npar])
            
    ! Vector of data observations read from file
    open(030,file='Sara/dobs_2.txt',status='old')
    if (allocated(tmp1Dr)) deallocate(tmp1Dr)
    allocate(tmp1Dr(nobs))
    do i = 1, nobs
        read(030,*) tmp1Dr(i)
    end do
    close(030)
    dobs = tmp1Dr
    
    ! Physical values for each of the categories 
    do i = 1, sV+1
        cat(i) = i-1
    end do
    
    ! Inverse covariance matrix of the data observations
    lambda = 1
    invCov = 0.00_dp
    do i = 1, size(dobs)
        invCov(i,i) = 1.00_dp / lambda**2
    end do    
    
    ! Specify the inverse problem 
    invProb%G = Gmat
    invProb%dobs = dobs
    invProb%cat = cat
    invProb%invCov = invCov
    
    ! --------------------- FM PARAMETERS -----------------
           
    ! Neighborhood parameters
    N = 1                                                           ! rectangular neighborhood 
    N(1,1,1) = 0; N(1,nN,1) = 0; N(nN,1,1) = 0; N(nN,nN,1) = 0      ! leaving out the corners
    Nmask%n = nN
    Nmask%m = mN
    Nmask%p = pN
    Nmask%nc = (nN+1)/2
    Nmask%mc = (mN+1)/2
    Nmask%pc = (pN+1)/2
    Nmask%mat = N

    ! Dimensions of the perturbation domain
    Dblock%n = nD
    Dblock%m = mD
    Dblock%p = pD
        
    ! Algorithm parameters
    options % runs = 1
    options % multigrid = 1
    options % condopt = .False.
    
    ! ---- TEST OF FMM ------ 
    
    ! Number of alpha parameters to solve for
    aval = 1
    allocate(avec(aval), t0(aval))

    ! Set alpha and minimum temperature
    avec = 1e-2
    t0 = 1e5
    
    ! Run for each alpha
    do i = 1, aval
        
        alpha = avec(i)     
        options % t0 = t0(i)
        options % maxIter = iter        
        options % tmin = 1d-5
        
        call FMM(Z0,Ztrain,Zcond,Nmask,sV,invProb,alpha,Dblock,options,solutions,Ps)
    end do
    
end program test_FMM