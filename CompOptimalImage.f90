!     
! File:   CompZtest.f03
! Author: Dr. Jan Frydendall
!
! Created on March 24, 2011, 2:18 PM
!
subroutine CompOptimalImage(Z0,Zcond,Ttrain,Htrain,nodes,sV,sN,invProb,alpha,Dblock,options,Zopt,Hopt,Ps)
    use precision
    use Interfaces, only : InferTree, Tree2Hist, CompObjFun, CopyTree, getNewIt, SimNewImage, UpdateSA, DeallocateTree
    use TypeDef
    use utils
    implicit none
    
    ! Inputs
    integer, intent(in), dimension(:,:,:)               :: Z0
    logical, intent(in), dimension(:,:,:)               :: Zcond
    type(streenode), intent(inout)                      :: Ttrain
    real(dp), intent(inout), dimension(:,:)             :: Htrain
    integer, intent(in), dimension(:,:)                 :: nodes
    integer, intent(in)                                 :: sV, sN
    type(inverseproblem), intent(in)                    :: invProb
    real(dp), intent(in)                                :: alpha
    type(DomainMask), intent(in)                        :: Dblock
    type(option), intent(in)                            :: options
    
    ! Outputs
    integer, intent(out), allocatable, dimension(:,:,:) :: Zopt
    real(dp), intent(out), allocatable, dimension(:,:)  :: Hopt
    real(dp), intent(out), allocatable, dimension(:,:)  :: Ps
    
    ! Local variables
    integer                                             :: Ntest
    integer                                             :: nZ, mZ, pZ, i, j, k
    integer, allocatable, dimension(:,:,:)              :: Ztest, Znew
    logical, allocatable, dimension(:,:,:)              :: Zex, Zexnew, Zexopt
    type(streenode)                                     :: Ttest, Tnew, Topt
    real(dp), allocatable, dimension(:,:)               :: Htest, Hnew
    real(dp), dimension(1:2)                            :: P0, Ptest, Pnew, Popt, deltaP 
    integer                                             :: it, maxiter, suggested, accepted
    real(dp)                                            :: t0, tmin, T, rate, p
    logical                                             :: newImage
    
    ! Dimensions
    nZ = size(Z0,1);    mZ = size(Z0,2);    pZ = size(Z0,3)
    
    ! Temperature and stopping criteria
    t0      = options%t0
    tmin    = options%tmin
    maxIter = floor(options % maxIter)
    rate    = (tmin/t0)**(1.0_dp/options % maxIter)
      
    ! Compute initial tree, histogram and misfit function value
    allocate(Ztest(nZ,mZ,pZ), source = Z0)
    call InferTree(Ztest, nodes, Ttrain, sV, sN, Ttest, Zex)
    Ntest = floor(sum(Ttest%repl))
    call Tree2Hist(sV, sN, Ttest, Htest)
    call CompObjFun(Htest, Htrain, Ntest, Ztest, invProb, alpha, Ptest)
    
    ! Initialization
    Tnew%depth = 0; allocate(Tnew%repl(sV+1), source = 0.0_dp); nullify(Tnew%next)  
    Topt%depth = 0; allocate(Topt%repl(sV+1), source = 0.0_dp); nullify(Topt%next)    
    allocate(Ps(maxIter+1, 2), source = -1.0_dp); Ps(1,:) = Ptest / [alpha, 1.00_dp]
    
    ! Initial optimal values
    allocate(Zopt(nZ,mZ,pZ), source = Ztest)
    allocate(Zexopt(nZ,mZ,pZ), source = Zex)
    call CopyTree(Ttest, sV, sN, Topt)
    allocate(Hopt(sV+1,size(Htest,2)), source = Htest)
    Popt = Ptest
        
    ! Initial algorithm parameters
    T = t0; it = 0; suggested = 0; accepted = 0
            
    ! For each iteration...
    do while (it < maxIter)
        
        ! Get new iteration
        call getNewIt(Zex, nodes, sN, i, j, k)
        
        ! Make movie
        if (mod(it, 50) == 0) then
            write(400,*)Ztest
        end if     
        
        ! Simulate the neighboring image
        call SimNewImage(Ztest,Zcond,Zex,Ttest,Htest,Ttrain,i,j,k,nodes,sV,sN,Dblock,Znew,Zexnew,Tnew,Hnew,newImage)      
        
        ! If a new image has been suggested
        if (newImage) then
            
            ! Continue after this iteration
            it = it + 1
            
            ! Compute new objective value
            call CompObjFun(Hnew, Htrain, Ntest, Znew, invProb, alpha, Pnew)            
            deltaP = Pnew - Ptest  
            
            ! If new solution is better ...
            if (sum(deltaP) < 0.0_dp) then

                ! Take the step
                call UpdateSA(sV, sN, Znew, Zexnew, Tnew, Hnew, Pnew, Ztest, Zex, Ttest, Htest, Ptest)
                Ps(it+1,:) = Pnew / [alpha, 1.00_dp]

                ! In case the new solution is the best so far
                if (sum(Ptest) < sum(Popt)) then
                    call UpdateSA(sV, sN, Ztest, Zex, Ttest, Htest, Ptest, Zopt, Zexopt, Topt, Hopt, Popt)
                end if

                ! ... if not better, choose it anyway?
            else
                suggested = suggested + 1

                ! Random number between 0 and 1
                call random_number(p)

                ! Transition probability
                if (p < exp(-sum(deltaP) / T)) then

                    ! Accept the step
                    call UpdateSA(sV, sN, Znew, Zexnew, Tnew, Hnew, Pnew, Ztest, Zex, Ttest, Htest, Ptest)
                    Ps(it+1,:) = Pnew / [alpha, 1.00_dp]
                    accepted = accepted + 1
                    
                end if
            end if
                         
            ! Update the temperature
            T = rate * T        
            
            ! Display the number of iterations run
            if ( mod(real(it,kind=dp)/real(maxIter,kind=dp)*100.0_dp, 5.0_dp) == 0) then
                print *, int(real(it,kind=dp)/real(maxIter,kind=dp)*100.0_dp), '%'
            end if  
        end if
 
        ! Deallocate the extra tree
        call DeallocateTree(sV, Tnew)
        if (it < maxIter) then
            Tnew%depth = 0; allocate(Tnew%repl(1:sV+1), source = 0.0_dp); nullify(Tnew%next)   
        end if
    end do
    print *, 'Accept ratio: ', real(accepted,kind=dp)/real(suggested,kind=dp), ' (out of ', suggested, ' suggested models)'      
   
end subroutine CompOptimalImage