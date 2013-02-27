!     
! File:   FMM.f90
! Author: Dr. Jan Frydendall < jf @ imm.dtu.dk >
!
! Created on May 4, 2011, 9:02 AM
!
subroutine FMM(Z0,Ztrain,Zcond,Nmask,sV,invProb,alpha,Dblock,options,solutions,Ps)
    use Interfaces, only : InferTrainTree, Tree2Hist, CompOptimalImage
    use TypeDef
    use precision
    use utils
    implicit none
    integer, intent(in), dimension(:,:,:)               :: Z0, Ztrain
    integer, intent(in)                                 :: sV
    logical, intent(in), dimension(:,:,:)               :: Zcond
    type(NeighborMask), intent(inout)                   :: Nmask
    type(inverseproblem), intent(in)                    :: invProb
    real(dp), intent(in)                                :: alpha
    type(DomainMask), intent(inout)                     :: Dblock
    type(option), intent(inout)                         :: options
    type(solution), intent(out)                         :: solutions
    real(dp), intent(out),allocatable, dimension(:,:)   :: Ps
        
    ! Local Variables
    integer                                             :: nc, mc, pc, nN, mN, pN
    integer                                             :: ncd, mcd, pcd, nD, mD, pD
    integer                                             :: nZ, mZ, pZ, nZt, mZt, pZt
    integer                                             :: sN, i, j, k, ijk, c
    integer, allocatable, dimension(:,:,:)              :: N, Zopt
    integer, allocatable, dimension(:)                  :: D
    integer, allocatable, dimension(:,:)                :: nodes, dnodes
    type(streenode)                                     :: Ttrain
    real(dp), allocatable, dimension(:,:)               :: H0, Htrain, Hopt
    integer, allocatable, dimension(:,:,:)              :: ZcondInt
      
    ! Image dimensions
    nZ  = size(Z0,1);       mZ  = size(Z0,2);       pZ  = size(Z0,3)
    nZt = size(Ztrain,1);   mZt = size(Ztrain,2);   pZt = size(Ztrain,3)
  
    ! Neighborhood mask
    if (allocated(Nmask%mat)) then
        
        ! Dimensions if array N is provided
        N = Nmask%mat
        nN = size(N,1); mN = size(N,2); pN = size(N,3)
        nc = Nmask%nc;  mc = Nmask%mc;  pc = Nmask%pc
        
        if (nc < 1 .or. nN < nc .or. mc < 1 .or. mN < mc .or. pc < 1 .or. pN < pc) then
            print *, 'ERROR in FMM: Invalid neighborhood center coordinates'
            stop
        end if        
    else
        
        ! Dimensions provided
        nN = Nmask%n;  mN = Nmask%m;  pN = Nmask%p

        if (floor(real(nN+1)/2.00_dp) .ne. ceiling(real(nN+1)/2.00_dp) .or. &
            floor(real(mN+1)/2.00_dp) .ne. ceiling(real(mN+1)/2.00_dp) .or. &
            floor(real(pN+1)/2.00_dp) .ne. ceiling(real(pN+1)/2.00_dp)) then
            print *, 'ERROR in FMM: Invalid neighborhood dimensions (they must be odd)'
            stop
        end if

        ! Center coordinates    
        nc = (nN+1)/2; mc = (mN+1)/2; pc = (pN+1)/2
        allocate(N(nN,mN,pN), source = 1)        
    end if
        
    ! Binary mask
    N(nc,mc,pc) = 0
    sN = sum(ravel(N));
    Nmask%mat = N
    Nmask%nc  = nc
    Nmask%mc  = mc
    Nmask%pc  = pc
    Nmask%n   = nN
    Nmask%m   = mN
    Nmask%p   = pN
        
    ! Dimension of the perturbation domain
    nD  = Dblock%n;
    mD  = Dblock%m;
    pD  = Dblock%p
    if (floor(real(nD+1)/2.00_dp) .ne. ceiling(real(nD+1)/2.00_dp) .or. &
        floor(real(mD+1)/2.00_dp) .ne. ceiling(real(mD+1)/2.00_dp) .or. &
        floor(real(pD+1)/2.00_dp) .ne. ceiling(real(pD+1)/2.00_dp)) then
        print *, 'ERROR in FMM: Invalid perturbation domain dimensions (they must be odd)'
        stop
    end if
    
    if (nD < 1 .or. nZ < nD .or. mD < 1 .or. mZ < mD .or. pD < 1 .or. pZ < pD) then
        print *, 'ERROR in FMM: Invalid perturbation domain dimensions (bigger than the image)'
        stop
    end if   
    
    ! Distances taking into account the correlations length
    allocate(D(sN), source = 0)       
    ijk = 0
    do k = 1, pN 
        do j = 1, mN 
            do i = 1, nN      
                if (N(i,j,k)==1) then
                    ijk = ijk + 1
                    D(ijk) = nN*mN*abs(pc-k) + nN*pN*abs(mc-j) + mN*pN*abs(nc-i)
                end if
            end do
        end do
    end do  
    
    ! Perturbation domain
    ncd = (nD+1)/2
    mcd = (mD+1)/2
    pcd = (pD+1)/2
    Dblock%mat = D
    Dblock%nc = ncd
    Dblock%mc = mcd
    Dblock%pc = pcd
       
    ! Compute neighborhood coordinates
    allocate(nodes(sN,3),source = 0)
    ijk = 0
    do k = 1, pN
        do j = 1, mN 
            do i = 1, nN
                if (N(i,j,k)==1) then
                    ijk = ijk + 1
                    nodes(ijk, :) = [i-nc, j-mc, k-pc]
                end if
            end do
        end do
    end do 
    Nmask%nodes = nodes
    
    ! Compute block coordinates
    allocate(dnodes(nD*mD*pD,3),source = 0)
    ijk = 0
    do k = -pcd+1, pD-pcd
        do j = -mcd+1, mD-mcd
            do i = -ncd+1, nD-ncd
                ijk = ijk + 1
                dnodes(ijk,:) = [i, j, k]
            end do
        end do
    end do    
    Dblock%nodes = dnodes
    
    ! Compute training search tree
    Ttrain%depth = 0
    allocate(Ttrain%repl(sV+1), source = 0.0_dp)
    nullify(Ttrain%next)   
    call InferTrainTree(ZTrain, nodes, sV, sN, Ttrain)    
    call Tree2Hist(sV, sN, Ttrain, Htrain)
    print *, 'Training search tree computed'
    
    call allocAs(Ztrain,solutions%Ztrain)
    solutions%sV = sV
       
    ! Save location of conditioned data
    if (allocated(Zcondint)) deallocate(Zcondint)
    allocate(Zcondint(nZ,mZ,pZ), source = 0)
    where (Zcond)
        Zcondint = 1
    end where
    
    ! Adjust algorithm paramters
    options%maxIter = options%maxIter * real((nZ * mZ * pZ - count(Zcond)), kind = dp)                  
    
    ! Compute optimal solution
    print *, 'Initiate optimization'
    call CompOptimalImage(Z0, Zcond, Ttrain, Htrain, nodes, sV, sN, invProb, alpha, Dblock, options, Zopt, Hopt, Ps)
    
    ! Save the computed histogram
    call allocAs(Z0, solutions%Z0)
    call allocAs(Zcondint, solutions%Zcond)
    call allocAs(Htrain, solutions%Htrain)
    call allocAs(Zopt, solutions%Zopt)
    call allocAs(Hopt, solutions%Hopt)
    
    ! Initial 
    write(100,*)solutions%Ztrain
    write(110,*)solutions%Htrain
    write(320,*)solutions%Z0
    
    ! Optimal solution
    write(300,*)solutions%Zopt
    write(350,*)solutions%Hopt
    write(370,*)Ps(:,1)
    write(371,*)Ps(:,2)

end subroutine FMM