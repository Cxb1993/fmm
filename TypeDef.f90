!     
! File:   TypeDef.f90
! Author: Dr. Jan Frydendall
!
! Created on May 12, 2011, 10:07 AM
!

MODULE TypeDef
    use precision
    
    type inverseproblem
        
        real(dp), allocatable, dimension(:,:)   :: G
        real(dp), allocatable, dimension(:)     :: dobs
        real(dp), allocatable, dimension(:,:)   :: invCov
        real(dp), allocatable, dimension(:)     :: cat
    end type inverseproblem
    
    type streenode
        real(dp), allocatable, dimension(:)     :: repl 
        integer                                 :: depth
        type(streenode), dimension(:), pointer  :: next
    end type streenode
        
    type info
        real(dp)                                :: PL
        real(dp)                                :: PP
        real(dp)                                :: steps
        real(dp)                                :: improved
        real(dp)                                :: accepted
        real(dp)                                :: suggest
    end type info
    
    type solution
        integer                                 :: sV
        integer, allocatable, dimension(:,:,:)  :: Zopt, Ztrain, Z0, Zcond
        real(dp), allocatable, dimension(:,:)   :: Htrain, H0, Hopt
        real(dp), allocatable, dimension(:)     :: Popt
        type(info), allocatable, dimension(:,:) :: infoH
    end type solution

    type option
        real(dp)                                :: t0, tmin, maxIter
        integer                                 :: runs, multigrid
        logical                                 :: condopt
    end type option
    
    type NeighborMask
        integer, allocatable, dimension(:,:,:)  :: mat
        integer                                 :: nc, mc, pc, n, m, p
        integer, allocatable, dimension(:,:)    :: nodes
    end type NeighborMask
   
    type DomainMask
        integer, allocatable, dimension(:)      :: mat
        integer                                 :: nc, mc, pc, n, m, p
        integer, allocatable, dimension(:,:)    :: nodes
    end type DomainMask    

END MODULE TypeDef
