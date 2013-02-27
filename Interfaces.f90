!     
! File:   Interfaces.f90
! Author: Dr. Jan Frydendall
!
! Created on May 12, 2011, 8:56 AM
!

MODULE Interfaces
    interface                      

        subroutine FMM(Z0,Ztrain,Zcond,Nmask,sV,invProb,alpha,Dblock,options,solutions,Ps)
            use TypeDef
            use precision
            use utils
            implicit none
            
            ! Inputs
            integer, intent(in), dimension(:,:,:)               :: Z0, Ztrain
            integer, intent(in)                                 :: sV
            logical, intent(in), dimension(:,:,:)               :: Zcond
            type(NeighborMask), intent(inout)                   :: Nmask
            type(inverseproblem), intent(in)                    :: invProb
            real(dp), intent(in)                                :: alpha
            type(DomainMask), intent(inout)                     :: Dblock
            
            ! Outputs
            type(option), intent(inout)                         :: options
            type(solution), intent(out)                         :: solutions
            real(dp), intent(out),allocatable, dimension(:,:)   :: Ps
        end subroutine FMM

        subroutine CompOptimalImage(Z0,Zcond,Ttrain,Htrain,nodes,sV,sN,invProb,alpha,Dblock,options,Zopt,Hopt,Ps)
            use precision
            use TypeDef
            use utils
            implicit none
            
            ! Inputs
            integer, intent(in), dimension(:,:,:)               :: Z0
            logical, intent(in), dimension(:,:,:)               :: Zcond
            type(streenode), intent(in)                         :: Ttrain
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
        end subroutine CompOptimalImage
                              
        subroutine SimNewImage(Ztest,Zcond,Zex,Ttest,Htest,Ttrain,i,j,k,nodes,sV,sN,Dblock,Znew,Zexnew,Tnew,Hnew,newImage)
            use precision
            use utils
            use TypeDef
            implicit none

            ! Inputs
            integer, intent(in), dimension(:,:,:)               :: Ztest
            logical, intent(in), dimension(:,:,:)               :: Zcond, Zex
            integer, intent(in)                                 :: i, j, k, sV, sN
            integer, intent(in), dimension(:,:)                 :: nodes
            type(streenode), intent(inout)                      :: Ttest, Ttrain
            real(dp), intent(in), dimension(:,:)                :: Htest
            type(DomainMask), intent(in)                        :: Dblock
            
            ! Outputs
            integer, intent(out), allocatable, dimension(:,:,:) :: Znew
            logical, intent(out), allocatable, dimension(:,:,:) :: Zexnew
            type(streenode), intent(out)                        :: Tnew
            real(dp), intent(out), allocatable, dimension(:,:)  :: Hnew
            logical, intent(out)                                :: newImage
        end subroutine SimNewImage
        
        subroutine SimVoxel(Zvec, Ttrain, sV, sN, D, hc)
            use precision
            use utils
            use TypeDef
            implicit none

            ! Input
            integer, intent(in), dimension(:)                   :: Zvec, D
            type(streenode), intent(in)                         :: Ttrain
            integer, intent(in)                                 :: sV, sN

            ! Output
            real(dp), intent(out), allocatable, dimension(:)    :: hc
        end subroutine SimVoxel
        
        
        ! ----------------- OBJECTIVE FUNCTIONS ----------------------
        
        subroutine CompObjFun(H,Htrain,N,Z,invProb,alpha,P)
            use TypeDef
            use precision
            use utils
            implicit none
            ! Input
            real(dp), intent(in), dimension(:,:)           :: H, Htrain   
            integer, intent(in)                            :: N
            integer, intent(in), dimension(:,:,:)          :: Z
            type(inverseproblem), intent(in)               :: invProb
            real(dp), intent(in)                           :: alpha
            ! Output
            real(dp), intent(out), dimension(1:2)          :: P
        end subroutine CompObjFun
   
        function CompChiDist(H,Htrain,N) result(X)
            use precision
            use utils
            implicit none
            real(dp), intent(in), dimension(:)             :: H, Htrain
            integer                                        :: N
            real(dp)                                       :: X
        end function
        
        function CompDataFit(Z,invProb) result(L)
            use TypeDef
            use precision
            use utils
            implicit none

            ! Input
            integer, intent(in), dimension(:,:,:)           :: Z
            type(inverseproblem), intent(in)                :: invProb
            real(dp)                                        :: L
        end function
        
        
        ! ----------------- AUXILIARY FUNCTIONS ----------------------
        
        subroutine getNeighborhood(Z, i, j, k, nodes, sN, Zvec, innervoxel)
            use TypeDef
            use utils
            implicit none

            ! Input
            integer, intent(in), dimension(:,:,:)           :: Z
            integer, intent(in)                             :: i, j, k, sN
            integer, intent(in), dimension(:,:)             :: nodes

            ! Output
            integer, intent(out), dimension(sN)             :: Zvec
            logical                                         :: innervoxel
        end subroutine getNeighborhood  
        
        subroutine getNewIt(Zex, nodes, sN, i, j, k)
            use precision
            use utils
            implicit none

            ! Inputs
            logical, intent(in), dimension(:,:,:)   :: Zex
            integer, intent(in), dimension(:,:)     :: nodes
            integer, intent(in)                     :: sN

            ! Outputs
            integer, intent(out)                    :: i, j, k
        end subroutine getNewIt
               
        subroutine UpdateSA(sV, sN, Zin, Zexin, Tin, Hin, Pin, Zout, Zexout, Tout, Hout, Pout)
            use precision
            use utils
            use TypeDef
            implicit none

            ! Input
            integer, intent(in)                         :: sV, sN
            integer, intent(in), dimension(:,:,:)       :: Zin
            logical, intent(in), dimension(:,:,:)       :: Zexin
            type(streenode), intent(inout)              :: Tin
            real(dp), intent(in), dimension(:,:)        :: Hin
            real(dp), intent(in), dimension(1:2)        :: Pin

            ! Output
            integer, intent(out), dimension(:,:,:)      :: Zout
            logical, intent(out), dimension(:,:,:)      :: Zexout
            type(streenode), intent(out)                :: Tout
            real(dp), intent(out), dimension(:,:)       :: Hout
            real(dp), intent(out), dimension(1:2)       :: Pout
        end subroutine UpdateSA  
        
        recursive subroutine Tree2Hist(sV, sN, T, H)
            use precision
            use TypeDef
            implicit none
            
            ! Input
            integer, intent(in)                                 :: sV, sN

            ! Input/Output
            type(streenode), intent(in)                         :: T
            real(dp), intent(out), allocatable, dimension(:,:)  :: H
        end subroutine Tree2Hist
        
        
        
        
        
        
        
        ! ----------------- TREE MANIPULATING FUNCTIONS ----------------------
        
        subroutine InferTrainTree(Ztrain, nodes, sV, sN, Ttrain)
            use precision
            use TypeDef
            implicit none
            ! InferTrainTree: Infers the search tree for a training image

            ! Input
            integer, intent(in), dimension(:,:,:)   :: Ztrain
            integer, intent(in), dimension(:,:)     :: nodes
            integer, intent(in)                     :: sV, sN

            ! Output
            type(streenode), intent(inout)          :: Ttrain
        end subroutine InferTrainTree
                     
        recursive subroutine UpdateTrainTree(Train, Zvec, cv, sV, sN)
            use precision
            use TypeDef
            implicit none
            ! UpdateTrainTree: Updates the search tree for a training image

            ! Input
            integer, intent(in), dimension(:)       :: Zvec
            integer, intent(in)                     :: cv, sV, sN
             
            ! Input/Output
            type(streenode), intent(inout)          :: Train
        end subroutine UpdateTrainTree
        
        
        subroutine InferTree(Z, nodes, Ttrain, sV, sN, T, Zex)
            use utils
            use precision
            use TypeDef
            implicit none
            ! InferTrainTree: Infers the search tree for a training image

            ! Input
            integer, intent(in), dimension(:,:,:)   :: Z
            integer, intent(in), dimension(:,:)     :: nodes
            type(streenode), intent(in)             :: Ttrain
            integer, intent(in)                     :: sV, sN

            ! Output
            type(streenode), intent(out)            :: T
            logical, intent(out), allocatable, dimension(:,:,:) :: Zex
        end subroutine InferTree
        
        subroutine wrapUpdateTree(T, Ttrain, Zvec, cv, sV, sN, exist)
            use TypeDef
            implicit none
            ! InferTrainTree: Infers the search tree for a training image

            ! Input
            type(streenode), intent(in)             :: Ttrain
            integer, intent(in), dimension(:)       :: Zvec
            integer, intent(in)                     :: cv, sV, sN

            ! Input/Output
            type(streenode), intent(inout)          :: T
            logical, intent(out)                    :: exist
        end subroutine wrapUpdateTree
            
        recursive subroutine UpdateTree(T, Ttrain, Zvec, cv, sV, sN, exist)
            use utils
            use precision
            use TypeDef
            implicit none
                        
            ! Input
            integer, intent(in), dimension(:)       :: Zvec
            integer, intent(in)                     :: cv, sV, sN

            ! Input/Output
            type(streenode), intent(inout)          :: T
            type(streenode), intent(inout)          :: Ttrain
            logical, intent(inout)                  :: exist
        end subroutine UpdateTree
        
        subroutine wrapUpdateTreeBoundary(T, Ttrain, cpdf, Zvec, cv, sV, sN, exist)
            use precision
            use TypeDef
            implicit none

            ! Input
            type(streenode), intent(in)             :: Ttrain
            real(dp), intent(in)                    :: cpdf
            integer, intent(in), dimension(:)       :: Zvec
            integer, intent(in)                     :: cv, sV, sN

            ! Input/Output
            type(streenode), intent(inout)          :: T
            logical, intent(out)                    :: exist
        end subroutine wrapUpdateTreeBoundary        
        
        recursive subroutine UpdateTreeBoundary(T, Ttrain, cpdfold, Zvec, cv, sV, sN, exist)
            use utils
            use precision
            use TypeDef
            implicit none

            ! Input
            integer, intent(inout), dimension(:)    :: Zvec
            integer, intent(in)                     :: cv, sV, sN
            type(streenode), intent(in)             :: Ttrain
            
            ! Input/Output
            type(streenode), intent(inout)          :: T
            real(dp), intent(inout)                 :: cpdfold
            logical, intent(inout)                  :: exist
        end subroutine UpdateTreeBoundary
       
        
        
        
        
        subroutine GrowTree(Z, zold, Zex, i, j, k, nodes, sV, sN, Ttrain, T)
            use TypeDef
            use utils
            implicit none

            ! Input
            integer, intent(in), dimension(:,:,:)   :: Z
            integer, intent(in)                     :: zold, i, j, k, sV, sN
            integer, intent(in), dimension(:,:)     :: nodes
            type(streenode), intent(in)             :: Ttrain

            ! Input/Output
            logical, intent(inout), allocatable, dimension(:,:,:) :: Zex
            type(streenode), intent(inout)          :: T
        end subroutine GrowTree
        
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
        end subroutine getCPDF
                  
        
        
        
        
        
        
        recursive subroutine CenterCount(T, sN, znew, zold, Zvec, exist)
            use TypeDef
            implicit none
            ! CenterCount: Updates the tree wrt. the center voxel when it is flipped

            ! Input
            integer, intent(in)                 :: sN, znew, zold
            integer, intent(in), dimension(:)   :: Zvec

            ! Input/Output
            type(streenode), intent(inout)      :: T
            logical, intent(inout)              :: exist
        end subroutine CenterCount
        
        recursive subroutine AddCount(T, sN, level, Zvec, zcen, exist)
            use precision
            use TypeDef
            implicit none

            ! Input
            integer, intent(in)                 :: sN, level, zcen
            integer, intent(in), dimension(:)   :: Zvec

            ! Input/Output
            type(streenode), intent(inout)      :: T
            logical, intent(inout)              :: exist
        end subroutine AddCount
        
        recursive subroutine SubtractCount(T, sN, level, Zvec, zcen)
            use precision
            use TypeDef
            implicit none

            ! Input
            integer, intent(in)                 :: sN, level, zcen
            integer, intent(in), dimension(:)   :: Zvec

            ! Input/Output
            type(streenode), intent(inout)      :: T
        end subroutine SubtractCount
        
        subroutine ExtendTree(sV, T) 
            use precision
            use TypeDef
            implicit none
            ! ExtendTree: Extends the search tree 'T' by adding a new node
            
            ! Input
            integer, intent(in)                 :: sV

            ! Input/Output
            type(streenode), intent(inout)      :: T
        end subroutine ExtendTree
       
        recursive subroutine ShapeTree(T, Ttrain, sV, sN)
            use precision
            use TypeDef
            implicit none

            ! Input
            type(streenode), intent(in)     :: Ttrain
            integer, intent(in)             :: sV, sN

            ! Input/Output
            type(streenode), intent(inout)  :: T
        end subroutine ShapeTree
        
        recursive subroutine CopyTree(Told, sV, sN, Tnew)
            use TypeDef
            use precision
            implicit none

            ! Input
            type(streenode), intent(in)         :: Told
            integer                             :: sV, sN    

            ! Output
            type(streenode), intent(inout)      :: Tnew
        end subroutine CopyTree
        
        recursive subroutine DeallocateTree(sV, T) 
            use TypeDef
            implicit none
            ! DeallocateTree: Deallocates recursively the search tree 'stree'
            ! starting from leaves to root.
            !
            ! Input
            integer, intent(in)                 :: sV

            ! Input/Output
            type(streenode), intent(inout)      :: T
        end subroutine DeallocateTree
        
    end interface
END MODULE Interfaces
