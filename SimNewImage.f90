!     
! File:   NewImage.f03
! Author: M.Sc. Katrine Lange
!
! Created on August 5, 2011, 13:01 PM
!
subroutine SimNewImage(Ztest,Zcond,Zex,Ttest,Htest,Ttrain,i,j,k,nodes,sV,sN,Dblock,Znew,Zexnew,Tnew,Hnew,newImage)
    use precision
    use utils
    use Interfaces, only : getNeighborhood, SimVoxel, GrowTree, Tree2Hist
    use TypeDef
    implicit none
    
    ! Input
    integer, intent(in), dimension(:,:,:)               :: Ztest
    logical, intent(in), dimension(:,:,:)               :: Zcond, Zex
    integer, intent(in)                                 :: i, j, k, sV, sN
    integer, intent(in), dimension(:,:)                 :: nodes
    type(DomainMask), intent(in)                        :: Dblock
    type(streenode), intent(inout)                      :: Ttest, Ttrain    
    real(dp), intent(in), dimension(:,:)                :: Htest
          
    ! Output
    integer, intent(out), allocatable, dimension(:,:,:) :: Znew
    logical, intent(out), allocatable, dimension(:,:,:) :: Zexnew
    type(streenode), intent(out)                        :: Tnew
    real(dp), intent(out), allocatable, dimension(:,:)  :: Hnew
    logical, intent(out)                                :: newImage
    
    ! Local variables
    integer                                             :: nZ, mZ, pZ, ii, jj, kk
    integer                                             :: sD, v, w
    integer                                             :: rindex
    real(dp)                                            :: rand, ri
    real(dp), allocatable, dimension(:)                 :: randvec, cpdf, hc
    integer, allocatable, dimension(:)                  :: ivec 
    integer, dimension(sN)                              :: Zvec, D
    integer, allocatable, dimension(:,:,:)              :: Ztmp
    integer, allocatable, dimension(:,:)                :: dnodes, dindex, dindextmp
    logical                                             :: innervoxel
    
    ! Dimensions
    nZ = size(Ztest, 1); mZ = size(Ztest, 2); pZ = size(Ztest, 3)
    dnodes = Dblock%nodes; D = Dblock%mat
    sD = size(dnodes,1)
    
    ! Initialization
    call allocAs(Ztest, Ztmp)
    
    ! Part of the image to be simulated
    allocate(dindextmp(sD,3), source = 0)
    w = 0
    do v = 1, sD
        ii = i + dnodes(v,1); jj = j + dnodes(v,2); kk = k + dnodes(v,3) 
        
        if (ii >= 1 .and. ii <= nZ .and. &
            jj >= 1 .and. jj <= mZ .and. &
            kk >= 1 .and. kk <= pZ) then
            
            w = w + 1
            dindextmp(w,:) = [ii, jj, kk]
            Ztmp(ii,jj,kk) = -1
        end if
    end do    
    allocate(dindex(w,3), source = dindextmp(1:w,:))    
    
    ! Set random path for simulating the voxels    
    allocate(randvec(w), ivec(w))
    ivec = [(v, v = 1, w)]
    call random_number(randvec)    
    call disort(randvec, ivec, w, 2)    
    dindex = dindex(ivec,:)
    
    ! Simulate each voxel at a time
    do v = 1, w
        ii = dindex(v,1); jj = dindex(v,2); kk = dindex(v,3)
    
        ! If the voxel is not conditioned upon
        if (.NOT. Zcond(ii,jj,kk)) then

            ! Extract neighborhood
            call getNeighborhood(Ztmp, ii, jj, kk, nodes, sN, Zvec, innervoxel)
                              
            ! Get the conditional PDF           
            call SimVoxel(Zvec, Ttrain, sV, sN, D, hc)

            ! Pick with this probability
            call random_number(rand)
            ri = hc(sV+1)*rand
            call argfind(ri < hc, rindex)              
            
            ! Change voxel
            Ztmp(ii,jj,kk) = rindex-1  
        else
            Ztmp(ii,jj,kk) = Ztest(ii,jj,kk)
        end if        
    end do
    
    ! -------------------- UPDATE THE HISTOGRAM OF EACH CHANGED VOXEL ---------------------
    
    ! Initialization
    call allocAs(Ztest, Znew)
    allocate(Zexnew(nZ,mZ,pZ), source = Zex)
    call CopyTree(Ttest, sV, sN, Tnew)
    
    ! Update tree for each voxel
    newImage = .false.
    do v = 1, w
                
        ii = dindex(v,1); jj = dindex(v,2); kk = dindex(v,3) 
        
        Znew(ii,jj,kk) = Ztmp(ii,jj,kk)

        ! Only necessary if voxel has changed
        if (Znew(ii,jj,kk) /= Ztest(ii,jj,kk)) then
            newImage = .true.
            
            ! Update the tree    
            call GrowTree(Znew, Ztest(ii,jj,kk), Zexnew, ii, jj, kk, nodes, sV, sN, Ttrain, Tnew)
        end if 
                
    end do
    
    ! Compute the histogram from the new tree
    if (newImage) then
        call Tree2Hist(sV, sN, Tnew, Hnew)
    else
        call allocAs(Htest, Hnew)
    end if
        
end subroutine SimNewImage