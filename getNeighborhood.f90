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
    
    ! Local variables
    integer                                         :: nZ, mZ, pZ, v, ii, jj, kk
        
    ! Initialize
    nZ = size(Z,1); mZ = size(Z,2); pZ = size(Z,3)    
    Zvec = -1
    innervoxel = .true.

    ! For each neighbor   
    do v = 1, sN
        
        ii = i + nodes(v,1)
        jj = j + nodes(v,2)
        kk = k + nodes(v,3)
    
        if (ii >= 1 .and. ii <= nZ .and. &
            jj >= 1 .and. jj <= mZ .and. &
            kk >= 1 .and. kk <= pZ) then
        
            Zvec(v) = Z(ii,jj,kk) 
        else
            innervoxel = .false.
        end if
    end do   
end subroutine getNeighborhood
