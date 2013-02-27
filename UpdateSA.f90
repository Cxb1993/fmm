subroutine UpdateSA(sV, sN, Zin, Zexin, Tin, Hin, Pin, Zout, Zexout, Tout, Hout, Pout)
    use Interfaces
    use precision
    use utils
    use TypeDef
    implicit none
    
    ! Input
    integer, intent(in)                                 :: sV, sN
    integer, intent(in), dimension(:,:,:)               :: Zin
    logical, intent(in), dimension(:,:,:)               :: Zexin
    type(streenode), intent(inout)                      :: Tin
    real(dp), intent(in), dimension(:,:)                :: Hin
    real(dp), intent(in), dimension(1:2)                :: Pin
    
    ! Output
    integer, intent(out), dimension(:,:,:), allocatable :: Zout
    logical, intent(out), dimension(:,:,:), allocatable :: Zexout
    type(streenode), intent(out)                        :: Tout
    real(dp), intent(out), dimension(:,:), allocatable  :: Hout
    real(dp), intent(out), dimension(1:2)               :: Pout
    
    ! Local variables
    integer                                             :: nZ, mZ, pZ, nh
    
    ! Dimensions
    nZ = size(Zin,1); mZ = size(Zin,2); pZ = size(Zin,3)
    nh = size(Hin,2)
    
    ! The image
    if (allocated(Zout)) deallocate(Zout)
    allocate(Zout(nZ,mZ,pz), source = Zin)
    
    if (allocated(Zexout)) deallocate(Zexout)
    allocate(Zexout(nZ,mZ,pZ), source = Zexin)
    
    ! The tree
    call CopyTree(Tin, sV, sN, Tout)
        
    ! The histogram
    if (allocated(Hout)) deallocate(Hout)
    allocate(Hout(sV+1,nh), source = Hin)
    
    ! The misfit function value
    Pout = Pin
    
end subroutine UpdateSA
