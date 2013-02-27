subroutine CompObjFun(H,Htrain,N,Z,invProb,alpha,P)
    use interfaces, only : CompChiDist, CompDataFit
    use TypeDef
    use precision
    use utils
    implicit none
    
    ! Input
    real(dp), intent(in), dimension(:,:)            :: H, Htrain   
    integer, intent(in)                             :: N
    integer, intent(in), dimension(:,:,:)           :: Z
    type(inverseproblem), intent(in)                :: invProb
    real(dp), intent(in)                            :: alpha
    
    ! Output
    real(dp), intent(out), dimension(1:2)           :: P
    
    ! Local Variables    
    real(dp)                                        :: X, L
      
    ! Compute Chi square distance
    X = CompChiDist(ravel(H), ravel(Htrain), N)        
   
    ! Compute data fit
    if (alpha > 0) then
        L = CompDataFit(Z, invProb) 
    else 
        L = 0.00_dp
    end if
        
    ! Objective function
    P(1) = alpha * X
    P(2) = L 
        
end subroutine CompObjFun