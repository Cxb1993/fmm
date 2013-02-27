!     
! File:   utils.f03
! Author: jf
!
! Created on March 29, 2011, 2:41 PM
!
module utils
    !----------------------------------------------------------------------
    ! Description
    !----------------------------------------------------------------------
    ! Various help functions to be used together with fortran 90 and up.
    ! This package is put together by Jan Frydendall and Carsten V\"olcker.
    !
    ! ---------------------------------------------------------------------
    ! Created by
    ! ---------------------------------------------------------------------
    ! Author
    !   Jan Frydendall, Carsten V\"olcker
    ! Address
    !   Department of Informatics and Mathematical Modelling
    !   Technical University of Denmark
    !   DTU - Building 321
    !   Richard Petersens Plads
    !   DK-2800 Lyngby
    ! Email
    !   {jf ,cv}(at)imm.dtu.dk
    ! Date
    !   December 15, 2010
    ! ---------------------------------------------------------------------
    ! Initialize
    ! ---------------------------------------------------------------------
    use precision
    interface find
        module procedure find1, find2, find3dp, find3int
    end interface
    interface arange
        module procedure arange1
    end interface
    interface flatten
        module procedure flatten_int, flatten_dp
    end interface
    interface allocAs
        module procedure allocAs3dp, allocAs3int, allocAs2int, allocAs2dp, allocAs1int, allocAs1dp
    end interface
    interface cumsum
        module procedure cumsum_dp, cumsum_int
    end interface
    interface mean
        module procedure mean_dp, mean_int
    end interface
    interface variance
        module procedure variance_dp, variance_int
    end interface
    interface ravel
        module procedure ravel1, ravel2, ravel1dp, ravel2dp, ravel1log, ravel2log
    end interface
    interface cumprod
       module procedure cumprod1,cumprod2
    end interface
    interface printarray
        module procedure printarray1int,printarray1dble,printarray1logic,printarray2int,printarray2dble,printarray2logic &
                        ,printarray3int
    end interface
    
    interface
        subroutine DSORT(DX,DY,N,KFLAG)
            use precision
            implicit none
            real(dp),dimension(:),intent(inout) :: DX
            real(dp),dimension(:),intent(inout) :: DY
            integer,intent(in)                  :: N,KFLAG
        end subroutine DSORT
    end interface

    interface
        subroutine ISORT(DX,DY,N,KFLAG)
            implicit none
            integer,dimension(:),intent(inout) :: DX
            integer,dimension(:),intent(inout) :: DY
            integer,intent(in)                 :: N,KFLAG
        end subroutine ISORT
    end interface

    interface
        subroutine DISORT(DX,DY,N,KFLAG)
            use precision
        implicit none
            real(dp),dimension(:),intent(inout) :: DX
            integer,dimension(:),intent(inout)  :: DY
            integer,intent(in)                  :: N,KFLAG
        end subroutine DISORT
    end interface
    
    contains
    
    subroutine multMatVec(A, x, b) 
        use precision
        implicit none
        ! input
        real(dp), intent(in), dimension(:,:)  :: A
        real(dp), intent(in), dimension(:)    :: x
        real(dp), intent(out), dimension(:)   :: b
        ! local variables
        integer                             :: n, m

        n = size(A,1)
        m = size(A,2)

        call DGEMV('N', n, m, 1.0_dp, A, n, x, 1, 0.0_dp, b, 1)    
    end subroutine multMatVec
       
    function cumsum_dp(x) result(c)
        ! Calculates the cummulative sum of a 1 dimensional array
        implicit none
        real(dp), intent(in), dimension(:) :: x
        real(dp), dimension(1:size(x))     :: c
        integer                            :: i
        
        c = x
        do i = 2, size(x)
            c(i) = c(i-1)+x(i)
        end do
    end function cumsum_dp
    
    function cumsum_int(x) result(c)
        ! Calculates the cummulative sum of a 1 dimensional array
        implicit none
        integer, intent(in), dimension(:)  :: x
        integer, dimension(1:size(x))      :: c
        integer                            :: i
        
        c = x
        do i = 2, size(x)
            c(i) = c(i-1)+x(i)
        end do
    end function cumsum_int
    

    real(dp) function mean_dp(x)
        ! Calculates the mean of an 1 dimensional array
        implicit none
        real(dp), intent(in), dimension(:) :: x
        mean_dp = sum(x)
        mean_dp = mean_dp/size(x)
    end function mean_dp

    integer function mean_int(x)
        ! Calculates the mean of an 1 dimensional array
        implicit none
        integer, intent(in), dimension(:) :: x
        mean_int = sum(x)
        mean_int = nint(mean_int/real(size(x),kind=dp))
    end function mean_int
    
    real(dp) function variance_dp(x)
        ! Calculates the mean of an 1 dimensional array
        implicit none
        real(dp), intent(in), dimension(:) :: x
        real(dp)                           :: mean
        variance_dp = mean_dp(x**2) - mean_dp(x)**2
    end function variance_dp

    real(dp) function variance_int(x)
        ! Calculates the mean of an 1 dimensional array
        implicit none
        integer, intent(in), dimension(:) :: x
        variance_int = mean_dp(real(x,kind = dp)**2) - real(mean_int(x)**2,kind = dp)
    end function variance_int
    
    subroutine argfind(mask, indx)
        ! ---------------------------------------------------------------------
        ! Description
        ! ---------------------------------------------------------------------
        ! Find indices of true elements in a one dimensional logical array,
        ! e.g. mask is a logical expresion like A > 3, where A is a one
        ! dimensional array.
        ! Input
        !   mask (logical, dimension(1:n)) : One dimensional logical array.
        ! Output
        !   indx (integer, dimension(1:n)) : Indices of true elements in mask.
        ! ---------------------------------------------------------------------
        ! Created by
        ! ---------------------------------------------------------------------
        ! Author
        !   Carsten V\"olcker, Jan Frydendall
        ! Address
        !   Department of Informatics and Mathematical Modelling
        !   Technical University of Denmark
        !   DTU - Building 321
        !   Richard Petersens Plads
        !   DK-2800 Lyngby
        ! Email
        !   {cv,jf}(at)imm.dtu.dk
        ! Date
        !   December 15, 2010
        ! ---------------------------------------------------------------------
        ! Initialize
        ! ---------------------------------------------------------------------
        implicit none
        ! Input:
        logical, intent(in), dimension(:) :: mask
        ! Output:
        integer, intent(out)              :: indx
        ! Local variables:
        integer :: i, m
        ! ---------------------------------------------------------------------
        ! Main
        ! ---------------------------------------------------------------------
        ! Count no. of elements in mask:
        m = size(mask) ! No. of elements in mask.
        ! Assume all elements in mask are true and allocate memory (row index
        ! i = indx(:,1), column index j = indx(:,2)):
        ! Extract indices of true elements by looping over elements of mask:
        do i = 1, m
            if (mask(i)) then
                indx = i
                return
            end if
        end do
        indx = 0
    end subroutine argfind
    
    subroutine find1(mask, indx)
        ! ---------------------------------------------------------------------
        ! Description
        ! ---------------------------------------------------------------------
        ! Find indices of true elements in a one dimensional logical array,
        ! e.g. mask is a logical expresion like A > 3, where A is a one
        ! dimensional array.
        ! Input
        !   mask (logical, dimension(1:n)) : One dimensional logical array.
        ! Output
        !   indx (integer, dimension(1:n)) : Indices of true elements in mask.
        ! ---------------------------------------------------------------------
        ! Created by
        ! ---------------------------------------------------------------------
        ! Author
        !   Carsten V\"olcker, Jan Frydendall
        ! Address
        !   Department of Informatics and Mathematical Modelling
        !   Technical University of Denmark
        !   DTU - Building 321
        !   Richard Petersens Plads
        !   DK-2800 Lyngby
        ! Email
        !   {cv,jf}(at)imm.dtu.dk
        ! Date
        !   December 15, 2010
        ! ---------------------------------------------------------------------
        ! Initialize
        ! ---------------------------------------------------------------------
        implicit none
        ! Input:
        logical, intent(in), dimension(:) :: mask
        ! Output:
        integer, intent(inout), dimension(:), allocatable :: indx
        ! Local variables:
        integer :: i, k, m
        integer, dimension(:), allocatable :: temp
        ! ---------------------------------------------------------------------
        ! Main
        ! ---------------------------------------------------------------------
        ! Count no. of elements in mask:
        m = size(mask) ! No. of elements in mask.
        ! Assume all elements in mask are true and allocate memory (row index
        ! i = indx(:,1), column index j = indx(:,2)):
        allocate(temp(1:m))
        ! Extract indices of true elements by looping over elements of mask:
        k = 0
        do i = 1, m
            if (mask(i)) then
                k = k + 1
                temp(k) = i
            end if
        end do
        ! Avoid memory leak by rellaocating indx for new size and new content,
        ! this is why indx is declared intent(inout) and not only intent(out):
        if (allocated(indx)) deallocate(indx)
        allocate(indx(1:k))
        indx = temp(1:k)
        deallocate(temp)
    end subroutine find1

    subroutine find2(mask, indx)
        ! ---------------------------------------------------------------------
        ! Description
        ! ---------------------------------------------------------------------
        ! Find row and column indices of true elements in a two dimensional
        ! logical array, e.g. mask is a logical expresion like A > 3, where A
        ! is a two dimensional array.
        ! Input
        !   mask (logical, dimension(1:n,1:2)) : Two dimensional logical array.
        ! Output
        !   indx (integer, dimension(1:n,1:2)) : Indices of true elements in
        !       mask, where indx(:,1) denotes row indices and indx(:,2) denotes
        !       column indices.
        ! ---------------------------------------------------------------------
        ! Created by
        ! ---------------------------------------------------------------------
        ! Author
        !   Carsten V\"olcker, Jan Frydendall
        ! Address
        !   Department of Informatics and Mathematical Modelling
        !   Technical University of Denmark
        !   DTU - Building 321
        !   Richard Petersens Plads
        !   DK-2800 Lyngby
        ! Email
        !   {cv,jf}(at)imm.dtu.dk
        ! Date
        !   December 15, 2010
        ! ---------------------------------------------------------------------
        ! Initialize
        ! ---------------------------------------------------------------------
        implicit none
        ! Input:
        logical, intent(in), dimension(:,:) :: mask
        ! Output:
        integer, intent(inout), dimension(:,:), allocatable :: indx
        ! Local variables:
        integer :: i, j, k, m, n
        integer, dimension(:,:), allocatable :: temp
        ! ---------------------------------------------------------------------
        ! Main
        ! ---------------------------------------------------------------------
        ! Count no. of rows and columns in mask:
        m = size(mask, 1) ! No. of rows in mask.
        n = size(mask, 2) ! No. of columns in mask.
        ! Assume all elements in mask are true and allocate memory (row index
        ! i = indx(:,1), column index j = indx(:,2)):
        allocate(temp(1:m * n, 1:2))
        ! Extract indices of true elements by looping over rows and and columns
        ! of mask:
        k = 0
        do j = 1, n
            do i = 1, m
                if (mask(i, j)) then
                    k = k + 1
                    temp(k,:) = [i, j]
                end if
            end do
        end do
        ! Avoid memory leak by rellaocating indx for new size and new content,
        ! this is why indx is declared intent(inout) and not only intent(out):
        if (allocated(indx)) deallocate(indx)
        allocate(indx(1:k, 1:2))
        indx = temp(1:k,:)
        deallocate(temp)
    end subroutine find2

    subroutine find3dp(a, mask, indx)
        ! ---------------------------------------------------------------------
        ! Description
        ! ---------------------------------------------------------------------
        ! Find indices of true elements in a one dimensional logical array,
        ! e.g. mask is a logical expresion like A > 3, where A is a one
        ! dimensional array.
        ! Input
        !   mask (logical, dimension(1:n)) : One dimensional logical array.
        ! Output
        !   indx (integer, dimension(1:n)) : array indx contains the indices of
        !                                    true elements in mask mapped trough array a
        ! ---------------------------------------------------------------------
        ! Created by
        ! ---------------------------------------------------------------------
        ! Author
        !   Carsten V\"olcker, Jan Frydendall
        ! Address
        !   Department of Informatics and Mathematical Modelling
        !   Technical University of Denmark
        !   DTU - Building 321
        !   Richard Petersens Plads
        !   DK-2800 Lyngby
        ! Email
        !   {cv,jf}(at)imm.dtu.dk
        ! Date
        !   December 15, 2010
        ! ---------------------------------------------------------------------
        ! Initialize
        ! ---------------------------------------------------------------------
        implicit none
        ! Input:
        real(dp), intent(in), dimension(:) :: a
        logical, intent(in), dimension(:) :: mask
        ! Output:
        integer, intent(inout), dimension(:), allocatable :: indx
        ! Local variables:
        integer :: i, k, m
        integer, dimension(:), allocatable :: temp
        ! ---------------------------------------------------------------------
        ! Main
        ! ---------------------------------------------------------------------
        ! Count no. of elements in mask:
        m = size(mask) ! No. of elements in mask.
        ! Assume all elements in mask are true and allocate memory (row index
        ! i = indx(:,1), column index j = indx(:,2)):
        allocate(temp(1:m))
        ! Extract indices of true elements by looping over elements of mask:
        k = 0
        do i = 1, m
            if (mask(i)) then
                k = k + 1
                temp(k) = i
            end if
        end do
        ! Avoid memory leak by rellaocating indx for new size and new content,
        ! this is why indx is declared intent(inout) and not only intent(out):
        if (allocated(indx)) deallocate(indx)
        allocate(indx(1:k))
        indx = a(temp(1:k))
        deallocate(temp)
    end subroutine find3dp

    subroutine find3int(a, mask, indx)
        ! ---------------------------------------------------------------------
        ! Description
        ! ---------------------------------------------------------------------
        ! Find indices of true elements in a one dimensional logical array,
        ! e.g. mask is a logical expresion like A > 3, where A is a one
        ! dimensional array.
        ! Input
        !   mask (logical, dimension(1:n)) : One dimensional logical array.
        ! Output
        !   indx (integer, dimension(1:n)) : array indx contains the indices of 
        !                                    true elements in mask mapped trough array a
        ! ---------------------------------------------------------------------
        ! Created by
        ! ---------------------------------------------------------------------
        ! Author
        !   Carsten V\"olcker, Jan Frydendall
        ! Address
        !   Department of Informatics and Mathematical Modelling
        !   Technical University of Denmark
        !   DTU - Building 321
        !   Richard Petersens Plads
        !   DK-2800 Lyngby
        ! Email
        !   {cv,jf}(at)imm.dtu.dk
        ! Date
        !   December 15, 2010
        ! ---------------------------------------------------------------------
        ! Initialize
        ! ---------------------------------------------------------------------
        implicit none
        ! Input:
        integer, intent(in), dimension(:) :: a
        logical, intent(in), dimension(:) :: mask
        ! Output:
        integer, intent(inout), dimension(:), allocatable :: indx
        ! Local variables:
        integer :: i, k, m
        integer, dimension(:), allocatable :: temp
        ! ---------------------------------------------------------------------
        ! Main
        ! ---------------------------------------------------------------------
        ! Count no. of elements in mask:
        m = size(mask) ! No. of elements in mask.
        ! Assume all elements in mask are true and allocate memory (row index
        ! i = indx(:,1), column index j = indx(:,2)):
        allocate(temp(1:m))
        ! Extract indices of true elements by looping over elements of mask:
        k = 0
        do i = 1, m
            if (mask(i)) then
                k = k + 1
                temp(k) = i
            end if
        end do
        ! Avoid memory leak by rellaocating indx for new size and new content,
        ! this is why indx is declared intent(inout) and not only intent(out):
        if (allocated(indx)) deallocate(indx)
        allocate(indx(1:k))
        indx = a(temp(1:k))
!        call find(a(temp(1:k))>0, indx)
        deallocate(temp)
    end subroutine find3int
    
    recursive function binarySearch_R(a, value) result (bsresult)
        real(dp), intent(in) :: a(:), value
        integer :: bsresult, mid
        mid = size(a)/2 + 1
        if (size(a) == 0) then
            bsresult = 0 ! not found
        else if (a(mid) > value) then
            bsresult = binarySearch_R(a(:mid - 1), value)
        else if (a(mid) < value) then
            bsresult = binarySearch_R(a(mid + 1:), value)
            if (bsresult /= 0) then
            bsresult = mid + bsresult
            end if
        else
            bsresult = mid ! SUCCESS!!
        end if
    end function binarySearch_R   

    subroutine arange1(lowerbound,upperbound,a)
        ! ---------------------------------------------------------------------
        ! Initialize
        ! ---------------------------------------------------------------------
        implicit none
        ! Input:
        integer,intent(in)                               :: lowerbound, upperbound
        ! Output:
        integer,intent(inout),dimension(:),allocatable     :: a
        ! Local variables:
        integer,dimension(:),allocatable :: move
        integer,dimension(:,:),allocatable :: temp
        integer :: i,k,m,n
        ! ---------------------------------------------------------------------
        ! Main
        ! ---------------------------------------------------------------------
        m = upperbound - lowerbound + 1
        if (allocated(a)) deallocate(a)
        allocate(a(1:m))
        do i = 1 , m
            a(i) = lowerbound + (i-1)
        end do
    end subroutine arange1

    function ravel1(array) result(a)
        ! ---------------------------------------------------------------------
        ! Initialize
        ! ---------------------------------------------------------------------
        implicit none
        ! Input:
        integer,intent(in), dimension(:,:)               :: array
        ! Output:
        integer,dimension(1:size(array))                 :: a
        ! ---------------------------------------------------------------------
        ! Main
        ! ---------------------------------------------------------------------
        a = reshape(array, shape=[size(array)])
    end function ravel1

    function ravel2(array) result(a)
        ! ---------------------------------------------------------------------
        ! Initialize
        ! ---------------------------------------------------------------------
        implicit none
        ! Input:
        integer,intent(in), dimension(:,:,:)             :: array
        ! Output:
        integer ,dimension(1:size(array))                :: a
        ! ---------------------------------------------------------------------
        ! Main
        ! ---------------------------------------------------------------------
        a = reshape(array, shape=[size(array)])
    end function ravel2

    function ravel1dp(array) result(a)
        ! ---------------------------------------------------------------------
        ! Initialize
        ! ---------------------------------------------------------------------
        implicit none
        ! Input:
        real(dp),intent(in), dimension(:,:)               :: array
        ! Output:
        real(dp),dimension(1:size(array))                 :: a
        ! ---------------------------------------------------------------------
        ! Main
        ! ---------------------------------------------------------------------
        a = reshape(array, shape=[size(array)])
    end function ravel1dp

    function ravel2dp(array) result(a)
        ! ---------------------------------------------------------------------
        ! Initialize
        ! ---------------------------------------------------------------------
        implicit none
        ! Input:
        real(dp),intent(in), dimension(:,:,:)             :: array
        ! Output:
        real(dp) ,dimension(1:size(array))                :: a
        ! ---------------------------------------------------------------------
        ! Main
        ! ---------------------------------------------------------------------
        a = reshape(array, shape=[size(array)])
    end function ravel2dp

    function ravel1log(array) result(a)
        ! ---------------------------------------------------------------------
        ! Initialize
        ! ---------------------------------------------------------------------
        implicit none
        ! Input:
        logical,intent(in), dimension(:,:)               :: array
        ! Output:
        logical,dimension(1:size(array))                 :: a
        ! ---------------------------------------------------------------------
        ! Main
        ! ---------------------------------------------------------------------
        a = reshape(array, shape=[size(array)])
    end function ravel1log

    function ravel2log(array) result(a)
        ! ---------------------------------------------------------------------
        ! Initialize
        ! ---------------------------------------------------------------------
        implicit none
        ! Input:
        logical,intent(in), dimension(:,:,:)             :: array
        ! Output:
        logical ,dimension(1:size(array))                :: a
        ! ---------------------------------------------------------------------
        ! Main
        ! ---------------------------------------------------------------------
        a = reshape(array, shape=[size(array)])
    end function ravel2log

    subroutine flatten_int(a,b,c)
        ! ---------------------------------------------------------------------
        ! Initialize
        ! ---------------------------------------------------------------------
        implicit none
        ! Input:
        integer,intent(in), dimension(:)             :: a,b
        ! Output:
        integer, intent(out), allocatable ,dimension(:)                :: c
        ! ---------------------------------------------------------------------
        ! Main
        ! ---------------------------------------------------------------------
        if (allocated(c)) deallocate(c)
        allocate(c(1:size(b)))
        c = a(b)
    end subroutine flatten_int

    subroutine flatten_dp(a,b,c)
        ! ---------------------------------------------------------------------
        ! Initialize
        ! ---------------------------------------------------------------------
        implicit none
        ! Input:
        real(dp),intent(in), dimension(:)             :: a
        integer,intent(in), dimension(:)              :: b
        ! Output:
        real(dp), intent(out), allocatable ,dimension(:)                :: c
        ! ---------------------------------------------------------------------
        ! Main
        ! ---------------------------------------------------------------------
        if (allocated(c)) deallocate(c)
        allocate(c(1:size(b)))
        c = a(b)
    end subroutine flatten_dp

    subroutine ind2sub(siz, idx, sub)
        ! ---------------------------------------------------------------------
        ! Initialize
        ! ---------------------------------------------------------------------
        implicit none
        ! Input:
        integer,intent(in), dimension(:)                 :: siz
        integer,intent(in)                               :: idx
        ! Output:
        integer,intent(inout),dimension(:),allocatable   :: sub
        ! Local variables:
        integer                                          :: i, n
        real(dp)                                         :: ndx, vi, vj
        real(dp), dimension(1:size(siz))                 :: siz1, k
        ! ---------------------------------------------------------------------
        ! Main
        ! ---------------------------------------------------------------------
        n = size(siz)
        siz1 = real(siz, kind=dp)
        ndx = real(idx,kind=dp)
        k = [1.0_dp, cumprod(siz1(1:n-1))]
        if (allocated(sub)) deallocate(sub)
        allocate(sub(1:n))
        do i = n, 1, -1
          vi = mod(ndx - 1.0_dp, k(i)) + 1.0_dp
          vj = (ndx - vi)/k(i) + 1.0_dp
          sub(i) = int(vj)
          ndx = vi
        end do
    end subroutine ind2sub

    function cumprod1(x)
        implicit none
        integer :: i, n
        integer, intent(in), dimension(:) :: x
        integer, allocatable, dimension(:) :: cumprod1
        n = size(x, 1)
        allocate(cumprod1(1:n))
        cumprod1 = 0
        cumprod1(1) = x(1)
        do i = 2, n
            cumprod1(i) = cumprod1(i - 1) * x(i)
        end do
    end function cumprod1

    function cumprod2(x)
        implicit none
        integer :: i, n
        real(dp), intent(in), dimension(:) :: x
        real(dp), allocatable, dimension(:) :: cumprod2
        n = size(x, 1)
        allocate(cumprod2(1:n))
        cumprod2 = 0.0_dp
        cumprod2(1) = x(1)
        do i = 2, n
            cumprod2(i) = cumprod2(i - 1) * x(i)
        end do
    end function cumprod2

    subroutine allocAs3dp(source, target)
        implicit none
        real(dp), intent(in), dimension(:,:,:)                 :: source
        real(dp), intent(inout), allocatable, dimension(:,:,:) :: target
        integer                                                :: m, n, o
        m = size(source,1);n = size(source,2);o = size(source,3)
        if (allocated(target)) deallocate(target)
        allocate(target(1:m,1:n,1:o),source = source)
    end subroutine allocAs3dp

   subroutine allocAs3int(source, target)
        implicit none
        integer, intent(in), dimension(:,:,:)                 :: source
        integer, intent(inout), allocatable, dimension(:,:,:) :: target
        integer                                               :: m, n, o
        m = size(source,1);n = size(source,2);o = size(source,3)
        if (allocated(target)) deallocate(target)
        allocate(target(1:m,1:n,1:o),source = source)
    end subroutine allocAs3int

    subroutine allocAs2int(source, target)
        implicit none
        integer, intent(in), dimension(:,:)                   :: source
        integer, intent(inout), allocatable, dimension(:,:)   :: target
        integer                                               :: m, n
        m = size(source,1);n = size(source,2)
        if (allocated(target)) deallocate(target)
        allocate(target(1:m,1:n),source = source)
    end subroutine allocAs2int

    subroutine allocAs2dp(source, target)
        implicit none
        real(dp), intent(in), dimension(:,:)                   :: source
        real(dp), intent(inout), allocatable, dimension(:,:)   :: target
        integer                                                :: m, n
        m = size(source,1);n = size(source,2)
        if (allocated(target)) deallocate(target)
        allocate(target(1:m,1:n),source = source)
    end subroutine allocAs2dp

    subroutine allocAs1int(source, target)
        implicit none
        integer, intent(in), dimension(:)                   :: source
        integer, intent(inout), allocatable, dimension(:)   :: target
        integer                                             :: m
        m = size(source,1)
        if (allocated(target)) deallocate(target)
        allocate(target(1:m),source = source)
    end subroutine allocAs1int

    subroutine allocAs1dp(source, target)
        implicit none
        real(dp), intent(in), dimension(:)                   :: source
        real(dp), intent(inout), allocatable, dimension(:)   :: target
        integer                                              :: m
        m = size(source,1)
        if (allocated(target)) deallocate(target)
        allocate(target(1:m),source = source)
    end subroutine allocAs1dp

subroutine printarray1int(a)
        ! ---------------------------------------------------------------------
        ! Initialize
        ! ---------------------------------------------------------------------
        implicit none
        ! Input:
        integer,intent(in),dimension(:) :: a
        integer :: i
        ! ---------------------------------------------------------------------
        ! Main
        ! ---------------------------------------------------------------------
        do i = 1,size(a)
            print*,a(i)
        end do
    end subroutine printarray1int

    subroutine printarray1dble(a)
        ! ---------------------------------------------------------------------
        ! Initialize
        ! ---------------------------------------------------------------------
        implicit none
        ! Input:
        real(dp),intent(in),dimension(:) :: a
        integer :: i
        ! ---------------------------------------------------------------------
        ! Main
        ! ---------------------------------------------------------------------
        do i = 1,size(a)
            print*,a(i)
        end do
    end subroutine printarray1dble

    subroutine printarray1logic(a)
        ! ---------------------------------------------------------------------
        ! Initialize
        ! ---------------------------------------------------------------------
        implicit none
        ! Input:
        logical,intent(in),dimension(:) :: a
        integer :: i
        ! ---------------------------------------------------------------------
        ! Main
        ! ---------------------------------------------------------------------
        do i = 1,size(a)
            print*,a(i)
        end do
    end subroutine printarray1logic

    subroutine printarray2int(a)
        ! ---------------------------------------------------------------------
        ! Initialize
        ! ---------------------------------------------------------------------
        implicit none
        ! Input:
        integer,intent(in),dimension(:,:) :: a
        integer :: i
        ! ---------------------------------------------------------------------
        ! Main
        ! ---------------------------------------------------------------------
        do i = 1,size(a,1)
            print*,a(i,:)
        end do
    end subroutine printarray2int

    subroutine printarray2dble(a)
        ! ---------------------------------------------------------------------
        ! Initialize
        ! ---------------------------------------------------------------------
        implicit none
        ! Input:
        real(dp),intent(in),dimension(:,:) :: a
        integer :: i
        ! ---------------------------------------------------------------------
        ! Main
        ! ---------------------------------------------------------------------
        do i = 1,size(a,1)
            print*,a(i,:)
        end do
    end subroutine printarray2dble

    subroutine printarray2logic(a)
        ! ---------------------------------------------------------------------
        ! Initialize
        ! ---------------------------------------------------------------------
        implicit none
        ! Input:
        logical,intent(in),dimension(:,:) :: a
        integer :: i
        ! ---------------------------------------------------------------------
        ! Main
        ! ---------------------------------------------------------------------
        do i = 1,size(a,1)
            print*,a(i,:)
        end do
    end subroutine printarray2logic

    subroutine printarray3int(a)
        ! ---------------------------------------------------------------------
        ! Initialize
        ! ---------------------------------------------------------------------
        implicit none
        ! Input:
        integer,intent(in),dimension(:,:,:) :: a
        integer :: i
        ! ---------------------------------------------------------------------
        ! Main
        ! ---------------------------------------------------------------------
        do i = 1,size(a,3)
            print *, 'a( : , : ,',i,')'
            call printarray2int(a(:,:,i))
        end do
    end subroutine printarray3int
end module utils