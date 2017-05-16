module modLinearAlgebra



  implicit NONE



  type QR_factorization 
     
     ! -------------------------------------------------------------------------------
     ! mat contains the QR factorization computed with LAPACK ZGEQRF
     ! tau contains scalar to define elementary reflectors 
     ! -------------------------------------------------------------------------------

     complex(8),dimension(:,:),allocatable :: mat

     complex(8),dimension(:),allocatable :: tau

     
  end type QR_factorization



  type SVD_factorization

     ! --------------------------------------------------------------------------------
     ! mat_shape is a two dimensional vector that indicates matrix shape
     ! Sigma contains the singular values
     ! U contains the left singular vectors ( one per column)
     ! VT contains the right singular vectors ( one per column )
     ! --------------------------------------------------------------------------------

     integer, dimension(2) :: mat_shape

     complex (8), dimension(:,:), allocatable :: U, Vt

     real ( 8 ), dimension( : ), allocatable :: Sigma

     complex (8), dimension(:,:), allocatable :: pseudo_inverse


  end type SVD_factorization    



contains



  subroutine createQR_factorization ( this, m, n )

    ! --------------------------------------------------------------------------------
    ! Allocates QR matrices and reflectors of a QR_factorization type. The matrix has
    ! m rows and n columns.
    ! --------------------------------------------------------------------------------

    type (QR_factorization) :: this

    integer :: m, n


    allocate ( this % mat ( m, n ) )

    allocate ( this % tau ( min (m, n) ) )


  end subroutine createQR_factorization



  subroutine createSVD_factorization ( this, mat, tol )

    ! --------------------------------------------------------------------------------
    ! Allocates matrices to store matrices U, Sigma and Vt and calls a subroutine to
    ! compute the SVD factorization of mat
    ! --------------------------------------------------------------------------------

    type ( SVD_factorization ) :: this

    complex (8), dimension (:,:) :: mat

    real (8) :: tol

    
    integer :: r


    this % mat_shape = shape ( mat )

    r = minval ( this % mat_shape )


    allocate ( this % U ( this % mat_shape (1) , this % mat_shape (1) ) )

    allocate ( this % Vt ( this % mat_shape (2), this % mat_shape (2) ) )
    
    allocate ( this % Sigma ( r ) )

    allocate ( this % pseudo_inverse ( this % mat_shape (2), this % mat_shape(1) ) )


    call computeSvd_factorization ( this, mat )
    
    call compute_Pseudoinverse ( this, tol )
    

  end subroutine createSVD_factorization


  
  subroutine destroyQR_factorization ( this )

    ! --------------------------------------------------------------------------------
    ! Dealloactes matrix and reflector from QR_factorization type
    ! --------------------------------------------------------------------------------

    type ( QR_factorization ) :: this

    
    deallocate ( this % mat )

    deallocate ( this % tau )


  end subroutine destroyQR_factorization



  subroutine destroySVD_factorization ( this )

    ! --------------------------------------------------------------------------------
    !
    ! Dealloactes matrix and reflector from QR_factorization type
    !
    ! --------------------------------------------------------------------------------
    
    type (SVD_factorization) :: this


    deallocate ( this % U, this % Vt, this % Sigma )

    deallocate ( this % pseudo_inverse )


  end subroutine destroySVD_factorization



  subroutine MatVecMultiply ( mat , vec , res)

    ! --------------------------------------------------------------------------------
    !
    ! Wrapper to compute res = Mat * vec
    !
    ! --------------------------------------------------------------------------------
   
    complex(8),dimension(:,:) :: mat

    complex(8),dimension(:) :: vec,res


    integer,dimension(2)::NM    
    

    NM =  shape ( mat )


    call ZGEMM( 'N' , 'N' , NM (1) , 1 , NM (2) , 1.0d0 , mat, NM(1), vec, NM (2), 0.0d0, res, NM (1) )

        
  end subroutine MatVecMultiply



  subroutine computeQR_Factorization ( mat, QR)

    ! --------------------------------------------------------------------------------
    !
    ! Wrapper to compute QR factorization of mat
    !
    ! --------------------------------------------------------------------------------

    complex(8),dimension(:,:),INTENT(IN) :: mat

    type(QR_Factorization) :: QR


    integer,dimension(2) ::  MN
    
    integer :: lwork, info

    complex(8),dimension(:),allocatable :: work


    MN = shape ( mat )

    lwork = MN (2)

    allocate ( work ( 1 : lwork ) )


    call ZGEQRF ( MN(1), MN(2), mat, MN(1), QR % tau, work, lwork, info)


    deallocate ( work )

    
  end subroutine computeQR_Factorization


  
  subroutine computeSVD_factorization ( this, mat )

    ! --------------------------------------------------------------------------------
    ! Wrapper to compute SVD factorization of mat using LAPACK
    ! --------------------------------------------------------------------------------

    complex(8),dimension(:,:) ::  mat

    type ( SVD_factorization ) :: this

    
    integer :: lda, ldu, ldvt

    integer :: lwork, info

    real(8), dimension(:), allocatable :: rwork

    complex(8), dimension(:), allocatable :: work


    lda  = this % mat_shape ( 1 )

    ldu  = this % mat_shape ( 1 )

    ldvt = this % mat_shape ( 2 )


    lwork = 2 * minval ( this % mat_shape ) + maxval ( this % mat_shape )

    allocate ( work ( lwork) )


    allocate ( rwork( 5 * minval ( this % mat_shape ) ) )
    


    call zgesvd ( 'A', 'A', &
         this % mat_shape (1), this % mat_shape (2), mat, lda, &
         this % Sigma, &
         this % U, ldu, &
         this % Vt, ldvt, &
         work, lwork, rwork, info )


    deallocate ( rwork, work )
    
    
  end subroutine computeSVD_factorization



  subroutine outerProduct ( u, v, alpha,  res )

    ! ----------------------------------------------------------------------
    !
    ! computes res = res + alpha * u * v'
    !
    ! ----------------------------------------------------------------------

    complex(8), dimension(:) :: u, v

    real(8) :: alpha

    complex(8), dimension(:,:) :: res


    integer :: m, n


    m = size ( u )

    n = size ( v )

    res = res + alpha * matmul ( reshape ( u, (/ m, 1 /) ), reshape ( v, (/ 1, n /) ) )


  end subroutine outerProduct



  subroutine compute_Pseudoinverse ( svd, tol )

    ! ----------------------------------------------------------------------
    !
    ! computes the Penrose pseudo-inverse of the matrix U* Sigma * Vt. 
    ! tol is used to establish an upper bound for the condition number
    ! of a truncated SVD. The upper bound is 1/tol
    !
    ! ----------------------------------------------------------------------

    type ( SVD_factorization ), target :: svd

    real(8) :: tol


    complex(8), dimension(:,:), pointer :: U_r, V_r

    real(8), dimension(:), pointer ::Sigma_r
    
    integer :: j, r


    r = 1

    do while ( svd % Sigma ( 1 ) * tol < svd % Sigma ( r ) )
       
       r = r + 1

    end do

    svd % pseudo_inverse = 0.0d0

    do j = 1, r-1

       call outerProduct ( &
            conjg ( svd % Vt ( j, : ) ), &
            conjg ( svd % U ( :, j ) ), &
            1 / svd % Sigma ( j ), &
            svd % pseudo_inverse )

    end do


  end subroutine compute_Pseudoinverse



  subroutine LeastSquareSolve_SVD ( svd, rhs )

    ! ----------------------------------------------------------------------
    ! 
    ! Computes Least Square solution using the Moore-Penrose pseudo-inverse
    ! and stores it in x
    !
    ! ----------------------------------------------------------------------
    
    type ( SVD_factorization ) :: svd

    complex(8), dimension(:) :: rhs


    rhs = matmul ( svd % pseudo_inverse, rhs )


  end subroutine LeastSquareSolve_SVD



  subroutine LeastSquareSolve ( QR_f, rhs )

    ! --------------------------------------------------------------------------------
    ! wrapper to compute the l.l.s problem
    ! min || A * x - rhs ||_2 using QR
    ! factorization.
    ! --------------------------------------------------------------------------------
    
    TYPE(QR_factorization) :: QR_f

    complex(8),dimension(:) :: rhs

    
    INTEGER :: M, N

    INTEGER :: lwork, info

    complex(8),dimension(:),allocatable :: work


    M = size ( rhs )

    N = 1

    lwork = N


    allocate ( work( lwork ) )


    call ZUNMQR ( 'L', 'C', M, N, M, QR_f % mat, M, QR_f % tau, rhs, M, work, lwork, info)

    call ZTRTRS ( 'U', 'N', 'N', M, N, QR_f % mat, M, rhs, M, info)


    deallocate ( work )


  end subroutine LeastSquareSolve



end module modLinearAlgebra
