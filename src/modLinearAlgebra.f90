MODULE modLinearAlgebra



  IMPLICIT NONE



  TYPE QR_factorization 
     ! mat contains the QR factorization computed with LAPACK ZGEQRF
     ! tau contains scalar to define elementary reflectors 


     COMPLEX(8),DIMENSION(:,:),ALLOCATABLE :: mat

     COMPLEX(8),DIMENSION(:),ALLOCATABLE :: tau

     
  END TYPE QR_factorization



CONTAINS



  SUBROUTINE LinearSolve ( mat, vec, res)
    ! wrapper to compute the solution of 
    ! mat * res = vec


    COMPLEX(8),DIMENSION(:,:) :: mat

    COMPLEX(8),DIMENSION(:) :: vec, res


    COMPLEX(8),DIMENSION(:,:),ALLOCATABLE :: mat_copy

    COMPLEX(8),DIMENSION(:),ALLOCATABLE :: vec_copy

    INTEGER,DIMENSION(:),ALLOCATABLE :: ipiv

    INTEGER,DIMENSION(2) :: mat_size

    INTEGER :: info,N

    
    N = SIZE ( vec )


    ALLOCATE ( vec_copy ( 0 : N-1 ) )

    ALLOCATE ( mat_copy ( 0 : N-1, 0 : N-1 ) )    

    ALLOCATE ( ipiv ( 0 : N-1 ) )


    mat_copy = mat

    vec_copy = vec


    CALL ZGESV ( N, 1, mat_copy, N, ipiv, vec_copy, N, info) 
    
    res = vec_copy


    DEALLOCATE ( vec_copy, mat_copy, ipiv )


  END SUBROUTINE LinearSolve



  SUBROUTINE MatVecMultiply ( mat , vec , res)
    ! Wrapper to compute res = Mat * vec

   
    COMPLEX(8),DIMENSION(:,:) :: mat

    COMPLEX(8),DIMENSION(:) :: vec,res


    INTEGER,DIMENSION(2)::NM
    

    NM = SHAPE ( mat )


    CALL ZGEMM( 'N' , 'N' , NM (1) , 1 , NM (2) , 1.0d0 , mat, NM(1), vec, NM (2), 0.0d0, res, NM (1) )

        
  END SUBROUTINE MatVecMultiply



  SUBROUTINE computeQR_Factorization ( mat, QR)
    ! Wrapper to compute QR factorization of mat
    

    COMPLEX(8),DIMENSION(:,:),INTENT(IN) :: mat

    TYPE(QR_Factorization) :: QR


    INTEGER,DIMENSION(2) ::  MN
    
    INTEGER :: lwork, info

    COMPLEX(8),DIMENSION(:),ALLOCATABLE :: work


    MN = SHAPE ( mat )

    lwork = MN (2)

    ALLOCATE ( work ( 1 : lwork ) )


    CALL ZGEQRF ( MN(1), MN(2), mat, MN(1), QR % tau, work, lwork, info)


    DEALLOCATE ( work )

    
  END SUBROUTINE computeQR_Factorization



  SUBROUTINE LeastSquareSolve ( QR_f, rhs)
    ! wrapper to compute the l.l.s problem
    ! min || A * x - rhs ||_2 using QR
    ! factorization.
    
    TYPE(QR_factorization) :: QR_f

    COMPLEX(8),DIMENSION(:) :: rhs

    
    INTEGER :: M, N

    INTEGER :: lwork, info

    COMPLEX(8),DIMENSION(:),ALLOCATABLE :: work


    M = SIZE ( rhs )

    N = 1

    lwork = N


    ALLOCATE ( work( lwork ) )


    CALL ZUNMQR ( 'L', 'C', M, N, M, QR_f % mat, M, QR_f % tau, rhs, M, work, lwork, info)

    CALL ZTRTRS ( 'U', 'N', 'N', M, N, QR_f % mat, M, rhs, M, info)


    DEALLOCATE ( work )


  END SUBROUTINE LeastSquareSolve



END MODULE modLinearAlgebra
