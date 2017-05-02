MODULE modInterpolationReferenceCell



  USE modMathConstants

  USE modLinearAlgebra



  IMPLICIT NONE


 
  TYPE InterpolationReferenceCell


     !Size of Cell
     REAL(8) :: L_x, L_y
     
     !Number of sources in first quadrant
     INTEGER :: N_src

     !Number of plane waves in first quadrant to interpolate inside cell
     INTEGER :: N_wave

     !Direction of plane waves
     REAL(8), DIMENSION(:), ALLOCATABLE :: u_x, u_y

     !Location of horizontal and vertical equivalent sources    
     REAL(8), DIMENSION(:), ALLOCATABLE :: src_x, src_y

     ! QR factorization of blocks for Least Square problem
     TYPE(QR_factorization) :: QR_A, QR_B, QR_C, QR_D


  END TYPE InterpolationReferenceCell




CONTAINS



  SUBROUTINE createInterpolationReferenceCell (this, N_s_h, N_s_v, N_w, L_x, L_y, k)
    ! Class constructor
    
    TYPE(InterpolationReferenceCell) :: this

    INTEGER :: N_s_h, N_s_v, N_w

    REAL(8) :: L_x, L_y, k


    REAL,ALLOCATABLE,DIMENSION(:) :: disc

    INTEGER :: j


    this % L_x = L_x

    this % L_y = L_y


    this % N_src = N_s_h + N_s_v 

    this % N_wave = N_w

    
    ! Construct source points

    ALLOCATE ( this % src_x ( 1 : this % N_src ), this % src_y ( 1 : this % N_src ) )

    
    this % src_x ( 1 : N_s_h ) = L_x * (/ ( ( 2.0d0 * j + 1.0d0 ) / ( 4.0d0 * N_s_h ) - 0.5d0, j = N_s_h, 2 * N_s_h - 1 ) /)

    this % src_y ( 1 : N_s_h ) = L_y * (/ ( 0.5d0, j = 1, N_s_h ) /)


    this % src_x ( N_s_h + 1 : N_s_h + N_s_v ) = L_x * (/ ( 0.5d0, j = 1, N_s_v ) /)
    
    this % src_y ( N_s_h + 1 : N_s_h + N_s_v ) = L_y * (/ ( ( 2.0d0 * j + 1.0d0 ) / ( 4.0d0 * N_s_v ) - 0.5d0, j = 2 * N_s_v - 1, N_s_v, -1 ) /)


    ! Construct direction of plane waves

    ALLOCATE ( this % u_x ( 1 : N_w ), this % u_y ( 1 : N_w ) )

    
    this % u_x = (/ ( cos ( pi / 2.0d0 * ( 2.0d0 * j + 1.0d0) / ( 2.0d0 * N_w ) ), j =0, N_w - 1 ) /)

    this % u_y = (/ ( sin ( pi / 2.0d0 * ( 2.0d0 * j + 1.0d0) / ( 2.0d0 * N_w ) ), j = 0, N_w - 1 ) /)


    ! Construct QR factorizations of interaction matrices

    ALLOCATE ( this % QR_A % mat ( 1 : N_s_h + N_s_v , 1 : N_w ) )

    ALLOCATE ( this % QR_B % mat ( 1 : N_s_h + N_s_v , 1 : N_w ) )

    ALLOCATE ( this % QR_C % mat ( 1 : N_s_h + N_s_v , 1 : N_w ) )

    ALLOCATE ( this % QR_D % mat ( 1 : N_s_h + N_s_v , 1 : N_w ) )


    ALLOCATE ( this % QR_A % tau ( 1 : MIN ( N_s_h + N_s_v , N_w ) ) )

    ALLOCATE ( this % QR_B % tau ( 1 : MIN ( N_s_h + N_s_v , N_w ) ) )

    ALLOCATE ( this % QR_C % tau ( 1 : MIN ( N_s_h + N_s_v , N_w ) ) )

    ALLOCATE ( this % QR_D % tau ( 1 : MIN ( N_s_h + N_s_v , N_w ) ) )


    CALL createMatrices (this, k)

    
  END SUBROUTINE createInterpolationReferenceCell


  
  SUBROUTINE createMatrices (this, k)
    ! Computes the QR factorization of each block to solve
    ! the least square problem to interpolate the field at 
    ! the boundary inside the cell.

    ! Read notes for a better comprehension of how the
    ! blocks are constructed.
    
    TYPE(InterpolationReferenceCell),TARGET :: this

    REAL(8) :: k

    CHARACTER :: type_flag


    INTEGER :: j

    REAL(8),DIMENSION(:,:),ALLOCATABLE :: P_x, P_y

    REAL(8),DIMENSION(:,:,:),ALLOCATABLE :: Q_x, Q_y

    COMPLEX(8),DIMENSION(:,:,:),ALLOCATABLE :: E


    ALLOCATE ( P_x ( 1 : this % N_src, 1 : this % N_wave) )

    ALLOCATE ( P_y ( 1 : this % N_src, 1 : this % N_wave) )


    ALLOCATE ( Q_x ( 1 : this % N_src, 1 : this % N_wave, 1:4 ) )

    ALLOCATE ( Q_y ( 1 : this % N_src, 1 : this % N_wave, 1:4 ) )

    
    ALLOCATE ( E ( 1 : this % N_src, 1 : this % N_wave, 1:4 ) )


    !CREATE QR DECOMPOSITION FOR HORIZONTAL INTERACTION MATRICES

    Q_x (:,:,1) = SPREAD ( this % u_x, 1, this % N_src )

    Q_y (:,:,1) = SPREAD ( this % u_y, 1, this % N_src )


    Q_x (:,:,2) = -Q_x (:,:,1) 

    Q_y (:,:,2) =  Q_y (:,:,1) 


    Q_x (:,:,3) = -Q_x (:,:,1) 

    Q_y (:,:,3) = -Q_y (:,:,1) 


    Q_x (:,:,4) =  Q_x (:,:,1) 

    Q_y (:,:,4) = -Q_y (:,:,1) 


    P_x = SPREAD ( this % src_x, 2, this % N_wave ) 
       
    P_y = SPREAD ( this % src_y, 2, this % N_wave )
       

    !Compute Plane Waves Evaluated at Equiv. sources

    FORALL ( j = 1:4 ) E(:,:,j) = EXP ( I * k * ( P_x * Q_x (:,:,j) + P_y * Q_y (:,:,j) ) )


    this % QR_A % mat = E(:,:,1) + E(:,:,2) + E(:,:,3) + E(:,:,4) 

    this % QR_B % mat = E(:,:,1) - E(:,:,2) + E(:,:,3) - E(:,:,4)

    this % QR_C % mat = E(:,:,1) + E(:,:,2) - E(:,:,3) - E(:,:,4)

    this % QR_D % mat = E(:,:,1) - E(:,:,2) - E(:,:,3) + E(:,:,4)

    
    CALL computeQR_Factorization ( this % QR_A % mat, this % QR_A)
    
    CALL computeQR_Factorization ( this % QR_B % mat, this % QR_B)

    CALL computeQR_Factorization ( this % QR_C % mat, this % QR_C)
    
    CALL computeQR_Factorization ( this % QR_D % mat, this % QR_D)


    DEALLOCATE ( E )
    
    DEALLOCATE ( P_x, P_y, Q_x, Q_y )

        
  END SUBROUTINE createMatrices



  SUBROUTINE computePlaneWaveWeights (this, rhs, xi)
    !Solves a least square problem to obtain the amplitudes
    !of the plane waves to interpolate inside the cell.

    !rhs contains the field in the boundary of the cell
    !xi stores the amplitudes.


    TYPE(InterpolationReferenceCell),TARGET :: this

    COMPLEX(8),TARGET,DIMENSION(:) :: rhs
    
    COMPLEX(8),TARGET,DIMENSION(:) :: xi


    INTEGER,POINTER :: N_s, N_w

    COMPLEX(8),POINTER :: rhs_1(:), rhs_2(:), rhs_3(:), rhs_4(:)

    COMPLEX(8),POINTER :: xi_1(:), xi_2(:), xi_3(:), xi_4(:)

    COMPLEX(8),ALLOCATABLE,DIMENSION(:) :: b_hat_a, b_hat_b, b_hat_c, b_hat_d


    N_s => this % N_src

    N_w => this % N_wave


    ALLOCATE ( b_hat_a ( 1 : N_s ) )

    ALLOCATE ( b_hat_b ( 1 : N_s ) )

    ALLOCATE ( b_hat_c ( 1 : N_s ) )

    ALLOCATE ( b_hat_d ( 1 : N_s ) )


    rhs_1 ( 1 : N_s ) => rhs ( 0 * N_s + 1 : 1 * N_s )

    rhs_2 ( 1 : N_s ) => rhs ( 1 * N_s + 1 : 2 * N_s )

    rhs_3 ( 1 : N_s ) => rhs ( 2 * N_s + 1 : 3 * N_s )

    rhs_4 ( 1 : N_s ) => rhs ( 3 * N_s + 1 : 4 * N_s)


    xi_1 ( 1 : N_w ) => xi ( 0 * N_w + 1 : 1 * N_w )

    xi_2 ( 1 : N_w ) => xi ( 1 * N_w + 1 : 2 * N_w )

    xi_3 ( 1 : N_w ) => xi ( 2 * N_w + 1 : 3 * N_w )

    xi_4 ( 1 : N_w ) => xi ( 3 * N_w + 1 : 4 * N_w ) 


    b_hat_a = 0.5d0 * (rhs_1 + rhs_2 + rhs_3 + rhs_4)

    b_hat_b = 0.5d0 * (rhs_1 - rhs_2 + rhs_3 - rhs_4)

    b_hat_c = 0.5d0 * (rhs_1 + rhs_2 - rhs_3 - rhs_4)

    b_hat_d = 0.5d0 * (rhs_1 - rhs_2 - rhs_3 + rhs_4)



    CALL LeastSquareSolve ( this % QR_A, b_hat_a)

    CALL LeastSquareSolve ( this % QR_B, b_hat_b)

    CALL LeastSquareSolve ( this % QR_C, b_hat_c)

    CALL LeastSquareSolve ( this % QR_D, b_hat_d)


    xi_1 = 0.5d0 * ( b_hat_a ( 1 : N_w ) + b_hat_b ( 1 : N_w ) + b_hat_c ( 1 : N_w ) + b_hat_d ( 1 : N_w ) )

    xi_2 = 0.5d0 * ( b_hat_a ( 1 : N_w ) - b_hat_b ( 1 : N_w ) + b_hat_c ( 1 : N_w ) - b_hat_d ( 1 : N_w ) )

    xi_3 = 0.5d0 * ( b_hat_a ( 1 : N_w ) + b_hat_b ( 1 : N_w ) - b_hat_c ( 1 : N_w ) - b_hat_d ( 1 : N_w ) )

    xi_4 = 0.5d0 * ( b_hat_a ( 1 : N_w ) - b_hat_b ( 1 : N_w ) - b_hat_c ( 1 : N_w ) + b_hat_d ( 1 : N_w ) )



    NULLIFY ( rhs_1, rhs_2, rhs_3, rhs_4 )

    NULLIFY ( xi_1, xi_2, xi_3, xi_4 )

    NULLIFY ( N_s, N_w )

    DEALLOCATE ( b_hat_a, b_hat_b, b_hat_c, b_hat_d ) 


  END SUBROUTINE computePlaneWaveWeights



  SUBROUTINE destroyInterpolationReferenceCell (this)

    
    TYPE(InterpolationReferenceCell) :: this


    DEALLOCATE ( this % u_x, this % u_y, this % src_x, this % src_y)

    DEALLOCATE ( this % QR_A % mat, this % QR_B % mat, this % QR_C % mat, this % QR_D % mat )

    DEALLOCATE ( this % QR_A % tau, this % QR_B % tau, this % QR_C % tau, this % QR_D % tau )


  END SUBROUTINE destroyInterpolationReferenceCell



END MODULE modInterpolationReferenceCell
