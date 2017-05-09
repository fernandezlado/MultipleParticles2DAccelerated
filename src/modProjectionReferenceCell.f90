module modProjectionReferenceCell



  use modLinearAlgebra

  use modSpecialFunctions



  implicit none


 
  type ProjectionReferenceCell


     !Kind of Cell (Horizontal or Vertical)
     character :: cellKind

     !Size of Cell
     real(8) :: L_x, L_y
     
     !Number of sources, collocation points. N_dim=4 in 2D.
     integer :: N_src, N_coll, N_wave, N_dim

     !Location of collocation points
     real(8), dimension(:), allocatable :: coll_x, coll_y

     !location of horizontal and vertical equivalent sources    
     real(8), dimension(:), allocatable :: src_x, src_y

     ! qr factorization of blocks for least square problem
     type(QR_factorization) :: QR_A, QR_B, QR_C, QR_D


  end type ProjectionReferenceCell



contains



  subroutine createProjectionReferenceCell (this, N_s, N_c, L_x, L_y, k, typeFlag)
    ! Class constructor

    
    type(ProjectionReferenceCell) :: this

    integer :: N_s, N_c

    real(8) :: L_x, L_y, k

    character :: typeFlag


    real,allocatable,dimension(:) :: disc

    integer :: j


    this % cellKind = typeFlag


    this % L_x = L_x

    this % L_y = L_y


    this % N_src = N_s

    this % N_coll = N_c


    this % N_dim = 4

    
    ! Construct source and collocation points. Read notes to learn
    ! about the location of src and collocation points.

    allocate ( this % coll_x ( N_c ), this % coll_y ( N_c ) )

    allocate ( this % src_x ( N_s ), this % src_y ( N_s ) )

    
    if ( this % cellKind == 'H' ) then
       

       this % src_x = L_x * (/ ( ( 2.0d0 * j + 1.0d0 ) / (4.0d0 * N_s ) - 0.5d0, j = N_s, 2 * N_s - 1 ) /)

       this % src_y = L_y * (/ ( 0.5d0, j = N_s, 2 * N_s - 1 )  /)


       this % coll_x = 3.0d0 * L_x * (/ ( ( 2.0d0 * j + 1.0d0 ) / (4.0d0 * N_c ) - 0.5d0, j = N_c, 2 * N_c - 1 ) /)

       this % coll_y = 3.0d0 * L_y * (/ ( 0.5d0, j = N_c, 2 * N_c - 1 )  /)


    else if ( this % cellKind == 'V' ) then


       this % src_x = L_x * (/ ( 0.5d0, j = N_s, 2 * N_s - 1 )  /)
       
       this % src_y = L_y * (/ ( ( 2.0d0 * j + 1.0d0 ) / (4.0d0 * N_s ) - 0.5d0, j = 0,  N_s - 1 ) /)


       this % coll_x = 3.0d0 * L_x * (/ ( 0.5d0, j = N_c, 2 * N_c - 1 ) /)
    
       this % coll_y = 3.0d0 * L_y * (/ ( ( 2.0d0 * j + 1.0d0 ) / (4.0d0 * N_c ) - 0.5d0, j = 0, N_c - 1 ) /)


    end if
    

    !Construct QR factorizations of interaction matrices

    call createQR_factorization ( this % QR_A, N_c, 2 * N_s )

    call createQR_factorization ( this % QR_B, N_c, 2 * N_s )

    call createQR_factorization ( this % QR_C, N_c, 2 * N_s )

    call createQR_factorization ( this % QR_D, N_c, 2 * N_s )


    call createMatrices (this, k)

    
  end subroutine createProjectionReferenceCell


  
  subroutine createMatrices (this, k)
    ! Computes the QR factorization of the 4 blocks to
    ! solve Least Square problem to compute equivalent 
    ! sources. Read notes to learn about how to obtain
    ! these 4 blocks.
    

    type(ProjectionReferenceCell),target :: this

    real(8) :: k

    character :: type_flag


    real(8),dimension(:,:),allocatable :: P_x, P_y

    real(8),dimension(:,:,:),allocatable :: Q_x, Q_y

    complex(8),dimension(:,:),allocatable :: G_1, G_2, G_3, G_4

    real(8) :: d_x, d_y


    ALLOCATE ( P_x ( this % N_coll, this % N_src ) )

    ALLOCATE ( P_y ( this % N_coll, this % N_src ) )


    ALLOCATE ( Q_x ( this % N_coll, this % N_src, 1:4 ) )

    ALLOCATE ( Q_y ( this % N_coll, this % N_src, 1:4 ) )

    
    ALLOCATE ( G_1 ( this % N_coll, this % N_src ) )

    ALLOCATE ( G_2 ( this % N_coll, this % N_src ) )

    ALLOCATE ( G_3 ( this % N_coll, this % N_src ) )

    ALLOCATE ( G_4 ( this % N_coll, this % N_src ) )


    ! THIS IF CONTROLS THE LOCATION OF THE REFLECTIONS OF EQUIVALENT 
    ! SOURCES WHETER THE CELL ACCOUNTS FOR THE HORIZONTAL OF VERTICAL
    ! CASE.

    IF ( this % cellKind == 'H' ) THEN

       
       Q_x (:,:,1) = SPREAD ( this % src_x, 1, this % N_coll )

       Q_y (:,:,1) = SPREAD ( this % src_y, 1, this % N_coll )


       Q_x (:,:,2) = -Q_x (:,:,1) 

       Q_y (:,:,2) =  Q_y (:,:,1) 


       Q_x (:,:,3) = -Q_x (:,:,1) 

       Q_y (:,:,3) = -Q_y (:,:,1) 


       Q_x (:,:,4) =  Q_x (:,:,1) 

       Q_y (:,:,4) = -Q_y (:,:,1) 


       d_x = 0.0d0

       d_y = 1.0d0


    ELSE IF ( this % cellKind == 'V' ) THEN

       
       Q_x (:,:,1) = SPREAD ( this % src_x, 1, this % N_coll )

       Q_y (:,:,1) = SPREAD ( this % src_y, 1, this % N_coll )


       Q_x (:,:,2) = Q_x (:,:,1) 

       Q_y (:,:,2) = -Q_y (:,:,1) 


       Q_x (:,:,3) = -Q_x (:,:,1) 

       Q_y (:,:,3) = -Q_y (:,:,1) 


       Q_x (:,:,4) = -Q_x (:,:,1) 

       Q_y (:,:,4) = Q_y (:,:,1) 

       
       d_x = 1.0d0

       d_y = 0.0d0

       
    END IF


    P_x = SPREAD ( this % coll_x, 2, this % N_src ) 
       
    P_y = SPREAD ( this % coll_y, 2, this % N_src )
       

    ! MONOPOLE TERMS

    G_1 = G_0 ( k, P_x - Q_x (:,:,1), P_y - Q_y (:,:,1) )

    G_2 = G_0 ( k, P_x - Q_x (:,:,2), P_y - Q_y (:,:,2) )

    G_3 = G_0 ( k, P_x - Q_x (:,:,3), P_y - Q_y (:,:,3) )

    G_4 = G_0 ( k, P_x - Q_x (:,:,4), P_y - Q_y (:,:,4) )

   
    this % QR_A % mat ( :, 1 : this % N_src ) = G_1 + G_2 + G_3 + G_4

    this % QR_B % mat ( :, 1 : this % N_src ) = G_1 - G_2 + G_3 - G_4

    this % QR_C % mat ( :, 1 : this % N_src ) = G_1 + G_2 - G_3 - G_4
 
    this % QR_D % mat ( :, 1 : this % N_src ) = G_1 - G_2 - G_3 + G_4


    ! DIPOLE TERMS

    G_1 = DG_0 ( k, P_x - Q_x (:,:,1), P_y - Q_y (:,:,1), d_x, d_y )

    G_2 = DG_0 ( k, P_x - Q_x (:,:,2), P_y - Q_y (:,:,2), d_x, d_y )
 
    G_3 = DG_0 ( k, P_x - Q_x (:,:,3), P_y - Q_y (:,:,3), -d_x, -d_y )

    G_4 = DG_0 ( k, P_x - Q_x (:,:,4), P_y - Q_y (:,:,4), -d_x, -d_y )


    this % QR_A % mat ( :, this % N_src + 1 : 2 * this % N_src ) = G_1 + G_2 + G_3 + G_4

    this % QR_B % mat ( :, this % N_src + 1 : 2 * this % N_src ) = G_1 - G_2 + G_3 - G_4

    this % QR_C % mat ( :, this % N_src + 1 : 2 * this % N_src ) = G_1 + G_2 - G_3 - G_4

    this % QR_D % mat ( :, this % N_src + 1 : 2 * this % N_src ) = G_1 - G_2 - G_3 + G_4


    !Compute QR factorization

    CALL computeQR_Factorization ( this % QR_A % mat, this % QR_A)
    
    CALL computeQR_Factorization ( this % QR_B % mat, this % QR_B)

    CALL computeQR_Factorization ( this % QR_C % mat, this % QR_C)
    
    CALL computeQR_Factorization ( this % QR_D % mat, this % QR_D)


    DEALLOCATE ( G_1, G_2, G_3, G_4 )
    
    DEALLOCATE ( P_x, P_y, Q_x, Q_y )

        
  END SUBROUTINE createMatrices



  SUBROUTINE computeEquivalentSourceAmplitude (this, rhs, x)


    TYPE(ProjectionReferenceCell),TARGET :: this

    COMPLEX(8),TARGET,DIMENSION(:) :: rhs
    
    COMPLEX(8),TARGET,DIMENSION(:) :: x


    INTEGER,POINTER :: N_s, N_c

    COMPLEX(8),POINTER :: rhs_1(:), rhs_2(:), rhs_3(:), rhs_4(:)

    COMPLEX(8),POINTER :: x_1(:), x_2(:), x_3(:), x_4(:)

    COMPLEX(8),ALLOCATABLE,DIMENSION(:) :: b_hat_a, b_hat_b, b_hat_c, b_hat_d


    N_s => this % N_src

    N_c => this % N_coll

    
    ALLOCATE ( b_hat_a ( 1 : N_c ) )

    ALLOCATE ( b_hat_b ( 1 : N_c ) )

    ALLOCATE ( b_hat_c ( 1 : N_c ) )

    ALLOCATE ( b_hat_d ( 1 : N_c ) )

    
    rhs_1 ( 1 : N_c ) => rhs ( 0 * N_c + 1 : 1 * N_c )

    rhs_2 ( 1 : N_c ) => rhs ( 1 * N_c + 1 : 2 * N_c )

    rhs_3 ( 1 : N_c ) => rhs ( 2 * N_c + 1 : 3 * N_c )

    rhs_4 ( 1 : N_c ) => rhs ( 3 * N_c + 1 : 4 * N_c )
    
    
    x_1 ( 1 : 2 * N_s ) => x ( 0 * N_s + 1 : 2 * N_s )

    x_2 ( 1 : 2 * N_s ) => x ( 2 * N_s + 1 : 4 * N_s )

    x_3 ( 1 : 2 * N_s ) => x ( 4 * N_s + 1 : 6 * N_s )

    x_4 ( 1 : 2 * N_s ) => x ( 6 * N_s + 1 : 8 * N_s )


    b_hat_a = 0.5d0 * ( rhs_1 + rhs_2 + rhs_3 + rhs_4 )
    
    b_hat_b = 0.5d0 * ( rhs_1 - rhs_2 + rhs_3 - rhs_4 )

    b_hat_c = 0.5d0 * ( rhs_1 + rhs_2 - rhs_3 - rhs_4 )

    b_hat_d = 0.5d0 * ( rhs_1 - rhs_2 - rhs_3 + rhs_4 )


    
    CALL LeastSquareSolve ( this % QR_A, b_hat_a)

    CALL LeastSquareSolve ( this % QR_B, b_hat_b)

    CALL LeastSquareSolve ( this % QR_C, b_hat_c)

    CALL LeastSquareSolve ( this % QR_D, b_hat_d)



    x_1 = 0.5d0 * ( b_hat_a ( 1 : 2 * N_s ) + b_hat_b ( 1 : 2 * N_s ) + b_hat_c ( 1 : 2 * N_s ) + b_hat_d ( 1 : 2 * N_s ) )

    x_2 = 0.5d0 * ( b_hat_a ( 1 : 2 * N_s ) - b_hat_b ( 1 : 2 * N_s ) + b_hat_c ( 1 : 2 * N_s ) - b_hat_d ( 1 : 2 * N_s ) )

    x_3 = 0.5d0 * ( b_hat_a ( 1 : 2 * N_s ) + b_hat_b ( 1 : 2 * N_s ) - b_hat_c ( 1 : 2 * N_s ) - b_hat_d ( 1 : 2 * N_s ) )

    x_4 = 0.5d0 * ( b_hat_a ( 1 : 2 * N_s ) - b_hat_b ( 1 : 2 * N_s ) - b_hat_c ( 1 : 2 * N_s ) + b_hat_d ( 1 : 2 * N_s ) )


    nullify ( rhs_1, rhs_2, rhs_3, rhs_4 )

    nullify ( x_1, x_2, x_3, x_4 )

    nullify ( N_s, N_c )


  end subroutine computeEquivalentSourceAmplitude



  subroutine destroyProjectionReferenceCell (this)

    
    type(ProjectionReferenceCell)::this


    deallocate ( this % coll_x, this % coll_y, this % src_x, this % src_y)

    call destroyQR_factorization ( this % QR_A )

    call destroyQR_factorization ( this % QR_B ) 

    call destroyQR_factorization ( this % QR_C )

    call destroyQR_factorization ( this % QR_D )


  end subroutine destroyProjectionReferenceCell



end module modProjectionReferenceCell
