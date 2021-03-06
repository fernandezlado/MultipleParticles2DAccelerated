MODULE modForwardMap



  USE modParameters

  USE modFFT

  USE modCell

  USE modProjectionReferenceCell



  IMPLICIT NONE



  TYPE (FFT_2D) :: fft_hor_local, fft_ver_local

  TYPE (FFT_2D) :: fft_hor_global, fft_ver_global

  
  TYPE (ProjectionReferenceCell) :: proj_cell_hor, proj_cell_ver

  TYPE (InterpolationReferenceCell) :: interp_cell

  COMPLEX(8),DIMENSION(:),ALLOCATABLE :: x_hor, x_ver, xi

  
  TYPE(Cell), DIMENSION(:,:), ALLOCATABLE :: cellArray

  TYPE(Obstacle), DIMENSION(:,:), ALLOCATABLE :: obstacleArray


  COMPLEX(8),DIMENSION(:,:,:),ALLOCATABLE :: psi, res


  COMPLEX(8),DIMENSION(:,:),ALLOCATABLE :: mon_equiv_hor, mon_equiv_ver

  COMPLEX(8),DIMENSION(:,:),ALLOCATABLE :: dip_equiv_hor, dip_equiv_ver

  COMPLEX(8),DIMENSION(:,:),ALLOCATABLE :: local_hor, local_ver

  
  COMPLEX(8),DIMENSION(:,:),ALLOCATABLE :: mon_hor_ker_global_fft, dip_hor_ker_global_fft

  COMPLEX(8),DIMENSION(:,:),ALLOCATABLE :: mon_ver_ker_global_fft, dip_ver_ker_global_fft

  COMPLEX(8),DIMENSION(:,:),ALLOCATABLE :: mon_hor_ker_local_fft, dip_hor_ker_local_fft

  COMPLEX(8),DIMENSION(:,:),ALLOCATABLE :: mon_ver_ker_local_fft, dip_ver_ker_local_fft


  COMPLEX(8),DIMENSION(:,:),ALLOCATABLE,TARGET :: field_equiv_hor, field_equiv_ver

  COMPLEX(8),DIMENSION(:,:),ALLOCATABLE :: aux_global_hor, aux_global_ver

  COMPLEX(8),DIMENSION(:,:),ALLOCATABLE :: aux_local_hor, aux_local_ver

  COMPLEX(8),DIMENSION(:),ALLOCATABLE,TARGET :: field_permuted

  COMPLEX(8),DIMENSION(:),POINTER :: field_per_1, field_per_2, field_per_3, field_per_4



  PRIVATE :: proj_cell_hor, proj_cell_ver, interp_cell

  PRIVATE :: cellArray, obstacleArray

  PRIVATE :: mon_equiv_hor, mon_equiv_ver

  PRIVATE :: dip_equiv_hor, dip_equiv_ver

  PRIVATE :: mon_hor_ker_global_fft, dip_hor_ker_global_fft

  PRIVATE :: mon_ver_ker_global_fft, dip_ver_ker_global_fft

  PRIVATE :: x_hor, x_ver, xi



CONTAINS


  
  SUBROUTINE initForwardMap ( phy, geo, alg )


    TYPE (phy_Parameters) :: phy
    
    TYPE (geo_Parameters) :: geo

    TYPE (alg_Parameters) :: alg


    INTEGER :: N, M, N_s
    

    ! --------------------------------------------------------------------------------
    !
    ! Create array to store densities and the result of the Forward Map
    !
    ! --------------------------------------------------------------------------------


    ALLOCATE ( psi ( alg % N_dis, geo % N_row, geo % N_col ) )

    ALLOCATE ( res ( alg % N_dis, geo % N_row, geo % N_col) )
    
    psi = 1.0d0

    res = 1.0d0
    

    ! --------------------------------------------------------------------------------
    !
    ! Create Projection and Interpolation reference Cells and x_hor, x_ver, xi that
    ! contain the solutions of the least squares problem to project and interpolate
    ! the field with the equivalent sources.
    !
    ! --------------------------------------------------------------------------------


    CALL createProjectionReferenceCell ( &
         proj_cell_hor, &
         alg % N_src_hor, &
         alg % N_coll_hor, &
         geo % L_x, &
         geo % L_y, &
         phy % k, & 
         'H' )


    CALL createProjectionReferenceCell ( &
         proj_cell_ver, &
         alg % N_src_ver, &
         alg % N_coll_ver, &
         geo % L_x, &
         geo % L_y, &
         phy % k, &
         'V' )


    CALL createInterpolationReferenceCell ( &
         interp_cell, &
         alg % N_src_hor, &
         alg % N_src_ver, &
         alg % N_wave, &
         geo % L_x, &
         geo % L_y, &
         phy % k )


    ! Create x_hor, x_ver: solutions of least square problem to obtain equiv. sources

    ALLOCATE ( x_hor ( 2 * 4 * alg % N_src_hor ) )

    ALLOCATE ( x_ver ( 2 * 4 * alg % N_src_ver ) )

    ! Create xi : solution of least square problem to expand field as sum of plane waves

    ALLOCATE ( xi ( 4 * alg % N_wave) )


    ! --------------------------------------------------------------------------------
    !
    ! Create Cell and Obstacles arrays
    !
    ! --------------------------------------------------------------------------------


    ALLOCATE ( cellArray ( 0 : geo % N_row + 1, 0 : geo % N_col + 1 ) )

    ALLOCATE ( obstacleArray ( 1 : geo % N_row, 1 : geo % N_col ) )    
    

    CALL createArrayGeometry ( phy, geo, alg )


    ! --------------------------------------------------------------------------------
    !
    ! Create FFTs handlers, compute FFT of global and local Kernels
    ! Create equiv. sources matrices, field at equiv. sources and
    ! auxiliary field evaluation. See ForwardMap below 
    !
    ! --------------------------------------------------------------------------------
    
    ! --------------------------------------------------------------------------------
    ! Horizontal case
    ! --------------------------------------------------------------------------------

    ! Global 

    !! Size of matrices

    N_s = 2 * alg % N_src_hor

    M   = N_s * ( geo % N_col + 2 )

    N   = geo % N_row + 3

    !! 2D FFT Handler
    
    CALL createFFT_2D ( fft_hor_global, 2 * M, 2 * N )

    !! FFTs of Green function and derivative
    
    ALLOCATE ( mon_hor_ker_global_fft ( 2 * M, 2 * N ) )
  
    ALLOCATE ( dip_hor_ker_global_fft ( 2 * M, 2 * N ) ) 
    

    CALL computeKernelFFT ( &
         phy % k, &
         M, N, &
         (/ geo % L_x / N_s , geo % L_y /), &
         mon_hor_ker_global_fft, &
         'M', &
         fft_hor_global )

    CALL computeKernelFFT ( &
         phy % k, &
         M, N, &
         (/ geo % L_x / N_s , geo % L_y /), &
         dip_hor_ker_global_fft, &
         'D', &
         fft_hor_global )


    !! Monopole and dipole equiv. sources

    ALLOCATE ( mon_equiv_hor ( 2 * M, 2 * N ) )

    ALLOCATE ( dip_equiv_hor ( 2 * M, 2 * N ) )
    

    !! Total field at equiv. sources and aux. matrix to compute
    !! Monopole/Dipole field

    ALLOCATE ( field_equiv_hor ( 2 * M, 2 * N ) )

    ALLOCATE ( aux_global_hor ( 2 * M, 2 * N ) )

    ! Local

    !! Size of matrices

    N_s = 2 * alg % N_src_hor

    M   = N_s

    N   = 2

    !! 2D FFT Handler
    
    CALL createFFT_2D ( fft_hor_local, 2 * M, 2 * N )

    !! FFTs of Green function and derivative
    
    ALLOCATE ( mon_hor_ker_local_fft ( 2 * M, 2 * N ) )
  
    ALLOCATE ( dip_hor_ker_local_fft ( 2 * M, 2 * N ) ) 
    

    CALL computeKernelFFT ( &
         phy % k, &
         M, N, &
         (/ geo % L_x / N_s , geo % L_y /), &
         mon_hor_ker_local_fft, &
         'M', &
         fft_hor_local )

    CALL computeKernelFFT ( &
         phy % k, &
         M, N, &
         (/ geo % L_x / N_s , geo % L_y /), &
         dip_hor_ker_local_fft, &
         'D', &
         fft_hor_local )


    ALLOCATE ( aux_local_hor ( 2 * M, 2 * N ) )

    ALLOCATE ( local_hor ( 2 * M, 2 * N) )


    ! --------------------------------------------------------------------------------
    ! Vertical case
    ! --------------------------------------------------------------------------------

    ! Global 

    !! Size of matrices

    N_s = 2 * alg % N_src_ver

    M = N_s * ( geo % N_row + 2 )

    N = geo % N_col + 3

    !! 2D FFT Handler

    CALL createFFT_2D ( fft_ver_global, 2 * M, 2 * N )

    !! FFTs of Green function and derivative

    ALLOCATE ( mon_ver_ker_global_fft ( 2 * M, 2 * N ) )
  
    ALLOCATE ( dip_ver_ker_global_fft ( 2 * M, 2 * N ) ) 
    

    CALL computeKernelFFT( phy % k, &
         M, N, &
         (/ geo % L_y / N_s , geo % L_x /), &
         mon_ver_ker_global_fft, &
         'M', &
         fft_ver_global )

    CALL computeKernelFFT( phy % k, &
         M, N, &
         (/ geo % L_y / N_s , geo % L_x /), &
         dip_ver_ker_global_fft, &
         'D', &
         fft_ver_global )


    !! Monopole and dipole equiv. sources

    ALLOCATE ( mon_equiv_ver ( 2 * M, 2 * N ) )

    ALLOCATE ( dip_equiv_ver ( 2 * M, 2 * N ) )


    !! Total field at equiv. sources and aux. matrix to compute
    !! Monopole/Dipole field

    ALLOCATE ( field_equiv_ver ( 2 * M, 2 * N ) )

    ALLOCATE ( aux_global_ver ( 2 * M, 2 * N ) )

    ! Local

    !! Size of matrices

    N_s = 2 * alg % N_src_ver

    M   = N_s

    N   = 2 

    !! 2D FFT Handler
    
    CALL createFFT_2D ( fft_ver_local, 2 * M, 2 * N )

    !! FFTs of Green function and derivative
    
    ALLOCATE ( mon_ver_ker_local_fft ( 2 * M, 2 * N ) )
  
    ALLOCATE ( dip_ver_ker_local_fft ( 2 * M, 2 * N ) ) 
    

    CALL computeKernelFFT ( &
         phy % k, &
         M, N, &
         (/ geo % L_y / N_s , geo % L_x /), &
         mon_ver_ker_local_fft, &
         'M', &
         fft_ver_local )

    CALL computeKernelFFT ( &
         phy % k, &
         M, N, &
         (/ geo % L_y / N_s , geo % L_x /), &
         dip_ver_ker_local_fft, &
         'D', &
         fft_ver_local )


    ALLOCATE ( aux_local_ver ( 2 * M, 2 * N ) )

    ALLOCATE ( local_ver ( 2 * M, 2 * N) )


    ! --------------------------------------------------------------------------------
    !
    ! Create field_permuted: array of size 4 * (N_src_hor + N_src_ver) containing the 
    ! field at the equivalent sources of a cell ordered in quadrants to apply 
    ! reduction of LeastSquare problem.
    !
    ! --------------------------------------------------------------------------------


    M = ( alg % N_src_hor + alg % N_src_ver )

    ALLOCATE ( field_permuted ( 4 * M ) )


    ! Create pointers to slices of field_permuted. Each pointer gives the slice that
    ! contains the field in the i-th quadrant


    field_per_1 => field_permuted ( 0 * M + 1 : 1 * M )

    field_per_2 => field_permuted ( 1 * M + 1 : 2 * M )

    field_per_3 => field_permuted ( 2 * M + 1 : 3 * M )

    field_per_4 => field_permuted ( 3 * M + 1 : 4 * M )


  END SUBROUTINE initForwardMap



  SUBROUTINE createArrayGeometry ( phy, geo, alg ) 


    TYPE (phy_Parameters) :: phy
    
    TYPE (geo_Parameters) :: geo

    TYPE (alg_Parameters) :: alg


    INTEGER :: n, m, l


    DO m = 1, geo % N_col

       DO n = 1, geo % N_row          

          CALL createObstacle ( &
               obstacleArray ( n, m ), &
               0, &
               (/ geo % radius, geo % L_x * m, geo % L_y * n /), &
               alg % N_dis ) 

          CALL createCell ( &
               cellArray ( n, m ), &
               geo % L_x * m, &
               geo % L_y * n, &
               obstacleArray ( n, m ), &
               psi (:, n, m), &
               proj_cell_hor, &
               proj_cell_ver, &
               interp_cell, &
               phy % k )

       END DO

    END DO


     DO l = 1, geo % N_def

        n = geo % defects_location ( l, 1 )

        m = geo % defects_location ( l, 2 )

        NULLIFY ( cellArray ( n, m ) % innerObstacle )

     END DO


  END SUBROUTINE createArrayGeometry


  
  SUBROUTINE computeKernelFFT ( k, N_row, N_col, step, mat_out, pole_flag, fft_handler)

    ! orient_flag = 'H' or 'V' for horizontal or vertical cell
    ! pole_flag = 'M' or 'D' for monopole or dipole kernel
    

    REAL(8) :: k

    INTEGER :: N_row, N_col

    REAL(8),DIMENSION(:) :: step

    COMPLEX(8),DIMENSION(:,:) :: mat_out

    CHARACTER :: pole_flag
    
    TYPE (FFT_2D) :: fft_handler

    
    INTEGER :: n, m

    REAL(8),DIMENSION(:),ALLOCATABLE :: src_1, src_2

    COMPLEX(8),DIMENSION(:,:),ALLOCATABLE :: kernel


    ALLOCATE ( kernel ( 2 * N_row, 2 * N_col ) )

    ALLOCATE ( src_1 ( 2 * N_row ) )

    ALLOCATE ( src_2 ( 2 * N_col ) )


    src_1 ( 0 * N_row + 1 : 1 * N_row ) = step(1) * (/ ( n, n = 0, N_row - 1 ) /)
   
    src_2 ( 0 * N_col + 1 : 1 * N_col ) = step(2) * (/ ( m, m = 0, N_col - 1 ) /)


    src_1 ( 1 * N_row + 1 : 2 * N_row ) = step(1) * (/ ( n, n = - N_row, -1 ) /)

    src_2 ( 1 * N_col + 1 : 2 * N_col ) = step(2) * (/ ( m, m = - N_col, -1 ) /)


    IF ( pole_flag == 'M' ) THEN


          FORALL ( n = 1 : 2 * N_row, m = 1 : 2 * N_col )

             kernel(n,m) = G_0 ( k, src_1(n), src_2(m) )
             
          END FORALL


    ELSE

       
       FORALL ( n = 1 : 2 * N_row, m = 1 : 2 * N_col )

          kernel(n,m) = DG_0 ( k, src_1(n), src_2(m), 0.0d0, 1.0d0 )
             
       END FORALL


    END IF


    ! Singular term src_1 = 0, src_2 = 0

    kernel ( 1, 1 ) = 0
    

    ! Compute FFT
 
    CALL directFFT ( fft_handler, kernel )

    mat_out = fft_handler % fft_mat


    DEALLOCATE ( src_1, src_2, kernel )


  END SUBROUTINE computeKernelFFT


  SUBROUTINE StoreEquivalentSourceAmplitude ( x, weights, N_src, type_flag, pole_flag)

    ! Stores the values of the equivalent sources in weights.
    ! mon_weights,dip_weights contain the weights of all the equivalent sources
    ! to perform, afterwards, the FFT. 

    ! Read notes for a better comprehension of the order in which the weights 
    ! are stored


    COMPLEX(8),DIMENSION(:),TARGET :: x

    COMPLEX(8),DIMENSION(:,:),TARGET :: weights

    INTEGER :: N_src

    CHARACTER :: type_flag, pole_flag


    REAL(8) :: sgn

    COMPLEX(8),DIMENSION(:),POINTER :: x_1, x_2, x_3, x_4

    COMPLEX(8),DIMENSION(:),POINTER :: panel_1, panel_2, panel_3, panel_4


    IF ( pole_flag == 'M' ) THEN
       

       x_1 => x ( 0 * N_src + 1 : 1 * N_src )

       x_4 => x ( 6 * N_src + 1 : 7 * N_src )

       x_2 => x ( 3 * N_src : 2 * N_src + 1 : -1 )

       x_3 => x ( 5 * N_src : 4 * N_src + 1 : -1 )

       sgn = 1.0d0

    ELSE 


       x_1 => x ( 1 * N_src + 1 : 2 * N_src )

       x_4 => x ( 7 * N_src + 1 : 8 * N_src )

       x_2 => x ( 4 * N_src : 3 * N_src + 1 : -1 )

       x_3 => x ( 6 * N_src : 5 * N_src + 1 : -1 )

       sgn = -1.0d0

    END IF


    IF ( type_flag == 'H' ) THEN
       

       panel_1 => weights ( 1 * N_src + 1 : 2 * N_src, 2 )

       panel_2 => weights ( 0 * N_src + 1 : 1 * N_src, 2 )

       panel_3 => weights ( 0 * N_src + 1 : 1 * N_src, 1 )

       panel_4 => weights ( 1 * N_src + 1 : 2 * N_src, 1 )

              
    ELSE IF ( type_flag == 'V' ) THEN


       panel_1 => weights ( 0 * N_src + 1 : 1 * N_src, 2 )

       panel_2 => weights ( 1 * N_src + 1 : 2 * N_src, 2 )

       panel_3 => weights ( 1 * N_src + 1 : 2 * N_src, 1 )

       panel_4 => weights ( 0 * N_src + 1 : 1 * N_src, 1 )


    END IF


    panel_1 = panel_1 + x_1

    panel_2 = panel_2 + x_2

    panel_3 = panel_3 + sgn * x_3

    panel_4 = panel_4 + sgn * x_4


    NULLIFY ( x_1, x_2, x_3, x_4 )

    NULLIFY ( panel_1, panel_2, panel_3, panel_4 )


  END SUBROUTINE StoreEquivalentSourceAmplitude



  SUBROUTINE computeEquivalentField_FFT ( equiv, ker_fft, field, fft_handle )


    COMPLEX(8),DIMENSION(:,:) :: equiv

    COMPLEX(8),DIMENSION(:,:) :: ker_fft

    COMPLEX(8),DIMENSION(:,:) :: field

    TYPE(FFT_2D) :: fft_handle


    CALL directFFT ( fft_handle, equiv )

    field = ker_fft * fft_handle % fft_mat
   
    CALL inverseFFT ( fft_handle, field )

    field = fft_handle % fft_mat


  END SUBROUTINE computeEquivalentField_FFT


  
  SUBROUTINE PermuteCellField (field_hor, field_ver, N_h, N_v)


    COMPLEX(8),DIMENSION(:,:) :: field_hor, field_ver

    INTEGER :: N_h, N_v


    field_per_1 ( 1 : N_h ) = field_hor ( 1 * N_h + 1 : 2 * N_h, 2 )

    field_per_2 ( 1 : N_h ) = field_hor ( 0 * N_h + 1 : 1 * N_h, 2 )

    field_per_3 ( 1 : N_h ) = field_hor ( 0 * N_h + 1 : 1 * N_h, 1 )

    field_per_4 ( 1 : N_h ) = field_hor ( 1 * N_h + 1 : 2 * N_h, 1 )


    field_per_1 ( N_h + 1 : N_h + N_v ) = field_ver ( 1 * N_v + 1 : 2 * N_v, 2 )

    field_per_2 ( N_h + 1 : N_h + N_v ) = field_ver ( 1 * N_v + 1 : 2 * N_v, 1 )

    field_per_3 ( N_h + 1 : N_h + N_v ) = field_ver ( 0 * N_v + 1 : 1 * N_v, 1 )

    field_per_4 ( N_h + 1 : N_h + N_v ) = field_ver ( 0 * N_v + 1 : 1 * N_v, 2 )


    ! Invert order

    field_per_2 ( 1 : N_h ) = field_per_2 ( N_h : 1 : -1 )

    field_per_3 ( 1 : N_h ) = field_per_3 ( N_h : 1 : -1 )


    field_per_1 ( N_h + 1 : N_h + N_v ) = field_per_1 ( N_h + N_v :  N_h + 1 : - 1 )

    field_per_2 ( N_h + 1 : N_h + N_v ) = field_per_2 ( N_h + N_v :  N_h + 1 : - 1 )


  END SUBROUTINE PermuteCellField


  SUBROUTINE ObstacleToCellMap ( cell_mn, refCell_hor, refCell_ver, x_h, x_v )

    ! Computes the horizontal and vertical equivalent sources of the mn-th cell
    ! To represent the field of the inner obstacle
    

    TYPE(Cell) :: cell_mn

    TYPE(ProjectionReferenceCell) :: refCell_hor, refCell_ver

    COMPLEX(8),DIMENSION(:) :: x_h, x_v

    
    COMPLEX(8),DIMENSION(:),ALLOCATABLE :: field_at_hor_coll, field_at_ver_coll
    
    
    ALLOCATE ( field_at_hor_coll ( 1 : 4 * refCell_hor % N_coll ) )

    ALLOCATE ( field_at_ver_coll ( 1 : 4 * refCell_ver % N_coll ) )


    CALL computeFieldInCollocationPoints ( cell_mn, field_at_hor_coll, field_at_ver_coll )

    CALL computeEquivalentSourceAmplitude ( refCell_hor, field_at_hor_coll, x_h)

    CALL computeEquivalentSourceAmplitude ( refCell_ver, field_at_ver_coll, x_v)

    
    DEALLOCATE ( field_at_hor_coll, field_at_ver_coll )


  END SUBROUTINE ObstacleToCellMap



  SUBROUTINE CellToObstacleMap ( cell_mn, field, res )

    
    TYPE (Cell) :: cell_mn

    COMPLEX(8),DIMENSION(:) :: field

    COMPLEX(8),DIMENSION(:) :: res

    
    CALL computePlaneWaveWeights ( interp_cell, field, xi )

    CALL MatVecMultiply ( cell_mn % cellToObsMatrix, xi, res )


  END SUBROUTINE CellToObstacleMap



  SUBROUTINE ForwardMap (phy, geo, alg)

    
    TYPE ( phy_Parameters ) :: phy

    TYPE ( geo_Parameters ) :: geo

    TYPE ( alg_Parameters ) :: alg
    

    COMPLEX(8),DIMENSION(:,:),POINTER :: field_hor, field_ver

    INTEGER :: n, m, Left_ind, Right_ind, N_s


    ! --------------------------------------------------------------------------------
    !      
    ! Initialize equivalent sources
    !
    ! --------------------------------------------------------------------------------

    mon_equiv_hor = 0.0d0
    
    dip_equiv_hor = 0.0d0
    

    mon_equiv_ver = 0.0d0 

    dip_equiv_ver = 0.0d0


    ! --------------------------------------------------------------------------------
    !
    ! Loop over all Cells and obtain equivalent sources
    !
    ! --------------------------------------------------------------------------------


    DO m = 1, geo % N_row

       DO n = 1, geo % N_col

          
          IF ( ASSOCIATED ( cellArray ( n, m ) % innerObstacle ) ) THEN
          

             CALL ObstacleToCellMap ( &
                  cellArray ( n, m ), &
                  proj_cell_hor, &
                  proj_cell_ver, x_hor, &
                  x_ver )        

             ! Store Cell's equivalent horizontal sources
             
             Left_ind  = ( 2 * alg % N_src_hor ) * n + 1

             Right_ind = ( 2 * alg % N_src_hor ) * ( n + 1 )


             CALL StoreEquivalentSourceAmplitude ( &
                  x_hor, &
                  cellArray ( n, m ) % mon_equiv_hor, &
                  alg % N_src_hor, &
                  'H', &
                  'M' )
             
             CALL StoreEquivalentSourceAmplitude ( &
                  x_hor, &
                  cellArray ( n, m ) % dip_equiv_hor &
                  , alg % N_src_hor, &
                  'H', &
                  'D' )


             !! Add to global sources

             mon_equiv_hor ( Left_ind : Right_ind, m + 1 : m + 2 ) = &
                  mon_equiv_hor ( Left_ind : Right_ind, m + 1 : m + 2 ) + &
                  cellArray ( n, m ) % mon_equiv_hor

             dip_equiv_hor ( Left_ind : Right_ind, m + 1 : m + 2 ) = &
                  dip_equiv_hor ( Left_ind : Right_ind, m + 1 : m + 2 ) + &
                  cellArray ( n, m ) % dip_equiv_hor

             
          
             ! Store Cell's equivalent vertical sources


             Left_ind  = ( 2 * alg % N_src_ver ) * m + 1

             Right_ind = ( 2 * alg % N_src_ver ) * ( m + 1 )


             CALL StoreEquivalentSourceAmplitude ( &
                  x_ver, &
                  cellArray ( n, m ) % mon_equiv_ver, &
                  alg % N_src_ver, &
                  'V', &
                  'M' )
             
             CALL StoreEquivalentSourceAmplitude ( &
                  x_ver, &
                  cellArray ( n, m ) % dip_equiv_ver, &
                  alg % N_src_ver, &
                  'V', &
                  'D' )

             !! Add to global sources

             mon_equiv_ver ( Left_ind : Right_ind, n + 1 : n + 2 ) = &
                  mon_equiv_ver ( Left_ind : Right_ind, n + 1 : n + 2 ) + &
                  cellArray ( n, m ) % mon_equiv_ver

             dip_equiv_ver ( Left_ind : Right_ind, n + 1 : n + 2 ) = &
                  dip_equiv_ver ( Left_ind : Right_ind, n + 1 : n + 2 ) + &
                  cellArray ( n, m ) % dip_equiv_ver


          END IF


       END DO

    END DO


    ! --------------------------------------------------------------------------------
    !
    ! Compute field at horizontal equivalent sources using FFT
    !
    ! --------------------------------------------------------------------------------

    ! Monopole contributions

    CALL computeEquivalentField_FFT ( &
         mon_equiv_hor, &
         mon_hor_ker_global_fft, &
         aux_global_hor, &
         fft_hor_global )

    field_equiv_hor = aux_global_hor


    ! Dipole contributions

    CALL computeEquivalentField_FFT ( &
         dip_equiv_hor, &
         dip_hor_ker_global_fft, &
         aux_global_hor, &
         fft_hor_global )

    field_equiv_hor = field_equiv_hor + aux_global_hor


    ! --------------------------------------------------------------------------------
    !
    ! Compute field at vertical equivalent sources using a global FFT
    !
    ! --------------------------------------------------------------------------------

    ! Monopole contributions

    CALL computeEquivalentField_FFT ( &
         mon_equiv_ver, &
         mon_ver_ker_global_fft, &
         aux_global_ver, &
         fft_ver_global )

    field_equiv_ver = aux_global_ver


    ! Dipole contributions

    CALL computeEquivalentField_FFT ( &
         dip_equiv_ver, &
         dip_ver_ker_global_fft, &
         aux_global_ver, &
         fft_ver_global )

    field_equiv_ver = field_equiv_ver + aux_global_ver


    ! --------------------------------------------------------------------------------
    !
    ! Iterate through all cells and substract local (incorrect) contributions
    !
    ! --------------------------------------------------------------------------------


    DO m = 1, geo % N_row


       DO n = 1, geo % N_col
          
          ! Horizontal sources

          N_s = 2 * alg % N_src_hor

          Left_ind  = N_s * n + 1
          
          Right_ind = N_s * ( n + 1 )

          !! Monopoles 

          local_hor ( 1 : N_s, 1 : 2 ) = cellArray ( n, m ) % mon_equiv_hor

          CALL computeEquivalentField_FFT ( &
               local_hor, &
               mon_hor_ker_local_fft, &
               aux_local_hor, &
               fft_hor_local )
          
          field_equiv_hor ( Left_ind : Right_ind, m + 1 : m + 2 ) = &
               field_equiv_hor ( Left_ind : Right_ind, m + 1 : m + 2 ) - &
               aux_local_hor ( 1 : N_s, 1 : 2)
          
          !! Dipoles

          local_hor ( 1 : N_s, 1 : 2 ) = cellArray ( n, m ) % dip_equiv_hor

          CALL computeEquivalentField_FFT ( &
               local_hor, &
               dip_hor_ker_local_fft, &
               aux_local_hor, &
               fft_hor_local )

          field_equiv_hor ( Left_ind : Right_ind, m + 1 : m + 2 ) = &
               field_equiv_hor ( Left_ind : Right_ind, m + 1 : m + 2 ) - &
               aux_local_hor ( 1 : N_s, 1 : 2)


          ! Vertical sources

          N_s = 2 * alg % N_src_ver

          Left_ind  = N_s * m + 1
          
          Right_ind = N_s * ( m + 1 ) 

          !! Monopoles 

          local_ver ( 1 : N_s, 1 : 2 ) = cellArray ( n, m ) % mon_equiv_ver

          CALL computeEquivalentField_FFT ( &
               local_ver, &
               mon_ver_ker_local_fft, &
               aux_local_ver, &
               fft_ver_local )


          field_equiv_ver ( Left_ind : Right_ind, n + 1 : n + 2 ) = &
               field_equiv_ver ( Left_ind : Right_ind, n + 1 : n + 2 ) - &
               aux_local_ver ( 1 : N_s, 1 : 2)

          
          !! Dipoles

          local_ver ( 1 : N_s, 1 : 2 ) = cellArray ( n, m ) % dip_equiv_ver

          CALL computeEquivalentField_FFT ( &
               local_ver, &
               dip_ver_ker_local_fft, &
               aux_local_ver, &
               fft_ver_local )

          field_equiv_ver ( Left_ind : Right_ind, n + 1 : n + 2 ) = &
               field_equiv_ver ( Left_ind : Right_ind, n + 1 : n + 2 ) - &
               aux_local_ver( 1 : N_s, 1 : 2)

          
       END DO


    END DO


    ! --------------------------------------------------------------------------------
    !
    ! Iterate through all cells and obtain outer field at the inner obstacle.
    !
    ! --------------------------------------------------------------------------------

    ! DO m = 1, geo % N_row


    !    DO n = 1, geo % N_col

          
    !       IF ( ASSOCIATED ( cellArray ( n, m ) % innerObstacle ) ) THEN

             
    !          field_hor => field_equiv_hor ( ( 2 * alg % N_src_hor ) * n + 1 : ( 2 * alg % N_src_hor ) * (n+1), m+1 : m+2 )

    !          field_ver => field_equiv_ver ( ( 2 * alg % N_src_ver ) * m + 1 : ( 2 * alg % N_src_ver ) * (m+1), n+1 : n+2 )
             
    !          CALL PermuteCellField ( field_hor, field_ver, alg % N_src_hor, alg % N_src_ver )

    !          CALL CellToObstacleMap ( cellArray (m,n) , field_permuted, res ( :, m, n ) )
             
    !       END IF


    !    END DO


    ! END DO
    

    
  END SUBROUTINE ForwardMap



!   SUBROUTINE TestEquivalentSources ( phy, geo, alg ) 


!     TYPE ( phy_Parameters ) :: phy

!     TYPE ( geo_Parameters ) :: geo

!     TYPE ( alg_Parameters ) :: alg


!     COMPLEX (8) :: field_1(1)

!     COMPLEX(8),DIMENSION(1) :: field_2
    
!     REAL(8) :: tar_x, tar_y

!     REAL(8),DIMENSION(:),ALLOCATABLE :: src_x, src_y

!     INTEGER :: n, m, n_0, m_0
    
!     COMPLEX(8), DIMENSION(:,:), ALLOCATABLE :: matrix


!     ALLOCATE ( matrix ( 1, alg % N_dis ) )
    
!     ALLOCATE ( src_y ( 2 * alg % N_src_ver * ( geo % N_row + 2 ) ) )

!     ALLOCATE ( src_x ( geo % N_col + 3 ) )

    
!     src_y = geo % L_y * (/ ( ( 2.0d0 * n + 1.0d0 ) / ( 4 * alg % N_src_ver ) - 0.5d0 , n = 0 , 2 * alg % N_src_ver * ( geo % N_row + 2 ) - 1 ) /)

!     src_x = geo % L_x * (/ ( m, m = 0, geo % N_col + 2 ) /) - geo % L_x / 2 


!     n_0 = 1

!     m_0 = 1

!     tar_x = src_x (m_0)

!     tar_y = src_y (n_0)
    

!     field_1 = 0.0d0


!     DO m = 1 , geo % N_col + 3 

!        DO n = 1 , 2 * alg % N_src_ver * ( geo % N_row + 2 )
        
!            IF ( n .NE. n_0 .OR. m .NE. m_0 ) THEN
              
!               field_1 = field_1 + G_0 ( phy % k, tar_x - src_x (m), tar_y - src_y (n) ) * mon_equiv_ver ( n, m )

!               field_1 = field_1 + DG_0 ( phy % k, tar_x - src_x (m), tar_y - src_y (n), 1.0d0, 0.0d0 ) * dip_equiv_ver (n, m)

!            END IF

!         END DO

!      END DO

!    write(*,*) field_1




!   DO m = 1, geo % N_row

!      DO n = 1, geo % N_col

!           if (ASSOCIATED ( cellArray (m, n) % innerObstacle ) ) THEN
          
!              CALL createEvalFieldAtPointsMatrix ( matrix, cellArray ( m, n ) % innerObstacle, (/tar_x/), (/tar_y/), phy % k )
          
!              CALL MatVecMultiply ( matrix, cellArray ( m, n ) % psi, field_2)

!              field_1 = field_1 - field_2

!           END IF

!        END DO

!     END DO

! !    write(*,*) field_1

!    write(*,*) abs ( field_1 )


!   END SUBROUTINE TestEquivalentSources




!   SUBROUTINE TestEquivalentSources ( phy, geo, alg ) 


!     TYPE ( phy_Parameters ) :: phy

!     TYPE ( geo_Parameters ) :: geo

!     TYPE ( alg_Parameters ) :: alg


!     COMPLEX (8) :: field_1(1)

!     COMPLEX(8),DIMENSION(1) :: field_2
    
!     REAL(8) :: tar_x, tar_y

!     REAL(8),DIMENSION(:),ALLOCATABLE :: src_x, src_y

!     INTEGER :: n, m, n_0, m_0, N_s

!     INTEGER :: n_cell, m_cell, n_point, m_point
    
!     COMPLEX(8), DIMENSION(:,:), ALLOCATABLE :: matrix


!     N_s = 2 * alg % N_src_hor

!     ALLOCATE ( matrix ( 1, alg % N_dis ) )
    
!     ALLOCATE ( src_x ( N_s * ( geo % N_col + 2 ) ) )

!     ALLOCATE ( src_y ( geo % N_row + 3 ) )

    
!     src_x = geo % L_x * (/ ( ( 2.0d0 * n + 1.0d0 ) / ( 2 * N_s ) - 0.5d0 , n = 0 , 2 * N_s * ( geo % N_col + 2 ) - 1 ) /)

!     src_y = geo % L_y * (/ ( m, m = 0, geo % N_row + 2 ) /) - geo % L_y / 2 


!     n_cell = 1

!     m_cell = 5

!     n_point = 3

!     m_point = 2

!     n_0 = n_cell * N_s + n_point

!     m_0 = m_cell + m_point

    
!     tar_x = src_x (n_0)

!     tar_y = src_y (m_0)
    

!     field_1 = 0.0d0


!     ! GLOBAL CONVOLUTION
!     DO n = 1 , N_s * ( geo % N_col + 2 )

!        DO m = 1 , geo % N_row + 3 

!           IF ( n .NE. n_0 .OR. m .NE. m_0 ) THEN

!              field_1 = field_1 + G_0 ( phy % k, tar_x - src_x (n), tar_y - src_y (m) ) * mon_equiv_hor ( n, m )

!              field_1 = field_1 + DG_0 ( phy % k, tar_x - src_x (n), tar_y - src_y (m), 0.0d0, 1.0d0 ) * dip_equiv_hor (n, m)

!           END IF

!        END DO

!     END DO

!     ! LOCAL CONVOLUTION
!     DO n = 1,  N_s

!        DO m = 1, 2

!           IF ( n .NE. n_point .OR. m .NE. m_point ) THEN

!              field_1 = field_1 - G_0 ( phy % k, tar_x - src_x (n + n_cell * N_s), tar_y - src_y (m+m_cell) ) * cellArray ( n_cell, m_cell ) % mon_equiv_hor ( n, m )

!              field_1 = field_1 - DG_0 ( phy % k, tar_x - src_x (n + n_cell * N_s), tar_y - src_y (m+m_cell), 0.0d0, 1.0d0 ) * cellArray ( n_cell, m_cell ) % dip_equiv_hor (n, m)

!           END IF

!        END DO

!     END DO

!     write(*,*) abs(field_equiv_hor (n_0, m_0)- field_1)

! !    field_1 = 0.0d0

!     DO n = 1, geo % N_col

!        DO m = 1, geo % N_row

!           if ( n == n_cell .AND. m == m_cell ) THEN

!           ELSE

!              if (ASSOCIATED ( cellArray (m, n) % innerObstacle ) ) THEN

!                 CALL createEvalFieldAtPointsMatrix ( matrix, cellArray ( m, n ) % innerObstacle, (/tar_x/), (/tar_y/), phy % k )
          
!                 CALL MatVecMultiply ( matrix, cellArray ( m, n ) % psi, field_2)

!                 field_1 = field_1 - field_2

!              END IF

!           END IF

!        END DO

!     END DO

!     write(*,*) abs(field_1)


!  END SUBROUTINE TestEquivalentSources 



  SUBROUTINE TestEquivalentSources ( phy, geo, alg ) 


    TYPE ( phy_Parameters ) :: phy

    TYPE ( geo_Parameters ) :: geo

    TYPE ( alg_Parameters ) :: alg


    COMPLEX (8) :: field_1(1)

    COMPLEX(8),DIMENSION(1) :: field_2
    
    REAL(8) :: tar_x, tar_y

    REAL(8),DIMENSION(:),ALLOCATABLE :: src_x, src_y

    INTEGER :: n, m, n_0, m_0, N_s

    INTEGER :: n_cell, m_cell, n_point, m_point
    
    COMPLEX(8), DIMENSION(:,:), ALLOCATABLE :: matrix


    N_s = 2 * alg % N_src_ver

    ALLOCATE ( matrix ( 1, alg % N_dis ) )
    
    ALLOCATE ( src_y ( N_s * ( geo % N_row + 2 ) ) )

    ALLOCATE ( src_x ( geo % N_col + 3 ) )

    
    src_y = geo % L_y * (/ ( ( 2.0d0 * n + 1.0d0 ) / ( 2 * N_s ) - 0.5d0 , n = 0 , N_s * ( geo % N_row + 2 ) - 1 ) /)

    src_x = geo % L_x * (/ ( m, m = 0, geo % N_col + 2 ) /) - geo % L_x / 2 


    n_cell = 5

    m_cell = 5

    n_point = 2

    m_point = 1

    n_0 = n_cell * N_s + n_point

    m_0 = m_cell + m_point

    
    tar_x = src_x (m_0)

    tar_y = src_y (n_0)
    

    field_1 = 0.0d0


    ! GLOBAL CONVOLUTION


    DO m = 1 , geo % N_col + 3 

       DO n = 1 , N_s * ( geo % N_row + 2 )
          
          IF ( n .NE. n_0 .OR. m .NE. m_0 ) THEN

             field_1 = field_1 + G_0 ( phy % k, tar_x - src_x (m), tar_y - src_y (n) ) * mon_equiv_ver ( n, m )

             field_1 = field_1 + DG_0 ( phy % k, tar_x - src_x (m), tar_y - src_y (n), 1.0d0, 0.0d0 ) * dip_equiv_ver (n, m)

          END IF

       END DO

    END DO
    write(*,*) field_1
    ! LOCAL CONVOLUTION
    DO n = 1,  N_s

       DO m = 1, 2

          IF ( n .NE. n_point .OR. m .NE. m_point ) THEN

             field_1 = field_1 - G_0 ( phy % k, tar_x - src_x (m + m_cell ), tar_y - src_y (n + n_cell* N_s) ) * cellArray ( m_cell, n_cell ) % mon_equiv_ver ( n, m )

             field_1 = field_1 - DG_0 ( phy % k, tar_x - src_x (m + m_cell ), tar_y - src_y (n + n_cell * N_s), 1.0d0, 0.0d0 ) * cellArray ( m_cell, n_cell ) % dip_equiv_ver (n, m)

          END IF

       END DO

    END DO


   write(*,*) abs ( field_equiv_ver (n_0, m_0) - field_1 )

!    field_1 = 0.0d0

    DO n = 1, geo % N_col

       DO m = 1, geo % N_row

          if ( n == m_cell .AND. m == n_cell ) THEN

          ELSE

             if (ASSOCIATED ( cellArray (m, n) % innerObstacle ) ) THEN

                CALL createEvalFieldAtPointsMatrix ( matrix, cellArray ( m, n ) % innerObstacle, (/tar_x/), (/tar_y/), phy % k )
          
                CALL MatVecMultiply ( matrix, cellArray ( m, n ) % psi, field_2)

                field_1 = field_1 - field_2

             END IF

          END IF

       END DO

    END DO

    write(*,*) abs(field_1)


 END SUBROUTINE TestEquivalentSources 



 !  SUBROUTINE TestEquivalentSources ( phy, geo, alg ) 


 !    TYPE ( phy_Parameters ) :: phy

 !    TYPE ( geo_Parameters ) :: geo

 !    TYPE ( alg_Parameters ) :: alg


 !    COMPLEX(8),DIMENSION(:),ALLOCATABLE :: field_1, field_2

    
 !    REAL(8),DIMENSION(:),ALLOCATABLE :: tar_x, tar_y

 !    INTEGER :: n, m, n_0, m_0, l, Left_ind, Right_ind, row, col
    
 !    COMPLEX(8), DIMENSION(:,:), ALLOCATABLE :: matrix

 !    ALLOCATE ( field_1 ( 2 * alg % N_src_ver ) )

 !    ALLOCATE ( field_2 ( 2 * alg % N_src_ver ) )
    
 !    ALLOCATE ( matrix ( 1, alg % N_dis ) )
    
 !    ALLOCATE ( tar_y ( 2 * alg % N_src_ver ) )

 !    ALLOCATE ( tar_x ( 2 * alg % N_src_ver ) )

    
 !    row = 1

 !    col = 1


 !    tar_y = geo % L_y * (/ ( ( 2.0d0 * n + 1.0d0 ) / ( 4 * alg % N_src_ver ) - 0.5d0 , n = 2 * alg % N_src_ver * row + 1, 2 * alg % N_src_ver * ( row + 1 ) ) /)

 !    tar_x = geo % L_x * col - geo % L_x / 2 



 !    field_1 = 0.0d0

 !    DO n = 1, geo % N_col

 !       DO m = 1, geo % N_row

 !          if ( n == col .AND. m == row ) THEN

 !          ELSE

 !             if (ASSOCIATED ( cellArray (m, n) % innerObstacle ) ) THEN

 !                CALL createEvalFieldAtPointsMatrix ( matrix, cellArray ( m, n ) % innerObstacle, tar_x, tar_y,  phy % k )

 !                CALL MatVecMultiply ( matrix, cellArray ( m, n ) % psi, field_2)

 !                field_1 = field_1 + field_2

 !             END IF

 !          END IF

 !       END DO

 !    END DO



 !    Left_ind  = 2 * alg % N_src_ver * row + 1
          
 !    Right_ind = 2 * alg % N_src_ver * (row + 1)

 !    write(*,*) abs (field_equiv_ver ( Left_ind : Right_ind, col + 1 : col + 1 ) - RESHAPE ( field_1, (/ size(field_1), 1 /) ) )



 ! END SUBROUTINE TestEquivalentSources 


  
  SUBROUTINE destroyForwardMap ()

    
    DEALLOCATE ( cellArray, obstacleArray ) 

    DEALLOCATE ( psi, res )


  END SUBROUTINE destroyForwardMap



END MODULE



