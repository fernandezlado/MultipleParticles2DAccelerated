module modForwardMap



  use modParameters

  use modFFT

  use modCell

  use modProjectionReferenceCell



  implicit none


  
  type (ProjectionReferenceCell) :: proj_cell_hor, proj_cell_ver

  type (InterpolationReferenceCell) :: interp_cell

  complex(8),dimension(:),allocatable :: field_hor_coll, field_ver_coll

  complex(8),dimension(:),allocatable :: x_hor, x_ver, xi

 
  type(Cell), dimension(:,:), allocatable :: cellArray

  type(Obstacle), dimension(:,:), target, allocatable :: obstacleArray


  complex(8),dimension(:,:,:),allocatable :: psi, res



  complex(8),dimension(:,:),allocatable :: mon_equiv_hor, mon_equiv_ver

  complex(8),dimension(:,:),allocatable :: dip_equiv_hor, dip_equiv_ver

  complex(8),dimension(:,:),allocatable :: local_hor, local_ver

  
  type (FFT_2D) :: fft_hor_local, fft_ver_local

  type (FFT_2D) :: fft_hor_global, fft_ver_global

  complex(8),dimension(:,:),allocatable :: mon_hor_ker_global_fft

  complex(8),dimension(:,:),allocatable :: dip_hor_ker_global_fft

  complex(8),dimension(:,:),allocatable :: mon_ver_ker_global_fft
  
  complex(8),dimension(:,:),allocatable :: dip_ver_ker_global_fft

  complex(8),dimension(:,:),allocatable :: mon_hor_ker_local_fft

  complex(8),dimension(:,:),allocatable :: dip_hor_ker_local_fft

  complex(8),dimension(:,:),allocatable :: mon_ver_ker_local_fft
  
  complex(8),dimension(:,:),allocatable :: dip_ver_ker_local_fft


  complex(8),dimension(:,:),allocatable,target :: field_equiv_hor, field_equiv_ver

  complex(8),dimension(:,:),allocatable :: aux_global_hor, aux_global_ver

  complex(8),dimension(:,:),allocatable :: aux_local_hor, aux_local_ver

  complex(8),dimension(:),allocatable,target :: field_permuted

  complex(8),dimension(:),pointer :: field_per_1, field_per_2, field_per_3, field_per_4



!  private :: proj_cell_hor, proj_cell_ver, interp_cell

  private :: cellArray, obstacleArray

  private :: mon_equiv_hor, mon_equiv_ver

  private :: dip_equiv_hor, dip_equiv_ver

  private :: mon_hor_ker_global_fft, dip_hor_ker_global_fft

  private :: mon_ver_ker_global_fft, dip_ver_ker_global_fft

  private :: x_hor, x_ver, xi



contains


  
  subroutine initForwardMap ( phy, geo, alg )


    type (phy_Parameters) :: phy
    
    type (geo_Parameters) :: geo

    type (alg_Parameters) :: alg


    integer :: N, M, N_s
    

    ! --------------------------------------------------------------------------------
    !
    ! Create array to store densities and the result of the Forward Map
    !
    ! --------------------------------------------------------------------------------


    allocate ( psi ( alg % N_dis, geo % N_row, geo % N_col ) )

    allocate ( res ( alg % N_dis, geo % N_row, geo % N_col) )
    

    psi = 1.0d0

    res = 1.0d0
    

    ! --------------------------------------------------------------------------------
    !
    ! Create Projection and Interpolation reference Cells and x_hor, x_ver, xi that
    ! contain the solutions of the least squares problem to project and interpolate
    ! the field with the equivalent sources.
    !
    ! --------------------------------------------------------------------------------


    call createProjectionReferenceCell ( &
         proj_cell_hor, &
         alg % N_src_hor, &
         alg % N_coll_hor, &
         geo % L_x, &
         geo % L_y, &
         phy % k, & 
         'H' )


    call createProjectionReferenceCell ( &
         proj_cell_ver, &
         alg % N_src_ver, &
         alg % N_coll_ver, &
         geo % L_x, &
         geo % L_y, &
         phy % k, &
         'V' )


    call createInterpolationReferenceCell ( &
         interp_cell, &
         alg % N_src_hor, &
         alg % N_src_ver, &
         alg % N_wave, &
         geo % L_x, &
         geo % L_y, &
         phy % k )


    ! Create x_hor, x_ver: solutions of least square problem to obtain equiv. sources

    allocate ( field_hor_coll ( 1 : 4 * proj_cell_hor % N_coll ) )

    allocate ( field_ver_coll ( 1 : 4 * proj_cell_ver % N_coll ) )

    allocate ( x_hor ( 2 * 4 * alg % N_src_hor ) )

    allocate ( x_ver ( 2 * 4 * alg % N_src_ver ) )

    ! Create xi : solution of least square problem to expand field as sum 
    ! of plane waves

    allocate ( xi ( 4 * alg % N_wave) )

    ! --------------------------------------------------------------------------------
    !
    ! Create Cell and Obstacles arrays
    !
    ! --------------------------------------------------------------------------------


    allocate ( cellArray ( 0 : geo % N_row + 1, 0 : geo % N_col + 1 ) )

    allocate ( obstacleArray ( 1 : geo % N_row, 1 : geo % N_col ) )    
    

    call createArrayGeometry ( phy, geo, alg )


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
    
    call createFFT_2D ( fft_hor_global, 2 * M, 2 * N )

    !! FFTs of Green function and derivative
    
    allocate ( mon_hor_ker_global_fft ( 2 * M, 2 * N ) )
  
    allocate ( dip_hor_ker_global_fft ( 2 * M, 2 * N ) ) 
    

    call computeKernelFFT ( &
         phy % k, &
         M, N, &
         (/ geo % L_x / N_s , geo % L_y /), &
         mon_hor_ker_global_fft, &
         'M', &
         fft_hor_global )

    call computeKernelFFT ( &
         phy % k, &
         M, N, &
         (/ geo % L_x / N_s , geo % L_y /), &
         dip_hor_ker_global_fft, &
         'D', &
         fft_hor_global )


    !! Monopole and dipole equiv. sources

    allocate ( mon_equiv_hor ( 2 * M, 2 * N ) )

    allocate ( dip_equiv_hor ( 2 * M, 2 * N ) )
    

    !! Total field at equiv. sources and aux. matrix to compute
    !! Monopole/Dipole field

    allocate ( field_equiv_hor ( 2 * M, 2 * N ) )

    allocate ( aux_global_hor ( 2 * M, 2 * N ) )

    ! Local

    !! Size of matrices

    N_s = 2 * alg % N_src_hor

    M   = 3 * N_s

    N   = 4

    !! 2D FFT Handler
    
    call createFFT_2D ( fft_hor_local, 2 * M, 2 * N )

    !! FFTs of Green function and derivative
    
    allocate ( mon_hor_ker_local_fft ( 2 * M, 2 * N ) )
  
    allocate ( dip_hor_ker_local_fft ( 2 * M, 2 * N ) ) 
    

    call computeKernelFFT ( &
         phy % k, &
         M, N, &
         (/ geo % L_x / N_s , geo % L_y /), &
         mon_hor_ker_local_fft, &
         'M', &
         fft_hor_local )

    call computeKernelFFT ( &
         phy % k, &
         M, N, &
         (/ geo % L_x / N_s , geo % L_y /), &
         dip_hor_ker_local_fft, &
         'D', &
         fft_hor_local )


    allocate ( aux_local_hor ( 2 * M, 2 * N ) )

    allocate ( local_hor ( 2 * M, 2 * N) )


    ! --------------------------------------------------------------------------------
    ! Vertical case
    ! --------------------------------------------------------------------------------

    ! Global 

    !! Size of matrices

    N_s = 2 * alg % N_src_ver

    M = N_s * ( geo % N_row + 2 )

    N = geo % N_col + 3

    !! 2D FFT Handler

    call createFFT_2D ( fft_ver_global, 2 * M, 2 * N )

    !! FFTs of Green function and derivative

    allocate ( mon_ver_ker_global_fft ( 2 * M, 2 * N ) )
  
    allocate ( dip_ver_ker_global_fft ( 2 * M, 2 * N ) ) 
    

    call computeKernelFFT( phy % k, &
         M, N, &
         (/ geo % L_y / N_s , geo % L_x /), &
         mon_ver_ker_global_fft, &
         'M', &
         fft_ver_global )

    call computeKernelFFT( phy % k, &
         M, N, &
         (/ geo % L_y / N_s , geo % L_x /), &
         dip_ver_ker_global_fft, &
         'D', &
         fft_ver_global )


    !! Monopole and dipole equiv. sources

    allocate ( mon_equiv_ver ( 2 * M, 2 * N ) )

    allocate ( dip_equiv_ver ( 2 * M, 2 * N ) )


    !! Total field at equiv. sources and aux. matrix to compute
    !! Monopole/Dipole field

    allocate ( field_equiv_ver ( 2 * M, 2 * N ) )

    allocate ( aux_global_ver ( 2 * M, 2 * N ) )

    ! Local

    !! Size of matrices

    N_s = 2 * alg % N_src_ver

    M   = 3 * N_s

    N   = 4

    !! 2D FFT Handler
    
    call createFFT_2D ( fft_ver_local, 2 * M, 2 * N )

    !! FFTs of Green function and derivative
    
    allocate ( mon_ver_ker_local_fft ( 2 * M, 2 * N ) )
  
    allocate ( dip_ver_ker_local_fft ( 2 * M, 2 * N ) ) 
    

    call computeKernelFFT ( &
         phy % k, &
         M, N, &
         (/ geo % L_y / N_s , geo % L_x /), &
         mon_ver_ker_local_fft, &
         'M', &
         fft_ver_local )

    call computeKernelFFT ( &
         phy % k, &
         M, N, &
         (/ geo % L_y / N_s , geo % L_x /), &
         dip_ver_ker_local_fft, &
         'D', &
         fft_ver_local )


    allocate ( aux_local_ver ( 2 * M, 2 * N ) )

    allocate ( local_ver ( 2 * M, 2 * N) )


    ! --------------------------------------------------------------------------------
    !
    ! Create field_permuted: array of size 4 * (N_src_hor + N_src_ver) containing the 
    ! field at the equivalent sources of a cell ordered in quadrants to apply 
    ! reduction of LeastSquare problem.
    !
    ! --------------------------------------------------------------------------------


    M = ( alg % N_src_hor + alg % N_src_ver )

    allocate ( field_permuted ( 4 * M ) )


    ! Create pointers to slices of field_permuted. Each pointer gives the slice that
    ! contains the field in the i-th quadrant


    field_per_1 => field_permuted ( 0 * M + 1 : 1 * M )

    field_per_2 => field_permuted ( 1 * M + 1 : 2 * M )

    field_per_3 => field_permuted ( 2 * M + 1 : 3 * M )

    field_per_4 => field_permuted ( 3 * M + 1 : 4 * M )


    write(*,*) "Forward map initiated"


  end subroutine initForwardMap



  subroutine createArrayGeometry ( phy, geo, alg ) 


    type (phy_Parameters) :: phy
    
    type (geo_Parameters) :: geo

    type (alg_Parameters) :: alg

    type (Obstacle), pointer :: obs_pointer
    
    integer :: n, m, l


    do n = 0, geo % N_col + 1

       do m = 0, geo % N_row + 1 


          if ( m >= 1 .and. m <= geo % N_row .and. n >= 1 .and. n <= geo % N_col ) then

             call createObstacle ( &
                  obstacleArray ( m, n ), &
                  1, &
                  (/ 0.5*geo % radius, geo % L_x * n, geo % L_y * m /), &
                  alg % N_dis )
 
             obs_pointer => obstacleArray ( m, n )

          else 
          
             obs_pointer => null()

          end if
       
          call createCell ( &
                  cellArray ( m, n ), &
                  geo % L_x * n, &
                  geo % L_y * m, &
                  obs_pointer, &
                  psi (:, m, n), &
                  proj_cell_hor, &
                  proj_cell_ver, &
                  interp_cell, &
                  phy % k )


       end do

    end do


    do l = 1, geo % N_def

       m = geo % defects_location ( l, 1 )
        
       n = geo % defects_location ( l, 2 )

       nullify ( cellArray ( m, n ) % innerObstacle )

    end do


  end subroutine createArrayGeometry


  
  subroutine computeKernelFFT ( k, N_row, N_col, step, mat_out, pole_flag, fft_handler)

    ! orient_flag = 'H' or 'V' for horizontal or vertical cell
    ! pole_flag = 'M' or 'D' for monopole or dipole kernel
    

    real(8) :: k

    integer :: N_row, N_col

    real(8),dimension(:) :: step

    complex(8),dimension(:,:) :: mat_out

    CHARACTER :: pole_flag
    
    type (FFT_2D) :: fft_handler

    
    integer :: n, m

    real(8),dimension(:),allocatable :: src_1, src_2

    complex(8),dimension(:,:),allocatable :: kernel


    allocate ( kernel ( 2 * N_row, 2 * N_col ) )

    allocate ( src_1 ( 2 * N_row ) )

    allocate ( src_2 ( 2 * N_col ) )


    src_1 ( 0 * N_row + 1 : 1 * N_row ) = step(1) * (/ ( n, n = 0, N_row - 1 ) /)
   
    src_2 ( 0 * N_col + 1 : 1 * N_col ) = step(2) * (/ ( m, m = 0, N_col - 1 ) /)


    src_1 ( 1 * N_row + 1 : 2 * N_row ) = step(1) * (/ ( n, n = - N_row, -1 ) /)

    src_2 ( 1 * N_col + 1 : 2 * N_col ) = step(2) * (/ ( m, m = - N_col, -1 ) /)


    if ( pole_flag == 'M' ) then


          forall ( n = 1 : 2 * N_row, m = 1 : 2 * N_col )

             kernel(n,m) = G_0 ( k, src_1(n), src_2(m) )
             
          end forall


    else

       
       forall ( n = 1 : 2 * N_row, m = 1 : 2 * N_col )

          kernel(n,m) = DG_0 ( k, src_1(n), src_2(m), 0.0d0, 1.0d0 )
             
       end forall


    end if


    ! Singular term src_1 = 0, src_2 = 0

    kernel ( 1, 1 ) = 0
    

    ! Compute FFT
 
    call directFFT ( fft_handler, kernel )

    mat_out = fft_handler % fft_mat


    deallocate ( src_1, src_2, kernel )


  end subroutine computeKernelFFT



  subroutine StoreEquivalentSourceAmplitude ( x, weights, N_src, type_flag, pole_flag)

    ! Stores the values of the equivalent sources in weights.
    ! mon_weights,dip_weights contain the weights of all the equivalent sources
    ! to perform, afterwards, the FFT. 

    ! Read notes for a better comprehension of the order in which the weights 
    ! are stored


    complex(8),dimension(:),target :: x

    complex(8),dimension(:,:),target :: weights

    integer :: N_src

    CHARACTER :: type_flag, pole_flag


    real(8) :: sgn

    complex(8),dimension(:),pointer :: x_1, x_2, x_3, x_4

    complex(8),dimension(:),pointer :: panel_1, panel_2, panel_3, panel_4


    if ( pole_flag == 'M' ) then
       

       x_1 => x ( 0 * N_src + 1 : 1 * N_src )

       x_4 => x ( 6 * N_src + 1 : 7 * N_src )

       x_2 => x ( 3 * N_src : 2 * N_src + 1 : -1 )

       x_3 => x ( 5 * N_src : 4 * N_src + 1 : -1 )

       sgn = 1.0d0

    else 


       x_1 => x ( 1 * N_src + 1 : 2 * N_src )

       x_4 => x ( 7 * N_src + 1 : 8 * N_src )

       x_2 => x ( 4 * N_src : 3 * N_src + 1 : -1 )

       x_3 => x ( 6 * N_src : 5 * N_src + 1 : -1 )

       sgn = -1.0d0

    end if


    if ( type_flag == 'H' ) then
       

       panel_1 => weights ( 1 * N_src + 1 : 2 * N_src, 2 )

       panel_2 => weights ( 0 * N_src + 1 : 1 * N_src, 2 )

       panel_3 => weights ( 0 * N_src + 1 : 1 * N_src, 1 )

       panel_4 => weights ( 1 * N_src + 1 : 2 * N_src, 1 )

              
    else if ( type_flag == 'V' ) then


       panel_1 => weights ( 0 * N_src + 1 : 1 * N_src, 2 )

       panel_2 => weights ( 1 * N_src + 1 : 2 * N_src, 2 )

       panel_3 => weights ( 1 * N_src + 1 : 2 * N_src, 1 )

       panel_4 => weights ( 0 * N_src + 1 : 1 * N_src, 1 )


    end if


    panel_1 = panel_1 + x_1

    panel_2 = panel_2 + x_2

    panel_3 = panel_3 + sgn * x_3

    panel_4 = panel_4 + sgn * x_4


    nullify ( x_1, x_2, x_3, x_4 )

    nullify ( panel_1, panel_2, panel_3, panel_4 )


  end subroutine StoreEquivalentSourceAmplitude



  subroutine computeEquivalentField_FFT ( equiv, ker_fft, field, fft_handle )


    complex(8),dimension(:,:) :: equiv

    complex(8),dimension(:,:) :: ker_fft

    complex(8),dimension(:,:) :: field

    type(FFT_2D) :: fft_handle


    call directFFT ( fft_handle, equiv )

    field = ker_fft * fft_handle % fft_mat
   
    call inverseFFT ( fft_handle, field )

    field = fft_handle % fft_mat


  end subroutine computeEquivalentField_FFT


  
  subroutine PermuteCellField (field_hor, field_ver, N_h, N_v)


    complex(8),dimension(:,:) :: field_hor, field_ver

    integer :: N_h, N_v


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


  end subroutine PermuteCellField



  subroutine ObstacleToCellMap ( cell_mn, refCell_hor, refCell_ver, x_h, x_v )

    ! Computes the horizontal and vertical equivalent sources of the mn-th cell
    ! To represent the field of the inner obstacle
    

    type(Cell) :: cell_mn

    type(ProjectionReferenceCell) :: refCell_hor, refCell_ver

    complex(8),dimension(:) :: x_h, x_v


    call computeFieldInCollocationPoints ( cell_mn, field_hor_coll, field_ver_coll )

    call computeEquivalentSourceAmplitude ( refCell_hor, field_hor_coll, x_h)

    call computeEquivalentSourceAmplitude ( refCell_ver, field_ver_coll, x_v)


  end subroutine ObstacleToCellMap



  subroutine CellToObstacleMap ( cell_mn, field, res )

    
    type (Cell) :: cell_mn

    complex(8),dimension(:) :: field

    complex(8),dimension(:) :: res


    call computePlaneWaveWeights ( interp_cell, field, xi )

    call MatVecMultiply ( cell_mn % cellToObsMatrix, xi, res )


  end subroutine CellToObstacleMap



  subroutine ForwardMap (phy, geo, alg)

    
    type ( phy_Parameters ) :: phy

    type ( geo_Parameters ) :: geo

    type ( alg_Parameters ) :: alg
    

    complex(8),dimension(:,:),pointer :: field_hor, field_ver

    integer :: n, m, Left_ind, Right_ind, N_s, N_h, N_v


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


    do n = 1, geo % N_col

       do m = 1, geo % N_row


          if ( associated ( cellArray ( m, n ) % innerObstacle ) ) then


             call ObstacleToCellMap ( &
                  cellArray ( m, n ), &
                  proj_cell_hor, &
                  proj_cell_ver, &
                  x_hor, &
                  x_ver )        

             ! Store Cell's equivalent horizontal sources            

             call StoreEquivalentSourceAmplitude ( &
                  x_hor, &
                  cellArray ( m, n ) % mon_equiv_hor, &
                  alg % N_src_hor, &
                  'H', &
                  'M' )
             
             call StoreEquivalentSourceAmplitude ( &
                  x_hor, &
                  cellArray ( m, n ) % dip_equiv_hor, &
                  alg % N_src_hor, &
                  'H', &
                  'D' )


             !! Add to global sources

             Left_ind  = ( 2 * alg % N_src_hor ) * n + 1

             Right_ind = ( 2 * alg % N_src_hor ) * ( n + 1 )


             mon_equiv_hor ( Left_ind : Right_ind, m + 1 : m + 2 ) = &
                  mon_equiv_hor ( Left_ind : Right_ind, m + 1 : m + 2 ) + &
                  cellArray ( m, n ) % mon_equiv_hor

             dip_equiv_hor ( Left_ind : Right_ind, m + 1 : m + 2 ) = &
                  dip_equiv_hor ( Left_ind : Right_ind, m + 1 : m + 2 ) + &
                  cellArray ( m, n ) % dip_equiv_hor

             
          
             ! Store Cell's equivalent vertical sources


             Left_ind  = ( 2 * alg % N_src_ver ) * m + 1

             Right_ind = ( 2 * alg % N_src_ver ) * ( m + 1 )


             call StoreEquivalentSourceAmplitude ( &
                  x_ver, &
                  cellArray ( m, n ) % mon_equiv_ver, &
                  alg % N_src_ver, &
                  'V', &
                  'M' )
             
             call StoreEquivalentSourceAmplitude ( &
                  x_ver, &
                  cellArray ( m, n ) % dip_equiv_ver, &
                  alg % N_src_ver, &
                  'V', &
                  'D' )

             !! Add to global sources

             mon_equiv_ver ( Left_ind : Right_ind, n + 1 : n + 2 ) = &
                  mon_equiv_ver ( Left_ind : Right_ind, n + 1 : n + 2 ) + &
                  cellArray ( m, n ) % mon_equiv_ver

             dip_equiv_ver ( Left_ind : Right_ind, n + 1 : n + 2 ) = &
                  dip_equiv_ver ( Left_ind : Right_ind, n + 1 : n + 2 ) + &
                  cellArray ( m, n ) % dip_equiv_ver


          end if


       end do

    end do


    ! --------------------------------------------------------------------------------
    !
    ! Compute field at horizontal equivalent sources using FFT
    !
    ! --------------------------------------------------------------------------------

    ! Monopole contributions

    call computeEquivalentField_FFT ( &
         mon_equiv_hor, &
         mon_hor_ker_global_fft, &
         aux_global_hor, &
         fft_hor_global )

    field_equiv_hor = aux_global_hor


    ! Dipole contributions

    call computeEquivalentField_FFT ( &
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

    call computeEquivalentField_FFT ( &
         mon_equiv_ver, &
         mon_ver_ker_global_fft, &
         aux_global_ver, &
         fft_ver_global )

    field_equiv_ver = aux_global_ver


    ! Dipole contributions

    call computeEquivalentField_FFT ( &
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


    do n = 1, geo % N_col
       
       do m = 1, geo % N_row
             
          ! Horizontal sources

          N_s = 2 * alg % N_src_hor

          Left_ind  = N_s * n + 1
          
          Right_ind = N_s * ( n + 1 )


          cellArray ( m, n ) % field_equiv_hor = &
               field_equiv_hor ( Left_ind : Right_ind, m + 1 : m + 2 )


          ! !! Monopoles 

          local_hor = 0.0d0

          do n_v = 1, 3

             do n_h = 1, 3

                local_hor  ( ( n_h - 1 ) * N_s + 1 : n_h * N_s, n_v : n_v + 1 ) = &
                     local_hor ( ( n_h - 1 ) * N_s + 1 : n_h * N_s, n_v : n_v + 1 ) + &
                     cellArray ( m + n_v - 2, n + n_h - 2 ) % mon_equiv_hor 

             end do

          end do 


          call computeEquivalentField_FFT ( &
               local_hor, &
               mon_hor_ker_local_fft, &
               aux_local_hor, &
               fft_hor_local )

          
          cellArray ( m, n ) % field_equiv_hor = &
               cellArray ( m, n ) % field_equiv_hor - &
               aux_local_hor ( N_s + 1 : 2 * N_s, 2 : 3)
          

          !! Dipoles

          local_hor = 0.0d0

          do n_v = 1, 3

             do n_h = 1, 3

                local_hor ( ( n_h - 1 ) * N_s + 1 : n_h * N_s, n_v : n_v + 1 ) = &
                     local_hor ( ( n_h - 1 ) * N_s + 1 : n_h * N_s, n_v : n_v + 1 ) + &
                     cellArray ( m + n_v - 2, n + n_h - 2 ) % dip_equiv_hor 

             end do

          end do 


          call computeEquivalentField_FFT ( &
               local_hor, &
               dip_hor_ker_local_fft, &
               aux_local_hor, &
               fft_hor_local )

          cellArray ( m, n ) % field_equiv_hor = & 
               cellArray ( m, n ) % field_equiv_hor - & 
               aux_local_hor ( N_s + 1 : 2 * N_s, 2 : 3)


          !! Vertical sources

          N_s = 2 * alg % N_src_ver

          Left_ind  = N_s * m + 1
          
          Right_ind = N_s * ( m + 1 ) 


          cellArray ( m, n ) % field_equiv_ver = &
               field_equiv_ver ( Left_ind : Right_ind, n + 1 : n + 2 )

          !! Monopoles 

          local_ver = 0.0d0

          do n_h = 1, 3

             do n_v = 1, 3

                local_ver ( ( n_v - 1 ) * N_s + 1 : n_v * N_s, n_h : n_h + 1 ) = &
                     local_ver ( ( n_v - 1 ) * N_s + 1 : n_v * N_s, n_h : n_h + 1 ) + &
                     cellArray ( m + n_v - 2, n + n_h - 2 ) % mon_equiv_ver

             end do

          end do


          call computeEquivalentField_FFT ( &
               local_ver, &
               mon_ver_ker_local_fft, &
               aux_local_ver, &
               fft_ver_local )


          cellArray ( m, n ) % field_equiv_ver      = &
               cellArray ( m, n ) % field_equiv_ver - &
               aux_local_ver ( N_s + 1 : 2 * N_s, 2 : 3 )

          
          !! Dipoles

          local_ver = 0.0d0
          
          do n_h = 1, 3

             do n_v = 1, 3

                local_ver ( ( n_v - 1 ) * N_s + 1 : n_v * N_s, n_h : n_h + 1 ) = &
                     local_ver ( ( n_v - 1 ) * N_s + 1 : n_v * N_s, n_h : n_h + 1 ) + &
                     cellArray ( m + n_v - 2, n + n_h - 2 ) % dip_equiv_ver

             end do

          end do


          call computeEquivalentField_FFT ( &
               local_ver, &
               dip_ver_ker_local_fft, &
               aux_local_ver, &
               fft_ver_local )

          cellArray ( m, n ) % field_equiv_ver      = &
               cellArray ( m, n ) % field_equiv_ver - &
               aux_local_ver ( N_s + 1 : 2 * N_s, 2 : 3 )

          
       end do


    end do


    ! --------------------------------------------------------------------------------
    !
    ! Iterate through all cells and obtain outer field at the inner obstacle.
    !
    ! --------------------------------------------------------------------------------


    do n = 1, geo % N_col
          
       do m = 1, geo % N_row

          
          if ( associated ( cellArray ( m, n ) % innerObstacle ) ) then

             
             call PermuteCellField ( &
                  cellArray ( m, n ) % field_equiv_hor, &
                  cellArray ( m, n ) % field_equiv_ver, &
                  alg % N_src_hor, &
                  alg % N_src_ver )
             
             
            call CellToObstacleMap ( &
                 cellArray (m,n), &
                 field_permuted, &
                 res ( :, m, n ) )

             
          end if


       end do


    end do
    

    
  end subroutine ForwardMap


  
  subroutine destroyForwardMap ()

    
    deallocate ( cellArray, obstacleArray ) 

    deallocate ( psi, res )


  end subroutine destroyForwardMap



end module



