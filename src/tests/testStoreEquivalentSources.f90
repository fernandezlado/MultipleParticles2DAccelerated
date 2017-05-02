PROGRAM MultipleParticle

  USE modCell

  USE modReferenceCell

  USE modFarInteractions

  USE modSpecialFunctions

  USE modObstacle

  USE modOperators

  USE modLinearAlgebra

  USE modForwardMap


  IMPLICIT NONE



  !REFERENCE CELL
  TYPE(ReferenceCell)::ref_cell_hor, ref_cell_ver

  !CELL CONTAINING OBSTACLE
  TYPE(Cell) :: cell_00

  !Photonic array parameters
  INTEGER::N_rows, N_cols, N_obs, n_x, n_y
  
  !Accelerator parameters
  INTEGER,PARAMETER :: N_src = 6, N_coll = 2 * N_src

  !PHYSICAL PARAMETERS
  REAL(8)::k, d_x, d_y, x_0, y_0

  !GEOMETRY PARAMETERS
  REAL(8)::L_x, L_y, R 

  !source and target obstacles 
  TYPE(Obstacle) :: obs_src, obs_tar
 
  !Discretization points
  INTEGER :: N_i

  COMPLEX(8),DIMENSION(:),ALLOCATABLE :: x_hor, x_ver
  
  INTEGER :: j,l


  COMPLEX(8),DIMENSION(:,:),ALLOCATABLE :: mon_weights_hor, dip_weights_hor

  COMPLEX(8),DIMENSION(:,:),ALLOCATABLE :: mon_weights_ver, dip_weights_ver


  REAL(8),DIMENSION(:),ALLOCATABLE :: src_x, src_y

  
  COMPLEX(8),DIMENSION(:),ALLOCATABLE :: field_1, field_2
  
  COMPLEX(8),DIMENSION(:,:),ALLOCATABLE :: mat


  k = 1.0d0

  L_x = 1.5d0

  L_y = 1.5d0

  
  N_rows = 1

  N_cols = 1

  N_i = 32
  
  N_x = 2
 
  N_y = 3

  ALLOCATE ( mon_weights_hor ( 0 : 2 * N_src * N_cols - 1, 0 : N_rows ) )

  ALLOCATE ( dip_weights_hor ( 0 : 2 * N_src * N_rows - 1, 0 : N_cols ) )

  ALLOCATE ( mon_weights_ver ( 0 : 2 * N_src * N_cols - 1, 0 : N_rows ) )

  ALLOCATE ( dip_weights_ver ( 0 : 2 * N_src * N_rows - 1, 0 : N_cols ) )


  ALLOCATE ( src_x ( 0 : 2 * N_src - 1 ) )

  ALLOCATE ( src_y ( 0 : 2 * N_src - 1 ) )

  ALLOCATE ( x_hor ( 0 : 4 * 2 * N_src -1 ) )

  ALLOCATE ( x_ver( 0 : 4 * 2 * N_src -1 ) )


  ALLOCATE ( field_1 ( 0 : N_i - 1 ) )

  ALLOCATE ( field_2 ( 0 : N_i - 1 ) )


  ALLOCATE ( mat ( 0 : N_i - 1, 0 : N_i - 1 ) )

  
  CALL createReferenceCell (ref_cell_hor, N_src, N_coll, L_x, L_y, k, 'H')

  CALL createReferenceCell (ref_cell_ver, N_src, N_coll, L_x, L_y, k, 'V')


  CALL createObstacle ( obs_src, 1, (/ 0.25d0, N_x * L_x, N_y * L_y /), N_i, 0, k )

  CALL createObstacle ( obs_tar, 0, (/ 5.0d0, N_x * L_x, N_y * L_y /), N_i, 0, k )


  CALL createCell ( cell_00, N_y, N_x, obs_src, ref_cell_hor, ref_cell_ver, k )


  mon_weights_hor = 0.0d0

  dip_weights_hor = 0.0d0


  CALL ObstacleToCellMap ( cell_00, ref_cell_hor, ref_cell_ver, x_hor, x_ver )

  CALL StoreEquivalentSourceAmplitude ( x_hor, mon_weights_hor, dip_weights_hor, 0, 0, N_src)

  CALL StoreEquivalentSourceAmplitude ( x_ver, mon_weights_ver, dip_weights_ver, 0, 0, N_src)

  !EVALUATE EXACT FIELD
  CALL createEvalFieldAtPointsMatrix (mat, obs_src, obs_tar % C_x, obs_tar % C_y, k)

  CALL matVecMultiply ( mat, obs_src % psi, field_1)


  !EVALUATE FIELD EQUIV SOURCES; HORIZONTAL CASE

  src_x = N_x * L_x + L_x * (/ ( j / (2.0d0 * N_src - 1) - 0.5d0, j = 0, 2 * N_src - 1 ) /)
 
  DO l = 0, N_i

     src_y =  N_y * L_y + (/ ( -L_y / 2.0d0, j = 0, 2 * N_src - 1 ) /)

     field_2 ( l ) = 0.0d0   

     field_2 ( l ) = field_2 (l) + SUM ( mon_weights_hor(:,0) * G_0 (k,obs_tar % C_x(l) - src_x, obs_tar % C_y(l) - src_y) )

     field_2 ( l ) = field_2 (l) + SUM ( dip_weights_hor(:,0) * DG_0 (k, obs_tar % C_x(l) - src_x, obs_tar % C_y(l) - src_y, 0.0d0, 1.0d0) )

     src_y = N_y * L_y + (/ ( L_y / 2.0d0, j = 0, 2 * N_src - 1 ) /)

     field_2 ( l ) = field_2 (l) + SUM ( mon_weights_hor(:,1) * G_0 (k,obs_tar % C_x(l) - src_x, obs_tar % C_y(l) - src_y) )

     field_2 ( l ) = field_2 (l) + SUM ( dip_weights_hor(:,1) * DG_0 (k, obs_tar % C_x(l) - src_x, obs_tar % C_y(l) - src_y, 0.0d0, 1.0d0) )


  END DO

  write(*,*) abs(field_1-field_2)


  !EVALUATE FIELD EQUIV SOURCES; VERTICAL CASE

  src_y = N_y * L_y + L_y * (/ ( (2 * N_src - 1 - j) / (2.0d0 * N_src - 1) - 0.5d0, j = 0, 2 * N_src - 1 ) /)
 

  DO l = 0, N_i

     src_x = N_x * L_x + (/ ( -L_x / 2.0d0, j = 0, 2 * N_src - 1 ) /)

     field_2 ( l ) = 0.0d0   

     field_2 ( l ) = field_2 (l) + SUM ( mon_weights_ver(:,0) * G_0 (k, obs_tar % C_x(l) - src_x, obs_tar % C_y(l) - src_y) )

     field_2 ( l ) = field_2 (l) + SUM ( dip_weights_ver(:,0) * DG_0 (k, obs_tar % C_x(l) - src_x, obs_tar % C_y(l) - src_y, 1.0d0, 0.0d0) )

     src_x = N_x * L_x + (/ ( L_x / 2.0d0, j = 0, 2 * N_src - 1 ) /)

     field_2 ( l ) = field_2 (l) + SUM ( mon_weights_ver(:,1) * G_0 (k,obs_tar % C_x(l) - src_x, obs_tar % C_y(l) - src_y) )

     field_2 ( l ) = field_2 (l) + SUM ( dip_weights_ver(:,1) * DG_0 (k, obs_tar % C_x(l) - src_x, obs_tar % C_y(l) - src_y, 1.0d0, 0.0d0) )


  END DO

  write(*,*) abs(field_1-field_2)


  CALL destroyReferenceCell ( ref_cell_hor )
 
  CALL destroyReferenceCell ( ref_cell_ver )

  CALL destroyCell ( cell_00 )

  CALL destroyObstacle ( obs_src )


END PROGRAM MultipleParticle
