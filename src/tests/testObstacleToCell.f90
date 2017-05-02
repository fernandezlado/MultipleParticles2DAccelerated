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
  INTEGER::N_row, N_cols, N_obs
  
  !Accelerator parameters
  INTEGER,PARAMETER :: N_src = 4, N_coll = 2 * N_src

  !PHYSICAL PARAMETERS
  REAL(8)::k, d_x, d_y, x_0, y_0

  !GEOMETRY PARAMETERS
  REAL(8)::L_x, L_y, R 

  !source and target obstacles 
  TYPE(Obstacle) :: obs_src, obs_tar
 
  !Field
  COMPLEX(8),DIMENSION(:),ALLOCATABLE :: field_1, field_2

  COMPLEX(8),DIMENSION(:),ALLOCATABLE :: x_hor, x_ver

  !Discretization points
  INTEGER :: N_i

  REAL(8) :: x_test, y_test
  
  INTEGER :: j,l

  !test
  COMPLEX(8),DIMENSION(:,:), ALLOCATABLE :: mat


  k = 1.0d0

  L_x = 1.0d0

  L_y = 1.0d0


  N_i = 32


  ALLOCATE ( x_hor ( 0 : 2 * 4 * N_src -1 ), x_ver ( 0 : 2 * 4 * N_src - 1 ) )
  
  ALLOCATE ( mat ( 0 : N_i - 1, 0 : N_i - 1) )

  ALLOCATE ( field_1 ( 0 : N_i - 1 ) )

  ALLOCATE ( field_2 ( 0 : N_i - 1 ) )



  CALL createObstacle ( obs_src, 0, (/ 0.5d0, 0.0d0, 0.0d0 /), N_i, 0, k )

  CALL createObstacle ( obs_tar, 0, (/ 3.0d0, 0.0d0, 0.0d0 /), N_i, 0, k )

  
  CALL createReferenceCell (ref_cell_hor, N_src, N_coll, L_x, L_y, k, 'H')

  CALL createReferenceCell (ref_cell_ver, N_src, N_coll, L_x, L_y, k, 'V')


  CALL createCell ( cell_00, 0, 0, obs_src, ref_cell_hor, ref_cell_ver, k )


  
  CALL ObstacleToCellMap ( cell_00, ref_cell_hor, ref_cell_ver, x_hor, x_ver )


  !Exact field evaluation

  CALL createEvalFieldAtPointsMatrix ( mat, obs_src, obs_tar % C_x, obs_tar % C_y, k )

  CALL MatVecMultiply ( mat, obs_src % psi, field_1 )

  !Equivalent sources field evaluation

  


  DO l = 0, N_i

     field_2 ( l ) = 0.0d0

     x_test = obs_tar % C_x ( l )

     y_test = obs_tar % C_y ( l )

  DO j = 0, N_src - 1

     ! monopole terms
     
     field_2 ( l )  = field_2 ( l )  + x_ver ( j ) * G_0 (k, x_test - ref_cell_ver % src_x ( j ) , y_test - ref_cell_ver % src_y ( j ) )

     field_2 ( l )  = field_2 ( l )  + x_ver ( 2 * N_src + j ) * G_0  (k, x_test - ref_cell_ver % src_x ( j ) , y_test + ref_cell_ver % src_y ( j ) )

     field_2 ( l )  = field_2 ( l )  + x_ver ( 4 * N_src + j ) * G_0  (k, x_test + ref_cell_ver % src_x ( j ) , y_test + ref_cell_ver % src_y ( j ) )

     field_2 ( l )  = field_2 ( l )  + x_ver ( 6 * N_src + j ) * G_0  (k, x_test + ref_cell_ver % src_x ( j ) , y_test - ref_cell_ver % src_y ( j ) )

     ! dipole terms

     field_2 ( l )  = field_2 ( l )  + x_ver ( N_src + j ) * DG_0 (k, x_test - ref_cell_ver % src_x ( j ), y_test - ref_cell_ver % src_y ( j ), 1.0d0, 0.0d0 ) 

     field_2 ( l )  = field_2 ( l )  + x_ver ( 3 * N_src + j ) * DG_0 (k, x_test - ref_cell_ver % src_x ( j ), y_test + ref_cell_ver % src_y ( j ), 1.0d0, 0.0d0 ) 

     field_2 ( l )  = field_2 ( l )  + x_ver ( 5 * N_src + j ) * DG_0 (k, x_test + ref_cell_ver % src_x ( j ), y_test + ref_cell_ver % src_y ( j ), -1.0d0, 0.0d0 ) 

     field_2 ( l )  = field_2 ( l )  + x_ver ( 7 * N_src + j ) * DG_0 (k, x_test + ref_cell_ver % src_x ( j ), y_test - ref_cell_ver % src_y ( j ), -1.0d0, 0.0d0 ) 
     

  END DO

  END DO

  
  write(*,*) abs(field_1-field_2)


  CALL destroyReferenceCell ( ref_cell_hor )
 
  CALL destroyReferenceCell ( ref_cell_ver )

  CALL destroyCell ( cell_00 )

  CALL destroyObstacle ( obs_tar )
  


END PROGRAM MultipleParticle

  ! DO j = 0, N_src - 1

  !    !monopole terms
     
  !    field = field + x ( j ) * G_0 (k, x_test - cell % src_x ( j ) , y_test - cell % src_y ( j ) )

  !    field = field + x ( 2 * N_src + j ) * G_0  (k, x_test + cell % src_x ( j ) , y_test - cell % src_y ( j ) )

  !    field = field + x ( 4 * N_src + j ) * G_0  (k, x_test + cell % src_x ( j ) , y_test + cell % src_y ( j ) )

  !    field = field + x ( 6 * N_src + j ) * G_0  (k, x_test - cell % src_x ( j ) , y_test + cell % src_y ( j ) )

  !    ! dipole terms

  !    field = field + x ( N_src + j ) * DG_0 (k, x_test - cell % src_x ( j ), y_test - cell % src_y ( j ), 0.0d0, 1.0d0 ) 

  !    field = field + x ( 3 * N_src + j ) * DG_0 (k, x_test + cell % src_x ( j ), y_test - cell % src_y ( j ), 0.0d0, 1.0d0 ) 

  !    field = field + x ( 5 * N_src + j ) * DG_0 (k, x_test + cell % src_x ( j ), y_test + cell % src_y ( j ), 0.0d0, -1.0d0 ) 

  !    field = field + x ( 7 * N_src + j ) * DG_0 (k, x_test - cell % src_x ( j ), y_test + cell % src_y ( j ), 0.0d0, -1.0d0 ) 
     

  ! END DO


  ! write(*,*) abs ( field -G_0(k,x_test-x_0,y_test-y_0)) 
