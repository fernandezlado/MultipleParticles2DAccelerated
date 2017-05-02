PROGRAM MultipleParticle


  USE modReferenceCell

  USE modFarInteractions

  USE modSpecialFunctions

  USE modObstacle

  USE modOperators

  USE modLinearAlgebra



  IMPLICIT NONE



  !REFERENCE CELL
  TYPE(ReferenceCell)::cell

  !Photonic array parameters
  INTEGER::N_row,N_cols,N_obs
  
  !Accelerator parameters
  INTEGER,PARAMETER :: N_src = 6, N_coll = 2 * N_src

  !dummy variables
  INTEGER::j,l,id

  !PHYSICAL PARAMETERS
  REAL(8)::k,d_x,d_y, x_0,y_0

  !GEOMETRY PARAMETERS
  REAL(8)::L_x,L_y,R 
 
  COMPLEX(8),DIMENSION(0:4 * N_coll-1) :: rhs
 
  COMPLEX(8),DIMENSION(0:2 * 4 * N_src -1):: x

  !TEST
  REAL(8) :: x_test, y_test

  COMPLEX(8) :: field


  k = 1.0d0

  L_x = 1.0d0

  L_y = 1.0d0

  x_0 = 0.0d0

  y_0 = 0.0d0


  CALL createReferenceCell (cell, N_src, N_coll, L_x, L_y, k, 'H')

  rhs ( 0 : N_coll - 1 ) = G_0 ( k, cell % coll_x - x_0, cell % coll_y - y_0)

  rhs ( N_coll : 2 * N_coll - 1 ) = G_0 ( k, -cell % coll_x - x_0, cell % coll_y - y_0)

  rhs ( 2 * N_coll : 3 * N_coll - 1 ) = G_0 ( k, -cell % coll_x - x_0, -cell % coll_y - y_0)

  rhs ( 3 * N_coll : 4 * N_coll - 1 ) = G_0 ( k, cell % coll_x - x_0, -cell % coll_y - y_0)


  CALL computeEquivalentSourceAmplitude (cell, rhs, x)

  x_test = 3.0d0 * L_x /2.0d0

  y_test = 3.0d0 * L_y / 2.0d0

  field = 0

  DO j = 0, N_src - 1

     !monopole terms
     
     field = field + x ( j ) * G_0 (k, x_test - cell % src_x ( j ) , y_test - cell % src_y ( j ) )

     field = field + x ( 2 * N_src + j ) * G_0  (k, x_test + cell % src_x ( j ) , y_test - cell % src_y ( j ) )

     field = field + x ( 4 * N_src + j ) * G_0  (k, x_test + cell % src_x ( j ) , y_test + cell % src_y ( j ) )

     field = field + x ( 6 * N_src + j ) * G_0  (k, x_test - cell % src_x ( j ) , y_test + cell % src_y ( j ) )

     ! dipole terms

     field = field + x ( N_src + j ) * DG_0 (k, x_test - cell % src_x ( j ), y_test - cell % src_y ( j ), 0.0d0, 1.0d0 ) 

     field = field + x ( 3 * N_src + j ) * DG_0 (k, x_test + cell % src_x ( j ), y_test - cell % src_y ( j ), 0.0d0, 1.0d0 ) 

     field = field + x ( 5 * N_src + j ) * DG_0 (k, x_test + cell % src_x ( j ), y_test + cell % src_y ( j ), 0.0d0, -1.0d0 ) 

     field = field + x ( 7 * N_src + j ) * DG_0 (k, x_test - cell % src_x ( j ), y_test + cell % src_y ( j ), 0.0d0, -1.0d0 ) 
     

  END DO


  write(*,*) abs ( field -G_0(k,x_test-x_0,y_test-y_0))

 
  CALL destroyReferenceCell (cell)
  

END PROGRAM MultipleParticle

