PROGRAM MultipleParticle

  USE modLinearAlgebra

  USE modObstacle
  
  USE modProjectionReferenceCell
  
  USE modCell

  USE modForwardMap
  
  USE modIO
 
  USE modParameters

  USE modInterpolationReferenceCell


  IMPLICIT NONE


  TYPE (Obstacle) :: obs

  TYPE (ProjectionReferenceCell) :: ref_cell_hor, ref_cell_ver

  TYPE (InterpolationReferenceCell) :: interp

  TYPE (Cell) :: cell_obj

  TYPE (phy_Parameters) :: phy

  TYPE (geo_Parameters) :: geo

  TYPE (alg_Parameters) :: alg


  COMPLEX(8), DIMENSION(:), ALLOCATABLE :: density
  
  COMPLEX(8), DIMENSION(:), ALLOCATABLE :: x_h, x_v

  COMPLEX(8), DIMENSION(:,:), ALLOCATABLE :: mon, dip

  COMPLEX(8), DIMENSION(:,:), ALLOCATABLE :: matrix


  REAL(8),DIMENSION(:),ALLOCATABLE :: src_x, src_y

  REAL(8) :: tar_x, tar_y

  COMPLEX(8) :: field_1
  
  complex(8),DIMENSION(1) :: field_2

  INTEGER :: n,m


  CALL loadData ( phy, geo, alg )

  ALLOCATE ( density ( alg % N_dis ) )

  ALLOCATE ( x_h ( 8 * alg % N_src_hor ) )

  ALLOCATE ( x_v ( 8 * alg % N_src_ver ) )

  ALLOCATE ( mon ( 2 * alg % N_src_hor, 2 ) )

  ALLOCATE ( dip ( 2 * alg % N_src_hor, 2 ) )

  ALLOCATE ( src_x ( 2 * alg % N_src_hor ) )

  ALLOCATE ( src_y ( 2 ) )

  ALLOCATE ( matrix ( 1, alg % N_dis ) )


  density = 1.0d0


  mon = 0.0d0

  dip = 0.0d0
  
  CALL createProjectionReferenceCell ( ref_cell_hor, alg % N_src_hor, alg % N_coll_hor, geo % L_x, geo % L_y, phy % k, 'H' )

  CALL createProjectionReferenceCell ( ref_cell_ver, alg % N_src_ver, alg % N_coll_ver, geo % L_x, geo % L_y, phy % k, 'V' )


  CALL createObstacle (obs, 0, (/ geo % radius, 0.0d0, 0.0d0 /), alg % N_dis )

  CALL createCell ( cell_obj, 0.0d0, 0.0d0, obs, density, ref_cell_hor, ref_cell_ver, interp , phy % k)

  CALL ObstacleToCellMap ( cell_obj, ref_cell_hor, ref_cell_ver, x_h, x_v )


  field_1 = 0.0d0

  tar_x = 10 * geo % L_x

  tar_y = 10 * geo % L_y

  ! DO n = 1, alg % N_src_hor

  !   field_1 = field_1 + x_h ( 0 * alg % N_src_hor + n ) * G_0 ( phy % k, tar_x - ref_cell_hor % src_x ( n ), tar_y - ref_cell_hor % src_y ( n ) )

  !   field_1 = field_1 + x_h ( 1 * alg % N_src_hor + n ) * DG_0 ( phy % k, tar_x - ref_cell_hor % src_x ( n ), tar_y - ref_cell_hor % src_y ( n ), 0.0d0, 1.0d0 )

  !   field_1 = field_1 + x_h ( 2 * alg % N_src_hor + n ) * G_0 ( phy % k, tar_x + ref_cell_hor % src_x ( n ), tar_y - ref_cell_hor % src_y ( n ) )

  !   field_1 = field_1 + x_h ( 3 * alg % N_src_hor + n ) * DG_0 ( phy % k, tar_x + ref_cell_hor % src_x ( n ), tar_y - ref_cell_hor % src_y ( n ), 0.0d0, 1.0d0 )

  !   field_1 = field_1 + x_h ( 4 * alg % N_src_hor + n ) * G_0 ( phy % k, tar_x + ref_cell_hor % src_x ( n ), tar_y + ref_cell_hor % src_y ( n ))

  !   field_1 = field_1 + x_h ( 5 * alg % N_src_hor + n ) * DG_0 ( phy % k, tar_x + ref_cell_hor % src_x ( n ), tar_y + ref_cell_hor % src_y ( n ), 0.0d0, -1.0d0 )

  !   field_1 = field_1 + x_h ( 6 * alg % N_src_hor + n ) * G_0 ( phy % k, tar_x - ref_cell_hor % src_x ( n ), tar_y + ref_cell_hor % src_y ( n ) )

  !   field_1 = field_1 + x_h ( 7 * alg % N_src_hor + n ) * DG_0 ( phy % k, tar_x - ref_cell_hor % src_x ( n ), tar_y + ref_cell_hor % src_y ( n ), 0.0d0, -1.0d0 )

  ! END DO


  CALL StoreEquivalentSourceAmplitude (x_h, mon, alg % N_src_hor, 'H', 'M')

  CALL StoreEquivalentSourceAmplitude (x_h, dip, alg % N_src_hor, 'H', 'D')


  src_x = geo % L_x * (/ ( ( 2.0d0 * m + 1.0d0 ) / ( 4 * alg % N_src_hor ) - 0.5d0, m = 0 , 2 * alg % N_src_hor - 1 ) /)

  src_y = (/ -geo % L_y/2, geo % L_y/2 /)

  tar_x = 4 * geo % L_x

  tar_y = 4 * geo % L_y

  field_1 = 0.0d0
  

  DO m = 1, 2

     DO n = 1, 2 * alg % N_src_hor

        field_1 = field_1 + G_0 ( phy % k, tar_x - src_x (n), tar_y - src_y (m) ) * mon ( n, m )

        field_1 = field_1 + DG_0 ( phy % k, tar_x - src_x (n), tar_y - src_y (m), 0.0d0, 1.0d0 ) * dip(n,m)
        
     END DO

  END DO

  write(*,*) x_h ( 5 * alg % N_src_hor + 1 : 6 * alg % N_src_hor )

  write(*,*) 

  write(*,*) dip ( 1 : alg % N_src_hor , 1 )

  CALL createEvalFieldAtPointsMatrix (matrix, obs, (/tar_x/), (/tar_y/), phy % k)

  Call MatVecMultiply ( matrix, density, field_2 )

  write(*,*) abs(field_1-field_2)/abs(field_2)
  
END PROGRAM MultipleParticle








