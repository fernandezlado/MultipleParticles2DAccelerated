PROGRAM MultipleParticle
  



  USE modProjectionReferenceCell

  USE modCell
  
  USE modInterpolationReferenceCell

  USE modFarInteractions

  USE modSpecialFunctions

  USE modObstacle

  USE modLinearAlgebra


  IMPLICIT NONE


  !REFERENCE CELL
  TYPE(InterpolationReferenceCell) :: ref_interp_cell

  TYPE(ProjectionReferenceCell) :: refCell_hor, refCell_ver

  !InterpolationReferenceCell paramters
  INTEGER,PARAMETER :: N_src_hor = 8, N_src_ver = 8, N_src = N_src_hor + N_src_ver

  INTEGER,PARAMETER :: N_wave = N_src

  REAL(8) :: L_x, L_y

  !PHYSICAL PARAMETERS
  REAL(8) :: k, x_test, y_test

  !source and target obstacles 
  TYPE(Obstacle),TARGET :: obs_src, obs_tar

  TYPE(Obstacle),POINTER :: obs_pointer
 
  !Discretization points
  INTEGER :: N_i

  REAL(8),DIMENSION(:),ALLOCATABLE :: tar_x, tar_y

  COMPLEX(8),DIMENSION(:),ALLOCATABLE :: rhs, x, field_1, field_2, density

  COMPLEX(8),DIMENSION(:,:),ALLOCATABLE :: mat
  
  !CELL containing target obstacle
  TYPE(Cell) :: cell_0
  

  k = 10.0d0

  L_x = 1.0d0

  L_y = 1.0d0

  
  N_i = 32
  
  
  ALLOCATE ( x ( 4 * N_wave) )

  ALLOCATE ( rhs ( 4 * N_src ) )


  ALLOCATE ( tar_x ( 4 * N_src ) )

  ALLOCATE ( tar_y ( 4 * N_src ) )



  ALLOCATE ( field_1 ( N_i ), field_2 ( N_i), density (N_i) )

  ALLOCATE ( mat ( 0 : 4 * N_src - 1 , 0 : N_i - 1 ) )
  

  CALL createObstacle ( obs_tar , 0, (/ 0.5d0, 0.0d0, 0.0d0 /), N_i )

  CALL createObstacle ( obs_src , 0, (/ 0.5d0, 2.0d0, 0.0d0 /), N_i )

  obs_pointer => obs_tar

  density = 1.0d0


  CALL createInterpolationReferenceCell ( &
       ref_interp_cell, &
       N_src_hor, &
       N_src_ver, &
       N_wave, &
       L_x, &
       L_y, &
       k)


  CALL createCell ( &
       cell_0, &
       0.0d0, &
       0.0d0, &
       obs_pointer, &
       density, &
       refCell_hor, &
       refCell_ver, &
       ref_interp_cell, &
       k )

  tar_x ( 0 * N_src + 1 : 1 * N_src ) = ref_interp_cell % src_x
 
  tar_y ( 0 * N_src + 1 : 1 * N_src ) = ref_interp_cell % src_y


  tar_x ( 1 * N_src + 1 : 2 * N_src ) =-ref_interp_cell % src_x
 
  tar_y ( 1 * N_src + 1 : 2 * N_src ) = ref_interp_cell % src_y


  tar_x ( 2 * N_src + 1 : 3 * N_src ) =-ref_interp_cell % src_x
 
  tar_y ( 2 * N_src + 1 : 3 * N_src ) =-ref_interp_cell % src_y


  tar_x ( 3 * N_src + 1 : 4 * N_src ) = ref_interp_cell % src_x
 
  tar_y ( 3 * N_src + 1 : 4 * N_src ) =-ref_interp_cell % src_y


  CALL createEvalFieldAtPointsMatrix ( mat, obs_src, tar_x, tar_y, k)

  CALL MatVecMultiply ( mat, density, rhs )
  
  
  CALL computePlaneWaveWeights ( ref_interp_cell, rhs, x)

!  write(*,*) x

  write (*,*) cell_0 % cellToObsMatrix(1,:)

  CALL MatVecMultiply ( cell_0 % cellToObsMatrix, x, field_1 )

  
  DEALLOCATE ( mat )

  ALLOCATE ( mat ( N_i, N_i ) )

  CALL createEvalFieldAtPointsMatrix ( mat, obs_src, obs_tar % C_x, obs_tar % C_y, k )

  CALL MatVecMultiply ( mat, density, field_2 )

  write(*,*) "-----------------------------------------------------------------------"

  write(*,*) "Absoulte error between field projected with plane waves and real field"

  write(*,*) "-----------------------------------------------------------------------"

  write(*,*) abs(field_1-field_2)

  write(*,*) "-----------------------------------------------------------------------"
  
  CALL destroyInterpolationReferenceCell ( ref_interp_cell )

  CALL destroyObstacle ( obs_src )


END PROGRAM MultipleParticle
