PROGRAM MultipleParticle


  USE modFarInteractions
  USE modSpecialFunctions
  USE modObstacle
  USE modOperators


  IMPLICIT NONE

  
  !Photonic array parameters
  INTEGER::N_row,N_cols,N_obs

  !dummy variables
  INTEGER::j,l,id

  !PHYSICAL PARAMETERS
  REAL(8)::k,d_x,d_y

  !GEOMETRY PARAMETERS
  REAL(8)::L_x,L_y,R

  !Algorithm PARAMETERS
  REAL(8)::eta
  INTEGER::N_i
  
  type(Obstacle)::eval_lines
  type(Obstacle),DIMENSION(:),ALLOCATABLE::obs_array
  type(intOperator)::int_op
  type(PairwiseFarInteraction)::far

  COMPLEX(8),DIMENSION(:,:),ALLOCATABLE::matrix
  COMPLEX(8),DIMENSION(:),ALLOCATABLE::RHS, sol,field,total

  !SAMPLING PARAMETERS
  INTEGER::N_x,N_y
  REAL(8)::X_min,X_max,Y_min,Y_max,Y

  k = 5.0d0
  d_x = 1.0d0
  d_y = 0.0d0

  eta = k

  L_x = 5.0d0
  L_y = 5.0d0
  R = 1.0d0
  

  N_i = 32

  N_row = 3
  N_cols = 3
  N_obs = N_row*N_cols

  ALLOCATE ( obs_array ( 0:N_obs-1) )
  ALLOCATE ( matrix ( 0:N_obs * N_i-1,0:N_obs*N_i-1 ) )
  ALLOCATE ( RHS ( 0:N_obs * N_i-1), sol ( 0:N_obs * N_i-1) )

  !CREATE GEOMETRY AND RHS
  DO j=0,N_row-1

     DO l=0,N_cols-1
        

        id = j * N_cols + l

        CALL createObstacle (obs_array(id), 0, (/ R, L_x * l, L_y * j/), N_i, id)
        RHS ( id * N_i : (id+1) * N_i - 1) = -planeWave (k, d_x, d_y, obs_array(id) % C_x, obs_array(id) % C_y )

     END DO

  END DO

  !CREATE MATRIX OF SELF AND FAR INTERACTIONS
  DO j=0,N_obs-1
     
     DO l=0,N_obs-1
        
        IF (j==l) THEN
           
           CALL createOperator ( int_op, obs_array(j), (/ .TRUE.,.TRUE.,.FALSE.,.FALSE. /),k,0.0d0)
           matrix ( j * N_i : (j+1) * N_i-1, j * N_i : (j+1) * N_i-1 ) = int_op % mat_DL - I * eta * int_op % mat_SL

           CALL destroyOperator ( int_op )

        ELSE
           
           CALL createPairwiseFarInteraction (far, obs_array(j), obs_array(l), k)
           
           matrix ( j * N_i : (j+1) * N_i-1, l * N_i : (l+1) * N_i-1 ) = far % mat_DL - I * eta * far % mat_SL

           CALL destroyPairwiseFarInteraction (far)
           
        END IF

     END DO

  END DO

  !ADD DIAGONAL TERM
  FORALL (j=0:N_i*N_obs-1)  matrix(j,j) = 0.5d0 + matrix(j,j)


  !SOLVE LINEAR SYSTEM
  CALL LinearSolve ( matrix, RHS, sol )
!  write(*,*) matrix

  !ASSIGN TO EACH OBSTACLE ITS DENSITY
  FORALL (j=0:N_obs-1) obs_array(j) % psi = sol ( j * N_i : (j+1) * N_i -1) 

!  CALL destroyPairwiseFarInteraction(far)
  
  !EVALUATE SCATTERED AND TOTAL FIELD

  OPEN(1,FILE="scattered_real.txt")
  OPEN(2,FILE="scattered_imag.txt")
  OPEN(3,FILE="total_real.txt")
  OPEN(4,FILE="total_imag.txt")
  
  X_min = -3 * L_x
  X_max = (N_cols+3) * L_x

  Y_min = -3 * L_y
  Y_max = (N_row+3) * L_y

  N_x = ( N_cols + 6 ) * N_i
  N_y = ( N_row + 6 ) * N_i

  ALLOCATE ( field(0:N_x-1) )
  ALLOCATE ( total(0:N_x-1) )
 
  DO j=0,N_y-1

     Y = Y_min + (Y_max-Y_min) / N_y * j

     field = 0.0d0

     CALL createObstacle (eval_lines, 2, (/ X_min, Y, X_max, Y /), N_x, 0)
     DO l=0,N_obs-1


        CALL createPairwiseFarInteraction (far, eval_lines, obs_array(l), k)

        field = field + MATMUL ( far%mat_DL, obs_array(l)%psi ) - I * eta * MATMUL(far%mat_SL, obs_array(l)%psi )

        CALL destroyPairwiseFarInteraction(far)

     END DO

     total=field+planeWave(k,d_x,d_y,eval_lines % C_x,Y)
     CALL destroyObstacle (eval_lines)     


     DO l=0,N_x-1
        WRITE(1,*) REAL(field(l))
        WRITE(2,*) AIMAG(field(l))
        WRITE(3,*) REAL(total(l))
        WRITE(4,*) AIMAG(field(l))
     END DO
  END DO

  CLOSE(1)
  CLOSE(2)
  WRITE(*,*) N_y
CONTAINS

  ELEMENTAL REAL(8) FUNCTION planeWave ( k, d_x, d_y, x, y)
    
    REAL(8),INTENT(IN)::k,x,y,d_x,d_y
    
    planeWave = EXP ( I * k * d_x * x + I * k * d_y * y )

  END FUNCTION planeWave

!  ALLOCATE ( matrix(0:N-1,0:N-1) )
!  ALLOCATE ( RHS(0:N-1) )
!  ALLOCATE ( field(0:N-1) )
!
!
!  k = 1.0d0
!  eta = k
!
!
!!CALL createObstacle ( obs_1, 0, (/ 1.0d0 , 0.0d0 , 0.0d0 /), N, 1)
!CALL createObstacle ( obs_2, 0, (/ 1.0d0, 5.0d0, 10.0d0 /) , N, 1 )
!
!CALL createOperator ( op_1, obs_2, (/ .TRUE., .TRUE., .FALSE., .FALSE. /), k, 0.0d0)
!
!CALL createPairwiseFarInteraction (far, obs_1, obs_2, k)
!
!
!
!matrix = op_1 % mat_DL - I * eta * op_1 % mat_SL
!FORALL (j=0:N-1) matrix(j,j) = 0.5 + matrix(j,j)
!RHS = DBESH0 ( k * SQRT ( (obs_2 % C_x-5.0d0) ** 2 + (obs_2 % C_y -10.0d0)** 2) )
!
!CALL LinearSolve(matrix, RHS, obs_2 % Psi)
!
!field = MATMUL(far %mat_DL, obs_2 % Psi) - I * eta * MATMUL(far %mat_SL, obs_2 % Psi)
!
!WRITE (*,*) MAXVAL ( abs (  field -  DBESH0 ( k * SQRT ( (obs_1 % C_x-5.0d0) ** 2 + (obs_1 % C_y-10.0d0) ** 2) ) ) )
!
!CALL destroyOperator(op_1)
!CALL destroyObstacle(obs_1)
!DEALLOCATE( matrix, rhs, field)
!
  
  
END PROGRAM MultipleParticle
