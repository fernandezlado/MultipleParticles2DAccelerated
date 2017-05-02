module modObstacle


  USE modMathConstants


  IMPLICIT none


  TYPE Obstacle
     

     INTEGER :: Num_dis


     REAL(8),DIMENSION(:),ALLOCATABLE :: t

     REAL(8),DIMENSION(:),ALLOCATABLE :: C_x

     REAL(8),DIMENSION(:),ALLOCATABLE :: C_y

     REAL(8),DIMENSION(:),ALLOCATABLE :: Cp_x

     REAL(8),DIMENSION(:),ALLOCATABLE :: Cp_y

     REAL(8),DIMENSION(:),ALLOCATABLE :: Cpp_x

     REAL(8),DIMENSION(:),ALLOCATABLE :: Cpp_y

          
  END TYPE Obstacle


CONTAINS


  SUBROUTINE createObstacle(this, shape, parameters, N)


    type(Obstacle) :: this

    INTEGER :: shape, id, N

    REAL(8),DIMENSION(:) :: parameters


    INTEGER :: j
    

    ALLOCATE ( this % t (0:N-1) )

    ALLOCATE ( this % C_x(0:N-1), this % C_y(0:N-1) )

    ALLOCATE ( this % Cp_x(0:N-1), this % Cp_y(0:N-1) )

    ALLOCATE ( this % Cpp_x(0:N-1), this % Cpp_y(0:N-1) )


    this % num_dis = N
    
    this % t =  (/ (2*pi / N * j,j=0,N-1) /) 


    SELECT CASE (shape)
       !DICTIONARY OF SHAPES
       !0: CIRCLE
       !1: KITE
       !2: ELLIPSE

       CASE (0)
          !PARAMETERS IS A THREE DIMENSIONAL REAL ARRAY
          !1ST COORDINATE: RADIUS
          !2ND COORDINATE: X-COORD OF CENTER
          !3RD COORDINATE: Y-COORD OF CENTER

          this % C_x = parameters(2) + parameters(1) * COS( this % t )

          this % C_y = parameters(3) + parameters(1) * SIN( this % t )


          this % Cp_x = - parameters(1) * SIN( this % t )

          this % Cp_y = parameters(1) * COS( this % t )


          this % Cpp_x = - parameters(1) * COS( this % t )

          this % Cpp_y = - parameters(1) * SIN( this % t )

          
       CASE(1)
          !PARAMETERS IS A THREE DIMENSIONAL REAL ARRAY
          !1ST COORDINATE: SCALING FACTOR, IF = 1 THEN KITE FROM COLTON-KRESS
          !2ND COORDINATE: X-COORD OF CENTER
          !3RD COORDINATE: Y-COORD OF CENTER

          this % C_x = parameters(2) + parameters(1) * ( COS( this % t ) + 0.65 * ( COS( 2 * this % t ) - 1 ) )

          this % C_y = parameters(3) + parameters(1) * 1.5 * SIN( this % t )


          this % Cp_x = - parameters(1) * ( SIN( this % t ) + 2 * 0.65 * SIN ( 2 * this % t) )

          this % Cp_y = parameters(1) * 1.5 * COS( this % t )


          this % Cpp_x = - parameters(1) * ( COS( this % t ) + 4 * 0.65 * COS ( 2 * this % t ) )


          this % Cpp_y = - parameters(1) * 1.5 * SIN( this % t )

       CASE(2)
          !PARAMETERS IS A FOUR DIMENSIONAL REAL ARRAY
          !1ST COORDINATE: X-COORD OF INITIAL POINT
          !2ND COORDINATE: Y-COORD OF INITIAL POINT
          !3RD COORDINATE: X-COORD OF END POINT
          !4TH COORDINATE: Y-COORD OF END POINT

          this % C_x = parameters(1) * ( 1 - this % t / (2*pi) ) + parameters(3) * this % t / (2*pi)
          this % C_y = parameters(2) * ( 1 - this % t / (2*pi) ) + parameters(4) *  this % t / (2*pi) 


          this % Cp_x = ( parameters(3) - parameters(1) ) / (2*pi)

          this % Cp_y = ( parameters(4) - parameters(2) ) / (2*pi)


          this % Cpp_x = 0.0d0

          this % Cpp_y = 0.0d0


    END SELECT


  END SUBROUTINE createObstacle
  

  SUBROUTINE destroyObstacle (this)
    
    TYPE(Obstacle) :: this

    
    DEALLOCATE ( this % t )

    DEALLOCATE ( this % C_x, this % C_y )

    DEALLOCATE ( this % Cp_x, this % Cp_y )

    DEALLOCATE ( this % Cpp_x, this % Cpp_y )


  END SUBROUTINE DESTROYOBSTACLE


end module modObstacle
