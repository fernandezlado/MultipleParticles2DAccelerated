MODULE modCell



  USE modLinearAlgebra

  USE modMathConstants

  USE modObstacle

  USE modProjectionReferenceCell

  USE modInterpolationReferenceCell

  USE modFarInteractions



  IMPLICIT NONE


  
  TYPE Cell


     REAL(8) :: center_x, center_y

     TYPE(Obstacle),POINTER :: innerObstacle

     COMPLEX(8),DIMENSION(:),POINTER :: psi

     COMPLEX(8),DIMENSION(:,:),ALLOCATABLE :: obsToCellMatrix_hor, obsToCellMatrix_ver
     
     COMPLEX(8),DIMENSION(:,:),ALLOCATABLE :: cellToObsMatrix
     
     COMPLEX(8),DIMENSION(:,:),ALLOCATABLE :: mon_equiv_hor, dip_equiv_hor

     COMPLEX(8),DIMENSION(:,:),ALLOCATABLE :: mon_equiv_ver, dip_equiv_ver

     COMPLEX(8),DIMENSION(:,:),ALLOCATABLE :: field_equiv_hor, field_equiv_ver


  END TYPE Cell



CONTAINS


    
  SUBROUTINE createCell (this, c_x, c_y, obs, density, refCell_hor, refCell_ver, interpRefCell, k )

    
    TYPE(Cell) :: this

    TYPE(InterpolationReferenceCell) :: interpRefCell

    REAL(8) :: c_x, c_y
    
    TYPE (Obstacle),POINTER :: obs

    COMPLEX(8),DIMENSION(:),TARGET :: density

    TYPE(ProjectionReferenceCell) :: refCell_hor, refCell_ver

    REAL(8) :: k


    this % center_x = c_x

    this % center_y = c_y
      
    this % psi => density


    ALLOCATE ( this % mon_equiv_hor ( 2 * refCell_hor % N_src, 2 ) )

    ALLOCATE ( this % dip_equiv_hor ( 2 * refCell_hor % N_src, 2 ) )

    ALLOCATE ( this % mon_equiv_ver ( 2 * refCell_ver % N_src, 2 ) )

    ALLOCATE ( this % dip_equiv_ver ( 2 * refCell_ver % N_src, 2 ) )

    ALLOCATE ( this % field_equiv_hor ( 2 * refCell_hor % N_src, 2 ) )

    ALLOCATE ( this % field_equiv_ver ( 2 * refCell_ver % N_src, 2 ) )


    this % mon_equiv_hor = 0.0d0

    this % dip_equiv_hor = 0.0d0

    this % mon_equiv_ver = 0.0d0

    this % dip_equiv_ver = 0.0d0

    this % field_equiv_hor = 0.0d0

    this % field_equiv_ver = 0.0d0


    this % innerObstacle => obs


    IF (ASSOCIATED ( this % innerObstacle ) .eqv. .TRUE. ) THEN


       ALLOCATE ( this % obsToCellMatrix_hor ( 4 * refCell_hor % N_coll,  obs % num_dis ) )

       ALLOCATE ( this % obsToCellMatrix_ver ( 4 * refCell_ver % N_coll,  obs % num_dis ) )

       ALLOCATE ( this % cellToObsMatrix ( obs % num_dis, 4 * interpRefCell % N_wave ) )


       CALL createObstacleToCellMatrix ( this, refCell_hor, k )

       CALL createObstacleToCellMatrix ( this, refCell_ver, k )

       CALL createCellToObstacleMatrix ( this, interpRefCell, k )


    END IF


  END SUBROUTINE createCell



  SUBROUTINE createObstacleToCellMatrix (this, refCell, k)


    TYPE(Cell) :: this
      
    TYPE(ProjectionReferenceCell) :: refCell
      
    REAL(8) :: k


    REAL(8),ALLOCATABLE,DIMENSION(:) :: tar_x, tar_y


    ALLOCATE ( tar_x ( 1 : 4 * refCell % N_coll ) )

    ALLOCATE ( tar_y ( 1 : 4 * refCell % N_coll ) )


    IF ( refCell % cellKind == 'H' ) THEN

       
       tar_x ( 0 * refCell % N_coll + 1 : 1 * refCell % N_coll ) = refCell % coll_x

       tar_y ( 0 * refCell % N_coll + 1 : 1 * refCell % N_coll ) = refCell % coll_y


       tar_x ( 1 * refCell % N_coll + 1 : 2 * refCell % N_coll ) =-refCell % coll_x

       tar_y ( 1 * refCell % N_coll + 1 : 2 * refCell % N_coll ) = refCell % coll_y


       tar_x ( 2 * refCell % N_coll + 1 : 3 * refCell % N_coll ) =-refCell % coll_x

       tar_y ( 2 * refCell % N_coll + 1 : 3 * refCell % N_coll ) =-refCell % coll_y


       tar_x ( 3 * refCell % N_coll + 1 : 4 * refCell % N_coll ) = refCell % coll_x

       tar_y ( 3 * refCell % N_coll + 1 : 4 * refCell % N_coll ) =-refCell % coll_y


       !translate reference Cell to the Cell centered at (m L_x, n L_y)

       tar_x = this % center_x + tar_x

       tar_y = this % center_y + tar_y


       CALL createEvalFieldAtPointsMatrix ( this % obsToCellMatrix_hor, this % innerObstacle, tar_x, tar_y, k)


    ELSE IF ( refCell % cellKind == 'V' ) THEN


       tar_x ( 0 * refCell % N_coll + 1 : 1 * refCell % N_coll ) = refCell % coll_x

       tar_y ( 0 * refCell % N_coll + 1 : 1 * refCell % N_coll ) = refCell % coll_y


       tar_x ( 1 * refCell % N_coll + 1 : 2 * refCell % N_coll ) = refCell % coll_x

       tar_y ( 1 * refCell % N_coll + 1 : 2 * refCell % N_coll ) =-refCell % coll_y


       tar_x ( 2 * refCell % N_coll + 1 : 3 * refCell % N_coll ) =-refCell % coll_x

       tar_y ( 2 * refCell % N_coll + 1 : 3 * refCell % N_coll ) =-refCell % coll_y


       tar_x ( 3 * refCell % N_coll + 1 : 4 * refCell % N_coll ) =-refCell % coll_x

       tar_y ( 3 * refCell % N_coll + 1 : 4 * refCell % N_coll ) = refCell % coll_y


       !translate reference Cell to the Cell centered at (m L_x, n L_y)

       tar_x = this % center_x + tar_x

       tar_y = this % center_y + tar_y


       CALL createEvalFieldAtPointsMatrix ( this % obsToCellMatrix_ver, this % innerObstacle, tar_x, tar_y, k)

    END IF

    
    DEALLOCATE ( tar_x, tar_y ) 


  END SUBROUTINE createObstacleToCellMatrix



  SUBROUTINE createCellToObstacleMatrix (this, refCell, k)


    TYPE(Cell) :: this

    TYPE(InterpolationReferenceCell) :: refCell

    REAL(8) :: k

    
    REAL(8),DIMENSION(:,:),ALLOCATABLE :: U_x, U_y, P_x, P_y

    INTEGER :: N_row, N_col, N_w

    
    N_row = this % innerObstacle % num_dis

    N_w = refCell % N_wave

    N_col = 4 * N_w


    ALLOCATE ( U_x ( 1 : N_row, 1 : N_w ), U_y ( 1 : N_row, 1 : N_w ) )

    ALLOCATE ( P_x ( 1 : N_row, 1 : N_w ), P_y ( 1 : N_row, 1 : N_w ) )


    U_x = SPREAD ( refCell % U_x, 1, N_row )

    U_y = SPREAD ( refCell % U_y, 1, N_row )


    P_x = SPREAD ( this % innerObstacle % C_x, 2, N_w )

    P_y = SPREAD ( this % innerObstacle % C_y, 2, N_w )


    this % cellToObsMatrix ( : , 0 * N_w + 1 : 1 * N_w ) = EXP ( I * k * ( P_x * U_x + P_y * U_y ) )

    this % cellToObsMatrix ( : , 1 * N_w + 1 : 2 * N_w ) = EXP ( I * k * (-P_x * U_x + P_y * U_y ) )

    this % cellToObsMatrix ( : , 2 * N_w + 1 : 3 * N_w ) = EXP ( I * k * (-P_x * U_x - P_y * U_y ) )

    this % cellToObsMatrix ( : , 3 * N_w + 1 : 4 * N_w ) = EXP ( I * k * ( P_x * U_x - P_y * U_y ) )


    DEALLOCATE ( U_x, U_y )

    DEALLOCATE ( P_x, P_y )


  END SUBROUTINE createCellToObstacleMatrix



  SUBROUTINE computeFieldInCollocationPoints ( this, hor, ver )


    TYPE(Cell) :: this

    COMPLEX(8),DIMENSION(:) :: hor, ver


    CALL MatVecMultiply ( this % obsToCellMatrix_hor, this % psi, hor )

    CALL MatVecMultiply ( this % obsToCellMatrix_ver, this % psi, ver )


  END SUBROUTINE computeFieldInCollocationPoints



  SUBROUTINE computeFieldAtInnerObstacle ( this, xi, field )


    TYPE(Cell) :: this
    
    COMPLEX(8),DIMENSION(:) :: xi, field


    CALL MatVecMultiply ( this % cellToObsMatrix, xi, field )


  END SUBROUTINE computeFieldAtInnerObstacle



  SUBROUTINE destroyCell (this)


    TYPE(Cell) :: this


    NULLIFY ( this % innerObstacle )

    DEALLOCATE ( this % obsToCellMatrix_hor, this % obsToCellMatrix_ver )

    DEALLOCATE ( this % cellToObsMatrix )


  END SUBROUTINE destroyCell



END MODULE modCell
