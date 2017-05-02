MODULE modFarInteractions

  USE modLinearAlgebra

  USE modSpecialFunctions
  USE modObstacle
  USE modMathConstants
  

  IMPLICIT NONE

  
!   TYPE PairwiseFarInteraction

!      TYPE(Obstacle)::target_obs,source_obs
!      COMPLEX(8),DIMENSION(:,:),ALLOCATABLE::mat_SL
!      COMPLEX(8),DIMENSION(:,:),ALLOCATABLE::mat_DL     
     
!   END TYPE PairwiseFarInteraction

  
!   TYPE FarInteractions

!      TYPE(Obstacle)::obs
!      TYPE(PairwiseFarInteraction),DIMENSION(:),ALLOCATABLE::neighbor_interactions
!      INTEGER::N_neigh
     
!   END TYPE FarInteractions

  
CONTAINS

  
!     SUBROUTINE createInteractionMatrices ( target_obs, source_obs, mat_SL, mat_DL , k)

!       TYPE(Obstacle)::target_obs,source_obs
!       COMPLEX(8),DIMENSION(:,:)::mat_SL,mat_DL
!       REAL(8)::k
    
!       REAL(8),DIMENSION(:,:),ALLOCATABLE::diff_X,diff_Y,dist
!       REAL(8),DIMENSION(:,:),ALLOCATABLE::X_p,Y_p,dS
    

!       ALLOCATE ( diff_X ( 0:target_obs % num_dis -1, 0:source_obs % num_dis -1 ) )
!       ALLOCATE ( diff_Y ( 0:target_obs % num_dis -1, 0:source_obs % num_dis -1 ) )
!       ALLOCATE ( dist ( 0:target_obs % num_dis -1, 0:source_obs % num_dis -1 ) )
      
!       ALLOCATE ( X_p ( 0:target_obs % num_dis -1, 0:source_obs % num_dis -1 ) )
!       ALLOCATE ( Y_p ( 0:target_obs % num_dis -1, 0:source_obs % num_dis -1 ) )
!       ALLOCATE ( dS ( 0:target_obs % num_dis -1, 0:source_obs % num_dis -1 ) )


!       diff_X = TRANSPOSE ( SPREAD ( target_obs % C_X , 1, source_obs % num_dis) ) &
!            -SPREAD ( source_obs % C_X , 1 , target_obs % num_dis)
    
!       diff_Y = TRANSPOSE ( SPREAD ( target_obs % C_Y , 1, source_obs % num_dis) ) &
!            -SPREAD ( source_obs % C_Y , 1 , target_obs % num_dis)
      
!       dist = SQRT ( diff_X**2 + diff_Y**2 )


!       X_p = SPREAD ( source_obs % Cp_x, 1, target_obs % num_dis)
      
!       Y_p = SPREAD ( source_obs % Cp_y, 1, target_obs % num_dis)

!       dS = SQRT ( X_p**2 + Y_p**2 )
      
!       mat_SL = ( 2*pi / source_obs % num_dis) *  I / 4.0d0 * DBESH0 ( k * dist ) * dS

!       mat_DL = ( 2*pi / source_obs % num_dis) * ( I * k) / 4.0d0 * DBESH1 ( k * dist) &
!             * ( diff_X * Y_p - diff_Y * X_p ) / dist


!       DEALLOCATE ( diff_X, diff_Y, dist )
!       DEALLOCATE ( X_p, Y_p, dS ) 


!    END SUBROUTINE createInteractionMatrices



!    SUBROUTINE createPairwiseFarInteraction(this,target_obs,source_obs,k)

!      TYPE(PairwiseFarInteraction)::this
!      TYPE(Obstacle)::target_obs,source_obs
!      REAL(8)::k


!      this % target_obs = target_obs
!      this % source_obs = source_obs


!      ALLOCATE ( this % mat_SL ( 0 : target_obs % num_dis - 1 , 0 : source_obs % num_dis - 1 ) )
!      ALLOCATE ( this % mat_DL ( 0 : target_obs % num_dis - 1 , 0 : source_obs % num_dis - 1 ) )


!      CALL createInteractionMatrices ( this % target_obs, this % source_obs, this % mat_SL, this % mat_DL ,k)


!    END SUBROUTINE createPairwiseFarInteraction


!    SUBROUTINE destroyPairwiseFarInteraction (this)

!      TYPE(PairwiseFarInteraction)::this

!      DEALLOCATE ( this % mat_DL, this % mat_SL)

!    END SUBROUTINE destroyPairwiseFarInteraction


!    SUBROUTINE createFarInteractions (this,self_obs,neighbor_obs,k)

!      TYPE(FarInteractions)::this
!      TYPE(Obstacle)::self_obs
!      TYPE(Obstacle),DIMENSION(:)::neighbor_obs
!      REAL(8)::k

!      INTEGER::j

!      ALLOCATE ( this % neighbor_interactions ( 0 : size (neighbor_obs) - 1 ) )

!      this % obs = self_obs

!      DO j = 0 , SIZE ( neighbor_obs )-1

!        CALL createPairwiseFarInteraction( this % neighbor_interactions ( j ), self_obs, neighbor_obs (j+1) , k )

!     END DO


!   END SUBROUTINE createFarInteractions

  
  SUBROUTINE createEvalFieldAtPointsMatrix (mat, obs, tar_x, tar_y, k)


    COMPLEX(8),DIMENSION(:,:) ::  mat

    TYPE(Obstacle) :: obs
    
    REAL(8),DIMENSION(:) :: tar_X, tar_Y

    REAL(8) ::  k


    REAL(8),DIMENSION(:,:),ALLOCATABLE :: diff_X, diff_Y, dist
    
    REAL(8),DIMENSION(:,:),ALLOCATABLE :: X_p, Y_p, dS

    COMPLEX(8),DIMENSION(:,:),ALLOCATABLE ::  mat_DL, mat_SL

    INTEGER :: N_tar

    
    N_tar = size ( tar_x )


    ALLOCATE ( diff_X ( 0 : N_tar - 1, 0 : obs % num_dis - 1 ) )

    ALLOCATE ( diff_Y ( 0 : N_tar - 1, 0 : obs % num_dis - 1 ) )

    ALLOCATE ( dist ( 0 : N_tar - 1, 0 : obs % num_dis - 1 ) )

    
    ALLOCATE ( X_p ( 0 : N_tar - 1, 0 : obs % num_dis - 1 ) )

    ALLOCATE ( Y_p ( 0 : N_tar - 1, 0 : obs % num_dis - 1 ) )

    ALLOCATE ( dS ( 0 : N_tar - 1, 0 : obs % num_dis - 1 ) )
    

    ALLOCATE ( mat_DL ( 0 : N_tar - 1, 0 : obs % num_dis - 1 ) )

    ALLOCATE ( mat_SL ( 0 : N_tar - 1, 0 : obs % num_dis - 1 ) )

      
    diff_X = SPREAD ( tar_x , 2 , obs % num_dis ) &
         - SPREAD ( obs % C_X , 1 , N_tar )
    
    diff_Y = SPREAD ( tar_y , 2 , obs % num_dis)  &
         -SPREAD ( obs % C_Y , 1 , N_tar )
    
    dist = SQRT ( diff_X**2 + diff_Y**2 )
    

    X_p = SPREAD ( obs % Cp_x, 1, N_tar )
      
    Y_p = SPREAD ( obs % Cp_y, 1, N_tar )
      
    dS = SQRT ( X_p**2 + Y_p**2 )
     
 
    mat_SL = ( 2*pi / obs % num_dis ) *  I / 4.0d0 * DBESH0 ( k * dist ) * dS

    mat_DL = ( 2*pi / obs % num_dis ) * ( I * k ) / 4.0d0 * DBESH1 ( k * dist ) &
         * ( diff_X * Y_p - diff_Y * X_p ) / dist
    

    mat = mat_DL - I * k * mat_SL


    DEALLOCATE ( mat_DL, mat_SL )

    DEALLOCATE ( diff_X, diff_Y, dist )

    DEALLOCATE ( X_p, Y_p, dS ) 

      
  END SUBROUTINE createEvalFieldAtPointsMatrix


  ! SUBROUTINE nonSelfInteractions (this, res, eta)
    

  !   TYPE(FarInteractions)::this
  !   COMPLEX(8),DIMENSION(:)::res
  !   REAL(8)::eta


  !   INTEGER::N_nei,j
  !   COMPLEX(8),DIMENSION(:),ALLOCATABLE::DL_field,SL_field

  !   ALLOCATE ( DL_FIELD (0: this % obs % num_dis ), SL_FIELD (0: this % obs % num_dis ) )

  !   N_nei = size ( this % neighbor_interactions )

  !   DO j=0,N_nei-1

  !      DL_field = 0.0d0
  !      CALL MatVecMultiply ( this % neighbor_interactions(j) % mat_DL, this % neighbor_interactions(j) % source_obs % psi , DL_field)

  !      SL_field = 0.0d0 
  !      CALL MatVecMultiply ( this % neighbor_interactions(j) % mat_SL, this % neighbor_interactions(j) % source_obs % psi , SL_field)

  !      res = res + DL_field - I * eta * SL_field

  !   END DO
    


  ! END SUBROUTINE nonSelfInteractions


END MODULE modFarInteractions
