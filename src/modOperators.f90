MODULE modOperators

  USE modLinearAlgebra
  USE modMathConstants
  USE modSpecialFunctions
  USE modObstacle

  
  IMPLICIT none


  TYPE intOperator
     
     TYPE(Obstacle)::obs

     COMPLEX(8),DIMENSION(:,:),ALLOCATABLE::mat_SL
     COMPLEX(8),DIMENSION(:,:),ALLOCATABLE::mat_DL
     
  END TYPE intOperator


CONTAINS  

  
  SUBROUTINE createMatrices(this,op_flags,k_ext,k_int)

    INTERFACE

       REAL(8) ELEMENTAL FUNCTION DBESJ0(t)
         REAL(8),INTENT(IN)::t
       END FUNCTION DBESJ0

       REAL(8) ELEMENTAL FUNCTION DBESJ1(t)
         REAL(8),INTENT(IN)::t
       END FUNCTION DBESJ1

    END INTERFACE
    
    TYPE(intOperator)::this
    LOGICAL,DIMENSION(4)::op_flags
    REAL(8)::k_ext,k_int

    REAL(8),DIMENSION(:,:),ALLOCATABLE::T,Tau
    REAL(8),DIMENSION(:,:),ALLOCATABLE::diff_X,diff_Y,dist
    REAL(8),DIMENSION(:,:),ALLOCATABLE::X_p,Y_p,dS
    REAL(8),DIMENSION(:,:),ALLOCATABLE::MK_Weights

    INTEGER::j,N_i

    COMPLEX(8),DIMENSION(:,:),ALLOCATABLE::K_s,K_r    


    N_i = this % obs % num_dis


    ALLOCATE ( T ( 0:N_i-1 , 0:N_i-1 ), Tau ( 0:N_i-1 , 0:N_i-1 ) )

    ALLOCATE ( diff_X ( 0:N_i-1 , 0:N_i-1 ) , diff_Y ( 0:N_i-1 , 0:N_i-1 ), dist ( 0:N_i-1 , 0:N_i-1 ) )

    ALLOCATE ( X_p ( 0:N_i-1 , 0:N_i-1 ) , Y_p ( 0:N_i-1 , 0:N_i-1 ), dS ( 0:N_i-1 , 0:N_i -1 ) )

    ALLOCATE ( K_s ( 0:N_i-1 , 0:N_i-1 ) , K_r ( 0:N_i-1 , 0:N_i-1 ) ) 

    ALLOCATE ( MK_Weights ( 0:N_i-1 , 0:N_i-1 ) )


    Tau = SPREAD ( this % obs % t, 1 , N_i) 
    T = TRANSPOSE ( Tau)
 
    
    diff_X = SPREAD ( this % obs % C_x , 1 , N_i )
    diff_X = TRANSPOSE ( diff_X ) - diff_X

    diff_Y = SPREAD ( this % obs % C_Y , 1 , N_i )
    diff_Y = TRANSPOSE ( diff_Y ) - diff_Y

    dist = SQRT ( diff_X**2 + diff_Y**2 )
    

    X_p = SPREAD ( this % obs % Cp_x , 1, N_i )
    Y_p = SPREAD ( this % obs % Cp_y , 1, N_i )

    dS = sqrt ( X_p ** 2 + Y_p ** 2)

    MK_Weights = -4.0d0 * Pi / (N_i**2) * COS (N_i/2 *(T-Tau) )
    
    DO j=1,N_i/2-1
       MK_Weights = MK_Weights - 4.0d0 * Pi / N_i *  COS ( j * (T-Tau) ) / j
    END DO


    !Single Layer Matrix
    IF ( op_flags(1) ) THEN


       ALLOCATE (this % mat_SL (0:N_i-1, 0:N_i-1 ) )
       
       !SINGULAR AND REGULAR KERNEL SPLITTING
       K_s = - 1 / (4*pi) * DBESJ0 (k_ext * dist) * dS


       K_r = I/4.0d0 * DBESH0 ( k_ext * dist ) * dS - log ( 4.0d0 * sin ( (T-Tau)/2.0d0 ) **2 ) * K_s

       !DIAGONAL TERMS
       FORALL (j=0:N_i-1) K_r(j,j) = ( I/4.0d0 - Gamma / (2.0d0*pi) - 1 / (2.0d0*pi) * log (k_ext/2 * dS(j,j) ) ) * dS(j,j)

       this % mat_SL = 2 * Pi / N_i * K_r + MK_Weights * K_s 

    END IF


    !Double Layer Matrix
    IF ( op_flags(2) ) THEN
        

       ALLOCATE (this % mat_DL (0:N_i-1, 0:N_i-1) )

       K_s = - k_ext / (4 * Pi) * ( diff_X * Y_p - diff_Y * X_p ) * DBESJ1 ( k_ext * dist) / dist

       K_r = I*k_ext/4.0d0 * ( diff_X * Y_p - diff_Y * X_p ) * DBESH1(k_ext*dist) / dist &
            - log ( 4.0d0 * sin ( (T-Tau)/2.0d0 ) **2 ) * K_s

       !DIAGONAL TERMS       
       FORALL (j=0:N_i-1) 
          K_s(j,j) = 0.0d0
          K_r(j,j) = - 1 / (4*pi) * (this % obs % Cp_x(j) * this % obs % Cpp_y(j) &
            - this % obs % Cp_y(j) * this % obs % Cpp_x(j) ) / dS(j,j)**2
       END FORALL

       this % mat_DL = 2 * Pi / N_i * K_r + MK_Weights * K_s


    END IF


    DEALLOCATE(T,Tau)
    DEALLOCATE(diff_X,diff_Y,dist)
    DEALLOCATE(X_p,Y_p,dS)
    DEALLOCATE(K_s,K_r)
    DEALLOCATE(MK_Weights)


  END SUBROUTINE createMatrices


  SUBROUTINE  createOperator (this, obs , op_flags , k_ext, k_int)

    TYPE(Obstacle)::obs
    TYPE(intOperator)::this
    LOGICAL,DIMENSION(4)::op_flags
    REAL(8)::k_ext,k_int

    this % obs = obs

    CALL createMatrices(this,op_flags,k_ext,k_int)
    
  END SUBROUTINE createOperator


  SUBROUTINE destroyOperator (this)
    
    type(intOperator)::this

    DEALLOCATE(this % mat_SL, this % mat_DL)

  END SUBROUTINE destroyOperator


  SUBROUTINE selfInteraction ( this, res, eta )
    
    type(intOperator)::this
    COMPLEX(8),DIMENSION(:)::res
    COMPLEX(8),DIMENSION(:),ALLOCATABLE::DL_self,SL_self
    REAL(8)::eta

    ALLOCATE ( DL_Self (0:size(res)-1), SL_Self(0:size(res) - 1) )

    
    DL_self = 0.0d0
    CALL MatVecMultiply (this % mat_DL, this % obs % psi, DL_Self)

    SL_self = 0.0d0
    CALL MatVecMultiply (this % mat_SL, this % obs % psi, SL_Self)


    res = DL_self - I * eta * SL_Self


    DEALLOCATE ( DL_Self, SL_Self )
    
  END SUBROUTINE selfInteraction

END MODULE modOperators
