MODULE modSpecialFunctions

  USE IFPORT
  USE modMathConstants

  IMPLICIT NONE


CONTAINS


  COMPLEX(8) ELEMENTAL FUNCTION  DBESH0 ( t )
    !COMPUTES 0-th ordern HANKEL FUNCTION OF THE FIRST KIND

    INTERFACE

       REAL(8) ELEMENTAL FUNCTION DBESJ0(t)
         REAL(8),INTENT(IN)::t
       END FUNCTION DBESJ0

       REAL(8) ELEMENTAL FUNCTION DBESY0(t)
         REAL(8),INTENT(IN)::t
       END FUNCTION DBESY0

    END INTERFACE
    

    REAL(8),INTENT(IN)::t


    DBESH0= DCMPLX (DBESJ0(t), DBESY0(t) )

  END FUNCTION DBESH0


  COMPLEX(8) ELEMENTAL FUNCTION  DBESH1 ( t )
    !COMPUTES 0-th ordern HANKEL FUNCTION OF THE FIRST KIND
    
    INTERFACE

       REAL(8) ELEMENTAL FUNCTION DBESJ1(t)
         REAL(8),INTENT(IN)::t
       END FUNCTION DBESJ1

       REAL(8) ELEMENTAL FUNCTION DBESY1(t)
         REAL(8),INTENT(IN)::t
       END FUNCTION DBESY1

    END INTERFACE
    

    REAL(8),INTENT(IN)::t


    DBESH1= DCMPLX( DBESJ1(t), DBESY1(t) )

  END FUNCTION DBESH1

  
  COMPLEX(8) ELEMENTAL FUNCTION G_0 (k,X,Y)
    !COMPUTES HELMHOLTZ GREEN FUNCTION AT X,Y

    REAL(8),INTENT(IN)::k,X,Y
    
    
    G_0 = I / 4.0d0 * DBESH0 ( k * SQRT ( X**2 + Y**2 ) )

  END FUNCTION G_0
    

  COMPLEX(8) ELEMENTAL FUNCTION DG_0 (k,X,Y,D_x,D_y)
    !COMPUTES (D_x,D_y) DIRECTIONAL DERIVATIVE OF HELMHOLTZ GREEN FUNCTION AT X,Y
    !ASSUMES D_X **2 + D_y ** 2=1

    REAL(8),INTENT(IN)::k,X,Y,D_x,D_y
    
    
    DG_0 = I / 4.0d0 * DBESH1 ( k * SQRT ( X**2 + Y**2 ) ) * ( X * D_x + Y * D_y ) / SQRT ( X**2 + Y**2 )

  END FUNCTION DG_0


END MODULE modSpecialFunctions
