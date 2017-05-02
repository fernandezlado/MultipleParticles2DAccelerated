MODULE modParameters



  IMPLICIT NONE


  
  TYPE phy_Parameters


     REAL(8) :: k

     REAL(8) :: angle
     
     
  end type phy_Parameters



  TYPE geo_Parameters


     REAL (8) :: radius

     REAL (8) :: L_x, L_y

     INTEGER :: N_row, N_col

     INTEGER :: N_def

     INTEGER,DIMENSION(:,:),ALLOCATABLE :: defects_location


  END TYPE geo_Parameters



  TYPE alg_Parameters

     
     INTEGER :: N_src_hor, N_src_ver

     INTEGER :: N_coll_hor, N_coll_ver

     INTEGER :: N_wave

     INTEGER :: N_dis


  END type alg_Parameters



END MODULE
