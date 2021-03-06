MODULE modIO



  USE modMathConstants

  USE modParameters



  IMPLICIT NONE



CONTAINS



  SUBROUTINE loadData (phy, geo, alg)

    
    TYPE(phy_Parameters) :: phy

    TYPE(geo_Parameters) :: geo

    TYPE(alg_Parameters) :: alg


    CHARACTER(200) :: line

    INTEGER :: IOstatus, Comments, j

    ! --------------------------------------------------------------------------------
    ! OPEN FILE
    ! --------------------------------------------------------------------------------

    OPEN (unit = 1, file = "parameters/parameters_1.txt", iostat = iostatus)


    IF ( iostatus .NE. 0 ) THEN
       
       WRITE (*,*) "Parameters file failed to open."

    ELSE 
       
       WRITE (*,*) "Parameters filed opened succesfully."

    END IF
    

    ! --------------------------------------------------------------------------------
    ! Physical Parameters Load
    ! --------------------------------------------------------------------------------

    IOSTATUS = 0

    COMMENTS = 0


    DO WHILE ( IOSTATUS>=0 .and. comments == 0)

   
       READ (1, FMT = '(A)',IOSTAT=IOSTATUS) line

       Comments =  index (line, '# Values')


    END DO


    READ(1,*) phy  % k, phy % angle
 
    phy % angle = Pi * phy % angle
    

    ! --------------------------------------------------------------------------------
    ! Array geometry Parameters Load
    ! --------------------------------------------------------------------------------


    IOSTATUS = 0

    COMMENTS = 0
    
    j = 0

    DO WHILE ( IOSTATUS>=0 .AND. COMMENTS == 0) 
    
       j = j + 1

       READ (1, FMT = '(A)',IOSTAT=IOSTATUS) line

       Comments =  index (line, '# Values')

    END DO


    READ (1, *) geo % radius

    READ (1, *) geo % L_x, geo % L_y

    READ (1, *) geo % N_row, geo % N_col

    READ (1, *) geo % N_def


    ALLOCATE ( geo % defects_location ( geo % N_def,2 ) )


    DO j = 1, geo % N_def

       
       READ (1,*) geo % defects_location ( j, 1 ), geo % defects_location ( j, 2 )


    END DO
    

    ! --------------------------------------------------------------------------------
    ! Algorithm Parameters Load
    ! --------------------------------------------------------------------------------


    IOSTATUS = 0

    COMMENTS = 0


    DO WHILE ( IOSTATUS>=0 .AND. COMMENTS == 0) 


       READ (1, FMT = '(A)',IOSTAT=IOSTATUS) line

       Comments =  index (line, '# Values')
    

    END DO


    READ(1,*) alg % N_src_hor, alg % N_src_ver

    READ(1,*) alg % N_dis


    alg % N_coll_hor = 2 * alg % N_src_hor

    alg % N_coll_ver = 2 * alg % N_src_ver

    alg % N_wave = ( alg % N_src_hor + alg % N_src_ver )


    ! --------------------------------------------------------------------------------


    CLOSE (1)

    
  END SUBROUTINE loadData



END MODULE modIO
