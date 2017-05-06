PROGRAM MultipleParticle



  USE modForwardMap

  USE modParameters
  
  USE modIO


  IMPLICIT NONE

  REAL(8) :: start, finish

  TYPE (phy_Parameters) :: phy

  TYPE (geo_Parameters) :: geo

  TYPE (alg_Parameters) :: alg


  CALL loadData ( phy, geo, alg )

  write(*,*) "------------------------------"
  
  call cpu_time ( start )

  CALL initForwardMap ( phy, geo, alg )

  call cpu_time ( finish )


  write(*,*) "------------------------------"

  write(*,*) "Array summary"

  write(*,*) "------------------------------"
  
  write(*,*) geo % N_row, " rows,", geo % N_col, " cols." 

  write(*,*) "------------------------------"

  write(*,*) "Time elapsed: ", finish-start

  write(*,*) "------------------------------"

  write(*,*) " Compute Forward map"

  write(*,*) "------------------------------"

  call cpu_time ( start )

  CALL forwardMap ( phy, geo, alg )
  
  call cpu_time ( finish )

  write(*,*) "Time elapsed: ", finish-start

  CALL TestEquivalentSources ( phy, geo, alg ) 

  call destroyForwardMap()
  
  DEALLOCATE ( geo % defects_location )


END PROGRAM MultipleParticle








