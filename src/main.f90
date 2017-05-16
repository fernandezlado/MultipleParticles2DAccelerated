program MultipleParticle


  use modLinearAlgebra
  
  use modForwardMap

  use modParameters
  
  use modIO


  implicit none

  real(8) :: start, finish

  type (phy_Parameters) :: phy

  type (geo_Parameters) :: geo

  type (alg_Parameters) :: alg


  call loadData ( phy, geo, alg )


  write (*,*) "------------------------------"
  
  call cpu_time ( start )

  call initForwardMap ( phy, geo, alg )

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

  call cpu_time ( start )

  call TestForwardMap ( phy, geo, alg )

  call cpu_time (finish) 

  write(*,*) "Time elapsed: ", finish-start

  call destroyForwardMap()
  
  deallocate ( geo % defects_location )


end program MultipleParticle








