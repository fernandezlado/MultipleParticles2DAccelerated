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

  CALL initForwardMap ( phy, geo, alg )

  CALL forwardMap ( phy, geo, alg )
  
  CALL TestEquivalentSources ( phy, geo, alg ) 

  call destroyForwardMap()
  
  DEALLOCATE ( geo % defects_location )


END PROGRAM MultipleParticle








