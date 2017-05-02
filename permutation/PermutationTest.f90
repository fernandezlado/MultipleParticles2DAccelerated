PROGRAM main

  IMPLICIT NONE

  TYPE IntegerPointer

     REAL,POINTER :: ptr
     
  END type IntegerPointer

  REAL(8)::start,finish
  
  INTEGER,PARAMETER :: n = 100

  REAL,TARGET,DIMENSION ( 4 * n ) :: storage
  REAL,DIMENSION( 4*n, 4*n) :: matrix

  INTEGER,DIMENSION ( 4 * n ) :: sigma

  TYPE (IntegerPointer),DIMENSION( 4 * n )  :: original, permuted

  INTEGER :: i,j,res


  FORALL ( i = 1:4*n ) original(i) % ptr => storage ( i )

  FORALL ( i = 0:4*n - 1 ) sigma (i+1) = n * mod ( i, 4 ) +  (i-mod ( i , 4 ) ) / 4 + 1

  FORALL ( i = 1:4*n ) permuted (i) % ptr => original ( sigma(i) ) % ptr

  CALL CPU_TIME (start)
  
  CALL random_number(matrix)

  do j=1,10000
     res=0
     CALL random_number(storage)
     do i=1,4*N
        res=res+permuted(i) % ptr
     end do
  END DO

  CALL CPU_TIME (finish)

  write(*,*) finish-start
  
  CALL CPU_TIME (start)

  do j=1,10000
     res=0
     CALL random_number(storage)
     do i=1,4*N
        res=res+storage(sigma(i))
     end do
  END DO
  
  CALL CPU_TIME (finish)

  write(*,*) finish-start


END PROGRAM main


