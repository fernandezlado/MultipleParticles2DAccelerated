MODULE modFFT



  USE MKL_DFTI



  IMPLICIT NONE

  

  TYPE FFT_2D


     TYPE(DFTI_DESCRIPTOR), POINTER :: desc_handle

     INTEGER :: status
     
     INTEGER :: N_x, N_y
     
     COMPLEX(8),DIMENSION(:,:),ALLOCATABLE :: fft_mat


  END TYPE FFT_2D



CONTAINS


  
  SUBROUTINE createFFT_2D ( this, N_x, N_y )

    
    TYPE(FFT_2D) :: this

    INTEGER :: N_x, N_y

    
    this % N_x = N_x

    this % N_y = N_y
    

    this % status = DftiCreateDescriptor ( this % desc_handle, DFTI_DOUBLE, DFTI_COMPLEX, 2, (/ N_x, N_y /) )

    this % status = DftiSetValue ( this % desc_handle, DFTI_BACKWARD_SCALE, 1.0d0 / ( N_x * N_y ) )

    this % status = DftiSetValue ( this % desc_handle, DFTI_PLACEMENT, DFTI_NOT_INPLACE )

    this % status = DftiCommitDescriptor ( this % desc_handle )


    ALLOCATE ( this % fft_mat ( N_x, N_y ) )

    
  END SUBROUTINE createFFT_2D
  
  

  SUBROUTINE directFFT ( this, x_in )

   
    TYPE(FFT_2D) :: this

    COMPLEX(8),DIMENSION(:,:) :: x_in
    
    
    this % status = DftiComputeForward ( this % desc_handle, x_in (:, 1), this % fft_mat ( :, 1 ) )


  END SUBROUTINE directFFT



  SUBROUTINE inverseFFT ( this, x_in)


    TYPE(FFT_2D) :: this

    COMPLEX(8),DIMENSION(:,:) :: x_in
    

    this % status = DftiComputeBackward ( this % desc_handle, x_in (:, 1), this % fft_mat ( :, 1 )  )


  END SUBROUTINE inverseFFT


  
  SUBROUTINE destroyFFT_2D ( this, x )

    
    TYPE(FFT_2D) :: this

    COMPLEX(8),DIMENSION(:) :: x

    
    this % status = DftiFreeDescriptor ( this % desc_handle )

    DEALLOCATE( this % fft_mat )


  END SUBROUTINE destroyFFT_2D


END MODULE modFFT
