module modInterpolationReferenceCell
!-------------------------------------------------------------------------------
! CALTECH, CMS, Oscar Bruno's Group
!-------------------------------------------------------------------------------
! modInterpolationReferenceCell.f90 - Module that defines the 
! InterpolationReferenceCellell type and its methods.
!
! DESCRIPTION: 
!> defines Type projectionReferenceCell. Implements the class constructor, 
!> a routine to compute the QR factorization of the LeastSquare problem and
!> a routine to obtain the equivalent sources amplitudes.
!
!> @author
!> Agustin G. Fernandez-Lado
!
! MODIFIED: 15 May 2017
!-------------------------------------------------------------------------------


  use modMathConstants

  use modLinearAlgebra



  implicit none


 
  type InterpolationReferenceCell


     !Size of Cell
     real(8) :: L_x, L_y
     
     !Number of sources in first quadrant
     integer :: N_src

     !Number of plane waves in first quadrant to interpolate inside cell
     integer :: N_wave

     !Direction of plane waves
     real(8), dimension(:), allocatable :: u_x, u_y

     !Location of horizontal and vertical equivalent sources    
     real(8), dimension(:), allocatable :: src_x, src_y

     ! QR factorization of blocks for Least Square problem
     type(QR_factorization) :: QR_A, QR_B, QR_C, QR_D

     ! Auxiliary vectors to solve Least Square problems
     complex(8),dimension(:),allocatable :: b_hat_a, b_hat_b, b_hat_c, b_hat_d


  end type InterpolationReferenceCell



contains



  subroutine createInterpolationReferenceCell (this, N_s_h, N_s_v, N_w, L_x, L_y, k)
    
    ! --------------------------------------------------------------------------------
    !
    ! Class constructor. N_s_h, N_s_v are the number equivalent sources in the
    ! first quadrant paralell to the x and y axis (h=horizontal, v=vertical)
    ! N_w is the number of plane waves (with directions in the first quarant) to
    ! perform the plane wave expansion. L_x, L_y are the lengths of the sides of
    ! the cell. k is the wavenumber
    !
    ! --------------------------------------------------------------------------------

    type(InterpolationReferenceCell) :: this

    integer :: N_s_h, N_s_v, N_w

    real(8) :: L_x, L_y, k


    real(8),allocatable,dimension(:) :: disc

    integer :: j


    this % L_x = L_x

    this % L_y = L_y


    this % N_src = N_s_h + N_s_v 

    this % N_wave = N_w

    ! --------------------------------------------------------------------------------
    !
    ! Construct source points
    !
    ! --------------------------------------------------------------------------------

    allocate ( this % src_x ( 1 : this % N_src ), this % src_y ( 1 : this % N_src ) )


    this % src_y ( 1 : N_s_h ) = L_y * (/ ( 0.5d0, j = 1, N_s_h ) /)    

    this % src_x ( 1 : N_s_h ) = L_x * &
         (/ &
         ( ( 2.0d0 * j + 1.0d0 ) / ( 4.0d0 * N_s_h ) - 0.5d0, &
         j = N_s_h, 2 * N_s_h - 1 ) &
         /)


    this % src_x ( N_s_h + 1 : N_s_h + N_s_v ) = L_x * (/ ( 0.5d0, j = 1, N_s_v ) /)
    
    this % src_y ( N_s_h + 1 : N_s_h + N_s_v ) = L_y * &
         (/ &
         ( ( 2.0d0 * j + 1.0d0 ) / ( 4.0d0 * N_s_v ) - 0.5d0, &
         j= 2 * N_s_v - 1, N_s_v, -1 ) &
         /)

    ! --------------------------------------------------------------------------------
    !
    ! Construct direction of plane waves
    !
    ! --------------------------------------------------------------------------------

    allocate ( this % u_x ( 1 : N_w ), this % u_y ( 1 : N_w ) )

    
    this % u_x = &
         (/ & 
         ( cos ( pi / 2.0d0 * ( 2.0d0 * j + 1.0d0) / ( 2.0d0 * N_w ) ), &
         j =0, N_w - 1 ) &
         /)

    this % u_y = &
         (/ &
         ( sin ( pi / 2.0d0 * ( 2.0d0 * j + 1.0d0) / ( 2.0d0 * N_w ) ), &
         j = 0, N_w - 1 ) &
         /)


    ! --------------------------------------------------------------------------------
    !
    ! Construct QR factorizations of interaction matrices
    !
    ! --------------------------------------------------------------------------------


    call createQR_factorization ( this % QR_A , N_s_h + N_s_v, N_w )

    call createQR_factorization ( this % QR_B, N_s_h + N_s_v, N_w )

    call createQR_factorization ( this % QR_C , N_s_h + N_s_v, N_w )

    call createQR_factorization ( this % QR_D , N_s_h + N_s_v, N_w )


    call createMatrices (this, k)


    ! --------------------------------------------------------------------------------
    ! 
    ! Allocate vectors to store RHS of the least square problems
    !
    ! --------------------------------------------------------------------------------


    allocate ( this % b_hat_a ( 1 : this % N_src ) )

    allocate ( this % b_hat_b ( 1 : this % N_src ) )

    allocate ( this % b_hat_c ( 1 : this % N_src ) )

    allocate ( this % b_hat_d ( 1 : this % N_src ) )

    
  end subroutine createInterpolationReferenceCell


  
  subroutine createMatrices (this, k)

    ! --------------------------------------------------------------------------------
    !
    ! Computes the QR factorization of each block to solve
    ! the least square problem to interpolate the field at 
    ! the boundary inside the cell.
    !
    ! Read notes for a better comprehension of how the
    ! blocks are constructed.
    !
    ! --------------------------------------------------------------------------------

    type(InterpolationReferenceCell),target :: this

    real(8) :: k

    character :: type_flag


    integer :: j

    
    real(8),dimension(:,:),allocatable :: Q_x, Q_y

    real(8),dimension(:,:,:),allocatable :: U_x, U_y

    complex(8),dimension(:,:,:),allocatable :: E


    allocate ( Q_x ( 1 : this % N_src, 1 : this % N_wave) )

    allocate ( Q_y ( 1 : this % N_src, 1 : this % N_wave) )


    allocate ( U_x ( 1 : this % N_src, 1 : this % N_wave, 1:4 ) )

    allocate ( U_y ( 1 : this % N_src, 1 : this % N_wave, 1:4 ) )

    
    allocate ( E ( 1 : this % N_src, 1 : this % N_wave, 1:4 ) )

    ! --------------------------------------------------------------------------------
    !
    ! Spread direction of plane waves and location of equivalent sources
    !
    ! --------------------------------------------------------------------------------

    U_x (:,:,1) = spread ( this % u_x, 1, this % N_src )

    U_y (:,:,1) = spread ( this % u_y, 1, this % N_src )


    U_x (:,:,2) = -U_x (:,:,1) 

    U_y (:,:,2) =  U_y (:,:,1) 


    U_x (:,:,3) = -U_x (:,:,1) 

    U_y (:,:,3) = -U_y (:,:,1) 


    U_x (:,:,4) =  U_x (:,:,1) 

    U_y (:,:,4) = -U_y (:,:,1) 


    Q_x = spread ( this % src_x, 2, this % N_wave ) 
       
    Q_y = spread ( this % src_y, 2, this % N_wave )

       
    ! --------------------------------------------------------------------------------
    ! 
    ! Compute Plane Waves Evaluated at Equiv. sources
    !
    ! --------------------------------------------------------------------------------

    forall ( j = 1:4 ) 

       E(:,:,j) = exp ( I * k * ( Q_x * U_x (:,:,j) + Q_y * U_y (:,:,j) ) )

    end forall

    ! --------------------------------------------------------------------------------
    !
    ! Construct matrices of eigenvalues
    !
    ! --------------------------------------------------------------------------------


    this % QR_A % mat = E(:,:,1) + E(:,:,2) + E(:,:,3) + E(:,:,4) 

    this % QR_B % mat = E(:,:,1) - E(:,:,2) + E(:,:,3) - E(:,:,4)

    this % QR_C % mat = E(:,:,1) + E(:,:,2) - E(:,:,3) - E(:,:,4)

    this % QR_D % mat = E(:,:,1) - E(:,:,2) - E(:,:,3) + E(:,:,4)

    ! --------------------------------------------------------------------------------
    !
    ! Compute QR factorization
    !
    ! --------------------------------------------------------------------------------
    
    call computeQR_Factorization ( this % QR_A % mat, this % QR_A)
    
    call computeQR_Factorization ( this % QR_B % mat, this % QR_B)

    call computeQR_Factorization ( this % QR_C % mat, this % QR_C)
    
    call computeQR_Factorization ( this % QR_D % mat, this % QR_D)


    deallocate ( E )
    
    deallocate ( Q_x, Q_y, U_x, U_y )

        
  end subroutine createMatrices



  subroutine computePlaneWaveWeights (this, rhs, xi)

    ! --------------------------------------------------------------------------------
    !
    ! Solves a least square problem to obtain the amplitudes
    ! of the plane waves to interpolate inside the cell.
    !
    ! rhs contains the field in the boundary of the cell
    ! xi stores the amplitudes.
    !
    ! --------------------------------------------------------------------------------

    type(InterpolationReferenceCell),target :: this

    complex(8),target,dimension(:) :: rhs
    
    complex(8),target,dimension(:) :: xi


    integer,pointer :: N_s, N_w

    complex(8),pointer :: rhs_1(:), rhs_2(:), rhs_3(:), rhs_4(:)

    complex(8),pointer :: xi_1(:), xi_2(:), xi_3(:), xi_4(:)


    N_s => this % N_src

    N_w => this % N_wave


    ! --------------------------------------------------------------------------------
    !
    ! Split RHS
    !
    ! --------------------------------------------------------------------------------

    rhs_1 ( 1 : N_s ) => rhs ( 0 * N_s + 1 : 1 * N_s )

    rhs_2 ( 1 : N_s ) => rhs ( 1 * N_s + 1 : 2 * N_s )

    rhs_3 ( 1 : N_s ) => rhs ( 2 * N_s + 1 : 3 * N_s )

    rhs_4 ( 1 : N_s ) => rhs ( 3 * N_s + 1 : 4 * N_s)


    ! --------------------------------------------------------------------------------
    !
    ! Split unknown
    !
    ! --------------------------------------------------------------------------------

    xi_1 ( 1 : N_w ) => xi ( 0 * N_w + 1 : 1 * N_w )

    xi_2 ( 1 : N_w ) => xi ( 1 * N_w + 1 : 2 * N_w )

    xi_3 ( 1 : N_w ) => xi ( 2 * N_w + 1 : 3 * N_w )

    xi_4 ( 1 : N_w ) => xi ( 3 * N_w + 1 : 4 * N_w ) 


    ! --------------------------------------------------------------------------------
    !
    ! Apply orthonormal transformation to RHS
    !
    ! --------------------------------------------------------------------------------

    this % b_hat_a = 0.5d0 * (rhs_1 + rhs_2 + rhs_3 + rhs_4)

    this % b_hat_b = 0.5d0 * (rhs_1 - rhs_2 + rhs_3 - rhs_4)

    this % b_hat_c = 0.5d0 * (rhs_1 + rhs_2 - rhs_3 - rhs_4)

    this % b_hat_d = 0.5d0 * (rhs_1 - rhs_2 - rhs_3 + rhs_4)


    ! --------------------------------------------------------------------------------
    !
    ! Solve Least Square Problems
    !
    ! --------------------------------------------------------------------------------

    call LeastSquareSolve ( this % QR_A, this % b_hat_a)

    call LeastSquareSolve ( this % QR_B, this % b_hat_b)

    call LeastSquareSolve ( this % QR_C, this % b_hat_c) 

    call LeastSquareSolve ( this % QR_D, this % b_hat_d)


    ! --------------------------------------------------------------------------------
    !
    ! Apply orthonormal transformation to solutions
    !
    ! --------------------------------------------------------------------------------

    xi_1 = 0.5d0 * &
         ( &
         this % b_hat_a ( 1 : N_w ) + &
         this % b_hat_b ( 1 : N_w ) + &
         this % b_hat_c ( 1 : N_w ) + &
         this % b_hat_d ( 1 : N_w )   &
         )

    xi_2 = 0.5d0 * &
         ( &
         this % b_hat_a ( 1 : N_w ) - &
         this % b_hat_b ( 1 : N_w ) + &
         this % b_hat_c ( 1 : N_w ) - &
         this % b_hat_d ( 1 : N_w )   &
         )

    xi_3 = 0.5d0 * &
         ( &
         this % b_hat_a ( 1 : N_w ) + &
         this % b_hat_b ( 1 : N_w ) - &
         this % b_hat_c ( 1 : N_w ) - &
         this % b_hat_d ( 1 : N_w )   &
         )

    xi_4 = 0.5d0 * &
         ( &
         this % b_hat_a ( 1 : N_w ) - &
         this % b_hat_b ( 1 : N_w ) - &
         this % b_hat_c ( 1 : N_w ) + &
         this % b_hat_d ( 1 : N_w )   &
         )



    nullify ( rhs_1, rhs_2, rhs_3, rhs_4 )

    nullify ( xi_1, xi_2, xi_3, xi_4 )

    nullify ( N_s, N_w )


  end subroutine computePlaneWaveWeights



  subroutine destroyInterpolationReferenceCell (this)

    
    type(InterpolationReferenceCell) :: this


    deallocate ( this % b_hat_a, this % b_hat_b, this % b_hat_c, this % b_hat_d )
    
    deallocate ( this % u_x, this % u_y, this % src_x, this % src_y)


    call destroyQR_factorization ( this % QR_A )

    call destroyQR_factorization ( this % QR_B ) 

    call destroyQR_factorization ( this % QR_C )

    call destroyQR_factorization ( this % QR_D )


  end subroutine destroyInterpolationReferenceCell



end module modInterpolationReferenceCell
