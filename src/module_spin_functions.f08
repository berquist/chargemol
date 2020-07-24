 !-------------------------------------------
 ! Module to calculate various spin related functions.
 ! Copyright (c) 2014, 2015, 2016, 2017 by Thomas A. Manz and Nidia Gabaldon Limas. Rights reserved.
 !-------------------------------------------

 MODULE module_spin_functions

 USE module_precision
 USE module_global_parameters

 IMPLICIT NONE

 !-----------------------------------------------------------------------------------
 !Public variables
 !-----------------------------------------------------------------------------------

 REAL(kind=dp), ALLOCATABLE, DIMENSION(:) :: Xi_lookup,inverse_Xi_lookup

 CONTAINS


 SUBROUTINE generate_spin_lookup_tables()
 !===================================================================================

 INTEGER :: i   
 ALLOCATE(Xi_lookup(num_lookup_points))
 ALLOCATE(inverse_Xi_lookup(num_lookup_points))
 Xi_lookup=0.0_dp
 inverse_Xi_lookup=0.0_dp
 DO i=1,num_lookup_points
   Xi_lookup(i) = calculate_Xi(real(i,dp)/num_lookup_points)
   inverse_Xi_lookup(i) = calculate_inverse_Xi(i*2.0_dp*pi/num_lookup_points)
 END DO
 
 END SUBROUTINE generate_spin_lookup_tables


 !-------------------------------------------
 !Calculating Xi.
 !-------------------------------------------
 ELEMENTAL FUNCTION calculate_Xi(tau) RESULT(Xi)

 REAL(kind=dp), INTENT(IN) :: tau
 REAL(kind=dp) :: Xi, Xi_derivative,inv_Xi, tau_estimate, Xi_estimate

 IF (tau < 0.0_dp) THEN
   Xi = 0.0_dp
 ELSE IF (tau > 1.0_dp) THEN
   Xi=2.0_dp*pi
 ELSE IF (tau < Xi_threshold) THEN
   Xi=4.0_dp*pi*tau/3.0_dp
 ELSE IF(tau > (1.0_dp - Xi_threshold) ) THEN
   Xi=2.0_dp*pi + (tau - 1.0_dp)*49.956_dp + (tau - 1.0_dp)**2*8370.6_dp
 ELSE
   Xi = pi*((1.0_dp/(tau**2) - 1.0_dp)*log((1.0_dp - tau)/(1.0_dp + tau)) + 2.0_dp/tau)
 END IF

 RETURN 

 END FUNCTION calculate_Xi


 !-------------------------------------------
 !Calculating Xi_derivative.
 !-------------------------------------------
 ELEMENTAL FUNCTION calculate_Xi_derivative(tau) RESULT(Xi_derivative)

 REAL(kind=dp), INTENT(IN) :: tau
 REAL(kind=dp) ::  Xi, Xi_derivative,inv_Xi, tau_estimate, Xi_estimate

 IF ( tau < Xi_threshold ) THEN
   Xi_derivative = 4.0_dp*pi/3.0_dp
 ELSE IF ( tau > ( 1.0_dp - Xi_threshold )) THEN
   Xi_derivative=50.0_dp
 ELSE
   Xi_derivative=2.0_dp * pi * (log (( 1.0_dp +tau ) / ( 1.0_dp -tau)) - 2.0_dp*tau)/(tau**3)
 END IF

 END FUNCTION calculate_Xi_derivative


 !-------------------------------------------
 !Calculating inverse_Xi.
 !-------------------------------------------
 ELEMENTAL FUNCTION calculate_inverse_Xi(Xi_value) RESULT(inv_Xi) 

 REAL(kind=dp), INTENT(IN) :: Xi_value
 REAL(kind=dp) ::  Xi_derivative,inv_Xi, tau_estimate, Xi_estimate
 INTEGER:: i

 IF (Xi_value < 0.0_dp) THEN
   inv_Xi=0.0_dp
 ELSE IF (Xi_value < Xi_zero_tolerance) THEN
   inv_Xi=3.0_dp * Xi_value/(4.0_dp *pi)
 ELSE IF (Xi_value > 2.0_dp*pi - Xi_zero_tolerance) THEN
   inv_Xi=1.0_dp
 ELSE
   tau_estimate= (Xi_value-0.036835_dp*Xi_value*Xi_value )/( 4.178319_dp-0.136129_dp*Xi_value+0.038148_dp*Xi_value*Xi_value)
   Xi_derivative = calculate_Xi_derivative(tau_estimate)
     do i=1,5
        Xi_estimate=calculate_Xi(tau_estimate)
        tau_estimate=tau_estimate+(Xi_value-Xi_estimate)/Xi_derivative
     END DO
     inv_Xi=tau_estimate
 END IF

 END FUNCTION calculate_inverse_Xi


 PURE FUNCTION fast_calculate_Xi(tau,Xi_lookup)
 !-----------------------------------------------------------------------------------
 
 INTEGER :: i  
 REAL(kind=dp),INTENT(IN) :: Xi_lookup(num_lookup_points)
 REAL(kind=dp),INTENT(IN) :: tau
 REAL(kind=dp) :: fast_calculate_Xi,threshold,f

 threshold = 1.0_dp/num_lookup_points
 IF (tau <= 0.0_dp) THEN
   fast_calculate_Xi = 0.0_dp
 ELSE IF (tau >= 1.0_dp) THEN
   fast_calculate_Xi = 2.0_dp*pi
 ELSE IF (tau < threshold) THEN
   fast_calculate_Xi = 4.0_dp*pi*tau/3.0_dp
 ELSE  
    i = floor(tau*num_lookup_points)
    f = tau*num_lookup_points - i
    fast_calculate_Xi = (1.0_dp-f)*Xi_lookup(i) + f*Xi_lookup(i+1)
 END IF

 END FUNCTION fast_calculate_Xi


 PURE FUNCTION fast_calculate_inverse_Xi(Xi_value,inverse_Xi_lookup)
 !-----------------------------------------------------------------------------------
 
 REAL(kind=dp),INTENT(IN) :: Xi_value,inverse_Xi_lookup(num_lookup_points)
 REAL(kind=dp) :: fast_calculate_inverse_Xi,threshold,scaled_Xi_value,f
 INTEGER :: i
 
 threshold = 1.0_dp/num_lookup_points
 scaled_Xi_value = Xi_value/(2.0_dp*pi)
 IF (scaled_Xi_value <= 0.0_dp) THEN
    fast_calculate_inverse_Xi = 0.0_dp
 ELSE IF (scaled_Xi_value >= 1.0_dp) THEN
    fast_calculate_inverse_Xi = 1.0_dp
 ELSE IF (scaled_Xi_value < threshold) THEN
    fast_calculate_inverse_Xi = Xi_value*0.75_dp/pi
 ELSE  
    i = floor(scaled_Xi_value*num_lookup_points)
    f = scaled_Xi_value*num_lookup_points - i
    fast_calculate_inverse_Xi = (1.0_dp-f)*inverse_Xi_lookup(i) + f*inverse_Xi_lookup(i+1)
 END IF

 END FUNCTION fast_calculate_inverse_Xi 
 

 PURE FUNCTION calculate_theta_scalar(a,b_projection,Xi_lookup)
 !-----------------------------------------------------------------------------------
 
 REAL(kind=dp),INTENT(IN) :: a,b_projection,Xi_lookup(num_lookup_points)
 REAL(kind=dp) :: calculate_theta_scalar,tau,Xi
 REAL(kind=dp), PARAMETER :: local_tolerance = 1.0e-10_dp 
 
 IF ((abs(b_projection) < local_tolerance) .or. (a < local_tolerance)) THEN
    calculate_theta_scalar = 0.0_dp
 ELSE
    tau = (abs(b_projection)/a)
    Xi = fast_calculate_Xi(tau,Xi_lookup)
    calculate_theta_scalar = b_projection*Xi/(2.0_dp*abs(b_projection))
 END IF
 
 END FUNCTION calculate_theta_scalar


 PURE FUNCTION calculate_theta_vector(a,b_vector,Xi_lookup)
 !-----------------------------------------------------------------------------------
 
 REAL(kind=dp),INTENT(IN) :: a,b_vector(3),Xi_lookup(num_lookup_points)
 REAL(kind=dp) :: calculate_theta_vector(3),mag_b_vector,tau,Xi,b(3)
 REAL(kind=dp), PARAMETER :: local_tolerance = 1.0e-10_dp 
 
 mag_b_vector = sqrt(dot_product(b_vector,b_vector))
 IF ((mag_b_vector > local_tolerance) .and. (a > local_tolerance)) THEN
   tau = (mag_b_vector/a)
   Xi = fast_calculate_Xi(tau,Xi_lookup)
   calculate_theta_vector = b_vector*Xi/(2.0_dp*mag_b_vector)
 ELSE
   calculate_theta_vector = 0.0_dp
 END IF   

 END FUNCTION calculate_theta_vector
 

 END MODULE module_spin_functions
