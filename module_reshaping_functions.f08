 MODULE module_reshaping_functions
 !======================================================================================
 ! Reshaping functions that makes sure that the atomic density decays monotonically and
 ! its buried tails do not become too contracted or too diffuse.
 ! Copyright (c) 2014, 2015, 2016, 2017 by Thomas A. Manz and Nidia Gabaldon Limas. Rights reserved.
 !======================================================================================
 USE module_precision

 IMPLICIT NONE
 
 CONTAINS

 PURE SUBROUTINE monotonic_decay_function(unreshaped_partial_density,nshells,radial_shell_integration_weight,&
 reshaped_partial_density,local_reshaping_iterations)

 !--------------------------------------------------------------------------------------------------------- 
 !Make the conditioned reference density monotonically decreasing with increasing radius
 !--------------------------------------------------------------------------------------------------------- 
 
 INTEGER,INTENT(IN)::nshells
 REAL(kind=dp),INTENT(IN),DIMENSION(:)::unreshaped_partial_density,radial_shell_integration_weight
 REAL(kind=dp),DIMENSION(:),INTENT(OUT)::reshaped_partial_density 
 INTEGER,INTENT(OUT)::local_reshaping_iterations
 INTEGER :: k,step,trial_num
 REAL(kind=dp)::temp,reshape_error,reshaping_error(5),add_coefficient(5),single_second_integration_sum,m_factor,&
 weighted_points,single_first_integration_sum
 REAL(kind=dp),PARAMETER::reshaping_convergence_tolerance=1.0e-10_dp
 REAL(kind=dp),ALLOCATABLE,DIMENSION(:)::sqrt_sigma

 !Monotonic decay with increasing radius
 ALLOCATE(sqrt_sigma(nshells))
 sqrt_sigma = 0.0_dp
 reshaped_partial_density = 0.0_dp
 reshaped_partial_density = unreshaped_partial_density
 !Calculate sqrt_sigma and weighted_points
 weighted_points = 0.0_dp
 single_first_integration_sum = 0.0_dp
 local_reshaping_iterations = 1
 DO k = 1,nshells
   sqrt_sigma(k)=sqrt(unreshaped_partial_density(k))
   weighted_points=weighted_points+sqrt_sigma(k)*radial_shell_integration_weight(k)
   single_first_integration_sum = single_first_integration_sum + unreshaped_partial_density(k)*radial_shell_integration_weight(k)
 END DO
 add_coefficient = 0.0_dp
 reshaping_error = 0.0_dp
 IF ((weighted_points<reshaping_convergence_tolerance) .or. (single_first_integration_sum<reshaping_convergence_tolerance)) THEN
   RETURN
 END IF
 reshaped_partial_density(1) = max(reshaped_partial_density(1),0.0_dp) 
 temp = reshaped_partial_density(1) 
 DO k = 2,nshells
    reshaped_partial_density(k) = max(reshaped_partial_density(k),0.0_dp) 
    IF ((radial_shell_integration_weight(k) > 0.0_dp) .and. (temp > 0.0_dp)) THEN
       reshaped_partial_density(k) = min(reshaped_partial_density(k),temp) 
    END IF
    temp = reshaped_partial_density(k) 
 END DO
 single_second_integration_sum = 0.0_dp
 DO k = 1,nshells
    single_second_integration_sum = single_second_integration_sum + reshaped_partial_density(k)*radial_shell_integration_weight(k)  
 END DO
 reshape_error = single_first_integration_sum - single_second_integration_sum
 reshaping_error(1) = reshape_error 
 IF (reshape_error < reshaping_convergence_tolerance) THEN
   RETURN
 END IF

 ! Set the upper bound
 DO step = 1,50
   add_coefficient(2) = 2.0_dp*add_coefficient(1) + reshape_error/weighted_points
   DO k = 1,nshells
     reshaped_partial_density(k)=unreshaped_partial_density(k) + add_coefficient(2)*sqrt_sigma(k)
   END DO
   reshaped_partial_density(1) = max(reshaped_partial_density(1),0.0_dp) 
   temp = reshaped_partial_density(1) 
   DO k = 2,nshells
      reshaped_partial_density(k) = max(reshaped_partial_density(k),0.0_dp) 
      IF ((radial_shell_integration_weight(k) > 0.0_dp) .and. (temp > 0.0_dp)) THEN
         reshaped_partial_density(k) = min(reshaped_partial_density(k),temp) 
      END IF
      temp = reshaped_partial_density(k) 
   END DO
   single_second_integration_sum = 0.0_dp
   DO k = 1,nshells
      single_second_integration_sum = single_second_integration_sum + reshaped_partial_density(k)*radial_shell_integration_weight(k)
   END DO
   reshape_error = single_first_integration_sum - single_second_integration_sum
   IF (reshape_error > 0.0_dp) THEN
      add_coefficient(1) = add_coefficient(2)
      reshaping_error(1) = reshape_error
   ELSE
      add_coefficient(5) = add_coefficient(2)
      reshaping_error(5) = reshape_error
      EXIT
   END IF
 END DO 

 DO trial_num = 1,50
   IF (abs(reshaping_error(5) - reshaping_error(1)) < reshaping_convergence_tolerance)  THEN
      EXIT
   END IF
   local_reshaping_iterations = trial_num 
   !Set the midpoint
   add_coefficient(3) = 0.5_dp*(add_coefficient(1) + add_coefficient(5))
   DO k = 1,nshells
     reshaped_partial_density(k)=unreshaped_partial_density(k) + add_coefficient(3)*sqrt_sigma(k) 
   END DO
   reshaped_partial_density(1) = max(reshaped_partial_density(1),0.0_dp) 
   temp = reshaped_partial_density(1) 
   DO k = 2,nshells
      reshaped_partial_density(k) = max(reshaped_partial_density(k),0.0_dp) 
      IF ((radial_shell_integration_weight(k) > 0.0_dp) .and. (temp > 0.0_dp)) THEN
         reshaped_partial_density(k) = min(reshaped_partial_density(k),temp) 
      END IF
      temp = reshaped_partial_density(k) 
   END DO
   single_second_integration_sum = 0.0_dp
   DO k = 1,nshells
      single_second_integration_sum = single_second_integration_sum + reshaped_partial_density(k)*radial_shell_integration_weight(k)
   END DO
   reshaping_error(3) = single_first_integration_sum - single_second_integration_sum
   IF (abs(reshaping_error(3)) < reshaping_convergence_tolerance) THEN
      EXIT
   END IF
   ! Fit the curve y = m_factor*x^2 + (1-m_factor)*x
   m_factor = 2.0_dp - 4.0_dp*(reshaping_error(1) - reshaping_error(3))/(reshaping_error(1) - reshaping_error(5)) 
   ! Compute the target add_coefficient
   IF (abs(m_factor) < reshaping_convergence_tolerance) THEN
      add_coefficient(4)=add_coefficient(1)+(add_coefficient(5)-add_coefficient(1))*(reshaping_error(1)/(reshaping_error(1)&
      -reshaping_error(5)))
   ELSE
     temp = reshaping_error(1)/(reshaping_error(1) - reshaping_error(5))
     add_coefficient(4) = add_coefficient(1) + (add_coefficient(5) - add_coefficient(1))*(m_factor-1.0_dp+sqrt((1.0_dp-m_factor&
     )**2 + 4.0_dp*m_factor*temp))/(2.0_dp*m_factor)
   END IF
   DO k = 1,nshells
     reshaped_partial_density(k)=unreshaped_partial_density(k) + add_coefficient(4)*sqrt_sigma(k) 
   END DO
   reshaped_partial_density(1) = max(reshaped_partial_density(1),0.0_dp) 
   temp = reshaped_partial_density(1) 
   DO k = 2,nshells
      reshaped_partial_density(k) = max(reshaped_partial_density(k),0.0_dp) 
      IF ((radial_shell_integration_weight(k) > 0.0_dp) .and. (temp > 0.0_dp)) THEN
         reshaped_partial_density(k) = min(reshaped_partial_density(k),temp) 
      END IF
      temp = reshaped_partial_density(k) 
   END DO
   single_second_integration_sum = 0.0_dp
   DO k = 1,nshells
      single_second_integration_sum=single_second_integration_sum+reshaped_partial_density(k)*radial_shell_integration_weight(k)
   END DO
   reshaping_error(4) = single_first_integration_sum - single_second_integration_sum
   IF (abs(reshaping_error(4)) < reshaping_convergence_tolerance) THEN
      EXIT
   END IF
   ! Compute the second target add_coefficient
   IF ((reshaping_error(3) > 0.0_dp) .AND. (reshaping_error(4) > 0.0_dp)) THEN
      IF (3.0_dp*abs(reshaping_error(5)) < abs(reshaping_error(4))) THEN
         add_coefficient(2)=add_coefficient(5)+(add_coefficient(5)-add_coefficient(4))*(2.0_dp*reshaping_error(5))/(&
         reshaping_error(4) - reshaping_error(5)) 
      ELSE IF (3.0_dp*abs(reshaping_error(4)) < abs(reshaping_error(5))) THEN
         add_coefficient(2)=add_coefficient(4)+(add_coefficient(5)-add_coefficient(4))*(2.0_dp*reshaping_error(4))/(&
         reshaping_error(4) - reshaping_error(5))
      ELSE
         add_coefficient(2) = 0.5_dp*(add_coefficient(4) + add_coefficient(5))
      END IF
   ELSE IF ((reshaping_error(3) > 0.0_dp) .AND. (reshaping_error(4) <= 0.0_dp)) THEN 
      IF (3.0_dp*abs(reshaping_error(3)) < abs(reshaping_error(4))) THEN
         add_coefficient(2)=add_coefficient(3)+(add_coefficient(4)-add_coefficient(3))*(2.0_dp*reshaping_error(3))/(&
         reshaping_error(3) - reshaping_error(4)) 
      ELSE IF (3.0_dp*abs(reshaping_error(4)) < abs(reshaping_error(3))) THEN
         add_coefficient(2)=add_coefficient(4)+(add_coefficient(4)-add_coefficient(3))*(2.0_dp*reshaping_error(4))/(&
         reshaping_error(3) - reshaping_error(4))
      ELSE
         add_coefficient(2) = 0.5_dp*(add_coefficient(4) + add_coefficient(3))
      END IF
   ELSE IF ((reshaping_error(3) <= 0.0_dp) .AND. (reshaping_error(4) > 0.0_dp)) THEN
      IF (3.0_dp*abs(reshaping_error(3)) < abs(reshaping_error(4))) THEN
         add_coefficient(2) = add_coefficient(3) + (add_coefficient(3) - add_coefficient(4))*(2.0_dp*reshaping_error(3))/&
         (reshaping_error(4) - reshaping_error(3)) 
      ELSE IF (3.0_dp*abs(reshaping_error(4)) < abs(reshaping_error(3))) THEN
         add_coefficient(2) = add_coefficient(4) + (add_coefficient(3) - add_coefficient(4))*(2.0_dp*reshaping_error(4))/(&
         reshaping_error(4) - reshaping_error(3))
      ELSE
         add_coefficient(2) = 0.5_dp*(add_coefficient(4) + add_coefficient(3))
      END IF        
   ELSE
      IF (3.0_dp*abs(reshaping_error(1)) < abs(reshaping_error(4))) THEN
         add_coefficient(2) = add_coefficient(1) + (add_coefficient(4) - add_coefficient(1))*(2.0_dp*reshaping_error(1))/(&
         reshaping_error(1) - reshaping_error(4)) 
      ELSE IF (3.0_dp*abs(reshaping_error(4)) < abs(reshaping_error(1))) THEN
         add_coefficient(2) = add_coefficient(4) + (add_coefficient(4) - add_coefficient(1))*(2.0_dp*reshaping_error(4))/(&
         reshaping_error(1) - reshaping_error(4))
      ELSE
         add_coefficient(2) = 0.5_dp*(add_coefficient(4) + add_coefficient(1))
      END IF           
   END IF

   DO k = 1,nshells
     reshaped_partial_density(k)=unreshaped_partial_density(k) + add_coefficient(2)*sqrt_sigma(k) 
   END DO
   reshaped_partial_density(1) = max(reshaped_partial_density(1),0.0_dp) 
   temp = reshaped_partial_density(1) 
   DO k = 2,nshells
      reshaped_partial_density(k) = max(reshaped_partial_density(k),0.0_dp) 
      IF ((radial_shell_integration_weight(k) > 0.0_dp) .and. (temp > 0.0_dp)) THEN
         reshaped_partial_density(k) = min(reshaped_partial_density(k),temp) 
      END IF
      temp = reshaped_partial_density(k) 
   END DO
   single_second_integration_sum = 0.0_dp
   DO k = 1,nshells
      single_second_integration_sum = single_second_integration_sum+reshaped_partial_density(k)*radial_shell_integration_weight(k)
   END DO
   reshaping_error(2) = single_first_integration_sum - single_second_integration_sum
   IF (abs(reshaping_error(2)) < reshaping_convergence_tolerance) THEN
      EXIT
   END IF
   IF (reshaping_error(2) > 0.0_dp) THEN
      add_coefficient(1) = add_coefficient(2)
      reshaping_error(1) = reshaping_error(2)
   ELSE
      add_coefficient(5) = add_coefficient(2)
      reshaping_error(5) = reshaping_error(2)
   END IF
   IF (reshaping_error(3) > 0.0_dp) THEN
     add_coefficient(1) = max(add_coefficient(3),add_coefficient(1))
     reshaping_error(1) = min(reshaping_error(3),reshaping_error(1))
   ELSE
     add_coefficient(5) = min(add_coefficient(3),add_coefficient(5))
     reshaping_error(5) = max(reshaping_error(3),reshaping_error(5))
   END IF
   IF (reshaping_error(4) > 0.0_dp) THEN
     add_coefficient(1) = max(add_coefficient(4),add_coefficient(1))
     reshaping_error(1) = min(reshaping_error(4),reshaping_error(1))
   ELSE
     add_coefficient(5) = min(add_coefficient(4),add_coefficient(5))
     reshaping_error(5) = max(reshaping_error(4),reshaping_error(5))
   END IF
 END DO


 DEALLOCATE(sqrt_sigma)


 END SUBROUTINE monotonic_decay_function
 !========================================================================================================================
 


 PURE SUBROUTINE tail_exponential_decay_function(unreshaped_partial_density,single_first_exp_const,single_second_exp_const,&
 nshells,radial_shell_integration_weight,reshaped_partial_density,local_reshaping_iterations)

 !--------------------------------------------------------------------------------------------------------- 
 ! Prevent the atoms from becoming too diffuse or too contracted
 !---------------------------------------------------------------------------------------------------------
 
 INTEGER,INTENT(IN)::nshells
 REAL(kind=dp),INTENT(IN),DIMENSION(:)::single_first_exp_const,single_second_exp_const
 REAL(kind=dp),INTENT(IN),DIMENSION(:)::unreshaped_partial_density,radial_shell_integration_weight
 REAL(kind=dp),DIMENSION(:),INTENT(OUT)::reshaped_partial_density
 INTEGER,INTENT(OUT)::local_reshaping_iterations
 INTEGER :: k,step,trial_num
 REAL(kind=dp)::temp,reshape_error,reshaping_error(5),add_coefficient(5),m_factor,&
 weighted_points,single_first_integration_sum,single_second_integration_sum
 REAL(kind=dp),PARAMETER::reshaping_convergence_tolerance=1.0e-10_dp
 REAL(kind=dp),ALLOCATABLE,DIMENSION(:)::sqrt_sigma

 ALLOCATE(sqrt_sigma(nshells))
 sqrt_sigma = 0.0_dp
 reshaped_partial_density = 0.0_dp
 reshaped_partial_density = unreshaped_partial_density
 !Calculate sqrt_sigma and weighted_points
 single_first_integration_sum=0.0_dp
 weighted_points=0.0_dp
 local_reshaping_iterations = 1
 DO k = 1,nshells
   sqrt_sigma(k)=sqrt(unreshaped_partial_density(k))
   single_first_integration_sum=single_first_integration_sum+unreshaped_partial_density(k)*radial_shell_integration_weight(k)
   weighted_points=weighted_points+sqrt_sigma(k)*radial_shell_integration_weight(k)
 END DO
 !Reshape the atomic weighting factors to not be too diffuse
 add_coefficient = 0.0_dp
 reshaping_error = 0.0_dp
 IF ((weighted_points<reshaping_convergence_tolerance) .or. (single_first_integration_sum<reshaping_convergence_tolerance)) THEN
   RETURN
 END IF
 reshaped_partial_density(1) = max(reshaped_partial_density(1),0.0_dp)
 temp = reshaped_partial_density(1)
 DO k = 2,nshells
   reshaped_partial_density(k) = max(reshaped_partial_density(k),0.0_dp)
   IF ((radial_shell_integration_weight(k) > 0.0_dp) .and. (temp > 0.0_dp)) THEN
     reshaped_partial_density(k) = min(reshaped_partial_density(k),temp*single_first_exp_const(k))
   END IF
   temp = reshaped_partial_density(k)
 END DO

 single_second_integration_sum = 0.0_dp
 DO k = 1,nshells
   single_second_integration_sum = single_second_integration_sum + reshaped_partial_density(k)*radial_shell_integration_weight(k)
 END DO
 reshape_error = single_first_integration_sum - single_second_integration_sum
 reshaping_error(1) = reshape_error 
 IF (reshape_error < reshaping_convergence_tolerance) THEN
   RETURN
 END IF
 ! Set the upper bound
 DO step = 1,50
   add_coefficient(2) = 2.0_dp*add_coefficient(1) + reshape_error/weighted_points
   DO k = 1,nshells
     reshaped_partial_density(k) = unreshaped_partial_density(k) + add_coefficient(2)*sqrt_sigma(k) 
   END DO
   reshaped_partial_density(1) = max(reshaped_partial_density(1),0.0_dp)
   temp = reshaped_partial_density(1)
   DO k = 2,nshells
     reshaped_partial_density(k) = max(reshaped_partial_density(k),0.0_dp)
     IF ((radial_shell_integration_weight(k) > 0.0_dp) .and. (temp > 0.0_dp)) THEN
       reshaped_partial_density(k) = min(reshaped_partial_density(k),temp*single_first_exp_const(k))
     END IF
     temp = reshaped_partial_density(k)
   END DO
   single_second_integration_sum = 0.0_dp
   DO k = 1,nshells
     single_second_integration_sum = single_second_integration_sum + reshaped_partial_density(k)*radial_shell_integration_weight(k)
   END DO
   reshape_error = single_first_integration_sum - single_second_integration_sum
   IF (reshape_error > 0.0_dp) THEN
     add_coefficient(1) = add_coefficient(2)
     reshaping_error(1) = reshape_error
   ELSE
     add_coefficient(5) = add_coefficient(2)
     reshaping_error(5) = reshape_error
     EXIT
   END IF
 END DO    
 DO trial_num = 1,50
   IF (abs(reshaping_error(5) - reshaping_error(1)) < reshaping_convergence_tolerance)  THEN
     EXIT
   END IF
   local_reshaping_iterations = trial_num
   !Set the midpoint
   add_coefficient(3) = 0.5_dp*(add_coefficient(1) + add_coefficient(5))
   DO k = 1,nshells
     reshaped_partial_density(k) = unreshaped_partial_density(k) + add_coefficient(3)*sqrt_sigma(k) 
   END DO
   reshaped_partial_density(1) = max(reshaped_partial_density(1),0.0_dp)
   temp = reshaped_partial_density(1)
   DO k = 2,nshells
     reshaped_partial_density(k) = max(reshaped_partial_density(k),0.0_dp)
     IF ((radial_shell_integration_weight(k) > 0.0_dp) .and. (temp > 0.0_dp)) THEN
       reshaped_partial_density(k) = min(reshaped_partial_density(k),temp*single_first_exp_const(k))
     END IF
       temp = reshaped_partial_density(k)
   END DO
   single_second_integration_sum = 0.0_dp
   DO k = 1,nshells
     single_second_integration_sum = single_second_integration_sum + reshaped_partial_density(k)*radial_shell_integration_weight(k)
   END DO
   reshaping_error(3) = single_first_integration_sum - single_second_integration_sum
   IF (abs(reshaping_error(3)) < reshaping_convergence_tolerance) THEN
     EXIT
   END IF
   ! Fit the curve y = m_factor*x^2 + (1-m_factor)*x
   m_factor = 2.0_dp - 4.0_dp*(reshaping_error(1) - reshaping_error(3))/(reshaping_error(1) - reshaping_error(5)) 
   ! Compute the target add_coefficient
   IF (abs(m_factor) < reshaping_convergence_tolerance) THEN
     add_coefficient(4)=add_coefficient(1)+(add_coefficient(5)-add_coefficient(1))*(reshaping_error(1)/(reshaping_error(1)&
     -reshaping_error(5)))
   ELSE
     temp = reshaping_error(1)/(reshaping_error(1) - reshaping_error(5))
     add_coefficient(4) = add_coefficient(1) + (add_coefficient(5) - add_coefficient(1))&
     *(m_factor - 1.0_dp + sqrt( (1.0_dp - m_factor)**2 + 4.0_dp*m_factor*temp))/(2.0_dp*m_factor)
   END IF
   DO k = 1,nshells
     reshaped_partial_density(k) = unreshaped_partial_density(k) + add_coefficient(4)*sqrt_sigma(k) 
   END DO
   reshaped_partial_density(1) = max(reshaped_partial_density(1),0.0_dp)
   temp = reshaped_partial_density(1)
   DO k = 2,nshells
     reshaped_partial_density(k) = max(reshaped_partial_density(k),0.0_dp)
     IF ((radial_shell_integration_weight(k) > 0.0_dp) .and. (temp > 0.0_dp)) THEN
       reshaped_partial_density(k) = min(reshaped_partial_density(k),temp*single_first_exp_const(k))
     END IF
     temp = reshaped_partial_density(k)
   END DO
   single_second_integration_sum = 0.0_dp
   DO k = 1,nshells
     single_second_integration_sum = single_second_integration_sum + reshaped_partial_density(k)*radial_shell_integration_weight(k)
   END DO
   reshaping_error(4) = single_first_integration_sum - single_second_integration_sum
   IF (abs(reshaping_error(4)) < reshaping_convergence_tolerance) THEN
     EXIT
   END IF
   ! Compute the second target add_coefficient
   IF ((reshaping_error(3) > 0.0_dp) .AND. (reshaping_error(4) > 0.0_dp)) THEN
     IF (3.0_dp*abs(reshaping_error(5)) < abs(reshaping_error(4))) THEN
       add_coefficient(2) = add_coefficient(5) + (add_coefficient(5) - add_coefficient(4))*(2.0_dp*reshaping_error(5))/(&
       reshaping_error(4) - reshaping_error(5)) 
     ELSE IF (3.0_dp*abs(reshaping_error(4)) < abs(reshaping_error(5))) THEN
       add_coefficient(2) = add_coefficient(4) + (add_coefficient(5) - add_coefficient(4))*(2.0_dp*reshaping_error(4))/(&
       reshaping_error(4) - reshaping_error(5))
     ELSE
       add_coefficient(2) = 0.5_dp*(add_coefficient(4) + add_coefficient(5))
     END IF
   ELSE IF ((reshaping_error(3) > 0.0_dp) .AND. (reshaping_error(4) <= 0.0_dp)) THEN 
     IF (3.0_dp*abs(reshaping_error(3)) < abs(reshaping_error(4))) THEN
       add_coefficient(2) = add_coefficient(3) + (add_coefficient(4) - add_coefficient(3))*(2.0_dp*reshaping_error(3))/(&
       reshaping_error(3) - reshaping_error(4)) 
     ELSE IF (3.0_dp*abs(reshaping_error(4)) < abs(reshaping_error(3))) THEN
       add_coefficient(2) = add_coefficient(4) + (add_coefficient(4) - add_coefficient(3))*(2.0_dp*reshaping_error(4))/(&
       reshaping_error(3) - reshaping_error(4))
     ELSE
       add_coefficient(2) = 0.5_dp*(add_coefficient(4) + add_coefficient(3))
     END IF
   ELSE IF ((reshaping_error(3) <= 0.0_dp) .AND. (reshaping_error(4) > 0.0_dp)) THEN
     IF (3.0_dp*abs(reshaping_error(3)) < abs(reshaping_error(4))) THEN
       add_coefficient(2) = add_coefficient(3) + (add_coefficient(3) - add_coefficient(4))*(2.0_dp*reshaping_error(3))/(&
       reshaping_error(4) - reshaping_error(3)) 
     ELSE IF (3.0_dp*abs(reshaping_error(4)) < abs(reshaping_error(3))) THEN
       add_coefficient(2) = add_coefficient(4) + (add_coefficient(3) - add_coefficient(4))*(2.0_dp*reshaping_error(4))/(&
       reshaping_error(4) - reshaping_error(3))
     ELSE
       add_coefficient(2) = 0.5_dp*(add_coefficient(4) + add_coefficient(3))
     END IF        
   ELSE
     IF (3.0_dp*abs(reshaping_error(1)) < abs(reshaping_error(4))) THEN
       add_coefficient(2) = add_coefficient(1) + (add_coefficient(4) - add_coefficient(1))*(2.0_dp*reshaping_error(1))/(&
       reshaping_error(1) - reshaping_error(4)) 
     ELSE IF (3.0_dp*abs(reshaping_error(4)) < abs(reshaping_error(1))) THEN
       add_coefficient(2) = add_coefficient(4) + (add_coefficient(4) - add_coefficient(1))*(2.0_dp*reshaping_error(4))/(&
       reshaping_error(1) - reshaping_error(4))
     ELSE
       add_coefficient(2) = 0.5_dp*(add_coefficient(4) + add_coefficient(1))
     END IF           
   END IF
   DO k = 1,nshells
     reshaped_partial_density(k) = unreshaped_partial_density(k) + add_coefficient(2)*sqrt_sigma(k) 
   END DO
   reshaped_partial_density(1) = max(reshaped_partial_density(1),0.0_dp)
   temp = reshaped_partial_density(1)
   DO k = 2,nshells
     reshaped_partial_density(k) = max(reshaped_partial_density(k),0.0_dp)
     IF ((radial_shell_integration_weight(k) > 0.0_dp) .and. (temp > 0.0_dp)) THEN
       reshaped_partial_density(k) = min(reshaped_partial_density(k),temp*single_first_exp_const(k))
     END IF
     temp = reshaped_partial_density(k)
   END DO
   single_second_integration_sum = 0.0_dp
   DO k = 1,nshells
     single_second_integration_sum = single_second_integration_sum + reshaped_partial_density(k)*radial_shell_integration_weight(k)
   END DO
   reshaping_error(2) = single_first_integration_sum - single_second_integration_sum
   IF (abs(reshaping_error(2)) < reshaping_convergence_tolerance) THEN
      EXIT
   END IF
   IF (reshaping_error(2) > 0.0_dp) THEN
     add_coefficient(1) = add_coefficient(2)
     reshaping_error(1) = reshaping_error(2)
   ELSE
     add_coefficient(5) = add_coefficient(2)
     reshaping_error(5) = reshaping_error(2)
   END IF
   IF (reshaping_error(3) > 0.0_dp) THEN
     add_coefficient(1) = max(add_coefficient(3),add_coefficient(1))
     reshaping_error(1) = min(reshaping_error(3),reshaping_error(1))
   ELSE
     add_coefficient(5) = min(add_coefficient(3),add_coefficient(5))
     reshaping_error(5) = max(reshaping_error(3),reshaping_error(5))
   END IF
   IF (reshaping_error(4) > 0.0_dp) THEN
     add_coefficient(1) = max(add_coefficient(4),add_coefficient(1))
     reshaping_error(1) = min(reshaping_error(4),reshaping_error(1))
   ELSE
     add_coefficient(5) = min(add_coefficient(4),add_coefficient(5))
     reshaping_error(5) = max(reshaping_error(4),reshaping_error(5))
   END IF
 END DO
 !Reshaping to prevent atoms from becoming too contracted    
 reshaped_partial_density(1) = max(reshaped_partial_density(1),0.0_dp)
 temp = reshaped_partial_density(1)
 DO k = 2,nshells
   reshaped_partial_density(k) = max(reshaped_partial_density(k),temp*single_second_exp_const(k))
   temp = reshaped_partial_density(k)
 END DO
 single_second_integration_sum = 0.0_dp
 DO k = 1,nshells
     single_second_integration_sum = single_second_integration_sum + reshaped_partial_density(k)*radial_shell_integration_weight(k)
 END DO
 IF (abs(single_second_integration_sum) > reshaping_convergence_tolerance) THEN
   DO k = 1,nshells
     reshaped_partial_density(k) = reshaped_partial_density(k)*single_first_integration_sum/single_second_integration_sum
   END DO 
 END IF
 
 END SUBROUTINE tail_exponential_decay_function

 !========================================================================================================================
 
 END MODULE module_reshaping_functions
  