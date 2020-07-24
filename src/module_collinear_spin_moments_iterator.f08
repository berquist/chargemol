 MODULE module_collinear_spin_moments_iterator
 ! Copyright (c) 2014, 2015, 2016, 2017 by Thomas A. Manz and Nidia Gabaldon Limas. Rights reserved.
 
 USE module_precision
 USE module_common_variable_declarations
 USE module_global_parameters
 USE module_spin_functions
 USE module_generate_atomic_spin_moment_file

 IMPLICIT NONE
 
 CONTAINS
 
 SUBROUTINE collinear_spin_moments_iterator()
 !===================================================================================
 
 INTEGER :: k,local_atomic_spin_sign
 REAL(kind=dp) :: a,b,local_atomic_spin,proportional_spin,L,temp_vector(3),distance,local_atomic_density,&
 proportional_theta_scalar,theta_scalar,mag_L_vector,mag_local_atomic_spin_vector,inv_Xi
 REAL(kind=dp), ALLOCATABLE, DIMENSION(:) :: old
 REAL(kind=dp), ALLOCATABLE, DIMENSION(:,:) :: sum_spin_density
 REAL(kind=dp), ALLOCATABLE, DIMENSION(:,:,:) :: Ypsilon_projection,Omega_projection,trial_spin_density,&
 sum_Ypsilon_projection,sum_spin_pseudodensity
 INTEGER :: clock_count,clock_count_rate,clock_count_max,clock_count2
 
 WRITE(output_FID,*)' '
 WRITE(output_FID,*)'******************** TIME CONTROL ********************'
 WRITE(output_FID,*)'Starting collinear_spin_moments_iterator'
 FLUSH(output_FID)
 CALL system_clock(clock_count,clock_count_rate,clock_count_max)
 WRITE(output_FID,*)'******************************************************'
 WRITE(output_FID,*)' '
 FLUSH(output_FID)
 
 ! Initialize some arrays
 ALLOCATE(spherical_average_atomic_spin(nshells,natoms))
 ALLOCATE(Ypsilon_projection(totnumA,totnumB,totnumC))
 ALLOCATE(Omega_projection(totnumA,totnumB,totnumC))
 ALLOCATE(spin_population(natoms))
 ALLOCATE(old(natoms))
 ALLOCATE(trial_spin_density(totnumA,totnumB,totnumC))
 ALLOCATE(sum_spin_density(nshells,natoms))
 ALLOCATE(sum_Ypsilon_projection(totnumA,totnumB,totnumC))
 ALLOCATE(sum_spin_pseudodensity(totnumA,totnumB,totnumC))
 a = 0.0_dp
 b = 0.0_dp
 Ypsilon_projection=0.0_dp
 Omega_projection=0.0_dp
 spin_population=0.0_dp
 trial_spin_density=0.0_dp
 spherical_average_atomic_spin=0.0_dp
 local_atomic_spin = 0.0_dp
 proportional_spin = 0.0_dp
 L = 0.0_dp
 tot_spin_moment = 0.0_dp
 ! Iterative solution for the atomic spin distributions
 WRITE(output_FID,*)'Iteratively solving for the atomic spin distributions:'
 FLUSH(output_FID)
 DO spin_iter = 1,500
   sum_points = 0
   sum_spin_density = 0.0_dp
   sum_Ypsilon_projection = 0.0_dp
   sum_spin_pseudodensity = 0.0_dp
!$omp parallel default(none) &
!$omp private(j,nc,kc,nb,kb,na,ka,temp_vector,distance,shell_index,local_atomic_density,proportional_spin,&
!$omp proportional_theta_scalar,local_atomic_spin,L,local_atomic_spin_sign,a,b,theta_scalar,mag_L_vector,&
!$omp mag_local_atomic_spin_vector,inv_Xi) &
!$omp shared(natoms,lower_nc,upper_nc,kCpoints,delta_nc,lower_nb,upper_nb,kBpoints,delta_nb,lower_na,upper_na,kApoints,&
!$omp delta_na,boundary,center_shift,total_pseudodensity,partial_density,spin_density,Xi_lookup,spin_iter,&
!$omp Ypsilon_projection,Omega_projection,trial_spin_density,inverse_Xi_lookup,sum_Ypsilon_projection,sum_spin_pseudodensity,&
!$omp spherical_average_density,spherical_average_atomic_spin,chunk_size,total_density) &
!$omp reduction(+:sum_points,sum_spin_density)
   DO j=1,natoms
!$omp do &
!$omp schedule(dynamic,chunk_size)
     DO nc = lower_nc(j),upper_nc(j)
       kc = kCpoints(j,delta_nc + nc + 1)
       DO nb = lower_nb(j),upper_nb(j)
         kb = kBpoints(j,delta_nb + nb + 1)
         DO na = lower_na(j),upper_na(j)
           ka = kApoints(j,delta_na + na + 1)  
           temp_vector(1) = na*boundary(1,1) + nb*boundary(2,1) + nc*boundary(3,1) - center_shift(1,j)
           temp_vector(2) = na*boundary(1,2) + nb*boundary(2,2) + nc*boundary(3,2) - center_shift(2,j)
           temp_vector(3) = na*boundary(1,3) + nb*boundary(2,3) + nc*boundary(3,3) - center_shift(3,j)
           distance = sqrt(temp_vector(1)*temp_vector(1) + temp_vector(2)*temp_vector(2) + temp_vector(3)*temp_vector(3))
           shell_index = ceiling(scalefactor*distance + zero_tolerance)
           IF((shell_index <= nshells) .and. (total_pseudodensity(ka,kb,kc) > zero_tolerance) .and. ((total_density(ka,kb,kc))&
           > zero_tolerance)) THEN
             sum_points(shell_index,j) = sum_points(shell_index,j) + 1
             local_atomic_density =  partial_density(shell_index,j)*(total_density(ka,kb,kc))/total_pseudodensity(ka,kb,kc)
             proportional_spin = partial_density(shell_index,j)*spin_density(ka,kb,kc)/total_pseudodensity(ka,kb,kc)
             proportional_theta_scalar = calculate_theta_scalar(local_atomic_density,proportional_spin,Xi_lookup)
             IF (spin_iter == 1) THEN
               local_atomic_spin = proportional_spin
               L = proportional_theta_scalar
               IF (local_atomic_spin < 0.0_dp) THEN
                 local_atomic_spin_sign = -1
               ELSE
                 local_atomic_spin_sign = 1
               END IF
             ELSE
               a = spherical_average_density(shell_index,j)
               b = spherical_average_atomic_spin(shell_index,j)
               theta_scalar = calculate_theta_scalar(a,b,Xi_lookup)
               L = Ypsilon_projection(ka,kb,kc) - Omega_projection(ka,kb,kc) + pi*(spin_density(ka,kb,kc) - &
               trial_spin_density(ka,kb,kc))/(total_density(ka,kb,kc)) + (1.0_dp-spin_ref_fraction)*theta_scalar + &
               spin_ref_fraction*proportional_theta_scalar
               mag_L_vector = abs(L)
               IF (L < 0.0_dp) THEN
                 local_atomic_spin_sign = -1
               ELSE
                 local_atomic_spin_sign = 1
               END IF
               IF (mag_L_vector > pi) THEN
                 mag_local_atomic_spin_vector = local_atomic_density
               ELSE
                 inv_Xi = fast_calculate_inverse_Xi(2.0_dp*mag_L_vector,inverse_Xi_lookup)
                 mag_local_atomic_spin_vector = local_atomic_density*inv_Xi
               END IF
                 local_atomic_spin = mag_local_atomic_spin_vector*local_atomic_spin_sign
             END IF 
             sum_spin_density(shell_index,j) = sum_spin_density(shell_index,j) + local_atomic_spin
             sum_Ypsilon_projection(ka,kb,kc) = sum_Ypsilon_projection(ka,kb,kc) + &
			 L*partial_density(shell_index,j)/total_pseudodensity(ka,kb,kc)
             sum_spin_pseudodensity(ka,kb,kc) = sum_spin_pseudodensity(ka,kb,kc) + local_atomic_spin
           END IF
         END DO    
       END DO
     END DO
!$omp end do
   END DO  
!$omp end parallel
   Ypsilon_projection = sum_Ypsilon_projection
   trial_spin_density = sum_spin_pseudodensity
   DO j = 1,natoms
     DO k = 1,nshells
       IF (sum_points(k,j) > 0) THEN
         spherical_average_atomic_spin(k,j) = sum_spin_density(k,j)/sum_points(k,j)
       ELSE 
         spherical_average_atomic_spin(k,j) = 0.0_dp
       END IF
     END DO
   END DO
   WRITE(output_FID,*) 'iteration: ',spin_iter
   FLUSH(output_FID)
   old=spin_population
   spin_population = sum(sum_spin_density,dim=1)*pixelvolume
   DO j=1,natoms
     spin_population(j) = spin_population(j) + occupancy_correction(2,j)
   END DO
   WRITE(output_FID,*) 'spin_population: ',spin_population
   FLUSH(output_FID)
   tot_spin_moment = sum(spin_population)
   WRITE(output_FID,*) 'tot_spin_moment: ',tot_spin_moment
   FLUSH(output_FID)
   !Compute occupations
   max_change=0.0_dp
   DO j=1,natoms
     max_change = max(max_change,abs(spin_population(j)-old(j)))
   END DO
   WRITE(output_FID,*) 'max_change: ',max_change
   WRITE(output_FID,*) ' '
   FLUSH(output_FID)   
   IF (((max_change < spin_convergence_tolerance) .and. (spin_iter > 6)) .or. ((max(maxval(abs(spin_population)),max_change) < &
   spin_convergence_tolerance) .and. (spin_iter > 1))) THEN
     CALL generate_atomic_spin_moment_file( )
     EXIT
   END IF
   Omega_projection = 0.0_dp
!$omp parallel default(none) &
!$omp private(nc,kc,nb,kb,na,ka,temp_vector,distance,shell_index,local_atomic_density,proportional_spin,a,&
!$omp b,theta_scalar,proportional_theta_scalar,j) &
!$omp shared(lower_nc,upper_nc,lower_nb,upper_nb,kCpoints,kBpoints,lower_na,upper_na,kApoints,delta_nc,delta_nb,delta_na,&
!$omp boundary,center_shift,total_pseudodensity,partial_density,spin_density,total_density,&
!$omp spherical_average_density,spherical_average_atomic_spin,Xi_lookup,Omega_projection,chunk_size,natoms)
   DO j=1,natoms
!$omp do schedule(dynamic,chunk_size)
     DO nc = lower_nc(j),upper_nc(j)
       kc = kCpoints(j,delta_nc + nc + 1)
       DO nb = lower_nb(j),upper_nb(j)
         kb = kBpoints(j,delta_nb + nb + 1)
         DO na = lower_na(j),upper_na(j)
           ka = kApoints(j,delta_na + na + 1)
           temp_vector(1) = na*boundary(1,1) + nb*boundary(2,1) + nc*boundary(3,1) - center_shift(1,j)
           temp_vector(2) = na*boundary(1,2) + nb*boundary(2,2) + nc*boundary(3,2) - center_shift(2,j)
           temp_vector(3) = na*boundary(1,3) + nb*boundary(2,3) + nc*boundary(3,3) - center_shift(3,j)
           distance = sqrt(temp_vector(1)*temp_vector(1) + temp_vector(2)*temp_vector(2) + temp_vector(3)*temp_vector(3))
           shell_index = ceiling(scalefactor*distance + zero_tolerance)
           IF ((shell_index <= nshells) .and. (total_pseudodensity(ka,kb,kc) > zero_tolerance) .and. &
           ((total_density(ka,kb,kc)) > zero_tolerance)) THEN
             local_atomic_density =  partial_density(shell_index,j)*(total_density(ka,kb,kc))/total_pseudodensity(ka,kb,kc)
             proportional_spin = partial_density(shell_index,j)*spin_density(ka,kb,kc)/total_pseudodensity(ka,kb,kc)
             a = spherical_average_density(shell_index,j)
             b = spherical_average_atomic_spin(shell_index,j)
             theta_scalar = calculate_theta_scalar(a,b,Xi_lookup)
             proportional_theta_scalar = calculate_theta_scalar(local_atomic_density,proportional_spin,Xi_lookup)
             Omega_projection(ka,kb,kc) = Omega_projection(ka,kb,kc) + ((1.0_dp-spin_ref_fraction)*theta_scalar + &
             spin_ref_fraction*proportional_theta_scalar)*partial_density(shell_index,j)/total_pseudodensity(ka,kb,kc)
           END IF
         END DO    
       END DO
     END DO
!$omp end do
   END DO
!$omp end parallel
 END DO

 WRITE(output_FID,*) 'Final spin populations: ',spin_population
 FLUSH(output_FID)
 tot_spin_moment = sum(spin_population)
 WRITE(output_FID,*)'The total spin moment of the unit cell is: ',tot_spin_moment
 WRITE(output_FID,*)' '
 FLUSH(output_FID)
 
 !Free some memory by clearing unneeded arrays
 DEALLOCATE(Ypsilon_projection)
 DEALLOCATE(Omega_projection)
 DEALLOCATE(sum_Ypsilon_projection)
 DEALLOCATE(trial_spin_density)
 DEALLOCATE(sum_spin_pseudodensity)
 
 WRITE(output_FID,*)' '
 WRITE(output_FID,*)'******************** TIME CONTROL ********************'
 CALL system_clock(clock_count2,clock_count_rate,clock_count_max)
 IF (clock_count_rate == 0) THEN
   WRITE(output_FID,*)'Processor clock not available, the time calculation wil not be performed.'
 ELSE
   seconds = modulo((clock_count2-clock_count),clock_count_max)/real(clock_count_rate)
   WRITE(output_FID,*)'Finished collinear_spin_moments_iterator in ',seconds, 'seconds'
   FLUSH(output_FID)
 END IF
 WRITE(output_FID,*)'******************************************************'
 WRITE(output_FID,*)' '
 FLUSH(output_FID)  
 
 END SUBROUTINE collinear_spin_moments_iterator
 
 END MODULE module_collinear_spin_moments_iterator


