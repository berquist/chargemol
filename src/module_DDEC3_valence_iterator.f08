 MODULE module_DDEC3_valence_iterator
 ! Copyright (c) 2014, 2015, 2016, 2017 by Thomas A. Manz and Nidia Gabaldon Limas. Rights reserved.

 USE module_precision
 USE module_global_parameters
 USE module_common_variable_declarations
 USE module_initialize_atomic_densities
 USE module_local_multipole_moment_analysis
 USE module_update_atomic_densities

 IMPLICIT NONE
 
 CONTAINS
 
 SUBROUTINE DDEC3_valence_iterator()
 !===================================================================================
 ! Solves the atomic partial charge distributions iteratively
 !===================================================================================
 
 REAL(kind=dp), ALLOCATABLE, DIMENSION(:) :: wA_renormalization,wA_renormalization_slope,old,first_integration_sum,&
 weighted_points,second_integration_sum
 REAL(kind=dp), ALLOCATABLE, DIMENSION (:,:) :: sum_PtoWref,sum_sqrt_Wref,&
 sum_wA_over_sqrt_Wref,old_spherical_average_density,avg_sqrt_Wref,avg_wA_over_sqrt_Wref,sigma
 REAL(kind=dp), ALLOCATABLE, DIMENSION (:,:,:) :: ref_pseudodensity,conditioned_ref_pseudodensity
 REAL(kind=dp) :: old_old_old_density_change,old_old_density_change,old_density_change,max_density_change,&
 nonoverlapping_atom_tolerance,local_atomic_density,temp_vector(3),&
 temporal_vector(8),ISA_part,IH_part,add_coefficient,temp,distance,&
 first_constraint_term,second_constraint_term,second_exp_const,first_exp_const
 INTEGER :: k,trial_num,reshaping_iterations
 INTEGER :: clock_count,clock_count_rate,clock_count_max,clock_count2
 
 WRITE(output_FID,*)' '
 WRITE(output_FID,*)'******************** TIME CONTROL ********************'
 WRITE(output_FID,*)'Starting valence_iterator'
 FLUSH(output_FID)
 CALL system_clock(clock_count,clock_count_rate,clock_count_max)
 WRITE(output_FID,*)'******************************************************'
 WRITE(output_FID,*)' '
 FLUSH(output_FID)
 
 WRITE(output_FID,*) 'Iteratively solving for the atomic partial charge distributions:'
 FLUSH(output_FID)
 reference_weighting = 3.0_dp/14.0_dp
 WRITE(output_FID,*)'Reference weighting: ',reference_weighting
 FLUSH(output_FID)
 ALLOCATE(wA_renormalization(natoms))
 wA_renormalization=1.0_dp
 old_old_old_density_change = 0.1_dp
 old_old_density_change = 0.1_dp
 old_density_change = 0.1_dp
 max_density_change = 0.1_dp
 nonoverlapping_atom_tolerance = charge_convergence_tolerance/(10.0_dp*wA_renormalization_max)
 ALLOCATE(total_density(totnumA,totnumB,totnumC))
 total_density = valence_density + core_density
 DEALLOCATE(valence_density)
 ALLOCATE(ref_pseudodensity(totnumA,totnumB,totnumC))
 ALLOCATE(conditioned_ref_pseudodensity(totnumA,totnumB,totnumC))
 ALLOCATE(sum_PtoWref(nshells,natoms))
 ALLOCATE(sum_sqrt_Wref(nshells,natoms))
 ALLOCATE(sum_wA_over_sqrt_Wref(nshells,natoms))
 ALLOCATE(wA_renormalization_slope(natoms))
 ALLOCATE(old_spherical_average_density(nshells,natoms))
 ALLOCATE(avg_sqrt_Wref(nshells,natoms))
 ALLOCATE(avg_wA_over_sqrt_Wref(nshells,natoms))
 ALLOCATE(old(natoms))
 ALLOCATE(first_integration_sum(natoms))
 ALLOCATE(weighted_points(natoms))
 ALLOCATE(sigma(nshells,natoms))
 ALLOCATE(second_integration_sum(natoms))
 ALLOCATE(net_atomic_charge(natoms))
 DO iter = 1,2000
   total_pseudodensity = 0.0_dp
   ref_pseudodensity = 0.0_dp
   conditioned_ref_pseudodensity = 0.0_dp 
   sum_points = 0
   sum_density = 0.0_dp
   sum_PtoWref = 0.0_dp
   sum_sqrt_Wref = 0.0_dp
   sum_wA_over_sqrt_Wref = 0.0_dp
   wA_renormalization_slope = 0.0_dp
!$omp parallel default(none) &
!$omp private(j,na,nb,nc,ka,kb,kc,temp_vector,distance,shell_index,local_atomic_density) &
!$omp shared(lower_nc,upper_nc,kCpoints,delta_nc,lower_nb,upper_nb,kBpoints,delta_nb,lower_na,upper_na,kApoints,delta_na,boundary,&
!$omp ref_pseudodensity,combined_oxidation_density,conditioned_ref_pseudodensity,conditioned_oxidation_density,total_pseudodensity,&
!$omp partial_density,natoms,pixelvolume,chunk_size,center_shift,total_density) &
!$omp reduction(+:sum_density,sum_points,sum_PtoWref,sum_sqrt_Wref,sum_wA_over_sqrt_Wref,wA_renormalization_slope)
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
           IF ((shell_index <= nshells)) THEN
             ref_pseudodensity(ka,kb,kc) = ref_pseudodensity(ka,kb,kc) + combined_oxidation_density(shell_index,j)
             conditioned_ref_pseudodensity(ka,kb,kc) = conditioned_ref_pseudodensity(ka,kb,kc) + &
             conditioned_oxidation_density(shell_index,j)
             total_pseudodensity(ka,kb,kc) = total_pseudodensity(ka,kb,kc) + partial_density(shell_index,j)
           END IF
         END DO    
       END DO
     END DO
!$omp end do
   END DO
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
           IF ((shell_index <= nshells) .and. (total_pseudodensity(ka,kb,kc) > zero_tolerance)) THEN 
             local_atomic_density =  partial_density(shell_index,j)*total_density(ka,kb,kc)/total_pseudodensity(ka,kb,kc)
             sum_density(shell_index,j) = sum_density(shell_index,j) + local_atomic_density
             sum_points(shell_index,j) = sum_points(shell_index,j) + 1
             IF (ref_pseudodensity(ka,kb,kc) > zero_tolerance**2) THEN
                sum_PtoWref(shell_index,j) = sum_PtoWref(shell_index,j) + total_density(ka,kb,kc)/&
                ref_pseudodensity(ka,kb,kc)
             ELSE
                sum_PtoWref(shell_index,j) = sum_PtoWref(shell_index,j) + 1.0_dp
             END IF
             sum_sqrt_Wref(shell_index,j) = sum_sqrt_Wref(shell_index,j) + sqrt&
             (conditioned_ref_pseudodensity(ka,kb,kc))
             IF (conditioned_ref_pseudodensity(ka,kb,kc) > zero_tolerance**2) THEN
                sum_wA_over_sqrt_Wref(shell_index,j) = sum_wA_over_sqrt_Wref(shell_index,j) + &
                conditioned_oxidation_density(shell_index,j)/sqrt(conditioned_ref_pseudodensity(ka,kb,kc))
             END IF
             wA_renormalization_slope(j) = wA_renormalization_slope(j) + pixelvolume*(1.0_dp - partial_density(shell_index,j)/&
             total_pseudodensity(ka,kb,kc))*local_atomic_density
           END IF
         END DO    
       END DO
     END DO
!$omp end do
   END DO
!$omp end parallel
   old_spherical_average_density = spherical_average_density
   DO j = 1,natoms
     DO k = 1,nshells
       IF (sum_points(k,j) > 0) THEN
         spherical_average_density(k,j) = sum_density(k,j)/sum_points(k,j)
         avg_PtoWref(k,j) = sum_PtoWref(k,j)/sum_points(k,j)
         avg_sqrt_Wref(k,j) = sum_sqrt_Wref(k,j)/sum_points(k,j)
         avg_wA_over_sqrt_Wref(k,j) = sum_wA_over_sqrt_Wref(k,j)/sum_points(k,j)
       ELSE 
         spherical_average_density(k,j) = 0.0_dp
         avg_PtoWref(k,j) = 1.0_dp
         avg_sqrt_Wref(k,j) = sqrt(zero_tolerance)
         avg_wA_over_sqrt_Wref(k,j) = sqrt(zero_tolerance)
       END IF
     END DO
   END DO
   WRITE(output_FID,*) "iteration= ",iter
   FLUSH(output_FID)
   ! Compute occupations and population changes
   old=valence_population
   old_old_old_change = old_old_change
   old_old_change = old_change
   old_change = max_change
   old_old_old_density_change = old_old_density_change
   old_old_density_change = old_density_change
   old_density_change = max_density_change
   valence_population = sum(sum_density,1)*pixelvolume
   DO j = 1,natoms
     valence_population(j) = valence_population(j) + occupancy_correction(1,j) - core_population(j)
   END DO
   normalization=nvalence/sum(valence_population)
   WRITE(output_FID,'(a,f13.9)')' Normalization: ', normalization
   FLUSH(output_FID)
   valence_population=normalization*valence_population
   net_atomic_charge = 0.0_dp
   DO j = 1,natoms
     net_atomic_charge(j) = atomic_number(j) - core_electrons(j) - valence_population(j)
   END DO
   WRITE(output_FID,*)'Net atomic charges for the current iteration: '
   WRITE(output_FID,'(10f14.6)')net_atomic_charge
   FLUSH(output_FID)
   change = valence_population - old
   max_change = maxval(abs(change))
   WRITE(output_FID,'(a,es13.6)')' Maximum change: ',max_change
   FLUSH(output_FID)
   max_density_change = maxval(abs(old_spherical_average_density - spherical_average_density))
   WRITE(output_FID,'(a,es13.6)')' Maximum density change: ',max_density_change
   FLUSH(output_FID)
   temporal_vector = [max_change,old_change,old_old_change,old_old_old_change,max_density_change,&
   old_density_change,old_old_density_change,old_old_old_density_change]
   IF  ((maxval(temporal_vector) < charge_convergence_tolerance) .and. (iter > 9)) THEN
     EXIT
   END IF
   IF (iter == 1) THEN
     WRITE(output_FID,*) 'Information for noniterative Hirshfeld method will be printed now.'
     FLUSH(output_FID)
     CALL local_multipole_moment_analysis()
     WRITE(output_FID,*) 'Hirshfeld analysis finished, calculation of iterative AIM will proceed.'
     FLUSH(output_FID)
   END IF
   !Update atomic densities
   IF (reference_weighting == 0.0_dp) THEN
     partial_density = spherical_average_density
   ELSE
     first_integration_sum = 0.0_dp
     weighted_points = 0.0_dp
     sigma = 0.0_dp
     DO j = 1,natoms
       DO k = 1,nshells
         ISA_part = spherical_average_density(k,j)**(1.0_dp-reference_weighting)
         conditioned_oxidation_density(k,j) = combined_oxidation_density(k,j)*avg_PtoWref(k,j)
         IH_part = conditioned_oxidation_density(k,j)**reference_weighting
         sigma(k,j) = ISA_part*IH_part
         first_integration_sum(j) = first_integration_sum(j) + sigma(k,j)*sum_points(k,j)
         weighted_points(j) = weighted_points(j) + sum_points(k,j)*sqrt(sigma(k,j))
       END DO
     END DO
     !Initialize the partial density array with an estimate
     partial_density = sigma
     reference_ion_charge = net_atomic_charge
     CALL update_atomic_densities()
     IF (iter > 3) THEN
       !Reshaping to prevent atoms from becoming too diffuse
       second_integration_sum = 0.0_dp
       reshaping_iterations = 1
       DO j = 1,natoms
         IF ((weighted_points(j) < zero_tolerance) .or. (first_integration_sum(j) < zero_tolerance)) THEN
           CYCLE
         END IF
         DO trial_num = 1,30
           partial_density(1,j) = max(partial_density(1,j),0.0_dp)
           temp = partial_density(1,j)
           DO k = 2,nshells
             first_constraint_term = minimum_buried_tail_exponent*(1.0_dp-(avg_wA_over_sqrt_Wref(k,j)/avg_sqrt_Wref(k,j))**2)
             !The decay exponent, convert to per radial shell
             first_exp_const = exp(-first_constraint_term*cutoff_radius*bohrperangstrom/(100.0_dp*nshells))
             partial_density(k,j) = max(partial_density(k,j),0.0_dp)
             IF ((sum_points(k,j) > 0) .and. (temp > 0.0_dp)) THEN
               partial_density(k,j) = min(partial_density(k,j),temp*first_exp_const)
             END IF
             temp = partial_density(k,j)
           END DO
           second_integration_sum(j) = 0.0_dp
           DO k = 1,nshells
             second_integration_sum(j) = second_integration_sum(j) + partial_density(k,j)*sum_points(k,j)
           END DO
           !Compute the constant to add to each point
           IF (weighted_points(j) > zero_tolerance) THEN
             add_coefficient = (first_integration_sum(j) - second_integration_sum(j))/weighted_points(j)
           ELSE
             add_coefficient = 0.0_dp
           END IF
           IF (abs(add_coefficient) > zero_tolerance) THEN
             reshaping_iterations = max(reshaping_iterations, trial_num)
             DO k = 1,nshells
               partial_density(k,j) = partial_density(k,j) + add_coefficient*sqrt(sigma(k,j))
             END DO
           ELSE
             EXIT
           END IF
         END DO
         !The following constraint makes sure the number of valence electrons assigned to each atom is non-negative
         IF (wA_renormalization_slope(j) > nonoverlapping_atom_tolerance) THEN
           wA_renormalization(j) = max(1.0_dp, (wA_renormalization(j) - 0.2_dp*valence_population(j)/&
           wA_renormalization_slope(j)))
         ELSE
           wA_renormalization(j) = 1.0_dp
         END IF
         DO k = 1,nshells
           partial_density(k,j) = partial_density(k,j)*wA_renormalization(j)
         END DO
       END DO
     WRITE(output_FID,'(A)')'Iterations to converge reshaping: '
     WRITE(output_FID,'(I5)') reshaping_iterations
     END IF
   END IF
 END DO
 IF (iter >= 2001) THEN
   WRITE(output_FID,*) 'Sorry, the calculation of atomic charges did not converge within 2000 iterations.'
   FLUSH(output_FID)
 END IF
 
 DEALLOCATE(conditioned_ref_pseudodensity)
 DEALLOCATE(ref_pseudodensity)
 
 WRITE(output_FID,*)' '
 WRITE(output_FID,*)'******************** TIME CONTROL ********************'
 CALL system_clock(clock_count2,clock_count_rate,clock_count_max)
 IF (clock_count_rate == 0) THEN
   WRITE(output_FID,*)'Processor clock not available, the time calculation wil not be performed.'
 ELSE
   seconds = modulo((clock_count2-clock_count),clock_count_max)/real(clock_count_rate)
   WRITE(output_FID,*)'Finished valence_iterator in ',seconds, 'seconds'
   FLUSH(output_FID)
 END IF
 WRITE(output_FID,*)'******************************************************'
 WRITE(output_FID,*)' '
 FLUSH(output_FID)
 
 END SUBROUTINE DDEC3_valence_iterator
 
 END MODULE module_DDEC3_valence_iterator
  