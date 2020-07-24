  MODULE module_DDEC6_valence_iterator
 ! Copyright (c) 2014, 2015, 2016, 2017 by Thomas A. Manz and Nidia Gabaldon Limas. Rights reserved.

 USE module_precision
 USE module_global_parameters
 USE module_common_variable_declarations
 USE module_initialize_atomic_densities
 USE module_local_multipole_moment_analysis
 USE module_update_atomic_densities
 USE module_compute_CM5
 USE module_reshaping_functions
!$ USE omp_lib

 IMPLICIT NONE
 
 CONTAINS
 
 SUBROUTINE DDEC6_valence_iterator()
 !===================================================================================
 ! Solves the atomic partial charge distributions iteratively
 !===================================================================================
 
 REAL(kind=dp), ALLOCATABLE, DIMENSION(:) :: wA_renormalization,wA_renormalization_slope,&
 conditioned_ref_population,target_conditioned_ref_population
 REAL(kind=dp), ALLOCATABLE, DIMENSION (:,:) :: sum_sqrt_Wref,sum_weighted_density,&
 sum_wA_over_sqrt_Wref,sigma,old_weighted_spherical_average_density,weighted_spherical_average_density,&
 first_exp_const,sum_weighted_points,second_exp_const,radial_shell_integration_weight
 REAL(kind=dp), ALLOCATABLE, DIMENSION (:,:,:) :: ref_pseudodensity,conditioned_ref_pseudodensity
 REAL(kind=dp) :: max_density_change,nonoverlapping_atom_tolerance,temp_vector(3),&
 temp,distance,first_constraint_term,second_constraint_term,local_atomic_density
 INTEGER :: k,reshaping_iterations,completed_steps,local_reshaping_iterations
 INTEGER :: clock_count,clock_count_rate,clock_count_max,clock_count2,clock_counta,clock_countb,clock_countc,clock_countd,&
 clock_counte
 LOGICAL :: update_kappa

 WRITE(output_FID,*)' '
 WRITE(output_FID,*)'******************** TIME CONTROL ********************'
 WRITE(output_FID,*)'Starting valence_iterator'
 FLUSH(output_FID)
 CALL system_clock(clock_count,clock_count_rate,clock_count_max)
 WRITE(output_FID,*)'******************************************************'
 WRITE(output_FID,*)' '
 FLUSH(output_FID)
 
 WRITE(output_FID,*) 'Iteratively solving for the atomic partial charge distributions:'
 ALLOCATE(wA_renormalization(natoms))
 wA_renormalization=1.0_dp
 max_density_change = 0.1_dp
 nonoverlapping_atom_tolerance = charge_convergence_tolerance/(10.0_dp*wA_renormalization_max)
 ALLOCATE(total_density(totnumA,totnumB,totnumC))
 total_density = valence_density + core_density
 DEALLOCATE(valence_density)
 ALLOCATE(tau(nshells,natoms))
 tau = 1.0_dp
 ALLOCATE(ref_pseudodensity(totnumA,totnumB,totnumC))
 ALLOCATE(sum_sqrt_Wref(nshells,natoms))
 ALLOCATE(sum_wA_over_sqrt_Wref(nshells,natoms))
 ALLOCATE(wA_renormalization_slope(natoms))
 ALLOCATE(sigma(nshells,natoms))
 sigma = 0.0_dp
 ALLOCATE(first_exp_const(nshells,natoms))
 first_exp_const = 0.0_dp
 ALLOCATE(Hirshfeld_net_atomic_charge(natoms))
 Hirshfeld_net_atomic_charge = 0.0_dp
 ALLOCATE(reference_ion_charge(natoms))
 reference_ion_charge = 0.0_dp
 ALLOCATE(localized_net_atomic_charge(natoms))
 localized_net_atomic_charge = 0.0_dp
 ALLOCATE(net_atomic_charge(natoms))
 net_atomic_charge = 0.0_dp
 ALLOCATE(old_net_atomic_charge(natoms))
 old_net_atomic_charge = 0.0_dp
 ALLOCATE(localized_pseudodensity(totnumA,totnumB,totnumC))
 ALLOCATE(localized_population(natoms))
 ALLOCATE(conditioned_ref_population(natoms))
 conditioned_ref_population = 0.0_dp
 ALLOCATE(target_conditioned_ref_population(natoms))
 target_conditioned_ref_population = 0.0_dp
 ALLOCATE(weighted_spherical_average_density(nshells,natoms))
 weighted_spherical_average_density = 0.0_dp
 ALLOCATE(old_weighted_spherical_average_density(nshells,natoms))
 ALLOCATE(sum_weighted_points(nshells,natoms))
 ALLOCATE(sum_weighted_density(nshells,natoms))
 ALLOCATE(second_exp_const(nshells,natoms))
 !The non-iterative part with exponential tail constraints
 DO iter=1,2
   ref_pseudodensity = 0.0_dp
   localized_pseudodensity = 0.0_dp
   valence_population = 0.0_dp 
   localized_population = 0.0_dp 
!$omp parallel default(none) &
!$omp private(j,na,nb,nc,ka,kb,kc,temp_vector,distance,shell_index) &
!$omp shared(lower_nc,upper_nc,kCpoints,delta_nc,lower_nb,upper_nb,kBpoints,delta_nb,lower_na,upper_na,kApoints,delta_na,&
!$omp boundary,natoms,chunk_size,center_shift,ref_pseudodensity,combined_oxidation_density,localized_pseudodensity,&
!$omp total_density,pixelvolume) &
!$omp reduction (+:valence_population,localized_population)
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
           IF (shell_index <= nshells) THEN
             ref_pseudodensity(ka,kb,kc) = ref_pseudodensity(ka,kb,kc) + combined_oxidation_density(shell_index,j)
             localized_pseudodensity(ka,kb,kc) = localized_pseudodensity(ka,kb,kc) + combined_oxidation_density(shell_index,j)&
             **localizing_power
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
           IF ((shell_index <= nshells) .and. (ref_pseudodensity(ka,kb,kc) > zero_tolerance) ) THEN
             valence_population(j) = valence_population(j)  + pixelvolume*combined_oxidation_density(shell_index,j)*&
             total_density(ka,kb,kc)/ref_pseudodensity(ka,kb,kc) 
             localized_population(j) = localized_population(j) + pixelvolume*(combined_oxidation_density(shell_index,j)**&
             localizing_power)*total_density(ka,kb,kc)/localized_pseudodensity(ka,kb,kc) 
           END IF
         END DO    
       END DO
     END DO
!$omp end do
   END DO
!$omp end parallel
   DO j=1,natoms
     valence_population(j) = valence_population(j) + occupancy_correction(1,j) - core_population(j) 
     localized_population(j) = localized_population(j) + occupancy_correction(1,j) - core_electrons(j) 
   END DO
   normalization=nvalence/sum(valence_population)
   WRITE(output_FID,*)'Normalization =',normalization
   valence_population=normalization*valence_population 
   localized_population=normalization*localized_population 
   DO j = 1,natoms
     net_atomic_charge(j) = atomic_number(j) - core_electrons(j) - valence_population(j) 
     localized_net_atomic_charge(j) = atomic_number(j) - core_electrons(j) - localized_population(j) 
   END DO
   IF (iter == 1) THEN
     Hirshfeld_net_atomic_charge = net_atomic_charge
     WRITE(output_FID,*) 'Information for noniterative Hirshfeld method will be printed now.'
     FLUSH(output_FID)
     partial_density = combined_oxidation_density
     total_pseudodensity = ref_pseudodensity
     CALL local_multipole_moment_analysis()
     WRITE(output_FID,*) 'Information for noniterative CM5 method will be printed now.'
     FLUSH(output_FID)
     CALL compute_CM5_NACs()
     WRITE(output_FID,*) 'Hirshfeld and CM5 analysis finished, calculation of iterative AIM will proceed.'
     FLUSH(output_FID)
   END IF   
   WRITE(output_FID,*) 'iter = ', iter
   WRITE(output_FID,*) 'Localized charges for the current iteration: '
   WRITE(output_FID,'(10f14.6)')localized_net_atomic_charge
   WRITE(output_FID,*) 'Net atomic charges for the current iteration: '
   WRITE(output_FID,'(10f14.6)') net_atomic_charge
   WRITE(output_FID,*) 'The updated reference ion charges will be: '
   reference_ion_charge=(2.0_dp*localized_net_atomic_charge+net_atomic_charge)/3.0_dp
   WRITE(output_FID,'(10f14.6)')reference_ion_charge
   CALL update_atomic_densities( )
 END DO
 
 !!!!!!!!!!!!!!!!!!!!Set timing to first 2 partitioning steps
 WRITE(output_FID,*)' '
 WRITE(output_FID,*)'******************************************************'
 CALL system_clock(clock_counta,clock_count_rate,clock_count_max)
 IF (clock_count_rate == 0) THEN
   WRITE(output_FID,*)'Processor clock not available, the time calculation wil not be performed.'
 ELSE
   seconds = modulo((clock_counta-clock_count),clock_count_max)/real(clock_count_rate)
   WRITE(output_FID,*)'Finished first 2 partitioning steps in ',seconds, 'seconds'
   FLUSH(output_FID)
 END IF
 WRITE(output_FID,*)'******************************************************'
 WRITE(output_FID,*)' '
 FLUSH(output_FID)
 
 !Compute the conditioned reference density
 ALLOCATE(sum_conditioned_density(nshells,natoms))
 sum_conditioned_density = 0.0_dp
 ref_pseudodensity = 0.0_dp 
 sum_points = 0 
 sum_conditioned_density = 0.0_dp
!$omp parallel default(none) &
!$omp private(j,na,nb,nc,ka,kb,kc,temp_vector,distance,shell_index) &
!$omp shared(lower_nc,upper_nc,kCpoints,delta_nc,lower_nb,upper_nb,kBpoints,delta_nb,lower_na,upper_na,kApoints,delta_na,&
!$omp boundary,natoms,chunk_size,center_shift,ref_pseudodensity,combined_oxidation_density,total_density) &
!$omp reduction(+:sum_conditioned_density,sum_points)
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
         IF (shell_index <= nshells) THEN
           ref_pseudodensity(ka,kb,kc) = ref_pseudodensity(ka,kb,kc) + combined_oxidation_density(shell_index,j) 
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
         IF ((shell_index <= nshells) .and. (ref_pseudodensity(ka,kb,kc) > zero_tolerance**2)) THEN
           sum_conditioned_density(shell_index,j)=sum_conditioned_density(shell_index,j)+combined_oxidation_density&
           (shell_index,j)*total_density(ka,kb,kc)/ref_pseudodensity(ka,kb,kc) 
           sum_points(shell_index,j) = sum_points(shell_index,j) + 1 
         END IF
       END DO    
     END DO
   END DO
!$omp end do
 END DO
!$omp end parallel
!$omp parallel do default(none)&
!$omp shared(natoms,sum_points,sigma,sum_conditioned_density) &
!$omp private(j,k) &
!$omp schedule(dynamic,1)
 DO j = 1,natoms
   DO k = 1,nshells
     IF (sum_points(k,j) > 0) THEN
       sigma(k,j) = sum_conditioned_density(k,j)/sum_points(k,j) 
     ELSE 
       sigma(k,j) = 0.0_dp
     END IF
   END DO
 END DO
!$omp end parallel do
 conditioned_ref_population = sum(sum_conditioned_density,1)*pixelvolume 
 valence_electrons=0.0_dp 
 DO j=1,natoms
   valence_electrons=valence_electrons+conditioned_ref_population(j)-core_population(j) + occupancy_correction(1,j)
 END DO
 normalization=nvalence/valence_electrons
 WRITE(output_FID,'(a, f14.6)') ' normalization= ',normalization
 FLUSH(output_FID)
!$omp parallel do &
!$omp default (none) &
!$omp private(j,k) &
!$omp shared (natoms,target_conditioned_ref_population,atomic_number,core_electrons,reference_ion_charge,normalization,&
!$omp core_population,occupancy_correction,conditioned_ref_population,sigma) &
!$omp schedule(dynamic,1)
 DO j = 1,natoms
   target_conditioned_ref_population(j) = (atomic_number(j) - core_electrons(j) - reference_ion_charge(j))/normalization + &
   core_population(j) - occupancy_correction(1,j) 
   DO k = 1,nshells
     !Normalize sigma to give the correct number of reference ion electrons
     IF (conditioned_ref_population(j) > charge_convergence_tolerance) THEN
       sigma(k,j) = max(sigma(k,j)*target_conditioned_ref_population(j)/conditioned_ref_population(j),0.0_dp)
     END IF
   END DO
 END DO  
!$omp end parallel do
 DEALLOCATE(ref_pseudodensity)
 DEALLOCATE(localized_pseudodensity)
 ALLOCATE(conditioned_ref_pseudodensity(totnumA,totnumB,totnumC))
 
 !!!!!!!!!!!!!!!!!!!!Set timing conditioned number of reference ion electrons
 WRITE(output_FID,*)' '
 WRITE(output_FID,*)'******************************************************'
 CALL system_clock(clock_countb,clock_count_rate,clock_count_max)
 IF (clock_counta == 0) THEN
   WRITE(output_FID,*)'Processor clock not available, the time calculation wil not be performed.'
 ELSE
   seconds = modulo((clock_countb-clock_counta),clock_count_max)/real(clock_count_rate)
   WRITE(output_FID,*)'Finished conditioned number of reference ion electrons in ',seconds, 'seconds'
   FLUSH(output_FID)
 END IF
 WRITE(output_FID,*)'******************************************************'
 WRITE(output_FID,*)' '
 FLUSH(output_FID)
 
 !Make the conditioned reference density monotonically decreasing with increasing radius
 ALLOCATE(radial_shell_integration_weight(nshells,natoms))
 radial_shell_integration_weight = 0.0_dp
 conditioned_oxidation_density = sigma
 reshaping_iterations = 1
!$omp parallel do default(none) &
!$omp private(j,k,local_reshaping_iterations) &
!$omp shared(natoms,sum_points,sigma,radial_shell_integration_weight,conditioned_oxidation_density,pixelvolume) &
!$omp reduction(max:reshaping_iterations) &
!$omp schedule(dynamic,1)
 DO j=1,natoms
   DO k=1,nshells
     radial_shell_integration_weight(k,j)=sum_points(k,j)*pixelvolume
   END DO
   CALL monotonic_decay_function(sigma(:,j),nshells,radial_shell_integration_weight(:,j),&
   conditioned_oxidation_density(:,j),local_reshaping_iterations)
   reshaping_iterations = max(reshaping_iterations,local_reshaping_iterations)
 END DO
!$omp end parallel do

 !!!!!!!!!!!!!!!!!!!!Set timing conditioned monotonic decay of reference density
 WRITE(output_FID,*)' '
 WRITE(output_FID,*)'******************************************************'
 CALL system_clock(clock_countc,clock_count_rate,clock_count_max)
 IF (clock_countb == 0) THEN
   WRITE(output_FID,*)'Processor clock not available, the time calculation wil not be performed.'
 ELSE
   seconds = modulo((clock_countc-clock_countb),clock_count_max)/real(clock_count_rate)
   WRITE(output_FID,*)'Finished conditioned monotonic decay of reference density in ',seconds, 'seconds'
   FLUSH(output_FID)
 END IF
 WRITE(output_FID,*)'******************************************************'
 WRITE(output_FID,*)' '
 FLUSH(output_FID)

 WRITE(output_FID,*) 'iter = 3'
 WRITE(output_FID,*)'Iterations to converge reshaping: '
 WRITE(output_FID,*) reshaping_iterations
 FLUSH(output_FID)
 !Now that the conditioned_oxidation_density has been computed, compute tau
 conditioned_ref_pseudodensity = 0.0_dp 
 sum_points = 0
 sum_sqrt_Wref = 0.0_dp 
 sum_wA_over_sqrt_Wref = 0.0_dp 
!$omp parallel default(none) &
!$omp private(j,na,nb,nc,ka,kb,kc,temp_vector,distance,shell_index) &
!$omp shared(lower_nc,upper_nc,kCpoints,delta_nc,lower_nb,upper_nb,kBpoints,delta_nb,lower_na,upper_na,kApoints,delta_na,&
!$omp boundary,natoms,chunk_size,center_shift,conditioned_ref_pseudodensity,conditioned_oxidation_density) &
!$omp reduction(+:sum_points,sum_sqrt_Wref,sum_wA_over_sqrt_Wref)
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
         IF (shell_index <= nshells) THEN
           conditioned_ref_pseudodensity(ka,kb,kc) = conditioned_ref_pseudodensity(ka,kb,kc) + conditioned_oxidation_density&
           (shell_index,j) 
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
         IF ((shell_index <= nshells) .and. (conditioned_ref_pseudodensity(ka,kb,kc) > zero_tolerance**2) ) THEN
           sum_points(shell_index,j) = sum_points(shell_index,j) + 1 
           sum_sqrt_Wref(shell_index,j) = sum_sqrt_Wref(shell_index,j) + sqrt(conditioned_ref_pseudodensity(ka,kb,kc)) 
           sum_wA_over_sqrt_Wref(shell_index,j) = sum_wA_over_sqrt_Wref(shell_index,j) + conditioned_oxidation_density&
          (shell_index,j)/sqrt(conditioned_ref_pseudodensity(ka,kb,kc)) 
         END IF
       END DO    
     END DO
   END DO
!$omp end do
 END DO
!$omp end parallel
!$omp parallel do &
!$omp default(none) &
!$omp private(j,k) &
!$omp shared(natoms,sum_points,tau,sum_wA_over_sqrt_Wref,sum_sqrt_Wref)
 DO j = 1,natoms
   DO k = 1,nshells       
     IF (sum_points(k,j) > 0) THEN
       tau(k,j) = sum_wA_over_sqrt_Wref(k,j)/sum_sqrt_Wref(k,j) 
     ELSE
       tau(k,j) = 1.0_dp
     END IF
   END DO
 END DO
!$omp end parallel do
 
 !!!!!!!!!!!!!!!!!!!!Set timing calculation of tau
 WRITE(output_FID,*)' '
 WRITE(output_FID,*)'******************************************************'
 CALL system_clock(clock_countd,clock_count_rate,clock_count_max)
 IF (clock_countc == 0) THEN
   WRITE(output_FID,*)'Processor clock not available, the time calculation wil not be performed.'
 ELSE
   seconds = modulo((clock_countd-clock_countc),clock_count_max)/real(clock_count_rate)
   WRITE(output_FID,*)'Finished tau calculation in ',seconds, 'seconds'
   FLUSH(output_FID)
 END IF
 WRITE(output_FID,*)'******************************************************'
 WRITE(output_FID,*)' '
 FLUSH(output_FID)
 
 
 WRITE(output_FID,*) 'Conditioned reference density computation complete.'
 FLUSH(output_FID)
!$omp parallel do &
!$omp default(none) &
!$omp private(j,k,first_constraint_term,second_constraint_term) &
!$omp shared(natoms,tau,first_exp_const,second_exp_const) &
!$omp schedule(dynamic,1)
 DO j = 1,natoms
   DO k = 1,nshells
      first_constraint_term = minimum_buried_tail_exponent*(1.0_dp - tau(k,j)**2)
      !The decay exponent, convert to per radial shell
      first_exp_const(k,j) = exp(-first_constraint_term/scalefactor)
      second_constraint_term = (1.0_dp-tau(k,j)**2)/maximum_buried_tail_exponent
      !The decay exponent, convert to per radial shell
      IF (second_constraint_term > 0.01_dp) THEN
         second_exp_const(k,j) = exp(-1.0_dp/(scalefactor*second_constraint_term))
      ELSE
         second_exp_const(k,j) = 0.0_dp
      END IF
   END DO
 END DO
!$omp end parallel do
 DEALLOCATE(conditioned_ref_pseudodensity)
 !Initialize the partial_density array
 partial_density = conditioned_oxidation_density
 
  !!!!!!!!!!!!!!!!!!!!Timing charge partitioning step 3
 WRITE(output_FID,*)' '
 WRITE(output_FID,*)'******************************************************'
 CALL system_clock(clock_countd,clock_count_rate,clock_count_max)
 IF (clock_counta == 0) THEN
   WRITE(output_FID,*)'Processor clock not available, the time calculation wil not be performed.'
 ELSE
   seconds = modulo((clock_countd-clock_counta),clock_count_max)/real(clock_count_rate)
   WRITE(output_FID,*)'Finished 3rd partitioning step in ',seconds, 'seconds'
   WRITE(output_FID,*)'This is the number of second to condition the reference densities plus the calculation of tau'
   FLUSH(output_FID)
 END IF
 WRITE(output_FID,*)'******************************************************'
 WRITE(output_FID,*)' '
 FLUSH(output_FID)
 
 
 !Perform the iterative update
 completed_steps = 3
 update_kappa = .false.
 DO iter = 4,500
   total_pseudodensity = 0.0_dp
   sum_points = 0
   sum_density = 0.0_dp
   sum_weighted_points = 0.0_dp
   sum_weighted_density = 0.0_dp
   wA_renormalization_slope = 0.0_dp  
!$omp parallel default(none) &
!$omp private(j,na,nb,nc,ka,kb,kc,temp_vector,distance,shell_index,local_atomic_density,temp) &
!$omp shared(lower_nc,upper_nc,kCpoints,delta_nc,lower_nb,upper_nb,kBpoints,delta_nb,lower_na,upper_na,kApoints,delta_na,boundary,&
!$omp total_pseudodensity,partial_density,natoms,pixelvolume,chunk_size,center_shift,total_density) & 
!$omp reduction(+:sum_density,sum_points,sum_weighted_points,sum_weighted_density)
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
             temp = partial_density(shell_index,j)/total_pseudodensity(ka,kb,kc)
             local_atomic_density =  temp*total_density(ka,kb,kc)
             sum_density(shell_index,j) = sum_density(shell_index,j) + local_atomic_density
             sum_weighted_density(shell_index,j) = sum_weighted_density(shell_index,j) + temp*local_atomic_density
             sum_points(shell_index,j) = sum_points(shell_index,j) + 1
             sum_weighted_points(shell_index,j) = sum_weighted_points(shell_index,j) + temp
           END IF
         END DO    
       END DO
     END DO
!$omp end do
   END DO
!$omp end parallel
   old_weighted_spherical_average_density = weighted_spherical_average_density
   DO j = 1,natoms
     DO k = 1,nshells
       IF (sum_points(k,j) > 0) THEN
         spherical_average_density(k,j) = sum_density(k,j)/sum_points(k,j)
         weighted_spherical_average_density(k,j)=(sum_density(k,j)-sum_weighted_density(k,j)+0.2_dp*sum_weighted_points(k,j)*&
         spherical_average_density(k,j))/(sum_points(k,j)-0.8_dp*sum_weighted_points(k,j))
         wA_renormalization_slope(j) = wA_renormalization_slope(j) + pixelvolume*(sum_density(k,j) - sum_weighted_density(k,j))
       ELSE 
         spherical_average_density(k,j) = 0.0_dp
         weighted_spherical_average_density(k,j) = 0.0_dp
       END IF
     END DO
   END DO
   WRITE(output_FID,*) "iteration= ",iter
   FLUSH(output_FID)
   ! Compute occupations and population changes
   old_net_atomic_charge = net_atomic_charge
   old_change = max_change
   valence_population = sum(sum_density,1)*pixelvolume
   DO j = 1,natoms
     valence_population(j) = valence_population(j) + occupancy_correction(1,j) - core_population(j)
   END DO
   normalization=nvalence/sum(valence_population)
   WRITE(output_FID,'(a,f13.9)')' Normalization: ', normalization
   FLUSH(output_FID)
   valence_population=normalization*valence_population
   DO j = 1,natoms
     net_atomic_charge(j) = atomic_number(j) - core_electrons(j) - valence_population(j)
   END DO
   WRITE(output_FID,*)'Net atomic charges for the current iteration: '
   WRITE(output_FID,'(10f14.6)')net_atomic_charge
   FLUSH(output_FID)
   max_change = maxval(abs(net_atomic_charge - old_net_atomic_charge))
   max_density_change = maxval(abs(old_weighted_spherical_average_density - weighted_spherical_average_density))
   WRITE(output_FID,'(a,es13.6)')' Max change: ',max_change
   WRITE(output_FID,'(a,es13.6)')' Maximum density change: ',max_density_change
   FLUSH(output_FID)
   ! Determine the value of update_kappa
   IF (iter == 4) THEN
     update_kappa = .false.
   ELSE
     DO j = 1,natoms
       IF (valence_population(j) < -charge_convergence_tolerance) THEN
          update_kappa = .true.
       END IF
     END DO
     IF ((update_kappa) .and. (max(max_change,old_change) < charge_convergence_tolerance)) THEN
        update_kappa = .false.
        wA_renormalization = 1.0_dp
     END IF
   END IF
   IF (update_kappa) THEN
     ! Update the wA_renormalization values
     ! The following constraint makes sure the number of valence electrons assigned to each atom is non-negative
     DO j = 1,natoms
       DO k = 1,nshells
         partial_density(k,j) = partial_density(k,j)/wA_renormalization(j)
       END DO   
       IF (wA_renormalization_slope(j) > nonoverlapping_atom_tolerance) THEN
         wA_renormalization(j) = max(1.0_dp, (wA_renormalization(j)*exp(-valence_population(j)/&
         wA_renormalization_slope(j))) )
       ELSE
         wA_renormalization(j) = 1.0_dp
       END IF
       DO k = 1,nshells
         partial_density(k,j) = partial_density(k,j)*wA_renormalization(j)
       END DO
     END DO      
   ELSE
     completed_steps = completed_steps + 1 
     IF  (completed_steps == 7) THEN
        EXIT
     END IF
     partial_density=weighted_spherical_average_density
     !Reshape the atomic weighting factors to not be too diffuse or too contracted
     reshaping_iterations = 1
!$omp parallel do default(none) &
!$omp private(j,k,local_reshaping_iterations) &
!$omp shared(natoms,radial_shell_integration_weight,sum_points,weighted_spherical_average_density,&
!$omp first_exp_const,second_exp_const,partial_density,pixelvolume) &
!$omp reduction(max:reshaping_iterations) &
!$omp schedule(dynamic,1)
     DO j=1,natoms
       DO k=1,nshells
         radial_shell_integration_weight(k,j) = sum_points(k,j)*pixelvolume
       END DO
       CALL tail_exponential_decay_function(weighted_spherical_average_density(:,j),first_exp_const(:,j),second_exp_const(:,j),&
       nshells,radial_shell_integration_weight(:,j),partial_density(:,j),local_reshaping_iterations)
       reshaping_iterations = max(reshaping_iterations,local_reshaping_iterations)
     END DO
!$omp end parallel do
     WRITE(output_FID,*)'Iterations to converge reshaping: '
     FLUSH(output_FID)
     WRITE(output_FID,*) reshaping_iterations
     FLUSH(output_FID)
   END IF
 END DO
 IF (iter >= 501) THEN
   WRITE(output_FID,*) 'Sorry, the calculation of atomic charges did not converge within 500 iterations.'
   FLUSH(output_FID)
 END IF

 !!!!!!!!!!!!!!!!!!!!Timing charge partitioning steps 4-7
 WRITE(output_FID,*)' '
 WRITE(output_FID,*)'******************************************************'
 CALL system_clock(clock_counte,clock_count_rate,clock_count_max)
 IF (clock_countd == 0) THEN
   WRITE(output_FID,*)'Processor clock not available, the time calculation wil not be performed.'
 ELSE
   seconds = modulo((clock_counte-clock_countd),clock_count_max)/real(clock_count_rate)
   WRITE(output_FID,*)'Finished charge partitioning steps 4-7 in ',seconds, 'seconds'
   FLUSH(output_FID)
 END IF
 WRITE(output_FID,*)'******************************************************'
 WRITE(output_FID,*)' '
 FLUSH(output_FID)

 
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
 
 END SUBROUTINE DDEC6_valence_iterator

 END MODULE module_DDEC6_valence_iterator
  