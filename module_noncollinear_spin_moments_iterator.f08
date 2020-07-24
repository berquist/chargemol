 MODULE module_noncollinear_spin_moments_iterator
 ! Copyright (c) 2014, 2015, 2016, 2017 by Thomas A. Manz and Nidia Gabaldon Limas. Rights reserved.
 
 USE module_precision
 USE module_common_variable_declarations
 USE module_global_parameters
 USE module_spin_functions
 USE module_generate_atomic_spin_moment_file

 IMPLICIT NONE
 
 CONTAINS
 
 SUBROUTINE noncollinear_spin_moments_iterator( )
 !===================================================================================
 
 INTEGER :: k,comp
 REAL(kind=dp) :: a,b(3),local_atomic_spin_vector(3),proportional_spin_vector(3),L_vector(3),&
 temp_vector(3),distance,local_atomic_density,proportional_theta_vector(3),theta_vector(3),mag_L_vector,&
 mag_local_atomic_spin_vector,inv_Xi,b_vector(3)
 REAL(kind=dp), ALLOCATABLE, DIMENSION(:,:) :: old
 REAL(kind=dp), ALLOCATABLE, DIMENSION(:,:,:) :: sum_spin_density_vector
 REAL(kind=dp), ALLOCATABLE, DIMENSION(:,:,:,:) :: sum_Ypsilon_vector,sum_spin_pseudodensity_vector,&
 Ypsilon_vector,Omega_vector,trial_spin_density_vector
 INTEGER :: clock_count,clock_count_rate,clock_count_max,clock_count2
 
 WRITE(output_FID,*)' '
 WRITE(output_FID,*)'******************** TIME CONTROL ********************'
 WRITE(output_FID,*)'Starting noncollinear_spin_moments_iterator'
 FLUSH(output_FID)
 CALL system_clock(clock_count,clock_count_rate,clock_count_max)
 WRITE(output_FID,*)'******************************************************'
 WRITE(output_FID,*)' '
 FLUSH(output_FID)
 
 !Initialize some arrays
 ALLOCATE(spherical_average_atomic_spin_vector(3,nshells,natoms))
 ALLOCATE(Ypsilon_vector(3,totnumA,totnumB,totnumC))
 ALLOCATE(Omega_vector(3,totnumA,totnumB,totnumC))
 ALLOCATE(spin_population_vector(3,natoms))
 ALLOCATE(old(3,natoms))
 ALLOCATE(trial_spin_density_vector(3,totnumA,totnumB,totnumC))
 ALLOCATE(sum_spin_density_vector(3,nshells,natoms))
 ALLOCATE(sum_Ypsilon_vector(3,totnumA,totnumB,totnumC))
 ALLOCATE(sum_spin_pseudodensity_vector(3,totnumA,totnumB,totnumC))
 a = 0.0_dp
 b_vector = 0.0_dp
 local_atomic_spin_vector = 0.0_dp
 proportional_spin_vector = 0.0_dp
 L_vector = 0.0_dp
 Ypsilon_vector = 0.0_dp
 Omega_vector = 0.0_dp
 spin_population_vector = 0.0_dp
 trial_spin_density_vector = 0.0_dp
 spherical_average_atomic_spin_vector = 0.0_dp

 ! Iterative solution for the atomic spin distributions
 WRITE(output_FID,*) 'Iteratively solving for the atomic spin distributions:'
 DO spin_iter = 1,500
   sum_points = 0
   sum_spin_density_vector = 0.0_dp
   sum_Ypsilon_vector = 0.0_dp
   sum_spin_pseudodensity_vector = 0.0_dp
!$omp parallel default(none) &
!$omp private(j,nc,kc,nb,kb,na,ka,temp_vector,distance,shell_index,local_atomic_density,proportional_spin_vector,&
!$omp proportional_theta_vector,local_atomic_spin_vector,L_vector,a,b_vector,theta_vector,mag_L_vector,&
!$omp mag_local_atomic_spin_vector,inv_Xi) &
!$omp shared(natoms,lower_nc,upper_nc,kCpoints,delta_nc,lower_nb,upper_nb,kBpoints,delta_nb,lower_na,upper_na,kApoints,&
!$omp delta_na,boundary,center_shift,total_pseudodensity,partial_density,spin_density_vector,Xi_lookup,total_density,&
!$omp spin_iter,Ypsilon_vector,Omega_vector,trial_spin_density_vector,inverse_Xi_lookup,sum_Ypsilon_vector,&
!$omp sum_spin_pseudodensity_vector,spherical_average_density,spherical_average_atomic_spin_vector,chunk_size) &
!$omp reduction(+:sum_points,sum_spin_density_vector)
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
           IF ((shell_index <= nshells) .and. (total_pseudodensity(ka,kb,kc) > zero_tolerance) .and. ((total_density(ka,kb,kc)) &
           > zero_tolerance)) THEN
             sum_points(shell_index,j) = sum_points(shell_index,j) + 1
             local_atomic_density =  partial_density(shell_index,j)*(total_density(ka,kb,kc))/total_pseudodensity(ka,kb,kc)
             DO comp=1,3 
               proportional_spin_vector(comp) = partial_density(shell_index,j)*spin_density_vector(comp,ka,kb,kc)/&
               total_pseudodensity(ka,kb,kc)
             END DO                        
             proportional_theta_vector = calculate_theta_vector(local_atomic_density,proportional_spin_vector,Xi_lookup)
             IF (spin_iter == 1) THEN
               local_atomic_spin_vector = proportional_spin_vector
               L_vector = proportional_theta_vector
             ELSE
               a = spherical_average_density(shell_index,j)
               DO comp=1,3 
                 b_vector(comp) = spherical_average_atomic_spin_vector(comp,shell_index,j)
               END DO
               theta_vector = calculate_theta_vector(a,b_vector,Xi_lookup)
               DO comp=1,3 
                 L_vector(comp) = Ypsilon_vector(comp,ka,kb,kc) - Omega_vector(comp,ka,kb,kc) + pi*(spin_density_vector&
                 (comp,ka,kb,kc) - trial_spin_density_vector(comp,ka,kb,kc))/(total_density(ka,kb,kc)) + &
                 (1.0_dp-spin_ref_fraction)*theta_vector(comp) + spin_ref_fraction*proportional_theta_vector(comp)
               END DO
               mag_L_vector = max(sqrt(dot_product(L_vector,L_vector)),zero_tolerance**2)
               IF (mag_L_vector > pi) THEN
                 mag_local_atomic_spin_vector = max(local_atomic_density,zero_tolerance**2)
               ELSE
                 inv_Xi = fast_calculate_inverse_Xi(2.0_dp*mag_L_vector,inverse_Xi_lookup)
                 mag_local_atomic_spin_vector = max(local_atomic_density*inv_Xi,zero_tolerance**2)
               END IF
                 DO comp=1,3
                 local_atomic_spin_vector(comp) = mag_local_atomic_spin_vector*L_vector(comp)/mag_L_vector
                 END DO
             END IF 
             DO comp=1,3 
               sum_spin_density_vector(comp,shell_index,j) = sum_spin_density_vector(comp,shell_index,j) + &
               local_atomic_spin_vector(comp)
               sum_Ypsilon_vector(comp,ka,kb,kc) = sum_Ypsilon_vector(comp,ka,kb,kc) + L_vector(comp)*partial_density&
               (shell_index,j)/total_pseudodensity(ka,kb,kc)
               sum_spin_pseudodensity_vector(comp,ka,kb,kc) = sum_spin_pseudodensity_vector(comp,ka,kb,kc)&
			   + local_atomic_spin_vector(comp)
             END DO
           END IF
         END DO    
       END DO
     END DO
!$omp end do
   END DO
!$omp end parallel
   Ypsilon_vector = sum_Ypsilon_vector
   trial_spin_density_vector = sum_spin_pseudodensity_vector
   DO j = 1,natoms
     DO k = 1,nshells
       IF (sum_points(k,j) > 0) THEN
         spherical_average_atomic_spin_vector(:,k,j) = sum_spin_density_vector(:,k,j)/sum_points(k,j)
       ELSE 
         spherical_average_atomic_spin_vector(1,k,j) = 0.0_dp
         spherical_average_atomic_spin_vector(2,k,j) = 0.0_dp
         spherical_average_atomic_spin_vector(3,k,j) = 0.0_dp
       END IF
     END DO
   END DO
   WRITE(output_FID,*) 'iteration: ',spin_iter
   FLUSH(output_FID)
   old=spin_population_vector
   spin_population_vector = 0.0_dp
   DO j = 1,natoms
     DO k = 1,nshells
       DO comp=1,3
         spin_population_vector(comp,j) = spin_population_vector(comp,j) + sum_spin_density_vector(comp,k,j)*pixelvolume
       END DO
     END DO
   END DO 
   tot_spin_moment_vector = 0.0_dp
   DO j = 1,natoms
     DO comp=1,3
       tot_spin_moment_vector(comp) = tot_spin_moment_vector(comp) + spin_population_vector(comp,j)
     END DO
   END DO
   WRITE(output_FID,*) 'spin_population_vector: '
   FLUSH(output_FID)
   DO i=1,natoms
     WRITE(output_FID,'(3f12.6)') (spin_population_vector(j,i),j=1,3)
   END DO
   FLUSH(output_FID)
   WRITE(output_FID,*) 'tot_spin_moment_vector: '
   WRITE(output_FID,*) tot_spin_moment_vector
   FLUSH(output_FID)
   !Compute occupations
   max_change=0.0_dp
   DO j=1,natoms
     DO comp = 1,3
       max_change = max(max_change,abs(spin_population_vector(comp,j)-old(comp,j)))
     END DO
   END DO
   WRITE(output_FID,*)'max_change: ',max_change
   WRITE(output_FID,*)' '
   FLUSH(output_FID)
   IF ((max_change < spin_convergence_tolerance) .and. (spin_iter > 6)) THEN
     CALL generate_atomic_spin_moment_file()
     EXIT
   END IF
   IF (spin_iter == 1) THEN
     WRITE(output_FID,*)'Initial spin populations based on proportional partitioning of the spin density are:'
     FLUSH(output_FID)
     DO i=1,natoms
       WRITE(output_FID,'(3f12.6)') (spin_population_vector(j,i),j=1,3)
     END DO
     FLUSH(output_FID)
     WRITE(output_FID,*) 'The total spin moment of the unit cell is:'
     WRITE(output_FID,*) tot_spin_moment_vector
     FLUSH(output_FID)
   END IF
   Omega_vector = 0.0_dp
!$omp parallel default(none) &
!$omp private(nc,kc,nb,kb,na,ka,temp_vector,distance,shell_index,local_atomic_density,proportional_spin_vector,a,&
!$omp b_vector,theta_vector,proportional_theta_vector,j,comp) &
!$omp shared(lower_nc,upper_nc,lower_nb,upper_nb,kCpoints,kBpoints,lower_na,upper_na,kApoints,delta_nc,delta_nb,delta_na,&
!$omp boundary,center_shift,total_pseudodensity,partial_density,spin_density_vector,total_density,&
!$omp spherical_average_density,spherical_average_atomic_spin_vector,Xi_lookup,Omega_vector,chunk_size,natoms)
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
           IF ((shell_index <= nshells) .and. (total_pseudodensity(ka,kb,kc) > zero_tolerance) .and. ((total_density(ka,kb,kc)) &
           > zero_tolerance)) THEN
             local_atomic_density =  partial_density(shell_index,j)*(total_density(ka,kb,kc))/total_pseudodensity(ka,kb,kc)
             DO comp=1,3
               proportional_spin_vector(comp) = partial_density(shell_index,j)*spin_density_vector(comp,ka,kb,kc)&
               /total_pseudodensity(ka,kb,kc)
             END DO                       
             a = spherical_average_density(shell_index,j)
             DO comp=1,3 
               b_vector(comp) = spherical_average_atomic_spin_vector(comp,shell_index,j)
             END DO
             theta_vector = calculate_theta_vector(a,b_vector,Xi_lookup)
             proportional_theta_vector = calculate_theta_vector(local_atomic_density,proportional_spin_vector,Xi_lookup)
             DO comp=1,3
               Omega_vector(comp,ka,kb,kc) = Omega_vector(comp,ka,kb,kc) + ((1.0_dp-spin_ref_fraction)&
               *theta_vector(comp) + spin_ref_fraction*proportional_theta_vector(comp))*partial_density(shell_index,j)&
			   /total_pseudodensity(ka,kb,kc)
             END DO
           END IF
         END DO    
       END DO
     END DO
!$omp end do
   END DO
!$omp end parallel
 END DO
 WRITE(output_FID,*) 'Final spin populations: '
 DO i=1,natoms
   WRITE(output_FID,'(3f12.6)') (spin_population_vector(j,i),j=1,3)
 END DO
 WRITE(output_FID,*) 'The total spin moment of the unit cell is '
 WRITE(output_FID,*) 'tot_spin_moment_vector: '
 WRITE(output_FID,*) tot_spin_moment_vector
 FLUSH(output_FID)
 !Free some memory by clearing unneeded arrays
 DEALLOCATE(Ypsilon_vector)
 DEALLOCATE(Omega_vector)
 DEALLOCATE(sum_Ypsilon_vector)
 DEALLOCATE(trial_spin_density_vector)
 DEALLOCATE(sum_spin_pseudodensity_vector)

 WRITE(output_FID,*)' '
 WRITE(output_FID,*)'******************** TIME CONTROL ********************'
 CALL system_clock(clock_count2,clock_count_rate,clock_count_max)
 IF (clock_count_rate == 0) THEN
   WRITE(output_FID,*)'Processor clock not available, the time calculation wil not be performed.'
 ELSE
   seconds = modulo((clock_count2-clock_count),clock_count_max)/real(clock_count_rate)
   WRITE(output_FID,*)'Finished noncollinear_spin_moments_iterator in ',seconds, 'seconds'
   FLUSH(output_FID)
 END IF
 WRITE(output_FID,*)'******************************************************'
 WRITE(output_FID,*)' '
 FLUSH(output_FID)  
 
 END SUBROUTINE noncollinear_spin_moments_iterator
 
 END MODULE module_noncollinear_spin_moments_iterator