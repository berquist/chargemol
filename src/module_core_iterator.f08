 MODULE module_core_iterator
 ! Copyright (c) 2014, 2015, 2016, 2017 by Thomas A. Manz and Nidia Gabaldon Limas. Rights reserved.

 USE module_precision
 USE module_global_parameters
 USE module_common_variable_declarations

 IMPLICIT NONE
 
 CONTAINS
 
 SUBROUTINE core_iterator()
 !===================================================================================
 
 REAL(kind=dp), ALLOCATABLE, DIMENSION(:) :: old,nA_correction,sum_cubes
 REAL(kind=dp) :: exp_const,distance,temp_vector(3),change,temp,local_core_density,max_nA_correction,temp_vector1(3)
 INTEGER :: iter,k,nthread
 INTEGER :: clock_count,clock_count_rate,clock_count_max,clock_count2
 REAL(kind=dp), ALLOCATABLE, DIMENSION(:,:,:) :: temp_core_density
 LOGICAL :: some_value
 
 WRITE(output_FID,*)' '
 WRITE(output_FID,*)'******************** TIME CONTROL ********************'
 WRITE(output_FID,*)'Starting core_iterator'
 FLUSH(output_FID)
 CALL system_clock(clock_count,clock_count_rate,clock_count_max)
 WRITE(output_FID,*)'******************************************************'
 WRITE(output_FID,*)' '
 FLUSH(output_FID) 
 
 !Allocate some arrays
 ALLOCATE(spherical_avg_core_density(nshells,natoms))
 spherical_avg_core_density = 0.0_dp
 ALLOCATE(sum_core_density(nshells,natoms))
 ALLOCATE(core_pseudodensity(totnumA,totnumB,totnumC))
 ALLOCATE(old(natoms))

 temp=max_core_electrons_per_pixel/pixelvolume
   DO j=1,natoms
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
               core_density(ka,kb,kc) = min(core_density(ka,kb,kc),temp)
             END IF
           END DO
         END DO
     END DO
   END DO
 
 !The decay exponent is 2 per bohr, convert to per radial shell
 exp_const = exp(-2.0_dp*cutoff_radius*bohrperangstrom/(100.0_dp*nshells))
 WRITE(output_FID,*)'Iteratively solving for the core charge distributions:'
 FLUSH(output_FID)
 
 DO iter = 1,200
   core_pseudodensity=0.0_dp
   sum_core_density = 0.0_dp 
   sum_points = 0
!$omp parallel default(none) &
!$omp private(j,na,nb,nc,ka,kb,kc,temp_vector,distance,shell_index) &
!$omp shared(natoms,core_electrons,lower_nc,upper_nc,kCpoints,delta_nc,lower_nb,upper_nb,kBpoints,lower_na,&
!$omp upper_na,kApoints,delta_na,boundary,center_shift,core_pseudodensity,core_density,delta_nb,partial_core_density,chunk_size) &
!$omp reduction(+:sum_core_density,sum_points)
   DO j=1,natoms
!$omp do schedule(dynamic,chunk_size)
     DO nc = lower_nc(j),upper_nc(j)
       IF (core_electrons(j) > zero_tolerance) THEN
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
               core_pseudodensity(ka,kb,kc) = core_pseudodensity(ka,kb,kc) + partial_core_density(shell_index,j)
             END IF
           END DO
         END DO
       END IF
     END DO
!$omp end do
   END DO
   DO j=1,natoms 
!$omp do schedule(dynamic,chunk_size)
     DO nc = lower_nc(j),upper_nc(j)
       IF (core_electrons(j) > zero_tolerance) THEN
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
             IF ((shell_index <= nshells) .and. (core_pseudodensity(ka,kb,kc) > zero_tolerance)) THEN
               sum_core_density(shell_index,j) = sum_core_density(shell_index,j) + &
               partial_core_density(shell_index,j)*core_density(ka,kb,kc)/core_pseudodensity(ka,kb,kc)
               sum_points(shell_index,j) = sum_points(shell_index,j) + 1
             END IF
           END DO    
         END DO   
       END IF
     END DO
!$omp end do
   END DO
!$omp end parallel
   DO j = 1,natoms
     DO k = 1,nshells
       IF (sum_points(k,j) > 0) THEN
         spherical_avg_core_density(k,j) = sum_core_density(k,j)/sum_points(k,j)
       ELSE 
         spherical_avg_core_density(k,j) = 0.0_dp
       END IF
     END DO
   END DO
   WRITE(output_FID,'(a,i3)')' Iteration= ',iter
   FLUSH(output_FID)
   old=core_population
   core_population=pixelvolume*(sum(sum_core_density,1))
   !Compute core occupations
   change = 0.0_dp
   DO j=1,natoms
     change = max(change,abs(core_population(j)-old(j)))
   END DO
   WRITE(output_FID,'(a,es14.6)')' change= ', change
   FLUSH(output_FID)
   IF ((change < charge_convergence_tolerance) .and. (iter > 5)) THEN
     EXIT
   END IF
   !Update atomic core density factors
   DO j = 1,natoms
     partial_core_density(1,j) = spherical_avg_core_density(1,j)
     temp = partial_core_density(1,j)
     DO k = 2,nshells
       IF ((sum_points(k,j) > 0) .and. (temp > 0.0_dp)) THEN
         partial_core_density(k,j) = max(min(spherical_avg_core_density(k,j),temp*exp_const), 0.0_dp)
       ELSE
         partial_core_density(k,j) = max(spherical_avg_core_density(k,j), 0.0_dp)
       END IF       
       temp = partial_core_density(k,j)   
     END DO
   END DO
 END DO
 !Correct the core density grid so it gives the proper number of core electrons for each atom
 WRITE(output_FID,*)'Correcting the core density grid'
 FLUSH(output_FID)
 !Initialize arrays
 ALLOCATE(nA_correction(natoms))
 ALLOCATE(sum_cubes(natoms))
 ALLOCATE(K_factor(natoms))
 ALLOCATE(max_density(natoms))
 ALLOCATE(temp_core_density(totnumA,totnumB,totnumC))
 !Compute the initial correction size 
 DO j=1,natoms
   nA_correction(j) = core_electrons(j) - core_population(j)
 END DO
 !Perform the correction
 DO iter = 1,40
   K_factor = 0.0_dp
   max_density = 0.0_dp
   sum_cubes = 0.0_dp
   temp_core_density = 0.0_dp
   local_core_density = 0.0_dp
!$omp parallel default(none) &
!$omp private(j,na,nb,nc,ka,kb,kc,temp_vector,distance,shell_index,local_core_density) &
!$omp shared(natoms,core_electrons,lower_nc,upper_nc,kCpoints,delta_nc,lower_nb,upper_nb,kBpoints,delta_nb,lower_na,upper_na,&
!$omp kApoints,delta_na,boundary,center_shift,core_pseudodensity,chunk_size,partial_core_density,core_density) &
!$omp reduction(+:sum_cubes) &
!$omp reduction(max:max_density)
   DO j=1,natoms
!$omp do schedule(dynamic,chunk_size)
     DO nc = lower_nc(j),upper_nc(j)
       IF (core_electrons(j) > zero_tolerance) THEN
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
             IF ((shell_index <= nshells) .and. (core_pseudodensity(ka,kb,kc) > zero_tolerance)) THEN
               local_core_density = partial_core_density(shell_index,j)*core_density(ka,kb,kc)/core_pseudodensity(ka,kb,kc)
               sum_cubes(j) = sum_cubes(j) + local_core_density**3
               max_density(j) = max(max_density(j),local_core_density)
             END IF                    
           END DO
         END DO
       END IF
     END DO 
!$omp end do
   END DO
!$omp end parallel
   DO j=1,natoms
     !Compute the adjustment constant
     IF ((sum_cubes(j) > zero_tolerance**3) .and. (max_density(j) > zero_tolerance)) THEN
       K_factor(j) = nA_correction(j)/(sum_cubes(j)*pixelvolume)
       K_factor(j) = min(K_factor(j), 0.25_dp/(max_density(j)**2))
     END IF
   END DO
!$omp parallel default(none) &
!$omp private(j,na,nb,nc,ka,kb,kc,temp_vector,distance,shell_index,local_core_density) &
!$omp shared(natoms,core_electrons,lower_nc,upper_nc,kCpoints,delta_nc,lower_nb,upper_nb,kBpoints,delta_nb,lower_na,upper_na,&
!$omp kApoints,delta_na,boundary,center_shift,core_pseudodensity,chunk_size,partial_core_density,core_density,&
!$omp K_factor,temp_core_density)
   DO j=1,natoms
!$omp do schedule(dynamic,chunk_size)
     DO nc = lower_nc(j),upper_nc(j)
       IF (core_electrons(j) > zero_tolerance) THEN
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
             IF ((shell_index <= nshells) .and. (core_pseudodensity(ka,kb,kc) > zero_tolerance)) THEN
               local_core_density = partial_core_density(shell_index,j)*core_density(ka,kb,kc)/core_pseudodensity(ka,kb,kc)
               IF (local_core_density > zero_tolerance**2) THEN
                 temp_core_density(ka,kb,kc) = temp_core_density(ka,kb,kc) + local_core_density/sqrt(1 - 2*K_factor(j)*&
                 local_core_density**2)
               END IF   
             END IF
           END DO
         END DO
       END IF
     END DO
!$omp end do
   END DO
!$omp end parallel
   core_density = temp_core_density
   sum_core_density = 0.0_dp
   sum_points = 0
!$omp parallel default(none) &
!$omp private(j,na,nb,nc,ka,kb,kc,temp_vector,distance,shell_index) &
!$omp shared(natoms,core_electrons,lower_nc,upper_nc,kCpoints,delta_nc,lower_nb,upper_nb,kBpoints,delta_nb,lower_na,upper_na,&
!$omp kApoints,delta_na,boundary,center_shift,core_pseudodensity,partial_core_density,chunk_size,core_density) &
!$omp reduction(+:sum_core_density,sum_points)
   DO j=1,natoms
!$omp do schedule(dynamic,chunk_size)
     DO nc = lower_nc(j),upper_nc(j)
       IF (core_electrons(j) > zero_tolerance) THEN
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
             IF ((shell_index <= nshells) .and. (core_pseudodensity(ka,kb,kc) > zero_tolerance)) THEN
               sum_core_density(shell_index,j) = sum_core_density(shell_index,j) + &
               partial_core_density(shell_index,j)*core_density(ka,kb,kc)/core_pseudodensity(ka,kb,kc)
               sum_points(shell_index,j) = sum_points(shell_index,j) + 1
             END IF
           END DO    
         END DO
       END IF
     END DO
!$omp end do
   END DO
!$omp end parallel
   DO j = 1,natoms
     DO k = 1,nshells
       IF (sum_points(k,j) > 0.0_dp) THEN
         spherical_avg_core_density(k,j) = sum_core_density(k,j)/sum_points(k,j)
       ELSE 
         spherical_avg_core_density(k,j) = 0.0_dp
       END IF
     END DO
   END DO
   WRITE(output_FID,*)'iteration= ',iter
   FLUSH(output_fid)
   core_population=pixelvolume*(sum(sum_core_density,1))
   !Compute the correction size
   max_nA_correction = 0.0_dp
   DO j=1,natoms
     nA_correction(j) = core_electrons(j) - core_population(j)
     max_nA_correction = max(max_nA_correction, abs(nA_correction(j)))
   END DO
   !Check for convergence
   IF (max_nA_correction < charge_convergence_tolerance) THEN
     WRITE(output_FID,*)'Core density grid corrected in the following number of iterations: '
     FLUSH(output_FID)
     WRITE(output_FID,*)'iter= ',iter
     FLUSH(output_FID)
     EXIT
   END IF  
   !Update atomic core density factors
   DO j = 1,natoms
     partial_core_density(1,j) = spherical_avg_core_density(1,j)
     temp = partial_core_density(1,j)
     DO k = 2,nshells
       IF ((sum_points(k,j) > 0) .and. (temp > 0.0_dp)) THEN
         partial_core_density(k,j) = max(min(spherical_avg_core_density(k,j),temp*exp_const), 0.0_dp)
       ELSE
         partial_core_density(k,j) = max(spherical_avg_core_density(k,j), 0.0_dp)
       END IF    
       temp = partial_core_density(k,j)
     END DO
   END DO
   core_pseudodensity = 0.0_dp
!$omp parallel default(none) &
!$omp private(na,nb,nc,ka,kb,kc,temp_vector,distance,shell_index,j) &
!$omp shared(lower_nc,upper_nc,kCpoints,delta_nc,lower_nb,upper_nb,kBpoints,delta_nb,lower_na,upper_na,kApoints,delta_na,boundary,&
!$omp center_shift,core_pseudodensity,partial_core_density,chunk_size,natoms,core_electrons) 
   DO j = 1,natoms
!$omp do schedule(dynamic,chunk_size)
     DO nc = lower_nc(j),upper_nc(j)
       IF (core_electrons(j) > zero_tolerance) THEN
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
               core_pseudodensity(ka,kb,kc) = core_pseudodensity(ka,kb,kc) + partial_core_density(shell_index,j)
             END IF
           END DO 
         END DO
       END IF
     END DO
!$omp end do
   END DO
!$omp end parallel

 END DO 

 DEALLOCATE(temp_core_density)
 
 WRITE(output_FID,*)' '
 WRITE(output_FID,*)'******************** TIME CONTROL ********************'
 CALL system_clock(clock_count2,clock_count_rate,clock_count_max)
 IF (clock_count_rate == 0) THEN
   WRITE(output_FID,*)'Processor clock not available, the time calculation wil not be performed.'
 ELSE
   seconds = modulo((clock_count2-clock_count),clock_count_max)/real(clock_count_rate)
   WRITE(output_FID,*)'Finished core_iterator in ',seconds, 'seconds'
   FLUSH(output_FID)
 END IF
 WRITE(output_FID,*)'******************************************************'
 WRITE(output_FID,*)' '
 FLUSH(output_FID)

 END SUBROUTINE core_iterator

 END MODULE module_core_iterator
 