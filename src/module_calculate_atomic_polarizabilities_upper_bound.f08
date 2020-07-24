 MODULE module_calculate_atomic_polarizabilities_upper_bound
 !======================================================================================
 ! Module that find the maximum polarizability value corresponging to a perfect conductor
 ! Copyright (c) 2014, 2015, 2016, 2017 by Thomas A. Manz and Nidia Gabaldon Limas. Rights reserved.
 !======================================================================================
  
 USE module_precision
 USE module_global_parameters
 USE module_common_variable_declarations
   
 IMPLICIT NONE
 
 CONTAINS
 
 SUBROUTINE calculate_atomic_polarizabilities_upper_bound()
 !====================================================================================
 ! Compute the positions of the charge centers in terms of na, nb, nc
 !====================================================================================
 
 INTEGER :: j
 REAL(kind=dp) :: wA_polarizability_volume,temp_vector(3),distance
 REAL(kind=dp),ALLOCATABLE,DIMENSION(:,:,:) :: W_polarizability_volume

 !Compute the atomic polarizability maximum value, which corresponds to perfect conduction
 ALLOCATE(W_polarizability_volume(totnumA,totnumB,totnumC))
 W_polarizability_volume = 0.0_dp
 wA_polarizability_volume = 0.0_dp
 ALLOCATE(atomic_polarizability_upper_bound(natoms))
 atomic_polarizability_upper_bound = 0.0_dp
!$omp parallel default(none) &
!$omp private(j,na,nb,nc,ka,kb,kc,temp_vector,distance,shell_index,wA_polarizability_volume) &
!$omp shared(lower_nc,upper_nc,kCpoints,delta_nc,lower_nb,upper_nb,kBpoints,delta_nb,lower_na,upper_na,kApoints,delta_na, &
!$omp natoms,boundary,chunk_size,center_shift,partial_density,atomic_number,W_polarizability_volume, &
!$omp CCSD_free_atom_polarizability_to_volume_ratio) &
!$omp reduction(+:atomic_polarizability_upper_bound)
 DO j=1,natoms
!$omp do schedule(dynamic,chunk_size) 
   DO na = lower_na(j),upper_na(j)
     ka = kApoints(j,delta_na + na + 1)
     DO nb = lower_nb(j),upper_nb(j)
       kb = kBpoints(j,delta_nb + nb + 1)
       DO nc = lower_nc(j),upper_nc(j)
         kc = kCpoints(j,delta_nc + nc + 1)
         temp_vector(1) = na*boundary(1,1) + nb*boundary(2,1) + nc*boundary(3,1) - center_shift(1,j)
         temp_vector(2) = na*boundary(1,2) + nb*boundary(2,2) + nc*boundary(3,2) - center_shift(2,j)
         temp_vector(3) = na*boundary(1,3) + nb*boundary(2,3) + nc*boundary(3,3) - center_shift(3,j)
         distance = sqrt(temp_vector(1)*temp_vector(1) + temp_vector(2)*temp_vector(2) + temp_vector(3)*temp_vector(3))
         shell_index = ceiling(scalefactor*distance + zero_tolerance)
         IF (shell_index <= nshells) THEN
           wA_polarizability_volume =  partial_density(shell_index,j)*distance**3*&
           CCSD_free_atom_polarizability_to_volume_ratio(atomic_number(j))
           W_polarizability_volume(ka,kb,kc) = W_polarizability_volume(ka,kb,kc) + wA_polarizability_volume
         END IF            
       END DO
     END DO
   END DO
!$omp end do
 END DO
 DO j=1,natoms
!$omp do schedule(dynamic,chunk_size) 
   DO na = lower_na(j),upper_na(j)
     ka = kApoints(j,delta_na + na + 1)
     DO nb = lower_nb(j),upper_nb(j)
       kb = kBpoints(j,delta_nb + nb + 1)
       DO nc = lower_nc(j),upper_nc(j)
         kc = kCpoints(j,delta_nc + nc + 1)
         temp_vector(1) = na*boundary(1,1) + nb*boundary(2,1) + nc*boundary(3,1) - center_shift(1,j)
         temp_vector(2) = na*boundary(1,2) + nb*boundary(2,2) + nc*boundary(3,2) - center_shift(2,j)
         temp_vector(3) = na*boundary(1,3) + nb*boundary(2,3) + nc*boundary(3,3) - center_shift(3,j)
         distance = sqrt(temp_vector(1)*temp_vector(1) + temp_vector(2)*temp_vector(2) + temp_vector(3)*temp_vector(3))
         shell_index = ceiling(scalefactor*distance + zero_tolerance)
         IF ((shell_index <= nshells) .and. (W_polarizability_volume(ka,kb,kc) > zero_tolerance)) THEN
           wA_polarizability_volume =  partial_density(shell_index,j)*distance**3*CCSD_free_atom_polarizability_to_volume_ratio&
           (atomic_number(j))
           atomic_polarizability_upper_bound(j) = atomic_polarizability_upper_bound(j) + wA_polarizability_volume/&
           W_polarizability_volume(ka,kb,kc)
         END IF
       END DO
     END DO
   END DO
!$omp end do
 END DO
!$omp end parallel
 !Now divide by 2pi which corresponds to the perfect conduction limit 
 atomic_polarizability_upper_bound = pixelvolume*atomic_polarizability_upper_bound/(2.0_dp*pi)

 DEALLOCATE(W_polarizability_volume)
 
 END SUBROUTINE calculate_atomic_polarizabilities_upper_bound
 
 END MODULE module_calculate_atomic_polarizabilities_upper_bound