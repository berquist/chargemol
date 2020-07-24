 MODULE module_integrate_bonding_terms
 ! Copyright (c) 2014, 2015, 2016, 2017 by Thomas A. Manz and Nidia Gabaldon Limas. Rights reserved.

 USE module_common_variable_declarations
 USE module_global_parameters

 IMPLICIT NONE
 
 CONTAINS
 
 SUBROUTINE compute_summed_contact_exchange()
 !=================================================================================== 
 !The summed contact exchange will be computed for each atom by appropriately integrating over the grid points
 !=================================================================================== 
 
 REAL(kind=dp) :: temp,temp_vector(3),distance
 INTEGER :: comp
 
 INTEGER :: clock_count,clock_count_rate,clock_count_max,clock_count2
 
 WRITE(output_FID,*)' '
 WRITE(output_FID,*)'******************** TIME CONTROL ********************'
 WRITE(output_FID,*)'Starting compute_summed_contact_exchange'
 FLUSH(output_FID)
 CALL system_clock(clock_count,clock_count_rate,clock_count_max)
 WRITE(output_FID,*)'******************************************************'
 WRITE(output_FID,*)' '
 FLUSH(output_FID) 
  
 ALLOCATE(summed_contact_exchange(natoms))
 summed_contact_exchange = 0.0_dp
!$omp parallel default(none) &
!$omp shared(natoms,lower_nc,upper_nc,kCpoints,delta_nc,lower_nb,upper_nb,kBpoints,delta_nb,lower_na,upper_na,kApoints,delta_na, &
!$omp dot_product_total_spherical_avg_atomic_exchange_vectors,num_nonzero_components,local_spherical_avg_atomic_exchange_vectors, &
!$omp total_local_spherical_avg_atomic_exchange_vectors,corrected_total_density,chunk_size,boundary,center_shift,pixelvolume) &
!$omp private(j,nc,nb,na,kc,kb,ka,temp_vector,distance,shell_index,temp,comp) &
!$omp reduction(+:summed_contact_exchange)
 DO j=1,natoms 
!$omp do &
!$omp schedule(static,chunk_size) 
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
         IF ((shell_index <= nshells) .and. (dot_product_total_spherical_avg_atomic_exchange_vectors(ka,kb,kc) &
         > zero_tolerance**2)) THEN
           temp = 0.0_dp
           DO comp=1,num_nonzero_components
             temp = temp + local_spherical_avg_atomic_exchange_vectors(comp,shell_index,j)*(&
             total_local_spherical_avg_atomic_exchange_vectors(comp,ka,kb,kc) - local_spherical_avg_atomic_exchange_vectors&
             (comp,shell_index,j))
           END DO
           summed_contact_exchange(j) = summed_contact_exchange(j) + 2.0_dp*pixelvolume*corrected_total_density&
           (ka,kb,kc)*temp/dot_product_total_spherical_avg_atomic_exchange_vectors(ka,kb,kc)
         END IF
       END DO    
     END DO
   END DO 
!$omp end do
 END DO
!$omp end parallel

 WRITE(output_FID,*)' '
 WRITE(output_FID,*)'******************** TIME CONTROL ********************'
 CALL system_clock(clock_count2,clock_count_rate,clock_count_max)
 IF (clock_count_rate == 0) THEN
   WRITE(output_FID,*)'Processor clock not available, the time calculation wil not be performed.'
 ELSE
   seconds = modulo((clock_count2-clock_count),clock_count_max)/real(clock_count_rate)
   WRITE(output_FID,*)'Finished compute_summed_contact_exchange in ',seconds, 'seconds'
   FLUSH(output_FID)
 END IF
 WRITE(output_FID,*)'******************************************************'
 WRITE(output_FID,*)' '
 FLUSH(output_FID) 

 END SUBROUTINE compute_summed_contact_exchange
 
 
 
 SUBROUTINE integrate_bonding_terms()
 !=================================================================================== 
 
 REAL(kind=dp) :: temp1,temp2,dist_i,dist_j,atomjimage_x,atomjimage_y,atomjimage_z,grid_x,grid_y,grid_z,&
 bond_pair_matrix_temp_12,bond_pair_matrix_temp_13,bond_pair_matrix_temp_14,bond_pair_matrix_temp_15
 INTEGER :: index_i,index_j,i,j,repeata,repeatb,repeatc,comp

 INTEGER :: clock_count,clock_count_rate,clock_count_max,clock_count2,nb_tot,nc_tot,nbc_tot,nbc
 
 WRITE(output_FID,*)' '
 WRITE(output_FID,*)'******************** TIME CONTROL ********************'
 WRITE(output_FID,*)'Starting integrate_bonding_terms'
 FLUSH(output_FID)
 CALL system_clock(clock_count,clock_count_rate,clock_count_max)
 WRITE(output_FID,*)'******************************************************'
 WRITE(output_FID,*)' '
 FLUSH(output_FID) 
 
temp1 = 0.0_dp
temp2 = 0.0_dp
index_i = 1
index_j = 1
grid_x = 0.0_dp
grid_y = 0.0_dp
grid_z = 0.0_dp
na = 0
nb = 0
nc = 0
dist_i = 0.0_dp
dist_j = 0.0_dp
i = 1
j = 1
atomjimage_x = 0.0_dp
atomjimage_y = 0.0_dp
atomjimage_z = 0.0_dp
repeata = 0
repeatb = 0
repeatc = 0
 !Below is a loop that integrates the contact exchange terms and updates the bond pair matrix
!$omp parallel default(none) &
!$omp shared(num_sym_unique_bond_pairs,bond_pair_matrix,periodicC,partial_density,&
!$omp periodicB,periodicA,totnumC,totnumB,totnumA,boundary,origin,coords,vector1,vector2,vector3,pixelvolume,&
!$omp corrected_total_density,num_nonzero_components,local_spherical_avg_atomic_exchange_vectors,total_pseudodensity,&
!$omp dot_product_total_spherical_avg_atomic_exchange_vectors,i,j,nb_tot,nc_tot,nbc_tot,&
!$omp total_local_spherical_avg_atomic_exchange_vectors,bond_pair_matrix_temp_12,bond_pair_matrix_temp_13,&
!$omp bond_pair_matrix_temp_14,bond_pair_matrix_temp_15) &
!$omp private(nbc,na,nb,nc,kc,kb,ka,grid_x,grid_y,grid_z,dist_i,repeata,repeatb,repeatc,atomjimage_x,atomjimage_y,atomjimage_z,&
!$omp dist_j,index_j,comp,index_i,current_pair,temp1,temp2) 
 DO current_pair = 1,num_sym_unique_bond_pairs
!$omp single
   bond_pair_matrix(12,current_pair) = 0.0_dp
   i = bond_pair_matrix(1,current_pair)
   j = bond_pair_matrix(2,current_pair)
   nb_tot = nint(bond_pair_matrix(9,current_pair) - bond_pair_matrix(8,current_pair)) + 1
   nc_tot = nint(bond_pair_matrix(11,current_pair)- bond_pair_matrix(10,current_pair)) + 1
   nbc_tot = nb_tot * nc_tot
   bond_pair_matrix_temp_12 = 0.0_dp
   bond_pair_matrix_temp_13 = 0.0_dp
   bond_pair_matrix_temp_14 = 0.0_dp
   bond_pair_matrix_temp_15 = 0.0_dp
!$omp end single
!$omp do &
!$omp reduction(+:bond_pair_matrix_temp_12,bond_pair_matrix_temp_13,bond_pair_matrix_temp_14,bond_pair_matrix_temp_15) &
!$omp schedule(static,1)
   DO nbc = 0,nbc_tot - 1
     nb = nint(modulo(nbc,nb_tot) + bond_pair_matrix(8,current_pair))
     nc = nint(modulo(nint((nbc-modulo(nbc,nb_tot))/real(nb_tot)),nc_tot) + bond_pair_matrix(10,current_pair))
     IF (periodicC) THEN
       kc = modulo(nc,totnumC) + 1
     ELSE
       kc = nc + 1
       IF ((kc < 1) .or. (kc > totnumC)) THEN
         CYCLE
       END IF
     END IF    
     IF (periodicB) THEN
       kb = modulo(nb,totnumB) + 1
     ELSE
       kb = nb + 1
       IF ((kb < 1) .or. (kb > totnumB)) THEN
         CYCLE
       END IF
     END IF 
     DO na = nint(bond_pair_matrix(6,current_pair)),nint(bond_pair_matrix(7,current_pair))
       IF (periodicA) THEN
         ka = modulo(na,totnumA) + 1
       ELSE
         ka = na + 1
         IF ((ka < 1) .or. (ka > totnumA)) THEN
           CYCLE
         END IF
       END IF
       IF ((dot_product_total_spherical_avg_atomic_exchange_vectors(ka,kb,kc) > zero_tolerance**2)) THEN
         grid_x = na*boundary(1,1) + nb*boundary(2,1) + nc*boundary(3,1) + origin(1)
         grid_y = na*boundary(1,2) + nb*boundary(2,2) + nc*boundary(3,2) + origin(2)
         grid_z = na*boundary(1,3) + nb*boundary(2,3) + nc*boundary(3,3) + origin(3)
         dist_i = sqrt( (grid_x - coords(1,i))**2 + (grid_y - coords(2,i))**2 + (grid_z - coords(3,i))**2 )
         index_i = ceiling(scalefactor*dist_i + zero_tolerance)
         IF (index_i > nshells) THEN
           CYCLE
         END IF
         repeata = nint(bond_pair_matrix(3,current_pair))
         repeatb = nint(bond_pair_matrix(4,current_pair))
         repeatc = nint(bond_pair_matrix(5,current_pair))
         atomjimage_x = coords(1,j) + repeata*vector1(1) + repeatb*vector2(1) + repeatc*vector3(1)
         atomjimage_y = coords(2,j) + repeata*vector1(2) + repeatb*vector2(2) + repeatc*vector3(2)
         atomjimage_z = coords(3,j) + repeata*vector1(3) + repeatb*vector2(3) + repeatc*vector3(3)
         dist_j = sqrt( (grid_x - atomjimage_x)**2 + (grid_y - atomjimage_y)**2 + (grid_z - atomjimage_z)**2 )
         index_j = ceiling(scalefactor*dist_j + zero_tolerance)
         IF (index_j > nshells) THEN
           CYCLE
         END IF
         temp1 = 0.0_dp
         temp2 = 0.0_dp 
         temp1 = 2.0_dp*pixelvolume*corrected_total_density(ka,kb,kc)
         DO comp = 1,num_nonzero_components
           temp2 = temp2 + local_spherical_avg_atomic_exchange_vectors(comp,index_i,i)*local_spherical_avg_atomic_exchange_vectors&
           (comp,index_j,j)
         END DO
         bond_pair_matrix_temp_12 = bond_pair_matrix_temp_12 + temp1*temp2&
         /dot_product_total_spherical_avg_atomic_exchange_vectors(ka,kb,kc)
         bond_pair_matrix_temp_13 = bond_pair_matrix_temp_13 + temp1*temp2&
         /total_local_spherical_avg_atomic_exchange_vectors(1,ka,kb,kc)**2
         bond_pair_matrix_temp_14 = bond_pair_matrix_temp_14 + temp1*partial_density(index_i,i)*partial_density(index_j,j)&
         /total_pseudodensity(ka,kb,kc)**2
         bond_pair_matrix_temp_15 = bond_pair_matrix_temp_15 + (10.0_dp/3.0_dp)*temp1&
         *(temp2/dot_product_total_spherical_avg_atomic_exchange_vectors(ka,kb,kc))**2 

!!! line removed         *((partial_density(index_i,i)+partial_density(index_j,j))/total_pseudodensity(ka,kb,kc))**2

       END IF    
     END DO
   END DO
!$omp end do
!$omp single
   bond_pair_matrix(12,current_pair) = max(bond_pair_matrix_temp_12, 0.0_dp)
   bond_pair_matrix(13,current_pair) = max(bond_pair_matrix_temp_13, 0.0_dp)
   bond_pair_matrix(14,current_pair) = max(bond_pair_matrix_temp_14, 0.0_dp)
   bond_pair_matrix(15,current_pair) = max(bond_pair_matrix_temp_15, 0.0_dp)
!$omp end single
 END DO
!$omp end parallel

 
 WRITE(output_FID,*)' '
 WRITE(output_FID,*)'******************** TIME CONTROL ********************'
 CALL system_clock(clock_count2,clock_count_rate,clock_count_max)
 IF (clock_count_rate == 0) THEN
   WRITE(output_FID,*)'Processor clock not available, the time calculation wil not be performed.'
 ELSE
   seconds = modulo((clock_count2-clock_count),clock_count_max)/real(clock_count_rate)
   WRITE(output_FID,*)'Finished integrate_bonding_terms in ',seconds, 'seconds'
   FLUSH(output_FID)
 END IF
 WRITE(output_FID,*)'******************************************************'
 WRITE(output_FID,*)' '
 FLUSH(output_FID) 

 DEALLOCATE(dot_product_total_spherical_avg_atomic_exchange_vectors)
 
END SUBROUTINE integrate_bonding_terms

END MODULE module_integrate_bonding_terms
