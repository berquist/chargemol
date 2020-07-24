 MODULE module_local_multipole_moment_analysis
 ! Copyright (c) 2014, 2015, 2016, 2017 by Thomas A. Manz and Nidia Gabaldon Limas. Rights reserved.
 
 USE module_precision
 USE module_global_parameters
 USE module_common_variable_declarations
 USE module_matrix_operations

 IMPLICIT NONE
 
 CONTAINS
 
 SUBROUTINE local_multipole_moment_analysis()
 !=====================================================================================
 ! Local multipole moment analysis.
 !=====================================================================================
 
 REAL(kind=dp),ALLOCATABLE,DIMENSION(:) :: sum_dipole_x,sum_dipole_y,sum_dipole_z,sum_quadrupole_xy,sum_quadrupole_xz,&
 sum_quadrupole_yz,sum_quadrupole_x2minusy2,sum_quadrupole_3z2minusr2

 REAL(kind=dp) :: local_core_density,local_valence_density,x,y,z
 REAL(kind=dp) :: temp,temp_vector(3),distance,local_atomic_density
 INTEGER :: clock_count,clock_count_rate,clock_count_max,clock_count2
 REAL(kind=dp) :: temp_quadrupole_matrix(3,3),quadrupole_eigenvals(3)
 
 WRITE(output_FID,*)' '
 WRITE(output_FID,*)'******************** TIME CONTROL ********************'
 WRITE(output_FID,*)'Starting local_multipole_moment_analysis'
 FLUSH(output_FID)
 CALL system_clock(clock_count,clock_count_rate,clock_count_max)
 WRITE(output_FID,*)'******************************************************'
 WRITE(output_FID,*)' '
 FLUSH(output_FID) 
 
 DO j=1,natoms
   final_result(j,1) = j
   final_result(j,2) = atomic_number(j)
   final_result(j,3) = coords(1,j)
   final_result(j,4) = coords(2,j)
   final_result(j,5) = coords(3,j)
 END DO
 final_result(:,6)= atomic_number(:) - core_electrons(:) - valence_population(:)
 !Compute the local dipole and quadrupole moments in atomic units
 temp = 0.0_dp
 ALLOCATE(sum_dipole_x(natoms))
 ALLOCATE(sum_dipole_y(natoms))
 ALLOCATE(sum_dipole_z(natoms))
 ALLOCATE(sum_quadrupole_xy(natoms))
 ALLOCATE(sum_quadrupole_xz(natoms))
 ALLOCATE(sum_quadrupole_yz(natoms))
 ALLOCATE(sum_quadrupole_x2minusy2(natoms))
 ALLOCATE(sum_quadrupole_3z2minusr2(natoms))
 sum_dipole_x = 0.0_dp
 sum_dipole_y = 0.0_dp
 sum_dipole_z = 0.0_dp
 sum_quadrupole_xy = 0.0_dp
 sum_quadrupole_xz = 0.0_dp
 sum_quadrupole_yz = 0.0_dp
 sum_quadrupole_x2minusy2 = 0.0_dp
 sum_quadrupole_3z2minusr2 = 0.0_dp
!$omp parallel default(none) &
!$omp private(j,na,nb,nc,ka,kb,kc,temp_vector,distance,shell_index,local_atomic_density,local_core_density,&
!$omp local_valence_density,temp,x,y,z) &
!$omp shared(natoms,lower_nc,upper_nc,kCpoints,delta_nc,lower_na,upper_na,kApoints,delta_na,final_result,occupancy_correction, &
!$omp lower_nb,upper_nb,kBpoints,delta_nb,boundary,center_shift,total_pseudodensity,partial_density,core_density,&
!$omp core_pseudodensity,partial_core_density,normalization,pixelvolume,chunk_size,total_density) &
!$omp reduction(+:sum_dipole_x,sum_dipole_y,sum_dipole_z,sum_quadrupole_xy,sum_quadrupole_xz,sum_quadrupole_yz,&
!$omp sum_quadrupole_x2minusy2,sum_quadrupole_3z2minusr2)
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
         IF ((shell_index <= nshells) .and. (total_pseudodensity(ka,kb,kc) > zero_tolerance)) THEN
           local_atomic_density =  partial_density(shell_index,j)*total_density(ka,kb,kc)/total_pseudodensity(ka,kb,kc)
           IF (core_pseudodensity(ka,kb,kc) >= zero_tolerance) THEN
             local_core_density = partial_core_density(shell_index,j)*core_density(ka,kb,kc)/core_pseudodensity(ka,kb,kc)
           ELSE
             local_core_density = 0.0_dp
           END IF
           local_valence_density = (local_atomic_density - local_core_density)
           temp = normalization*local_valence_density*pixelvolume
           x = temp_vector(1)
           y = temp_vector(2)
           z = temp_vector(3)
           sum_dipole_x(j) = sum_dipole_x(j) - x*temp
           sum_dipole_y(j) = sum_dipole_y(j) - y*temp
           sum_dipole_z(j) = sum_dipole_z(j) - z*temp
           sum_quadrupole_xy(j) = sum_quadrupole_xy(j) - x*y*temp
           sum_quadrupole_xz(j) = sum_quadrupole_xz(j) - x*z*temp
           sum_quadrupole_yz(j) = sum_quadrupole_yz(j) - y*z*temp
           sum_quadrupole_x2minusy2(j) = sum_quadrupole_x2minusy2(j) - (x*x - y*y)*temp
           sum_quadrupole_3z2minusr2(j) = sum_quadrupole_3z2minusr2(j) - (3.0_dp*z*z-(x**2 + y**2 + z**2))*temp
         END IF
       END DO    
     END DO
   END DO
!$omp end do
 END DO
!$omp end parallel
 DO j=1,natoms
   final_result(j,7) = sum_dipole_x(j) + occupancy_correction(3,j)
   final_result(j,8) = sum_dipole_y(j) + occupancy_correction(4,j)
   final_result(j,9) = sum_dipole_z(j) + occupancy_correction(5,j)
   final_result(j,10) = sqrt(final_result(j,7)**2 + final_result(j,8)**2 + final_result(j,9)**2)
   final_result(j,11) = sum_quadrupole_xy(j) + occupancy_correction(9,j)
   final_result(j,12) = sum_quadrupole_xz(j) + occupancy_correction(10,j)
   final_result(j,13) = sum_quadrupole_yz(j) + occupancy_correction(11,j)
   final_result(j,14) = sum_quadrupole_x2minusy2(j) + occupancy_correction(6,j) - occupancy_correction(7,j)
   final_result(j,15) = sum_quadrupole_3z2minusr2(j) + 2.0_dp*occupancy_correction(8,j) - occupancy_correction(6,j) - &
   occupancy_correction(7,j)
   temp_quadrupole_matrix(1,:) = [0.5_dp*(sum_quadrupole_x2minusy2(j)-sum_quadrupole_3z2minusr2(j)/3.0_dp),sum_quadrupole_xy(j),&
   sum_quadrupole_xz(j)]
   temp_quadrupole_matrix(2,:) = [sum_quadrupole_xy(j),0.5_dp*(-sum_quadrupole_x2minusy2(j)-sum_quadrupole_3z2minusr2(j)/3.0_dp),&
   sum_quadrupole_yz(j)]
   temp_quadrupole_matrix(3,:) = [sum_quadrupole_xz(j),sum_quadrupole_yz(j),sum_quadrupole_3z2minusr2(j)/3.0_dp]
   quadrupole_eigenvals = eig(temp_quadrupole_matrix)
   final_result(j,16) = quadrupole_eigenvals(1)
   final_result(j,17) = quadrupole_eigenvals(2)
   final_result(j,18) = quadrupole_eigenvals(3)
 END DO
 WRITE(output_FID,*) 'Multipole analysis for each of the expansion sites.'
 FLUSH(output_FID)
 WRITE(output_FID,*) 'XYZ coordinates, net charges, and multipoles are in atomic units. Dipoles and quadrupoles are for &
 &valence electrons.'
 WRITE(output_FID,*) 'center number, atomic number, x, y, z, net_charge, dipole_x, dipole_y, dipole_z, dipole_mag, Qxy, &
 &Qxz, Qyz, Q(x^2-y^2), Q(3z^2 - R^2),three eigenvalues of traceless quadrupole moment tensor'
 FLUSH(output_FID)
 DO j=1,natoms
   WRITE(output_FID,'(2i6,16f10.6)') int(final_result(j,1)),int(final_result(j,2)),(final_result(j,i),i=3,18)
 END DO
 FLUSH(output_FID)

 WRITE(output_FID,*)' '
 WRITE(output_FID,*)'******************** TIME CONTROL ********************'
 CALL system_clock(clock_count2,clock_count_rate,clock_count_max)
 IF (clock_count_rate == 0) THEN
   WRITE(output_FID,*)'Processor clock not available, the time calculation wil not be performed.'
 ELSE
   seconds = modulo((clock_count2-clock_count),clock_count_max)/real(clock_count_rate)
   WRITE(output_FID,*)'Finished local_multipole_moment_analysis in ',seconds, 'seconds'
   FLUSH(output_FID)
 END IF
 WRITE(output_FID,*)'******************************************************'
 WRITE(output_FID,*)' '
 FLUSH(output_FID) 
 
 END SUBROUTINE local_multipole_moment_analysis
 
 END MODULE module_local_multipole_moment_analysis