 MODULE module_check_noncollinear_XC_functional
 ! Copyright (c) 2014, 2015, 2016, 2017 by Thomas A. Manz and Nidia Gabaldon Limas. Rights reserved.

 USE module_precision
 USE module_global_parameters
 USE module_common_variable_declarations

 IMPLICIT NONE 

 CONTAINS
 
 SUBROUTINE check_noncollinear_XC_functional()
 !===================================================================================
 ! This program separately computes the integrated changes of the directional and 
 ! magnitude components to the spin magnetization density. These quantities can be
 ! used to determine whether a two or four-component non-collinear DFT functional 
 ! is needed.
 !===================================================================================
 
 REAL(kind=dp) :: spin_density_vector_magnitude,spin_density_vector_magnitude_a,spin_density_vector_magnitude_b,&
 spin_density_vector_magnitude_c,delta_magnitude_x,delta_magnitude_y,delta_magnitude_z,sq_gradient_mag,delta_mx_x,delta_mx_y,&
 delta_mx_z,sq_gradient_mx,delta_my_x,delta_my_y,delta_my_z,sq_gradient_my,delta_mz_x,delta_mz_z,sq_gradient_mz,&
 sum_sq_gradient_mag,sum_sq_gradient_tot,sum_sq_gradient_dir,delta_mz_y
 INTEGER :: clock_count,clock_count_rate,clock_count_max,clock_count2
 
 WRITE(output_FID,*)' '
 WRITE(output_FID,*)'******************** TIME CONTROL ********************'
 WRITE(output_FID,*)'Starting check_noncollinear_XC_functional'
 FLUSH(output_FID)
 CALL system_clock(clock_count,clock_count_rate,clock_count_max)
 WRITE(output_FID,*)'******************************************************'
 WRITE(output_FID,*)' '
 FLUSH(output_FID)
 
 sum_sq_gradient_mag = 0.0_dp
 sum_sq_gradient_tot = 0.0_dp
!$omp parallel do default(none) &
!$omp private(kc,kc_plus,kb,kb_plus,ka,ka_plus,spin_density_vector_magnitude,&
!$omp spin_density_vector_magnitude_a,spin_density_vector_magnitude_b,spin_density_vector_magnitude_c,delta_magnitude_x,&
!$omp delta_magnitude_y,delta_magnitude_z,sq_gradient_mag,delta_mx_x,delta_mx_y,delta_mx_z,sq_gradient_mx,delta_my_x,&
!$omp delta_my_y,delta_my_z,sq_gradient_my,delta_mz_x,delta_mz_y,delta_mz_z,sq_gradient_mz) &
!$omp shared(totnumC,periodicC,totnumB,periodicB,totnumA,periodicA,spin_density_vector,inv_boundary,pixelvolume) &
!$omp reduction(+:sum_sq_gradient_mag,sum_sq_gradient_tot) &
!$omp schedule(dynamic,chunk_size)
 DO kc = 1,totnumC
  IF ( .not. ((.not. periodicC) .and. (kc == totnumC)) ) THEN
   kc_plus = modulo(kc,totnumC) + 1
   DO kb = 1,totnumB
    IF ( .not. ((.not. periodicB) .and. (kb == totnumB)) ) THEN
     kb_plus = modulo(kb,totnumB) + 1
     DO ka = 1,totnumA
      IF ( .not. ((.not. periodicA) .and. (ka == totnumA)) ) THEN
       ka_plus = modulo(ka,totnumA) + 1
       !Spin density magnitude changes
       spin_density_vector_magnitude = sqrt(spin_density_vector(1,ka,kb,kc)**2 + spin_density_vector(2,ka,kb,kc)**2 + &
       spin_density_vector(3,ka,kb,kc)**2)
       spin_density_vector_magnitude_a = sqrt(spin_density_vector(1,ka_plus,kb,kc)**2 + spin_density_vector(2,ka_plus,kb,kc)**2 &
       + spin_density_vector(3,ka_plus,kb,kc)**2)
       spin_density_vector_magnitude_b = sqrt(spin_density_vector(1,ka,kb_plus,kc)**2 + spin_density_vector(2,ka,kb_plus,kc)**2 &
       + spin_density_vector(3,ka,kb_plus,kc)**2)
       spin_density_vector_magnitude_c = sqrt(spin_density_vector(1,ka,kb,kc_plus)**2 + spin_density_vector(2,ka,kb,kc_plus)**2 &
       + spin_density_vector(3,ka,kb,kc_plus)**2)
       delta_magnitude_x = (spin_density_vector_magnitude_a - spin_density_vector_magnitude)*inv_boundary(1,1) + &
       (spin_density_vector_magnitude_b - spin_density_vector_magnitude)*inv_boundary(2,1) + (spin_density_vector_magnitude_c &
       - spin_density_vector_magnitude)*inv_boundary(3,1)
       delta_magnitude_y = (spin_density_vector_magnitude_a - spin_density_vector_magnitude)*inv_boundary(1,2) + &
       (spin_density_vector_magnitude_b - spin_density_vector_magnitude)*inv_boundary(2,2) + (spin_density_vector_magnitude_c &
       - spin_density_vector_magnitude)*inv_boundary(3,2)
       delta_magnitude_z = (spin_density_vector_magnitude_a - spin_density_vector_magnitude)*inv_boundary(1,3) + &
       (spin_density_vector_magnitude_b - spin_density_vector_magnitude)*inv_boundary(2,3) + (spin_density_vector_magnitude_c &
       - spin_density_vector_magnitude)*inv_boundary(3,3)
       sq_gradient_mag = delta_magnitude_x**2 + delta_magnitude_y**2 + delta_magnitude_z**2
       !mx changes 
       delta_mx_x = (spin_density_vector(1,ka_plus,kb,kc) - spin_density_vector(1,ka,kb,kc))*inv_boundary(1,1) + &
       (spin_density_vector(1,ka,kb_plus,kc) - spin_density_vector(1,ka,kb,kc))*inv_boundary(2,1) + (spin_density_vector(1,ka,kb&
       ,kc_plus) - spin_density_vector(1,ka,kb,kc))*inv_boundary(3,1)
       delta_mx_y = (spin_density_vector(1,ka_plus,kb,kc) - spin_density_vector(1,ka,kb,kc))*inv_boundary(1,2) + &
       (spin_density_vector(1,ka,kb_plus,kc) - spin_density_vector(1,ka,kb,kc))*inv_boundary(2,2) + (spin_density_vector(1,ka,kb&
       ,kc_plus) - spin_density_vector(1,ka,kb,kc))*inv_boundary(3,2)
       delta_mx_z = (spin_density_vector(1,ka_plus,kb,kc) - spin_density_vector(1,ka,kb,kc))*inv_boundary(1,3) + &
       (spin_density_vector(1,ka,kb_plus,kc) - spin_density_vector(1,ka,kb,kc))*inv_boundary(2,3) + (spin_density_vector(1,ka,kb&
       ,kc_plus) - spin_density_vector(1,ka,kb,kc))*inv_boundary(3,3)   
       sq_gradient_mx = delta_mx_x**2 + delta_mx_y**2 + delta_mx_z**2
       !my changes 
       delta_my_x = (spin_density_vector(2,ka_plus,kb,kc) - spin_density_vector(2,ka,kb,kc))*inv_boundary(1,1) + &
       (spin_density_vector(2,ka,kb_plus,kc) - spin_density_vector(2,ka,kb,kc))*inv_boundary(2,1) + (spin_density_vector(2,ka,kb&
       ,kc_plus) - spin_density_vector(2,ka,kb,kc))*inv_boundary(3,1)
       delta_my_y = (spin_density_vector(2,ka_plus,kb,kc) - spin_density_vector(2,ka,kb,kc))*inv_boundary(1,2) + &
       (spin_density_vector(2,ka,kb_plus,kc) - spin_density_vector(2,ka,kb,kc))*inv_boundary(2,2) + (spin_density_vector(2,ka,kb&
       ,kc_plus) - spin_density_vector(2,ka,kb,kc))*inv_boundary(3,2)
       delta_my_z = (spin_density_vector(2,ka_plus,kb,kc) - spin_density_vector(2,ka,kb,kc))*inv_boundary(1,3) + &
       (spin_density_vector(2,ka,kb_plus,kc) - spin_density_vector(2,ka,kb,kc))*inv_boundary(2,3) + (spin_density_vector(2,ka,kb&
       ,kc_plus) - spin_density_vector(2,ka,kb,kc))*inv_boundary(3,3)
       sq_gradient_my = delta_my_x**2 + delta_my_y**2 + delta_my_z**2
       !mz changes 
       delta_mz_x = (spin_density_vector(3,ka_plus,kb,kc) - spin_density_vector(3,ka,kb,kc))*inv_boundary(1,1) + &
       (spin_density_vector(3,ka,kb_plus,kc) - spin_density_vector(3,ka,kb,kc))*inv_boundary(2,1) + (spin_density_vector(3,ka,kb&
       ,kc_plus) - spin_density_vector(3,ka,kb,kc))*inv_boundary(3,1)
       delta_mz_y = (spin_density_vector(3,ka_plus,kb,kc) - spin_density_vector(3,ka,kb,kc))*inv_boundary(1,2) + &
       (spin_density_vector(3,ka,kb_plus,kc) - spin_density_vector(3,ka,kb,kc))*inv_boundary(2,2) + (spin_density_vector(3,ka,kb&
       ,kc_plus) - spin_density_vector(3,ka,kb,kc))*inv_boundary(3,2)
       delta_mz_z = (spin_density_vector(3,ka_plus,kb,kc) - spin_density_vector(3,ka,kb,kc))*inv_boundary(1,3) + &
       (spin_density_vector(3,ka,kb_plus,kc) - spin_density_vector(3,ka,kb,kc))*inv_boundary(2,3) + (spin_density_vector(3,ka,kb&
       ,kc_plus) - spin_density_vector(3,ka,kb,kc))*inv_boundary(3,3)
       sq_gradient_mz = delta_mz_x**2 + delta_mz_y**2 + delta_mz_z**2
       IF (spin_density_vector_magnitude > zero_tolerance) THEN
         sum_sq_gradient_mag = sum_sq_gradient_mag + sq_gradient_mag*pixelvolume
         sum_sq_gradient_tot = sum_sq_gradient_tot + (sq_gradient_mx + sq_gradient_my + sq_gradient_mz)*pixelvolume
       END IF
      END IF
     END DO
    END IF
   END DO
  END IF     
 END DO      
!$omp end parallel do 
 WRITE(output_FID,*)'Changes in the magnitude of the spin magnetization density contribute an integrated weight of '
 WRITE(output_FID,*) sum_sq_gradient_mag
 WRITE(output_FID,*) 'The total changes in the spin magnetization density contribute an integrated weight of'
 WRITE(output_FID,*) sum_sq_gradient_tot
 WRITE(output_FID,*) 'Based on these calculations, directional changes of the spin magnetization density contribute an integrated &
 &weight of'
 sum_sq_gradient_dir = max((sum_sq_gradient_tot - sum_sq_gradient_mag), 0.0_dp)
 WRITE(output_FID,*) sum_sq_gradient_dir

 WRITE(output_FID,*)' '
 WRITE(output_FID,*)'******************** TIME CONTROL ********************'
 CALL system_clock(clock_count2,clock_count_rate,clock_count_max)
 IF (clock_count_rate == 0) THEN
   WRITE(output_FID,*)'Processor clock not available, the time calculation wil not be performed.'
 ELSE
   seconds = modulo((clock_count2-clock_count),clock_count_max)/real(clock_count_rate)
   WRITE(output_FID,*)'Finished check_noncollinear_XC_functional in ',seconds, 'seconds'
   FLUSH(output_FID)
 END IF
 WRITE(output_FID,*)'******************************************************'
 WRITE(output_FID,*)' '
 FLUSH(output_FID) 

 END SUBROUTINE check_noncollinear_XC_functional
 
 END MODULE module_check_noncollinear_XC_functional