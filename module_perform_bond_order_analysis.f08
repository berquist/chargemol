 MODULE module_perform_bond_order_analysis
 !=====================================================================================
 ! Copyright (c) 2014, 2015, 2016, 2017 by Thomas A. Manz and Nidia Gabaldon Limas. Rights reserved.
 
 USE module_precision
 USE module_global_parameters
 USE module_common_variable_declarations
 USE module_prepare_BO_density_grids
 USE module_compute_local_atomic_exchange_vectors
 USE module_initialize_bond_pair_matrix
 USE module_integrate_bonding_terms
 USE module_calculate_final_BOs

 IMPLICIT NONE
 
 CONTAINS
 
 SUBROUTINE perform_bond_order_analysis()
 !====================================================================================
 
 INTEGER :: clock_count,clock_count_rate,clock_count_max,clock_count2
 
 WRITE(output_FID,*)' '
 WRITE(output_FID,*)'******************** TIME CONTROL ********************'
 WRITE(output_FID,*)'Starting perform_BO_analysis'
 FLUSH(output_FID)
 CALL system_clock(clock_count,clock_count_rate,clock_count_max)
 WRITE(output_FID,*)'******************************************************'
 WRITE(output_FID,*)' '
 FLUSH(output_FID) 

 WRITE(output_FID,*)'Starting Effective Bond Order (BO) Analysis'
 FLUSH(output_FID)
 !Set the number of nonzero density and spin density components
 IF (spin_available) THEN
   IF (non_collinear) THEN
     num_nonzero_components = 4
   ELSE
     num_nonzero_components = 2
   END IF
 ELSE
   num_nonzero_components = 1
 END IF
 CALL prepare_BO_density_grids()
 CALL compute_local_atomic_exchange_vectors()
 CALL initialize_bond_pair_matrix()
 IF (flag == 1) THEN
   RETURN !  break
 END IF
 CALL compute_summed_contact_exchange()
 CALL integrate_bonding_terms()
 DEALLOCATE(total_local_spherical_avg_atomic_exchange_vectors)
 IF (flag == 1) THEN
   RETURN !  break
 END IF
 CALL calculate_final_BOs()
 WRITE(output_FID,*)'Bond order analysis complete.'
 FLUSH(output_FID)
 
 WRITE(output_FID,*)' '
 WRITE(output_FID,*)'******************** TIME CONTROL ********************'
 CALL system_clock(clock_count2,clock_count_rate,clock_count_max)
 IF (clock_count_rate == 0) THEN
   WRITE(output_FID,*)'Processor clock not available, the time calculation wil not be performed.'
 ELSE
   seconds = modulo((clock_count2-clock_count),clock_count_max)/real(clock_count_rate)
   WRITE(output_FID,*)'Finished perform_BO_analysis in ',seconds, 'seconds'
   FLUSH(output_FID)
 END IF
 WRITE(output_FID,*)'******************************************************'
 WRITE(output_FID,*)' '
 FLUSH(output_FID) 
 
 END SUBROUTINE perform_bond_order_analysis

 END MODULE module_perform_bond_order_analysis