 MODULE module_calculate_final_BOs
! Copyright (c) 2014, 2015, 2016, 2017 by Thomas A. Manz and Nidia Gabaldon Limas. Rights reserved.
 
 USE module_common_variable_declarations
 USE module_global_parameters
 USE module_generate_bond_order_file
 USE module_string_utilities

 IMPLICIT NONE

 CONTAINS
 
 SUBROUTINE calculate_final_BOs()
 !===================================================================================

 INTEGER :: count
 REAL(kind=dp),ALLOCATABLE,DIMENSION(:) :: check_summed_contact_exchange,SBO_correction,sum_squared_CE,bond_scaling_sum

 INTEGER :: clock_count,clock_count_rate,clock_count_max,clock_count2
 
 WRITE(output_FID,*)' '
 WRITE(output_FID,*)'******************** TIME CONTROL ********************'
 WRITE(output_FID,*)'Starting calculate_final_BOs'
 FLUSH(output_FID)
 CALL system_clock(clock_count,clock_count_rate,clock_count_max)
 WRITE(output_FID,*)'******************************************************'
 WRITE(output_FID,*)' '
 FLUSH(output_FID)  
 
 WRITE(output_FID,*)'Checking for consistency in the summed contact exchange'
 ALLOCATE(check_summed_contact_exchange(natoms))
 check_summed_contact_exchange = 0.0_dp
 ALLOCATE(sum_squared_CE(natoms))
 sum_squared_CE = 0.0_dp
 DO current_pair = 1,num_sym_unique_bond_pairs
   check_summed_contact_exchange(nint(bond_pair_matrix(1,current_pair))) = check_summed_contact_exchange(nint(&
   bond_pair_matrix(1,current_pair))) + bond_pair_matrix(12,current_pair)
   check_summed_contact_exchange(nint(bond_pair_matrix(2,current_pair))) = check_summed_contact_exchange(nint(&
   bond_pair_matrix(2,current_pair))) + bond_pair_matrix(12,current_pair)
   i = nint(bond_pair_matrix(1,current_pair))
   j = nint(bond_pair_matrix(2,current_pair))
   sum_squared_CE(i) = sum_squared_CE(i) + bond_pair_matrix(12,current_pair)**2
   sum_squared_CE(j) = sum_squared_CE(j) + bond_pair_matrix(12,current_pair)**2
 END DO
 MAD1 = maxval(abs(check_summed_contact_exchange - summed_contact_exchange))
 WRITE(output_FID,'(a,f12.6)')' The maximum error in the summed contact exchange is ',MAD1 
 FLUSH(output_FID)  
 WRITE(output_FID,*) ' '  
 WRITE(output_FID,*) 'The sum of contact exchanges for each atom are:' 
 DO j = 1,natoms
   WRITE(output_FID,'(f12.6)')(summed_contact_exchange(j))
 END DO 
 WRITE(output_FID,*) ' '
 FLUSH(output_FID)
 ALLOCATE(bond_scaling_sum(natoms))
 bond_scaling_sum = 0.0_dp
 DO current_pair = 1,num_sym_unique_bond_pairs
   i = nint(bond_pair_matrix(1,current_pair))
   j = nint(bond_pair_matrix(2,current_pair))
   IF ( (sum_squared_CE(i) > zero_tolerance) .and. (sum_squared_CE(j) > zero_tolerance) ) THEN
      bond_pair_matrix(16,current_pair) = 1.0_dp - (tanh( (summed_contact_exchange(i)**2/sum_squared_CE(i) +&
      summed_contact_exchange(j)**2/sum_squared_CE(j) - 2.0_dp)/26.0_dp))**2
      bond_pair_matrix(17,current_pair) = bond_pair_matrix(15,current_pair) + bond_pair_matrix(12,current_pair)**2/6.0_dp
      bond_pair_matrix(18,current_pair) = bond_pair_matrix(16,current_pair)*&
      min(bond_pair_matrix(17,current_pair),bond_pair_matrix(12,current_pair))
      bond_pair_matrix(19,current_pair) = bond_pair_matrix(12,current_pair) + bond_pair_matrix(18,current_pair)
      bond_scaling_sum(i) = bond_scaling_sum(i) + bond_pair_matrix(18,current_pair)
      bond_scaling_sum(j) = bond_scaling_sum(j) + bond_pair_matrix(18,current_pair)
   END IF
 END DO
 ALLOCATE(SBO_correction(natoms))
 SBO_correction = 0.0_dp
 count = 0
 DO current_pair = 1,num_sym_unique_bond_pairs
   i = nint(bond_pair_matrix(1,current_pair))
   j = nint(bond_pair_matrix(2,current_pair))
   IF ( (bond_scaling_sum(i) > zero_tolerance) .and. (bond_scaling_sum(j) > zero_tolerance) ) THEN
      bond_pair_matrix(20,current_pair) = bond_pair_matrix(12,current_pair) + bond_pair_matrix(18,current_pair)*&
      min(1.0_dp,(atomic_number(i) - final_result(i,6) - 0.5_dp*summed_contact_exchange(i))/bond_scaling_sum(i),&
      (atomic_number(j) - final_result(j,6) - 0.5_dp*summed_contact_exchange(j))/bond_scaling_sum(j))
      SBO_correction(i) = SBO_correction(i) + bond_pair_matrix(20,current_pair) - bond_pair_matrix(12,current_pair)
      SBO_correction(j) = SBO_correction(j) + bond_pair_matrix(20,current_pair) - bond_pair_matrix(12,current_pair)
   END IF
   IF (bond_pair_matrix(19,current_pair) - bond_pair_matrix(20,current_pair) > 0.1_dp*BO_print_cutoff) THEN
      count = count + 1
   END IF
 END DO
 ALLOCATE(atomic_summed_BO(natoms))
 atomic_summed_BO = 0.0_dp
 WRITE(output_FID,*) 'The sum of BOs for each atom are:'   
 DO j = 1,natoms
   atomic_summed_BO(j) = summed_contact_exchange(j) + SBO_correction(j)
   WRITE(output_FID,'(f12.6)')(atomic_summed_BO(j))
 END DO  
 WRITE(output_FID,*) ' '
 FLUSH(output_FID)
 WRITE(output_FID,'(A)')' The final bond pair matrix is'
 DO i=1,num_sym_unique_bond_pairs
   WRITE(output_FID,'(2i6,18f12.6)')int(bond_pair_matrix(1,i)),int(bond_pair_matrix(2,i)),(bond_pair_matrix(j,i),j=3,20)
 END DO
 WRITE(output_FID,'(A)')' The legend for the bond pair matrix follows:'
 WRITE(output_FID,'(A)')' Columns 1 and 2 = atom numbers'
 WRITE(output_FID,'(A)')' Columns 3 to 5 = repeata, repeatb, repeatc (i.e., periodic translation of second atom)'
 WRITE(output_FID,'(A)')' Columns 6 to 11 = corners of the shared parallepiped: min na, max na, min nb, max nb, min nc, max nc'
 WRITE(output_FID,'(A)')' Columns 12 =  contact exchange'
 WRITE(output_FID,'(A)')' Columns 13 =  a term for computing average spin polarization of bonding'
 WRITE(output_FID,'(A)')' Columns 14 =  overlap population'
 WRITE(output_FID,'(A)')' Columns 15 =  the integrated second order atomic exchange propensities for computing bond order'
 WRITE(output_FID,'(A)')' Columns 16 =  the coordination term constructed from tanh function'
 WRITE(output_FID,'(A)')' Columns 17 =  the pairwise term'
 WRITE(output_FID,'(A)')' Columns 18 =  expansion term combining coordination and pairwise terms'
 WRITE(output_FID,'(A)')' Columns 19 =  bond index before imposing the atom self-exchange constraint'
 WRITE(output_FID,'(A)')' Columns 20 =  final bond order'
 WRITE(output_FID,'(A)')' ' 
 WRITE(output_FID,'(A,I6,A)')' The atom self-exchange constraint affected (by at least 0.1_dp*BO_print_cutoff) ',&
 count,' bond orders.'
 WRITE(output_FID,'(A)')' '
 IF (StrLowCase(trim(adjustl(charge_type))) == 'ddec3') THEN
   OPEN(newunit = BO_FID, file = 'DDEC3_bond_orders.xyz', status = 'replace')
 ELSE
    OPEN(newunit = BO_FID, file = 'DDEC6_even_tempered_bond_orders.xyz', status = 'replace')
 END IF
 CALL  generate_bond_order_file()
 CLOSE(BO_FID)
 
 WRITE(output_FID,*)' '
 WRITE(output_FID,*)'******************** TIME CONTROL ********************'
 CALL system_clock(clock_count2,clock_count_rate,clock_count_max)
 IF (clock_count_rate == 0) THEN
   WRITE(output_FID,*)'Processor clock not available, the time calculation will not be performed.'
 ELSE
   seconds = modulo((clock_count2-clock_count),clock_count_max)/real(clock_count_rate)
   WRITE(output_FID,*)'Finished calculate_final_BOs in ',seconds, 'seconds'
   FLUSH(output_FID)
 END IF
 WRITE(output_FID,*)'******************************************************'
 WRITE(output_FID,*)' '
 FLUSH(output_FID) 
    
 END SUBROUTINE calculate_final_BOs
 
 END MODULE module_calculate_final_BOs
 