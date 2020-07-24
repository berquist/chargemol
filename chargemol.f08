 PROGRAM chargemol
!===================================================================================
!
! MAIN PROGRAM DDEC
! Copyright (c) 2014, 2015, 2016, 2017 by Thomas A. Manz and Nidia Gabaldon Limas. Rights reserved.
!
!===================================================================================

 USE module_precision
 USE module_global_parameters
 USE module_string_utilities
 USE module_common_variable_declarations
 USE module_read_wfx
 USE module_compute_atomic_radial_moments
 USE module_run_valence_core_densities
 USE module_core_iterator
 USE module_DDEC3_valence_iterator
 USE module_DDEC6_valence_iterator
 USE module_local_multipole_moment_analysis
 USE module_total_multipole_moment_analysis
 USE module_cloud_penetration
 USE module_generate_net_atomic_charge_file
 USE module_quote
 USE module_spin_functions
 USE module_collinear_spin_moments_iterator
 USE module_noncollinear_spin_moments_iterator
 USE module_perform_bond_order_analysis
 USE module_read_job_control
 USE module_calculate_atomic_polarizabilities_upper_bound
 USE module_generate_atomic_polarizability_upper_bound_file
 USE module_check_atomic_reference_polarizabilities_available
 USE module_print_overlap_populations
 USE module_print_atomic_densities_file
!$ USE omp_lib
 
 IMPLICIT NONE

 INTEGER :: k
 INTEGER :: clock_count,clock_count_rate,clock_count_max,clock_count2
 LOGICAl :: file_was_selected
 
 !Setting the default number of core electrons
 num_core=(/0, 0, 2, 2, 2, 2, 2, 2, 2, 2, 10, 10, 10, 10, 10, 10, 10, 10, 18, 18, 18, 18, 18, 18,&
 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36,&
 36, 36, 36, 36, 36, 54, 54, 54 ,54, 54, 54, 54, 54, 54, 54, 54, 54, 54, 54, 54, 54, 54, 54, 54, 54,&
 54, 54, 54, 54, 54, 54, 54, 54, 54, 54, 54, 54, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86,&
 86, 86, 86, 86, 86, 86, 86, 86, 86, 86/) 
 
 !Change directory
 IF(run_interactive) THEN
   DO k=1,10  
     PRINT*,'Enter the complete directory of the job.'
     PRINT*,'   Example: C:\Users\ngabaldo\DDEC_06_09_2014\H2'
     READ(*,*)job_directory
     WRITE(*,*) job_directory
     iostat_value = CHDIR(job_directory)
     IF (iostat_value == 0) EXIT        
   END DO  
 END IF     
 
 compute_BOs = .TRUE.
 charge_type = 'ddec6'
 print_atomic_densities = .FALSE.
 CALL read_job_control()
 
 WRITE(output_FID,*)' '
 WRITE(output_FID,*)'******************** TIME CONTROL ********************'
 WRITE(output_FID,*)'Starting ',version
 FLUSH(output_FID)
 CALL system_clock(clock_count,clock_count_rate,clock_count_max)
 WRITE(output_FID,*)'******************************************************'
 WRITE(output_FID,*)' '
 FLUSH(output_FID) 
 
 !Write date in output file
 CALL date_and_time(date,time,zone,values)
 CALL date_and_time(DATE=date,ZONE=zone)
 CALL date_and_time(TIME=time)
 CALL date_and_time(VALUES=values)
 WRITE(output_FID,*) 'Starting DDEC program'
 WRITE(output_FID,*) date(1:4),"/",date(5:6),"/",date(7:),"  ",time(1:2),":",time(3:4),":",time(5:6)
 WRITE(output_FID,*)'Copyright (c) 2014, 2015, 2016 by Thomas A. Manz and Nidia Gabaldon Limas. Rights reserved.'
 FLUSH(output_FID) 
 !$ number_of_threads = omp_get_max_threads()
 !Write OpenMP information
 run_parallel = .false.
!$ run_parallel = .true.
 WRITE(output_FID,*)' '
 WRITE(output_FID,*)'**************** THREAD INFORMATION ******************'
 IF (run_parallel) THEN
   WRITE(output_FID,*) 'Job running using OpenMP.'
   WRITE(output_FID,*) 'The number of parallel threads is: ',number_of_threads
 ELSE
   WRITE(output_FID,*) 'Job running using serial mode.'
 END IF
 WRITE(output_FID,*) ' '
 WRITE(output_FID,*)'******************************************************'
 WRITE(output_FID,*) ' '
 
 IF ((.not. periodicA).and.(periodicA .neqv. .FALSE.))THEN
   WRITE(output_FID,*)'An invalid option has been selected for the periodicity.'
   WRITE(output_FID,*)'Program will terminate'
   FLUSH(output_FID)
   STOP
 END IF
 IF ((.not. periodicB).and.(periodicB .neqv. .FALSE.))THEN
   WRITE(output_FID,*)'An invalid option has been selected for the periodicity.'
   WRITE(output_FID,*)'Program will terminate'
   FLUSH(output_FID)
   STOP
 END IF
 IF ((.not. periodicC).and.(periodicC .neqv. .FALSE.))THEN
   WRITE(output_FID,*)'An invalid option has been selected for the periodicity.'
   WRITE(output_FID,*)'Program will terminate'
   FLUSH(output_FID)
   STOP
 END IF
 IF (input_type == 'nan') THEN
   WRITE(output_FID,*) 'Program could not find an input file.'
   WRITE(output_FID,*) 'Possible solutions: '
   WRITE(output_FID,*) '1. Specify the name of input file on second line of job_control.txt.'
   WRITE(output_FID,*) 'NOTE: This step is optional for VASP and cube files, but required for wfx and xsf jobs'
   WRITE(output_FID,*) '2. Specify full directory path as a command line argument.'
   WRITE(output_FID,*) 'NOTES:'
   WRITE(output_FID,*) '(a) Linux directories should end in /.'
   WRITE(output_FID,*) '   Example: /home/ngabaldo/bin/atomic_densities/'
   WRITE(output_FID,*) '(b) Windows directories should end in \.'
   WRITE(output_FID,*) '   Example: C:\home\ngabaldo\bin\atomic_densities\'
   WRITE(output_FID,*) '(c) Some queuing systems may not recognize ~/ as your home directory. In this case, specify an explicit &
   &directory path.'
   WRITE(output_FID,*) '   Example: /home/ngabaldo/bin/atomic_densities/'
   FLUSH(output_FID)
   CLOSE(output_FID)
   STOP
 END IF
 
 ! Read the electron density, atomic coordinates, etc. and perform the analysis
 ALLOCATE(radial_shell_volume(nshells))
 DO k=1,nshells
   radial_shell_volume(k) = 4.0_dp*pi*(k**3 - (k-1)**3)/(3.0_dp*scalefactor**3)
 END DO

 CALL run_valence_core_densities()
 CALL core_iterator()
 
 IF(trim(adjustl(StrLowCase(charge_type))) == 'ddec3') THEN
   charge_type = 'ddec3'
   CALL DDEC3_valence_iterator()
 ELSE
   charge_type = 'ddec6'
   CALL DDEC6_valence_iterator()
 END IF
 
 CALL local_multipole_moment_analysis()
 DEALLOCATE(core_pseudodensity)
 DEALLOCATE(core_density)
 CAlL total_multipole_moment_analysis()
 CALL cloud_penetration()
 CALL generate_net_atomic_charge_file()
 IF (spin_available) THEN 
   CALL generate_spin_lookup_tables()
   IF (non_collinear) THEN
     CALL noncollinear_spin_moments_iterator( )
     DEALLOCATE(spin_density_vector)
   ELSE
     CALL collinear_spin_moments_iterator( )
     DEALLOCATE(spin_density)
   END IF    
 END IF

 IF (compute_BOs) THEN
   CALL perform_bond_order_analysis()
   CALL compute_atomic_radial_moments()
   DEALLOCATE(corrected_total_density)
   DEALLOCATE(total_pseudodensity)
   !CALL check_atomic_reference_polarizabilities_available()
   !IF  (atomic_reference_polarizabilities_available) THEN
   !   CALL calculate_atomic_polarizabilities_upper_bound() 
   !   CALL generate_atomic_polarizability_upper_bound_file()
   !END IF
   CALL print_overlap_populations()
 END IF

 IF (print_atomic_densities) THEN
   CALL print_atomic_densities_file()
 END IF
 
 CALL quote()
 CALL date_and_time(date,time,zone,values)
 CALL date_and_time(DATE=date,ZONE=zone)
 CALL date_and_time(TIME=time)
 CALL date_and_time(VALUES=values)
 WRITE(output_FID,*) date(1:4),"/",date(5:6),"/",date(7:),"  ",time(1:2),":",time(3:4),":",time(5:6)
 FLUSH(output_FID)
   
 WRITE(output_FID,*)' '
 WRITE(output_FID,*)'******************** TIME CONTROL ********************'
 CALL system_clock(clock_count2,clock_count_rate,clock_count_max)
 IF (clock_count_rate == 0) THEN
   WRITE(output_FID,*)'Processor clock not available, the time calculation wil not be performed.'
 ELSE
   seconds = modulo((clock_count2-clock_count),clock_count_max)/real(clock_count_rate)
   WRITE(output_FID,*)'Finished chargemol in ',seconds, 'seconds'
   FLUSH(output_FID)
 END IF
 WRITE(output_FID,*)'******************************************************'
 WRITE(output_FID,*)' '
 FLUSH(output_FID) 
 
 CLOSE(output_FID)
 
 END PROGRAM chargemol