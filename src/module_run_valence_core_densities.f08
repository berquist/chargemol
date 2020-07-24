 MODULE module_run_valence_core_densities
 ! Copyright (c) 2014, 2015, 2016, 2017 by Thomas A. Manz and Nidia Gabaldon Limas. Rights reserved.
 
 USE module_read_wfx
 USE module_format_vasp_densities
 USE module_gen_dens_grids_from_gaussian_basis_set_coeff
 USE module_check_noncollinear_XC_functional
 USE module_align_periodic_atoms
 USE module_format_total_cube_density
 USE module_format_valence_and_total_cube_densities
 USE module_format_valence_cube_density
 USE module_format_xsf_densities

 IMPLICIT NONE

 CONTAINS
 
 SUBROUTINE run_valence_core_densities()

 REAL(kind=dp) :: sum_valence_occupancy_correction,atomic_number_sum

 IF (iostat_value > 0) THEN
   WRITE(output_FID,*) 'Could not read job_control.txt file. Program will terminate.'
   FLUSH(output_FID)
   STOP
 ELSE
   WRITE(output_FID,*) atomic_densities_directory
 END IF
 IF (input_type == 'vtc') THEN
   CALL format_valence_and_total_cube_densities( )
 ELSE IF (input_type == 'vsp') THEN
   CALL format_vasp_densities( )
 ELSE IF (input_type == 'xsf') THEN
   CALL format_xsf_densities( )
 ELSE IF (input_type == 'val') THEN
   CALL format_valence_cube_density( )
 ELSE IF (input_type == 'wfx') THEN
   CALL read_wfx( )
   CALL generate_density_grids_from_gaussian_basis_set_coefficients( )
 ELSE IF (input_type == 'cub') THEN
   CALL format_total_cube_density( )
 ELSE
   WRITE(output_FID,*)'Fatal error, wrong input file type. Program will terminate.'
   FLUSH(output_FID) 
   STOP
 END IF
 IF (periodicA) THEN
   periA = 1
 ELSE
   periA = 0
 END IF
 IF (periodicB) THEN
   periB = 1
 ELSE
   periB = 0
 END IF 
 IF (periodicC) THEN
   periC = 1
 ELSE
   periC = 0   
 END IF 
 num_periodic_directions = periA + periB + periC
 ncore = 0.0_dp
 DO j=1,natoms
   ncore = ncore+core_electrons(j)
 END DO
 WRITE(output_FID,'(a,f13.4)')" ncore = ",ncore
 FLUSH(output_FID)
 atomic_number_sum = 0.0_dp
 DO j=1,natoms
   atomic_number_sum = atomic_number_sum + atomic_number(j)
 END DO
 nvalence = atomic_number_sum - netcharge - ncore
 WRITE(output_FID,'(a,f13.4)')" nvalence = ",nvalence
 FLUSH(output_FID)
 int_tolerance = max(integration_tolerance,integration_tolerance_percent*nvalence/100.0_dp)
 sum_valence_occupancy_correction = 0.0_dp
 DO j=1,natoms
    sum_valence_occupancy_correction = sum_valence_occupancy_correction + occupancy_correction(1,j)
 END DO
 checkme = abs(nvalence - sum(valence_density(:,:,:))*pixelvolume - sum_valence_occupancy_correction)
 WRITE(output_FID,'(a,es13.4)')" pixelvolume = ",pixelvolume
 WRITE(output_FID,'(a,es13.4)')" numerically integrated valence density = ",sum(valence_density(:,:,:))*pixelvolume
 WRITE(output_FID,'(a,es13.4)')" sum_valence_occupancy_correction = ",sum_valence_occupancy_correction
 WRITE(output_FID,'(a,es13.4)')" checkme = ",checkme
 FLUSH(output_FID)
 
 IF ((checkme > int_tolerance) .or. (abs(sum_negative_density) > int_tolerance)) THEN
   WRITE(output_FID,*)'The electrons are not properly accounted for.'
   FLUSH(output_FID)
   WRITE(output_FID,*)'Either the grid in your electron density input &
   &file is too coarse, you have specified the incorrect net charge in &
   &the chargemol_job.m file, or the number of core electrons for each &
   &atom has not been setup correctly.'
   FLUSH(output_FID)
   WRITE(output_FID,*)'Program will terminate'
   FLUSH(output_FID)
   STOP
 ELSE
   WRITE(output_FID,*)'The grid spacing is adequate and all electrons a&
   &re properly accounted for.'
   FLUSH(output_FID)
   WRITE(output_FID,*)'Calculation will proceed.'
   FLUSH(output_FID)
 END IF
 IF (spin_available) THEN 
   IF(non_collinear) THEN
     CALL check_noncollinear_XC_functional( )
   END IF
 END IF
 
 !Make sure the atoms for periodic directions are inside the unit cell
 CALL align_periodic_atoms( )
 
 END SUBROUTINE run_valence_core_densities
 
 END MODULE module_run_valence_core_densities