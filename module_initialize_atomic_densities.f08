 MODULE module_initialize_atomic_densities
 ! Copyright (c) 2014, 2015, 2016, 2017 by Thomas A. Manz and Nidia Gabaldon Limas. Rights reserved.
 
 USE module_precision
 USE module_global_parameters
 USE module_common_variable_declarations
 USE module_string_utilities
 USE module_oxidation_density

 IMPLICIT NONE
 
 CONTAINS
 
 SUBROUTINE initialize_atomic_densities()
 !===================================================================================
 
 CHARACTER (200) :: combinedstring,line_of_text
 INTEGER :: combinedstring_FID,io,k,atom_charge
 REAL(kind=dp),ALLOCATABLE,DIMENSION(:) :: temp_oxidation_density

 !Initialize some atomic density distributions 
 ALLOCATE(valence_population(natoms))
 ALLOCATE(core_population(natoms))
 ALLOCATE(spherical_average_density(nshells,natoms))
 !ALLOCATE(spherical_avg_ref_pseudodensity(natoms,nshells))
 ALLOCATE(avg_PtoWref(nshells,natoms))
 ALLOCATE(final_result(natoms,18))
 ALLOCATE(total_pseudodensity(totnumA,totnumB,totnumC))
 ALLOCATE(change(natoms))
 ALLOCATE(combined_oxidation_density(nshells,natoms))
 ALLOCATE(conditioned_oxidation_density(nshells,natoms))
 ALLOCATE(partial_core_density(nshells,natoms))
 ALLOCATE(neutral_density(nshells,natoms))
 ALLOCATE(temp_oxidation_density(natoms))
 avg_PtoWref = 1.0_dp
 max_change = 1.0_dp
 old_change = 1.0_dp
 old_old_change = 1.0_dp
 core_population=0.0_dp
 spherical_average_density=0.0_dp
 valence_population=0.0_dp
 DO j = 1,natoms
   atom_charge = 0
   temp_oxidation_density = oxidation_density(output_FID,atomic_densities_directory,density_set_prefix,atomic_number(j),&
   atom_charge,cutoff_radius,nshells,radial_shell_volume)
   neutral_density(:,j) = temp_oxidation_density(:)
 END DO  

 combined_oxidation_density=neutral_density
 partial_density=neutral_density
 conditioned_oxidation_density=neutral_density
 partial_core_density=neutral_density


 END SUBROUTINE initialize_atomic_densities
 
 END MODULE module_initialize_atomic_densities