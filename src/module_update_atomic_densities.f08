 MODULE module_update_atomic_densities
 ! Copyright (c) 2014, 2015, 2016, 2017 by Thomas A. Manz and Nidia Gabaldon Limas. Rights reserved.
 
 USE module_precision
 USE module_global_parameters
 USE module_common_variable_declarations
 USE module_oxidation_density

 IMPLICIT NONE
 
 CONTAINS
 
 SUBROUTINE update_atomic_densities()
 !===================================================================================
 
 REAL(kind=dp), ALLOCATABLE, DIMENSION (:,:) :: lower_oxidation_density,upper_oxidation_density
 REAL(kind=dp), ALLOCATABLE, DIMENSION (:) :: temp_oxidation_density
 REAL(kind=dp) :: f
 INTEGER :: k,atom_charge
 
 !Define boundaries for reference ions
 available_reference_ion_range(:,1) = [-2,-2,-2,-2,-2,-5,-4,-3,-2,-2,-2,-2,-2,-5,-4,-3,-3,-2,-2,-2,-2,-2,-2,-3,-4,-3,-2,&
 -2,-2,-2,-2,-5,-4,-3,-2,-2,-2,-2,-2,-2,-2,-3,-4,-3,-2,-2,-2,-2,-2,-5,-4,-3,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,&
 -2,-2,-2,-2,-2,-3,-4,-3,-4,-3,-2,-2,-2,-5,-4,-3,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2]
 available_reference_ion_range(:,2) = [1,2,2,3,4,5,6,3,2,2,2,3,4,5,6,7,8,2,2,3,4,5,6,7,8,7,6,5,5,3,4,5,6,7,8,3,2,3,4,5,6,7,8,9,7,&
 7,5,3,4,5,6,7,8,9,2,3,4,5,5,5,4,4,4,4,5,5,4,4,4,4,4,5,6,7,8,9,9,7,6,5,4,5,6,7,8,7,2,3,4,5,6,7,8,9,8,9,5,5,5,4,4,4,4,5,6,7,8,9,9]
 
 !Create the oxidation state reference densities
 ALLOCATE(lower_oxidation_density(nshells,natoms))
 ALLOCATE(upper_oxidation_density(nshells,natoms))
 ALLOCATE(temp_oxidation_density(nshells))
 lower_oxidation_density=0.0_dp
 upper_oxidation_density=0.0_dp
 temp_oxidation_density=0.0_dp
 !Read the lower and upper state reference densities
 DO j = 1,natoms
   IF ((reference_ion_charge(j) > available_reference_ion_range(atomic_number(j),1)) .and. (reference_ion_charge(j) < &
   available_reference_ion_range(atomic_number(j),2))) THEN
     atom_charge = ceiling(reference_ion_charge(j))
     temp_oxidation_density = oxidation_density(output_FID,atomic_densities_directory,density_set_prefix,atomic_number(j),&
     atom_charge,cutoff_radius,nshells,radial_shell_volume)
     lower_oxidation_density(:,j) = temp_oxidation_density(:)
     atom_charge = atom_charge - 1
     temp_oxidation_density = oxidation_density(output_FID,atomic_densities_directory,density_set_prefix,atomic_number(j),&
     atom_charge,cutoff_radius,nshells,radial_shell_volume)
     upper_oxidation_density(:,j) = temp_oxidation_density(:)
     f = ceiling(reference_ion_charge(j)) - reference_ion_charge(j)
     DO k = 1,nshells
       combined_oxidation_density(k,j) = upper_oxidation_density(k,j)*f + lower_oxidation_density(k,j)*(1.0_dp - f)
     END DO
   ELSE IF(reference_ion_charge(j) <= available_reference_ion_range(atomic_number(j),1)) THEN
     atom_charge = available_reference_ion_range(atomic_number(j),1)
     IF ( (atomic_number(j) - atom_charge) > zero_tolerance) THEN
       temp_oxidation_density = ((atomic_number(j) - reference_ion_charge(j))/(atomic_number(j) - atom_charge))*oxidation_density(&
       output_FID,atomic_densities_directory,density_set_prefix,atomic_number(j),atom_charge,cutoff_radius,nshells,&
       radial_shell_volume)
     ELSE
       temp_oxidation_density = 0.0_dp
     END IF
     combined_oxidation_density(:,j) = temp_oxidation_density(:)
   ELSE IF(reference_ion_charge(j) >= available_reference_ion_range(atomic_number(j),2)) THEN
     atom_charge = available_reference_ion_range(atomic_number(j),2)
     IF ( (atomic_number(j) - atom_charge) > zero_tolerance) THEN
       temp_oxidation_density = ((atomic_number(j) - reference_ion_charge(j))/(atomic_number(j) - atom_charge))*oxidation_density(&
       output_FID,atomic_densities_directory,density_set_prefix,atomic_number(j),atom_charge,cutoff_radius,nshells,&
       radial_shell_volume)
     ELSE
       temp_oxidation_density = 0.0_dp
     END IF
     combined_oxidation_density(:,j) = temp_oxidation_density(:)
   END IF   
 END DO   
   
 END SUBROUTINE update_atomic_densities
  
 END MODULE module_update_atomic_densities
