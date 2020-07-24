 MODULE module_check_atomic_reference_polarizabilities_available
 ! Copyright (c) 2014, 2015, 2016, 2017 by Thomas A. Manz and Nidia Gabaldon Limas. Rights reserved.
   
 USE module_precision
 USE module_common_variable_declarations
 USE module_global_parameters

 IMPLICIT NONE
 
 CONTAINS
 
 SUBROUTINE check_atomic_reference_polarizabilities_available()
 !===================================================================================
 !Making sure the reference densities are defined for atoms in the unit cell
 !===================================================================================
  
 INTEGER :: j

 CCSD_free_atom_polarizability_to_volume_ratio=[0.600000_dp,0.349889_dp,1.769052_dp,0.660831_dp,0.442599_dp,0.340437_dp,&
 0.279044_dp,0.234137_dp,0.199836_dp,0.174217_dp,1.423495_dp,0.687820_dp,0.492570_dp,0.369533_dp,0.296449_dp,0.257363_dp,&
 0.223031_dp,0.195699_dp,1.396956_dp,0.730674_dp,0.651409_dp,0.505595_dp,0.560241_dp,0.790767_dp,0.761859_dp,0.696586_dp,&
 0.683327_dp,0.673847_dp,0.492783_dp,0.484293_dp,0.443948_dp,0.344594_dp,0.279184_dp,0.244589_dp,0.214540_dp,0.189858_dp,&
 1.268292_dp,0.713673_dp,0.623889_dp,0.737771_dp,0.620816_dp,0.562386_dp,0.455743_dp,0.436590_dp,0.488037_dp,0.271461_dp,&
 0.486627_dp,0.390891_dp,0.402023_dp,0.320494_dp,0.263435_dp,0.234012_dp,0.206669_dp,0.184463_dp,1.192351_dp,0.713350_dp,&
 0.608058_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.485300_dp,&
 0.428533_dp,0.466757_dp,0.377559_dp,0.347274_dp,0.323348_dp,0.304531_dp,0.289975_dp,0.276755_dp,0.381399_dp,0.314427_dp,&
 0.261473_dp,0.234052_dp,0.207392_dp,0.185688_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,&
 0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp]

 atomic_reference_polarizabilities_available = .true.
 DO j=1,natoms
   IF (abs(CCSD_free_atom_polarizability_to_volume_ratio(atomic_number(j))) < zero_tolerance) THEN
     atomic_reference_polarizabilities_available = .false.
     WRITE(output_FID,*) 'There is no CCSD free atom polarizability to volume ratio reference for at least one element in &
     &your material.'
     WRITE(output_FID,*)'Program will skip the calculation of atomic polarizability upper bounds.'
     FLUSH(output_FID)
     RETURN
   END IF
 END DO
         
 END SUBROUTINE check_atomic_reference_polarizabilities_available 

 END MODULE module_check_atomic_reference_polarizabilities_available