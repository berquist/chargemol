 MODULE module_print_atomic_densities_file
 ! Copyright (c) 2014, 2015, 2016, 2017 by Thomas A. Manz and Nidia Gabaldon Limas. Rights reserved.
 
 USE module_precision
 USE module_global_parameters
 USE module_common_variable_declarations
 USE module_string_utilities
 USE module_atomic_number_to_symbol

 IMPLICIT NONE

 CONTAINS
 
 SUBROUTINE print_atomic_densities_file()
 !===================================================================================
 ! Print the spherically averaged density and the sum points for each atom
 !=================================================================================== 
 
 CHARACTER(41) :: filename 
 INTEGER :: FID
 REAL(kind=dp) :: V1(3),V2(3),V3(3),temp_array(nshells,natoms),r_inner,r_outer
 
 !Generate the output file
 filename = 'spherically_averaged_atomic_densities.txt'
 OPEN(NEWUNIT=FID, FILE=filename,STATUS='REPLACE')
 WRITE(FID,'(i5,a)')natoms,'  atoms'
 V1 = vector1*periA/bohrperangstrom
 V2 = vector2*periB/bohrperangstrom
 V3 = vector3*periC/bohrperangstrom
 IF ((periodicA .or. periodicB) .or. periodicC) THEN
   WRITE(FID,'(a,3f13.6,a,3f13.6,a,3f13.6,a)') '[{ ',V1(1),V1(2),V1(3),' }, { ',V2(1),V2(2),V2(3),' }, { ',V3(1),V3(2),V3(3),' }]'    
   FLUSH(FID)
 ELSE
   WRITE(FID,'(a)') 'Nonperiodic system'
   FLUSH(FID)
 END IF
 WRITE(FID,'(i5,a)') nshells,'  radial shells'
 WRITE(FID,'(a)')'Spherically averaged atomic densities in atomic units printed below.'
 DO i = 1,natoms
    WRITE(FID,'(10ES15.8)')(corrected_spherical_avg_density(j,i),j=1,nshells)
    WRITE(FID,*)' '
 END DO
 FLUSH(FID)
 WRITE(FID,*)' '
 WRITE(FID,'(a)')'Radial shell integration weights (summed grid point volumes) in atomic units printed below.'
 temp_array = pixelvolume*sum_points
 DO i = 1,natoms
    WRITE(FID,'(10ES15.8)')(temp_array(j,i),j=1,nshells)
    WRITE(FID,*)' '
 END DO
 WRITE(FID,*)' '
 FLUSH(FID)
 ALLOCATE(shell_analytic_radial_moments(nshells,6))
 DO i=1,nshells
    DO j=1,6
       r_inner = (i-1)/scalefactor 
       r_outer = i/scalefactor
       shell_analytic_radial_moments(i,j)=4.0_dp*pi*(r_outer**(1+j) - r_inner**(1+j))/(1.0_dp + j)
    END DO
 END DO
 WRITE(FID,'(a)') 'Average analytic radial shell moments (m = -1 to 4) in atomic units printed below.'
 DO j=1,6
    WRITE(FID,'(10ES15.8)')(shell_analytic_radial_moments(i,j),i=1,nshells)
    WRITE(FID,*)' '
 END DO  
 FLUSH(FID)

 CLOSE(FID)
  
 END SUBROUTINE print_atomic_densities_file
 END MODULE module_print_atomic_densities_file