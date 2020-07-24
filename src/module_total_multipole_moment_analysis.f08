 MODULE module_total_multipole_moment_analysis
 ! Copyright (c) 2014, 2015, 2016, 2017 by Thomas A. Manz and Nidia Gabaldon Limas. Rights reserved.
 
 USE module_precision
 USE module_global_parameters
 USE module_common_variable_declarations
 USE module_compute_center_of_mass
 USE module_matrix_operations

 IMPLICIT NONE
 
 CONTAINS
 
 SUBROUTINE total_multipole_moment_analysis()
 !=================================================================================== 
 ! The multipole moments are computed in two parts: an analytic part which contains 
 ! the atomic nuclei and the occupancy corrections and a numeric part which contains 
 ! the multipole moments of the valence density grid
 !===================================================================================  
 
 !REAL(kind=dp) :: partial_charge_dipole_x,partial_charge_dipole_y,partial_charge_dipole_z,partial_charge_quadrupole_xy,&
 !partial_charge_quadrupole_xz,partial_charge_quadrupole_yz,partial_charge_quadrupole_x2minusy2,&
 !partial_charge_quadrupole_3z2minusr2,dipole_x,dipole_y,dipole_z,quadrupole_xy,quadrupole_xz,&
 !quadrupole_yz,quadrupole_x2minusy2,quadrupole_3z2minusr2,dipole_magnitude,partial_charge_dipole_magnitude
 REAL(kind=dp) :: x,y,z,charge,temp_partial_charge_quadrupole_matrix(3,3),temp_quadrupole_matrix(3,3),&
 temp_full_quadrupole_matrix(3,3)
 
 IF (((.not. periodicA) .and. (.not. periodicB)) .and. (.not. periodicC)) THEN
   CALL compute_center_of_mass()
   WRITE(output_FID,*)'Since the system is nonperiodic, total multipole moment analysis will be performed.'
   FLUSH(output_FID)
   WRITE(output_FID,*)'Total multipole moments are computed using the center of mass as the origin.'
   FLUSH(output_FID)
   WRITE(output_FID,*)'Dipole and quadrupole moments using the net atomic charges.'
   FLUSH(output_FID)
   WRITE(output_FID,*)'This corresponds to truncatating the distributed multipole expansion at monopole order.'
   FLUSH(output_FID)
   partial_charge_dipole_x = 0.0_dp
   partial_charge_dipole_y = 0.0_dp
   partial_charge_dipole_z = 0.0_dp
   partial_charge_quadrupole_xy = 0.0_dp
   partial_charge_quadrupole_xz = 0.0_dp
   partial_charge_quadrupole_yz = 0.0_dp
   partial_charge_quadrupole_x2minusy2 = 0.0_dp
   partial_charge_quadrupole_3z2minusr2 = 0.0_dp
   DO j = 1,natoms
     x = coords(1,j) - center_of_mass(1)
     y = coords(2,j) - center_of_mass(2)
     z = coords(3,j) - center_of_mass(3)
     charge = final_result(j,6)
     partial_charge_dipole_x = partial_charge_dipole_x + charge*x
     partial_charge_dipole_y = partial_charge_dipole_y + charge*y
     partial_charge_dipole_z = partial_charge_dipole_z + charge*z
     partial_charge_quadrupole_xy = partial_charge_quadrupole_xy + charge*x*y
     partial_charge_quadrupole_xz = partial_charge_quadrupole_xz + charge*x*z
     partial_charge_quadrupole_yz = partial_charge_quadrupole_yz + charge*y*z
     partial_charge_quadrupole_x2minusy2 = partial_charge_quadrupole_x2minusy2 + charge*(x**2 - y**2)
     partial_charge_quadrupole_3z2minusr2 = partial_charge_quadrupole_3z2minusr2 + charge*(3*z*z - (x**2 + y**2 + z**2))
   END DO
   WRITE(output_FID,'(a,f13.6)')' partial_charge_dipole_x= ',partial_charge_dipole_x
   FLUSH(output_FID)
   WRITE(output_FID,'(a,f13.6)')' partial_charge_dipole_y= ',partial_charge_dipole_y
   FLUSH(output_FID)
   WRITE(output_FID,'(a,f13.6)')' partial_charge_dipole_z= ',partial_charge_dipole_z
   FLUSH(output_FID)
   partial_charge_dipole_magnitude = sqrt(partial_charge_dipole_x**2 + partial_charge_dipole_y**2 + partial_charge_dipole_z**2)
   WRITE(output_FID,'(a,f13.6)')' partial_charge_dipole_magnitude',partial_charge_dipole_magnitude
   FLUSH(output_FID)
   WRITE(output_FID,'(a,f13.6)')' partial_charge_quadrupole_xy= ',partial_charge_quadrupole_xy
   FLUSH(output_FID)
   WRITE(output_FID,'(a,f13.6)')' partial_charge_quadrupole_xz= ',partial_charge_quadrupole_xz
   FLUSH(output_FID)
   WRITE(output_FID,'(a,f13.6)')' partial_charge_quadrupole_yz= ',partial_charge_quadrupole_yz
   FLUSH(output_FID)
   WRITE(output_FID,'(a,f13.6)')' partial_charge_quadrupole_x2minusy2= ',partial_charge_quadrupole_x2minusy2
   FLUSH(output_FID)
   WRITE(output_FID,'(a,f13.6)')' partial_charge_quadrupole_3z2minusr2= ',partial_charge_quadrupole_3z2minusr2
   FLUSH(output_FID)
   temp_partial_charge_quadrupole_matrix(1,:)=[0.5_dp*(partial_charge_quadrupole_x2minusy2-partial_charge_quadrupole_3z2minusr2/&
   3.0_dp),partial_charge_quadrupole_xy,partial_charge_quadrupole_xz]
   temp_partial_charge_quadrupole_matrix(2,:)=[partial_charge_quadrupole_xy,0.5_dp*(-partial_charge_quadrupole_x2minusy2-&
   partial_charge_quadrupole_3z2minusr2/3.0_dp),partial_charge_quadrupole_yz]
   temp_partial_charge_quadrupole_matrix(3,:)=[partial_charge_quadrupole_xz,partial_charge_quadrupole_yz,&
   partial_charge_quadrupole_3z2minusr2/3.0_dp]
   partial_charge_quadrupole_eigenvals = eig(temp_partial_charge_quadrupole_matrix)
   WRITE(output_FID,'(a,f13.6)')' Dipole and quadrupole moments using the net atomic charges and the atomic dipoles.'
   FLUSH(output_FID)
   WRITE(output_FID,*)'This corresponds to truncatating the distributed multipole expansion at dipole order.'
   FLUSH(output_FID)
   dipole_x = partial_charge_dipole_x
   dipole_y = partial_charge_dipole_y
   dipole_z = partial_charge_dipole_z   
   quadrupole_xy = partial_charge_quadrupole_xy
   quadrupole_xz = partial_charge_quadrupole_xz
   quadrupole_yz = partial_charge_quadrupole_yz
   quadrupole_x2minusy2 = partial_charge_quadrupole_x2minusy2
   quadrupole_3z2minusr2 = partial_charge_quadrupole_3z2minusr2
   DO j = 1,natoms
     x = coords(1,j) - center_of_mass(1)
     y = coords(2,j) - center_of_mass(2)
     z = coords(3,j) - center_of_mass(3)
     dipole_x = dipole_x + final_result(j,7)
     dipole_y = dipole_y + final_result(j,8)
     dipole_z = dipole_z + final_result(j,9)       
     quadrupole_xy = quadrupole_xy + x*final_result(j,8) + y*final_result(j,7)
     quadrupole_xz = quadrupole_xz + x*final_result(j,9) + z*final_result(j,7)
     quadrupole_yz = quadrupole_yz + y*final_result(j,9) + z*final_result(j,8)
     quadrupole_x2minusy2 = quadrupole_x2minusy2 + 2.0_dp*(x*final_result(j,7) - y*final_result(j,8))
     quadrupole_3z2minusr2 = quadrupole_3z2minusr2 + 2.0_dp*(2.0_dp*z*final_result(j,9)-x*final_result(j,7)-y*final_result(j,8))
   END DO
   WRITE(output_FID,'(a,f13.6)')' dipole_x= ',dipole_x
   FLUSH(output_FID)
   WRITE(output_FID,'(a,f13.6)')' dipole_y= ',dipole_y
   FLUSH(output_FID)
   WRITE(output_FID,'(a,f13.6)')' dipole_z = ',dipole_z
   FLUSH(output_FID)
   dipole_magnitude = sqrt(dipole_x**2 + dipole_y**2 + dipole_z**2)
   WRITE(output_FID,'(a,f13.6)')' dipole_magnitude= ',dipole_magnitude
   FLUSH(output_FID)
   WRITE(output_FID,'(a,f13.6)')' quadrupole_xy= ',quadrupole_xy
   FLUSH(output_FID)
   WRITE(output_FID,'(a,f13.6)')' quadrupole_xz= ',quadrupole_xz
   FLUSH(output_FID)
   WRITE(output_FID,'(a,f13.6)')' quadrupole_yz= ',quadrupole_yz
   FLUSH(output_FID)
   WRITE(output_FID,'(a,f13.6)')' quadrupole_x2minusy2= ',quadrupole_x2minusy2
   FLUSH(output_FID)
   WRITE(output_FID,'(a,f13.6)')' quadrupole_3z2minusr2 = ',quadrupole_3z2minusr2
   FLUSH(output_FID)
   temp_quadrupole_matrix(1,:)=[0.5_dp*(quadrupole_x2minusy2-quadrupole_3z2minusr2/3.0_dp),quadrupole_xy,quadrupole_xz]
   temp_quadrupole_matrix(2,:)=[quadrupole_xy,0.5_dp*(-quadrupole_x2minusy2-quadrupole_3z2minusr2/3.0_dp),quadrupole_yz]
   temp_quadrupole_matrix(3,:)=[quadrupole_xz,quadrupole_yz,quadrupole_3z2minusr2/3.0_dp]
   quadrupole_eigenvals = eig(temp_quadrupole_matrix)
   WRITE(output_FID,*)'Dipole and quadrupole moments using the net atomic charges, atomic dipoles, and atomic &
   &quadrupoles.'
   FLUSH(output_FID)
   WRITE(output_FID,*)'This corresponds to truncatating the distributed multipole expansion at quadrupole order.'
   FLUSH(output_FID)
   WRITE(output_FID,'(a,f13.6)')' dipole_x= ',dipole_x
   FLUSH(output_FID)
   WRITE(output_FID,'(a,f13.6)')' dipole_y= ',dipole_y
   FLUSH(output_FID)
   WRITE(output_FID,'(a,f13.6)')' dipole_z= ',dipole_z
   FLUSH(output_FID)
   dipole_magnitude = sqrt(dipole_x**2 + dipole_y**2 + dipole_z**2)
   WRITE(output_FID,*) 'dipole_magnitude= ',dipole_magnitude
   FLUSH(output_FID)
   full_quadrupole_xy = quadrupole_xy
   full_quadrupole_xz = quadrupole_xz
   full_quadrupole_yz = quadrupole_yz
   full_quadrupole_x2minusy2 = quadrupole_x2minusy2
   full_quadrupole_3z2minusr2 = quadrupole_3z2minusr2    
   DO j = 1,natoms      
     full_quadrupole_xy = full_quadrupole_xy + final_result(j,11)
     full_quadrupole_xz = full_quadrupole_xz + final_result(j,12)
     full_quadrupole_yz = full_quadrupole_yz + final_result(j,13)
     full_quadrupole_x2minusy2 = full_quadrupole_x2minusy2 + final_result(j,14)
     full_quadrupole_3z2minusr2 = full_quadrupole_3z2minusr2 + final_result(j,15)     
   END DO  
   WRITE(output_FID,'(a,f13.6)')' full_quadrupole_xy= ',full_quadrupole_xy
   FLUSH(output_FID)
   WRITE(output_FID,'(a,f13.6)')' full_quadrupole_xz= ',full_quadrupole_xz
   FLUSH(output_FID)
   WRITE(output_FID,'(a,f13.6)')' full_quadrupole_yz= ',full_quadrupole_yz
   FLUSH(output_FID)
   WRITE(output_FID,'(a,f13.6)')' full_quadrupole_x2minusy2= ',full_quadrupole_x2minusy2
   FLUSH(output_FID)
   WRITE(output_FID,'(a,f13.6)')' full_quadrupole_3z2minusr2 = ',full_quadrupole_3z2minusr2
   FLUSH(output_FID)   
   temp_full_quadrupole_matrix(1,:)=[0.5_dp*(full_quadrupole_x2minusy2-full_quadrupole_3z2minusr2/3.0_dp),full_quadrupole_xy,&
   full_quadrupole_xz]
   temp_full_quadrupole_matrix(2,:)=[full_quadrupole_xy,0.5_dp*(-full_quadrupole_x2minusy2-full_quadrupole_3z2minusr2/3.0_dp),&
   full_quadrupole_yz]
   temp_full_quadrupole_matrix(3,:)=[full_quadrupole_xz,full_quadrupole_yz,full_quadrupole_3z2minusr2/3.0_dp]
   full_quadrupole_eigenvals = eig(temp_full_quadrupole_matrix)
 ELSE
   WRITE(output_FID,*)'Since the system is periodic, total multipole moment analysis will not be performed.'
 END IF
    
 END SUBROUTINE total_multipole_moment_analysis
 
 END MODULE module_total_multipole_moment_analysis