 MODULE module_generate_net_atomic_charge_file
 ! Copyright (c) 2014, 2015, 2016, 2017 by Thomas A. Manz and Nidia Gabaldon Limas. Rights reserved.
 
 USE module_precision
 USE module_global_parameters
 USE module_common_variable_declarations
 USE module_atomic_number_to_symbol
 USE module_string_utilities 

 IMPLICIT NONE
 
 CONTAINS
 SUBROUTINE generate_net_atomic_charge_file()
 !===================================================================================
 ! Generating output file.
 !===================================================================================
  
 CHARACTER(50) :: filename
 CHARACTER(2) :: atomic_symbol
 INTEGER :: charge_fid,z
 REAL(kind=dp) :: V1(3),V2(3),V3(3)
 
 V1 = vector1*periA/bohrperangstrom
 V2 = vector2*periB/bohrperangstrom
 V3 = vector3*periC/bohrperangstrom
  IF (StrLowCase(trim(adjustl(charge_type))) == 'ddec3') THEN
   filename = 'DDEC3_net_atomic_charges.xyz'
 ELSE
   filename = 'DDEC6_even_tempered_net_atomic_charges.xyz'
 END IF 
 OPEN(NEWUNIT=charge_fid, FILE=filename,STATUS='REPLACE')
 WRITE(charge_fid,'(i5)') natoms
 FLUSH(charge_fid)
 IF ((periodicA .or. periodicB) .or. periodicC) THEN
   WRITE(charge_fid,'(a,3f13.6,a,3f13.6,a,3f13.6,a)') 'jmolscript: load "" {1 1 1} spacegroup "x,y,z" unitcell [{ ',V1(1),V1(2),&
   V1(3),' }, { ',V2(1),V2(2),V2(3),' }, { ',V3(1),V3(2),V3(3),' }]'    
   FLUSH(charge_fid)
 ELSE
   WRITE(charge_fid,'(a)') 'Nonperiodic system'
   FLUSH(charge_fid)
 END IF
 DO i = 1,natoms
   z = nint(final_result(i,2))
   CALL atomic_number_to_symbol(z,output_fid,atomic_symbol)
   WRITE(charge_fid,'(2a,f13.6,a,f13.6,a,f13.6,a,f13.6)') atomic_symbol,'',final_result(i,3)/bohrperangstrom,'',final_result(i,4)&
   /bohrperangstrom,'',final_result(i,5)/bohrperangstrom,'',final_result(i,6)
   FLUSH(charge_fid)
 END DO
 
 !Print information at the bottom of the file
 WRITE(charge_fid,*)''
 WRITE(charge_fid,*) version
 FLUSH(charge_fid)
 WRITE(charge_fid,'(a)')'See ddec.sourceforge.net for latest version.'
 FLUSH(charge_fid)
 WRITE(charge_fid,*)''
 WRITE(charge_fid,'(a)')'Computational parameters:'
 FLUSH(charge_fid)
 IF (StrLowCase(trim(adjustl(charge_type))) == 'ddec3') THEN
    WRITE(charge_fid,'(a,f13.6)')'Reference_weighting = ',reference_weighting
    FLUSH(charge_fid)
 END IF
 WRITE(charge_fid,'(a,i6)')'Number of radial integration shells = ',nshells
 FLUSH(charge_fid)
 WRITE(charge_fid,'(a,i6)')'Cutoff radius (pm) = ',nint(cutoff_radius)
 FLUSH(charge_fid)
 WRITE(charge_fid,'(a,es13.4)')'Error in the integrated total number of electrons before renormalization (e) = ',checkme
 FLUSH(charge_fid)
 WRITE(charge_fid,'(a,es13.4)')'Charge convergence tolerance = ', charge_convergence_tolerance
 FLUSH(charge_fid)
 WRITE(charge_fid,'(a,i6)')'Minimum radius for electron cloud penetration fitting (pm) = ', nint(rmin_cloud_penetration)
 FLUSH(charge_fid)
 WRITE(charge_fid,'(a,f8.4)')'Minimum decay exponent for electron density of buried atom tails = ',minimum_buried_tail_exponent
 FLUSH(charge_fid)
 IF (StrLowCase(trim(adjustl(charge_type))) == 'ddec6') THEN
   WRITE(charge_fid,'(a,f8.4)') 'Maximum decay exponent for electron density of buried atom tails =',maximum_buried_tail_exponent
   FLUSH(charge_fid)
 END IF
 WRITE(charge_fid,'(a,i6)')'Number of iterations to convergence = ',iter
 FLUSH(charge_fid)
 WRITE(charge_fid,*)''
 WRITE(charge_fid,'(a)')'The following XYZ coordinates are in angstroms. The atomic dipoles and quadrupoles are in atomic units.'
 FLUSH(charge_fid)
 WRITE(charge_fid,'(a)')'atom number, atomic symbol, x, y, z, net_charge, dipole_x, dipole_y, dipole_z, dipole_mag, Qxy, Qxz, Qyz,&
 &Q(x^2-y^2), Q(3z^2 - R^2), three eigenvalues of traceless quadrupole moment tensor'
 FLUSH(charge_fid)
 DO i = 1,natoms
   z = nint(final_result(i,2))
   CALL atomic_number_to_symbol(z,output_fid,atomic_symbol)
   WRITE(charge_fid,'(i5,a,a,16(a,f13.6))') i,' ',atomic_symbol,'',final_result(i,3)/bohrperangstrom,'',final_result(i,4)/&
   bohrperangstrom,'',final_result(i,5)/bohrperangstrom,'',final_result(i,6),'',final_result(i,7),'',final_result(i,8),'',&
   final_result(i,9),'',final_result(i,10),'',final_result(i,11),'',final_result(i,12),'',final_result(i,13),'',final_result&
   (i,14),'',final_result(i,15),'',final_result(i,16),'',final_result(i,17),'',final_result(i,18)
   FLUSH(charge_fid)
 END DO
 WRITE(charge_fid,*)''
 WRITE(charge_fid,'(a)')'The sperically averaged electron density of each atom fit to a function of the form exp(a - br) for r >=&
 &rmin_cloud_penetration'
 FLUSH(charge_fid)
 WRITE(charge_fid,'(a)')'atom number,atomic symbol, x, y, z, a, b, Rsquared where a and b are in atomic units and Rsquared is the &
 &squared correlation coefficient'
 FLUSH(charge_fid)
 DO i = 1,natoms
   z = nint(final_result(i,2))
   CALL atomic_number_to_symbol(z,output_fid,atomic_symbol)
   WRITE(charge_fid,'(i5,a,a,6(a,f13.6))')i,' ',atomic_symbol,'',final_result(i,3)/bohrperangstrom,'',final_result(i,4)/&
   bohrperangstrom,'',final_result(i,5)/bohrperangstrom,'',fitted_tail_intercept(i),'',-fitted_tail_slope(i),'',&
   fitted_tail_Rsquared(i)
   FLUSH(charge_fid)
 END DO
 WRITE(charge_fid,*)''
 IF((.not. periodicA) .and. (.not. periodicB) .and. (.not. periodicC)) THEN
   WRITE(charge_fid,*)' '
   WRITE(charge_fid,*)'Since the system is non-periodic, the total dipole and quadrupole moments were computed relative to the &
   &system center of mass.'
   WRITE(charge_fid,*)'Multipole moments (in atomic units) using the net atomic charges:'
   WRITE(charge_fid,'(4(a,f9.4))')' Partial_charge_dipole_x = ',partial_charge_dipole_x,'   Partial_charge_dipole_y = ',parti&
   &al_charge_dipole_y, '   Partial_charge_dipole_z = ',partial_charge_dipole_z,'     Partial_charge_dipole_magnitude = ',&
   partial_charge_dipole_magnitude
   WRITE(charge_fid,'(5(a,f9.4))')' Partial_charge_quadrupole_xy = ',partial_charge_quadrupole_xy,'   Partial_charge_quadrupol&
   &e_xz = ',partial_charge_quadrupole_xz,'   Partial_charge_quadrupole_yz = ',partial_charge_quadrupole_yz,'   &
   &Partial_charge_quadrupole_x2minusy2 = ',partial_charge_quadrupole_x2minusy2,'   Partial_charge_quadrupole_3z2minusr2 = ',pa&
   &rtial_charge_quadrupole_3z2minusr2
   WRITE(charge_fid,'(a,3f9.4)')' Eigenvalues of traceless quadrupole moment tensor = ',partial_charge_quadrupole_eigenvals
   WRITE(charge_fid,'(a)')' Multipole moments (in atomic units) using the net atomic charges and atomic dipoles:'
   WRITE(charge_fid,'(4(a,f9.4))')' Dipole_x = ', dipole_x,'   Dipole_y = ', dipole_y,'   Dipole_z = ', dipole_z,'   Dipole_mag&
   &nitude = ',dipole_magnitude
   WRITE(charge_fid,'(5(a,f9.4))')' Quadrupole_xy = ',quadrupole_xy,'   Quadrupole_xz = ',quadrupole_xz,'   Quadrupole_yz = ',q&
   &uadrupole_yz,'   Quadrupole_x2minusy2 = ',quadrupole_x2minusy2,'   Quadrupole_3z2minusr2 = ', quadrupole_3z2minusr2
   WRITE(charge_fid,'(a,3f9.4)')' Eigenvalues of traceless quadrupole moment tensor = ',quadrupole_eigenvals
   WRITE(charge_fid,'(a)')' Multipole moments (in atomic units) using the net atomic charges, atomic dipoles, and atomic quadrupo&
   &les:'
   WRITE(charge_fid,'(4(a,f9.4))')' Dipole_x = ',dipole_x,'   Dipole_y = ',dipole_y,'   Dipole_z = ',dipole_z,'   Dipole_magnitud&
   &e = ',dipole_magnitude
   WRITE(charge_fid,'(5(a,f9.4))')' Full_quadrupole_xy = ',full_quadrupole_xy,'   Full_quadrupole_xz = ',full_quadrupole_xz,'   F&
   &ull_quadrupole_yz = ',full_quadrupole_yz,'   Full_quadrupole_x2minusy2 = ',full_quadrupole_x2minusy2,'   Full_quadrupole_3z2m&
   &inusr2 = ',full_quadrupole_3z2minusr2
   WRITE(charge_fid,'(a,3f9.4)')' Eigenvalues of traceless quadrupole moment tensor = ',full_quadrupole_eigenvals
   WRITE(charge_fid,*)' '
 END IF
 CALL date_and_time(date,time,zone,values)
 CALL date_and_time(DATE=date,ZONE=zone)
 CALL date_and_time(TIME=time)
 CALL date_and_time(VALUES=values)

 WRITE(charge_fid,*) date(1:4),"/",date(5:6),"/",date(7:),"  ",time(1:2),":",time(3:4),":",time(5:6)

 CLOSE (charge_fid)
 
 END SUBROUTINE generate_net_atomic_charge_file
 
 END MODULE module_generate_net_atomic_charge_file