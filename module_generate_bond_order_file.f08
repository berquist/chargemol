 MODULE module_generate_bond_order_file
 ! Copyright (c) 2014, 2015, 2016, 2017 by Thomas A. Manz and Nidia Gabaldon Limas. Rights reserved.
 
 USE module_common_variable_declarations
 USE module_global_parameters
 USE module_atomic_number_to_symbol
 USE module_atomic_number_to_symbol

 IMPLICIT NONE
 
 CONTAINS
 
 SUBROUTINE generate_bond_order_file()
 !====================================================================================
 
 CHARACTER(2) :: atomic_symbol
 INTEGER :: BO_column,other_atom,translation_x,translation_y,translation_z,z
 REAL(kind=dp) :: V1(3),V2(3),V3(3),temp_scalar
 
 
 !Generate a file summarizing the effective bond orders
 V1 = vector1*periA/bohrperangstrom
 V2 = vector2*periB/bohrperangstrom
 V3 = vector3*periC/bohrperangstrom
 WRITE(BO_FID,'(i5)') natoms
 IF ((periodicA) .or. (periodicB) .or. (periodicC)) THEN
   WRITE(BO_FID,'(a,3f13.6,a,3f13.6,a,3f13.6,a)') 'jmolscript: load "" {1 1 1} &
   &spacegroup "x,y,z" unitcell [{ ',V1(1),V1(2),V1(3),' }, { ',V2(1),V2(2),V2(3),' }, { ',V3(1),V3(2),V3(3),' }]'
   FLUSH(BO_FID)
 ELSE    
   WRITE(BO_FID,*)'nonperiodic system'
 END IF
 DO i = 1,natoms
   z = nint(final_result(i,2))
   CALL atomic_number_to_symbol(z,output_FID,atomic_symbol)
   WRITE(BO_FID,'(2a,f13.6,a,f13.6,a,f13.6,f15.6)')atomic_symbol,' ',final_result&
   (i,3)/bohrperangstrom,' ',final_result(i,4)/bohrperangstrom,' ',final_result(i,5)/bohrperangstrom,atomic_summed_BO(i)
 END DO
 !Print information at the bottom of the file
 WRITE(BO_FID,*)'   '
 WRITE(BO_FID,*)version
 WRITE(BO_FID,*)'See ddec.sourceforge.net for latest version.' 
 WRITE(BO_FID,*)'   '
 WRITE(BO_FID,*)'   '
 WRITE(BO_FID,*)'The sum of bond orders (SBO) for each atom in the unit &
 &cell are listed above.'
 WRITE(BO_FID,'(a,f12.6,a)')' All bond orders greater than ',BO_print_cutoff,' are printed below.'
 BO_column = 20 
 DO i = 1,natoms
   WRITE(BO_FID,*) ' '
   WRITE(BO_FID,*) ' '
   WRITE(BO_FID,*) '===============================================================&
   &===================='
   z = atomic_number(i)
   CALL atomic_number_to_symbol(z,output_FID,atomic_symbol)
   WRITE(BO_FID,'(a,i6,3a)')' Printing BOs for ATOM # ',i,' ( ',atomic_symbol,&
   ' ) in the reference unit cell.'
   WRITE(BO_FID,*)' '
   DO current_pair = 1,num_sym_unique_bond_pairs
     IF (bond_pair_matrix(BO_column,current_pair) < BO_print_cutoff) THEN
       CYCLE
     END IF
     IF (nint(bond_pair_matrix(1,current_pair)) == i) THEN
       other_atom = nint(bond_pair_matrix(2,current_pair))
       z = atomic_number(other_atom)
       CALL atomic_number_to_symbol(z,output_FID,atomic_symbol)
       translation_x = nint(bond_pair_matrix(3,current_pair))
       translation_y = nint(bond_pair_matrix(4,current_pair))
       translation_z = nint(bond_pair_matrix(5,current_pair))
       WRITE(BO_FID,'(3(a,i3),a,i5,3a,f10.4,a,f10.4)')' Bonded to the (',translation_x,', ',&
       translation_y,', ',translation_z,') translated image of atom number ',other_atom,' ( ',atomic_symbol, &
       ' ) with bond order = ', bond_pair_matrix(BO_column,current_pair), &
       '    The average spin polarization of this bonding = ',&
       sqrt(  max( (bond_pair_matrix(13,current_pair)/bond_pair_matrix(12,current_pair) - 1.00), 0.00_dp )  )
     END IF
     IF (nint(bond_pair_matrix(2,current_pair)) == i) THEN
       other_atom = nint(bond_pair_matrix(1,current_pair))
       z = atomic_number(other_atom)
       CALL atomic_number_to_symbol(z,output_FID,atomic_symbol)
       translation_x = -(nint(bond_pair_matrix(3,current_pair)))
       translation_y = -(nint(bond_pair_matrix(4,current_pair)))
       translation_z = -(nint(bond_pair_matrix(5,current_pair)))
       WRITE(BO_FID,'(3(a,i3),a,i5,3a,f10.4,a,f10.4)')' Bonded to the (',translation_x,', ',&
       translation_y,', ',translation_z,') translated image of atom number ',other_atom,' ( ',atomic_symbol, &
       ' ) with bond order = ', bond_pair_matrix(BO_column,current_pair), &
       '    The average spin polarization of this bonding = ',&
       sqrt(  max( (bond_pair_matrix(13,current_pair)/bond_pair_matrix(12,current_pair) - 1.00), 0.00_dp )  )
     END IF
   END DO
   WRITE(BO_FID,'(a,f13.6)')' The sum of bond orders for this atom is S&
   &BO =   ',atomic_summed_BO(i)
 END DO
 WRITE(BO_FID,*)'   '
 CALL date_and_time(date,time,zone,values)
 CALL date_and_time(DATE=date,ZONE=zone)
 CALL date_and_time(TIME=time)
 CALL date_and_time(VALUES=values)
 WRITE(BO_FID,*) date(1:4),"/",date(5:6),"/",date(7:),"  ",time(1:2),":",time(3:4),&
 ":",time(5:6)

 END SUBROUTINE generate_bond_order_file
 
 END MODULE module_generate_bond_order_file
