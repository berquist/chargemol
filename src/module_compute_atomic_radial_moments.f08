 MODULE module_compute_atomic_radial_moments
 ! Copyright (c) 2014, 2015, 2016, 2017 by Thomas A. Manz and Nidia Gabaldon Limas. Rights reserved.
 
 USE module_precision
 USE module_global_parameters
 USE module_common_variable_declarations
 USE module_string_utilities
 USE module_atomic_number_to_symbol
!$ USE omp_lib

 IMPLICIT NONE

 CONTAINS
 
 SUBROUTINE compute_atomic_radial_moments ()
 !===================================================================================
 ! Compute the Rcubed moment of each atom in the material
 !=================================================================================== 
 
 CHARACTER(2) :: atomic_symbol
 CHARACTER(41) :: filename
 INTEGER :: j,shell_index,z
 REAL(kind=dp) :: V1(3),V2(3),V3(3),distance,temp_vector(3),temp
 REAL(kind=dp),ALLOCATABLE,DIMENSION(:)::Rsquared_moment,Rcubed_moment,Rfourth_moment
 
 ALLOCATE(Rsquared_moment(natoms))
 ALLOCATE(Rcubed_moment(natoms))
 ALLOCATE(Rfourth_moment(natoms))
 Rsquared_moment = 0.0_dp
 Rcubed_moment = 0.0_dp
 Rfourth_moment = 0.0_dp
!$omp parallel default(none) &
!$omp private(j,na,nb,nc,ka,kb,kc,temp_vector,distance,shell_index,temp) &
!$omp shared(lower_nc,upper_nc,kCpoints,delta_nc,lower_nb,upper_nb,kBpoints,delta_nb,lower_na,upper_na,kApoints,delta_na,&
!$omp boundary,natoms,chunk_size,center_shift,partial_density,corrected_total_density,&
!$omp total_pseudodensity,pixelvolume) &
!$omp reduction (+:Rsquared_moment,Rcubed_moment,Rfourth_moment)
 DO j=1,natoms
!$omp do schedule(dynamic,chunk_size)
   DO na = lower_na(j),upper_na(j)
     ka = kApoints(j,delta_na + na + 1)
     DO nb = lower_nb(j),upper_nb(j)
       kb = kBpoints(j,delta_nb + nb + 1)
       DO nc = lower_nc(j),upper_nc(j)
         kc = kCpoints(j,delta_nc + nc + 1)
         temp_vector(1) = na*boundary(1,1) + nb*boundary(2,1) + nc*boundary(3,1) - center_shift(1,j)
         temp_vector(2) = na*boundary(1,2) + nb*boundary(2,2) + nc*boundary(3,2) - center_shift(2,j)
         temp_vector(3) = na*boundary(1,3) + nb*boundary(2,3) + nc*boundary(3,3) - center_shift(3,j)
         distance = sqrt(temp_vector(1)*temp_vector(1) + temp_vector(2)*temp_vector(2) + temp_vector(3)*temp_vector(3))
         shell_index = ceiling(scalefactor*distance + zero_tolerance)
         IF ((shell_index <= nshells) .and. (total_pseudodensity(ka,kb,kc) > zero_tolerance)) THEN
           temp = partial_density(shell_index,j)*pixelvolume*corrected_total_density(ka,kb,kc)/total_pseudodensity(ka,kb,kc)
           Rsquared_moment(j)=Rsquared_moment(j)+distance**2*temp
           Rcubed_moment(j)=Rcubed_moment(j)+distance**3*temp
           Rfourth_moment(j)=Rfourth_moment(j)+distance**4*temp
         END IF
       END DO    
     END DO
   END DO
 END DO
!$omp end parallel

 !Generate the file containing the Rsquared atomic moments
 WRITE(output_FID,*)'The computed Rsquared moments of the atoms (in bohr^2) are: '
 WRITE(output_FID,'(10f13.6)')Rsquared_moment
 FLUSH(output_FID)
 V1 = vector1*periA/bohrperangstrom
 V2 = vector2*periB/bohrperangstrom
 V3 = vector3*periC/bohrperangstrom
 filename = 'DDEC_atomic_Rsquared_moments.xyz'
 OPEN(NEWUNIT=radial_moment_FID, FILE=filename,STATUS='REPLACE')
 WRITE(radial_moment_FID,'(i5)') natoms
 FLUSH(radial_moment_FID)
 IF ((periodicA .or. periodicB) .or. periodicC) THEN
   WRITE(radial_moment_FID,'(a,3f13.6,a,3f13.6,a,3f13.6,a)') 'jmolscript: load "" {1 1 1} spacegroup "x,y,z" unitcell &
   &[{ ',V1(1),V1(2),V1(3),' }, { ',V2(1),V2(2),V2(3),' }, { ',V3(1),V3(2),V3(3),' }]'    
   FLUSH(radial_moment_FID)
 ELSE
   WRITE(radial_moment_FID,'(a)') 'Nonperiodic system'
   FLUSH(radial_moment_FID)
 END IF
 DO j = 1,natoms
   z = nint(final_result(j,2))
   CALL atomic_number_to_symbol(z,output_fid,atomic_symbol)
   WRITE(radial_moment_FID,'(2a,f13.6,a,f13.6,a,f13.6,a,f13.6)') atomic_symbol,' ',final_result(j,3)/bohrperangstrom,'',&
   final_result(j,4)/bohrperangstrom,'',final_result(j,5)/bohrperangstrom,'',Rsquared_moment(j)
   FLUSH(radial_moment_FID)
 END DO
 WRITE(radial_moment_FID,*)' '
 WRITE(radial_moment_FID,'(a)')'Same information as above printed with atom number.'
 DO j = 1,natoms
   z = nint(final_result(j,2))
   CALL atomic_number_to_symbol(z,output_fid,atomic_symbol)
   WRITE(radial_moment_FID,'(i4,3a,f13.6,a,f13.6,a,f13.6,a,f13.6)') j,' ',atomic_symbol,'',final_result(j,3)/bohrperangstrom,'',&
   final_result(j,4)/bohrperangstrom,'',final_result(j,5)/bohrperangstrom,'',Rsquared_moment(j)
   FLUSH(radial_moment_FID)
 END DO
 CLOSE(radial_moment_FID)

 !Generate the file containing the Rcubed atomic moments
 WRITE(output_FID,*)'The computed Rcubed moments of the atoms (in bohr^3) are: '
 WRITE(output_FID,'(10f13.6)')Rcubed_moment
 FLUSH(output_FID)
 V1 = vector1*periA/bohrperangstrom
 V2 = vector2*periB/bohrperangstrom
 V3 = vector3*periC/bohrperangstrom
 filename = 'DDEC_atomic_Rcubed_moments.xyz'
 OPEN(NEWUNIT=radial_moment_FID, FILE=filename,STATUS='REPLACE')
 WRITE(radial_moment_FID,'(i5)') natoms
 FLUSH(radial_moment_FID)
 IF ((periodicA .or. periodicB) .or. periodicC) THEN
   WRITE(radial_moment_FID,'(a,3f13.6,a,3f13.6,a,3f13.6,a)') 'jmolscript: load "" {1 1 1} spacegroup "x,y,z" unitcell &
   &[{ ',V1(1),V1(2),V1(3),' }, { ',V2(1),V2(2),V2(3),' }, { ',V3(1),V3(2),V3(3),' }]'    
   FLUSH(radial_moment_FID)
 ELSE
   WRITE(radial_moment_FID,'(a)') 'Nonperiodic system'
   FLUSH(radial_moment_FID)
 END IF
 DO j = 1,natoms
   z = nint(final_result(j,2))
   CALL atomic_number_to_symbol(z,output_fid,atomic_symbol)
   WRITE(radial_moment_FID,'(2a,f13.6,a,f13.6,a,f13.6,a,f13.6)') atomic_symbol,' ',final_result(j,3)/bohrperangstrom,'',&
   final_result(j,4)/bohrperangstrom,'',final_result(j,5)/bohrperangstrom,'',Rcubed_moment(j)
   FLUSH(radial_moment_FID)
 END DO
 WRITE(radial_moment_FID,*)' '
 WRITE(radial_moment_FID,'(a)')'Same information as above printed with atom number.'
 DO j = 1,natoms
   z = nint(final_result(j,2))
   CALL atomic_number_to_symbol(z,output_fid,atomic_symbol)
   WRITE(radial_moment_FID,'(i4,3a,f13.6,a,f13.6,a,f13.6,a,f13.6)') j,' ',atomic_symbol,'',final_result(j,3)/bohrperangstrom,'',&
   final_result(j,4)/bohrperangstrom,'',final_result(j,5)/bohrperangstrom,'',Rcubed_moment(j)
   FLUSH(radial_moment_FID)
 END DO
 CLOSE(radial_moment_FID)

 !Generate the file containing the Rfourth atomic moments
 WRITE(output_FID,*)'The computed Rfourth moments of the atoms (in bohr^4) are: '
 WRITE(output_FID,'(10f13.6)')Rfourth_moment
 FLUSH(output_FID)
 V1 = vector1*periA/bohrperangstrom
 V2 = vector2*periB/bohrperangstrom
 V3 = vector3*periC/bohrperangstrom
 filename = 'DDEC_atomic_Rfourth_moments.xyz'
 OPEN(NEWUNIT=radial_moment_FID, FILE=filename,STATUS='REPLACE')
 WRITE(radial_moment_FID,'(i5)') natoms
 FLUSH(radial_moment_FID)
 IF ((periodicA .or. periodicB) .or. periodicC) THEN
   WRITE(radial_moment_FID,'(a,3f13.6,a,3f13.6,a,3f13.6,a)') 'jmolscript: load "" {1 1 1} spacegroup "x,y,z" unitcell &
   &[{ ',V1(1),V1(2),V1(3),' }, { ',V2(1),V2(2),V2(3),' }, { ',V3(1),V3(2),V3(3),' }]'    
   FLUSH(radial_moment_FID)
 ELSE
   WRITE(radial_moment_FID,'(a)') 'Nonperiodic system'
   FLUSH(radial_moment_FID)
 END IF
 DO j = 1,natoms
   z = nint(final_result(j,2))
   CALL atomic_number_to_symbol(z,output_fid,atomic_symbol)
   WRITE(radial_moment_FID,'(2a,f13.6,a,f13.6,a,f13.6,a,f13.6)') atomic_symbol,' ',final_result(j,3)/bohrperangstrom,'',&
   final_result(j,4)/bohrperangstrom,'',final_result(j,5)/bohrperangstrom,'',Rfourth_moment(j)
   FLUSH(radial_moment_FID)
 END DO
 WRITE(radial_moment_FID,*)' '
 WRITE(radial_moment_FID,'(a)')'Same information as above printed with atom number.'
 DO j = 1,natoms
   z = nint(final_result(j,2))
   CALL atomic_number_to_symbol(z,output_fid,atomic_symbol)
   WRITE(radial_moment_FID,'(i4,3a,f13.6,a,f13.6,a,f13.6,a,f13.6)') j,' ',atomic_symbol,'',final_result(j,3)/bohrperangstrom,'',&
   final_result(j,4)/bohrperangstrom,'',final_result(j,5)/bohrperangstrom,'',Rfourth_moment(j)
   FLUSH(radial_moment_FID)
 END DO
 CLOSE(radial_moment_FID)
  
 END SUBROUTINE compute_atomic_radial_moments
 END MODULE module_compute_atomic_radial_moments