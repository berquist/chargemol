 MODULE module_format_total_cube_density
 ! Copyright (c) 2014, 2015, 2016, 2017 by Thomas A. Manz and Nidia Gabaldon Limas. Rights reserved.
 
 USE module_common_variable_declarations
 USE module_global_parameters
 USE module_charge_center_positions_and_parallelpiped
 USE module_check_grid_spacing
 USE module_initialize_atomic_densities
 USE module_matrix_operations
 USE module_read_spin_density_cube_files

 IMPLICIT NONE
 
 CONTAINS
 
 SUBROUTINE format_total_cube_density( )
 !===================================================================================
 
 CHARACTER(200) :: headerlines
 REAL(kind=dp) :: v
 REAL(kind=dp),ALLOCATABLE,DIMENSION(:,:,:) ::  raw_valence
 INTEGER :: M,k,clock_count,clock_count_rate,clock_count_max,clock_count2
 
 WRITE(output_FID,*)' '
 WRITE(output_FID,*)'******************** TIME CONTROL ********************'
 WRITE(output_FID,*)'Starting format_total_cube_density'
 FLUSH(output_FID)
 CALL system_clock(clock_count,clock_count_rate,clock_count_max)
 WRITE(output_FID,*)'******************************************************'
 WRITE(output_FID,*)' '
 FLUSH(output_FID) 

 
 !Read in the file containing the total density
 OPEN(NEWUNIT=total_density_FID, FILE='total_density.cube',POSITION='append',IOSTAT=iostat_value,STATUS='old')
 IF (iostat_value /= 0) THEN
   WRITE(output_FID,*)'Could not open total_density.cube file. Program will terminate'
   FLUSH(output_FID)
   STOP
 ELSE
   WRITE(output_FID,*)'inputfile = "total_density.cube"'
   FLUSH(output_FID)
 END IF
 BACKSPACE(total_density_FID)
 READ(total_density_FID,'(A)')headerlines
 IF (len(trim(adjustl(headerlines))) /= 0) THEN
   WRITE(total_density_FID,*)''
 END IF
 CLOSE(total_density_FID) 
 OPEN(NEWUNIT=total_density_FID, FILE='total_density.cube',IOSTAT=iostat_value,STATUS='old')
 READ(total_density_FID,*)headerlines
 READ(total_density_FID,*)headerlines
 READ(total_density_FID,*)natoms,(origin(i),i=1,3)
 READ(total_density_FID,*)totnumA,(boundary(1,i),i=1,3) 
 READ(total_density_FID,*)totnumB,(boundary(2,i),i=1,3)
 READ(total_density_FID,*)totnumC,(boundary(3,i),i=1,3)
 origin=origin*distance_scale
 boundary = boundary*distance_scale
 pixelvolume=determinant(boundary)
 parameters(1,:) = [real(natoms,dp),origin(1),origin(2),origin(3)]
 parameters(2,:) = [real(totnumA,dp),boundary(1,1),boundary(1,2),boundary(1,3)]
 parameters(3,:) = [real(totnumB,dp),boundary(2,1),boundary(2,2),boundary(2,3)]
 parameters(4,:) = [real(totnumC,dp),boundary(3,1),boundary(3,2),boundary(3,3)]
 WRITE(output_FID,*)'parameters'
 DO i=1,4
   WRITE(output_FID,*)(parameters(i,j),j=1,4)
   FLUSH(output_FID)
 END DO
 ALLOCATE(atomic_number(natoms))
 ALLOCATE(atomic_number2(natoms))
 ALLOCATE(effective_nuclear_charge(natoms))
 ALLOCATE(effective_nuclear_charge2(natoms))
 ALLOCATE(coords(3,natoms))  
 ALLOCATE(coords2(3,natoms))
 effective_nuclear_charge = 0.0_dp
 DO i=1,natoms
   READ(total_density_FID,*)atomic_number(i),effective_nuclear_charge(i),(coords(j,i),j=1,3)
 END DO
 coords = coords * distance_scale
 M = totnumA*totnumB*totnumC
 ALLOCATE(raw_valence(totnumC,totnumB,totnumA))
 DO i = 1, totnumA
   DO j = 1, totnumB
     READ(total_density_FID,*) (raw_valence(k,j,i),k=1,totnumC)
   END DO
 END DO
 CLOSE(total_density_FID)
 raw_valence = raw_valence*density_scale
 !Read in the file containing the total density
 flag = 0
 !Construct the atomic density array:
 vector1=parameters(2,1)*[parameters(2,2),parameters(2,3),parameters(2,4)]
 vector2=parameters(3,1)*[parameters(3,2),parameters(3,3),parameters(3,4)]
 vector3=parameters(4,1)*[parameters(4,2),parameters(4,3),parameters(4,4)]
 sum_negative_density = 0.0_dp
 ALLOCATE(valence_density(totnumA,totnumB,totnumC))
 valence_density=0.0_dp
!$omp parallel do default(none) &
!$omp private(ka,kb,kc,v) &
!$omp shared(M,totnumA,totnumB,totnumC,raw_valence,pixelvolume,valence_density) &
!$omp reduction(+:sum_negative_density) &
!$omp schedule(static,collapsed_chunk_size)
 DO ka = 1,totnumA
   DO kb = 1,totnumB
     DO kc = 1,totnumC
       v = max(raw_valence(kc,kb,ka),0.0_dp)
       IF (raw_valence(kc,kb,ka) < 0.0_dp) THEN
         sum_negative_density = sum_negative_density - pixelvolume*(v - raw_valence(kc,kb,ka))
       END IF
       valence_density(ka,kb,kc) = v
     END DO
   END DO
 END DO
!$omp end parallel do
 DEALLOCATE(raw_valence)
 WRITE(output_FID,*)'sum_negative_density= ',sum_negative_density
 FLUSH(output_FID)
 ALLOCATE(core_electrons(natoms))
 ALLOCATE(core_density(totnumA,totnumB,totnumC))
 core_electrons = 0
 core_density = 0.0_dp
 ALLOCATE(center_nabc(3,natoms))
 ALLOCATE(center_shift(3,natoms))
 center_nabc = 0.0_dp
 center_shift = 0.0_dp
 CALL charge_center_positions( )
 CALL parallelpiped( )
 CALL check_grid_spacing( )
 CALL read_spin_density_cube_files( )
 CALL initialize_atomic_densities( )
 !Valence occupancy corrections can be added here if desired
 ALLOCATE(occupancy_correction(11,natoms))
 occupancy_correction = 0.0_dp
  
 WRITE(output_FID,*)' '
 WRITE(output_FID,*)'******************** TIME CONTROL ********************'
 CALL system_clock(clock_count2,clock_count_rate,clock_count_max)
 IF (clock_count_rate == 0) THEN
   WRITE(output_FID,*)'Processor clock not available, the time calculation wil not be performed.'
 ELSE
   seconds = modulo((clock_count2-clock_count),clock_count_max)/real(clock_count_rate)
   WRITE(output_FID,*)'Finished format_total_cube_density in ',seconds, 'seconds'
   FLUSH(output_FID)
 END IF
 WRITE(output_FID,*)'******************************************************'
 WRITE(output_FID,*)' '
 FLUSH(output_FID)  

 END SUBROUTINE format_total_cube_density

 END MODULE module_format_total_cube_density