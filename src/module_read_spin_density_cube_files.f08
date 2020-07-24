 MODULE module_read_spin_density_cube_files
 !===================================================================================
 ! Module with subroutine to read gaussian cube densities type files
 ! Copyright (c) 2014, 2015, 2016, 2017 by Thomas A. Manz and Nidia Gabaldon Limas. Rights reserved.
 !===================================================================================

 USE module_precision
 USE module_global_parameters
 USE module_common_variable_declarations
!$ USE omp_lib

 IMPLICIT NONE
 
 CONTAINS
 
 SUBROUTINE read_spin_density_cube_files()
 !===================================================================================
 
 ! Read in the file containing the valence density
 
 INTEGER :: M,k,spin_flag,spin_FID,spin_x_FID,spin_y_FID,spin_z_FID,&
 iostat_value_x,iostat_value_y,iostat_value_z,start_index,end_index
 CHARACTER(200) :: headerlines
 REAL(kind=dp) :: parameters2(4,4),p,s
 REAL(kind=dp),ALLOCATABLE,DIMENSION(:,:,:) :: raw_spin,raw_spin_x,raw_spin_y,raw_spin_z

 !Open and read the spin density information
 !Read in the file containing the spin density
 spin_flag = 0
 OPEN(NEWUNIT=spin_FID,FILE='spin_density.cube',IOSTAT=iostat_value,STATUS='old')
 IF (iostat_value == 0) THEN
   WRITE(output_FID,*)'inputfile = "spin_density.cube"'
   FLUSH(output_FID)
   spin_available = .true.
   non_collinear = .false.
 ELSE
   OPEN(NEWUNIT=spin_x_FID,FILE='spin_density_x.cube',IOSTAT=iostat_value_x,STATUS='old')
   OPEN(NEWUNIT=spin_y_FID,FILE='spin_density_y.cube',IOSTAT=iostat_value_y,STATUS='old')
   OPEN(NEWUNIT=spin_z_FID,FILE='spin_density_z.cube',IOSTAT=iostat_value_z,STATUS='old')
   IF ((iostat_value_x == 0) .and. (iostat_value_y == 0) .and. (iostat_value_z == 0)) THEN
     spin_available = .true.
     non_collinear = .true.
   ELSE
     spin_available =  .false.
   END IF
 END IF
 IF (spin_available) THEN
   IF (.not. non_collinear) THEN
     parameters2 = 0.0_dp
     atomic_number2 = 0.0_dp
     coords2 = 0.0_dp
     effective_nuclear_charge2 = 0.0_dp
     READ(spin_FID,*)headerlines
     READ(spin_FID,*)headerlines
     DO i=1,4
       READ(spin_FID,*)(parameters2(i,j),j=1,4)
	   parameters2(i,2:4) = parameters2(i,2:4) * distance_scale
     END DO
     IF (sum(abs(parameters - parameters2)) > zero_tolerance) THEN
       spin_flag = 1
     END IF
       DO i=1,natoms
       READ(spin_FID,*)atomic_number2(i),effective_nuclear_charge2(i),(coords2(j,i),j=1,3)
     END DO
     coords2 = coords2 *distance_scale
     IF(sum(abs(coords - coords2)) + sum(abs(effective_nuclear_charge - effective_nuclear_charge2)) + sum(abs(atomic_number-&
     atomic_number2)) > zero_tolerance) THEN
       spin_flag = 1
       WRITE(output_FID,*)'The spin density file does not contain the same lattice vectors, grid spacing, or atom positions as &
       &the density files.'
       WRITE(output_FID,*)'Program will skip the spin moment analysis.'
     END IF  
     ALLOCATE(raw_spin(totnumC,totnumB,totnumA))
     DO i = 1, totnumA
       DO j = 1, totnumB
         READ(spin_FID,*) (raw_spin(k,j,i),k=1,totnumC)
       END DO
     END DO
     CLOSE(spin_FID)
     raw_spin = raw_spin*density_scale
     IF (spin_flag == 1) THEN
       spin_available = .false.
     ELSE
       ALLOCATE(spin_density(totnumA,totnumB,totnumC))
       spin_density = 0.0_dp
!$omp parallel do default(none) &
!$omp private(ka,kb,kc,p,s) &
!$omp shared(M,totnumA,totnumB,totnumC,valence_density,core_density,raw_spin,spin_density) &
!$omp schedule(dynamic,chunk_size)
       DO ka = 1,totnumA
         DO kb = 1,totnumB
           DO kc = 1,totnumC
             p = valence_density(ka,kb,kc) + core_density(ka,kb,kc)
             s = raw_spin(kc,kb,ka)
             IF ((abs(p) < zero_tolerance) .or. (abs(s) < zero_tolerance)) THEN
               spin_density(ka,kb,kc) = 0.0_dp
             ELSE
               spin_density(ka,kb,kc) = s*min(p/abs(s),1.0_dp)
             END IF
           END DO
         END DO
       END DO
!$omp end parallel do
       DEALLOCATE(raw_spin)
     END IF
   ELSE
     !Load the spin_density_x.cube file
     WRITE(output_FID,*)'inputfile = "spin_density_x.cube"'
     FLUSH(output_FID)
     parameters2 = 0.0_dp
     atomic_number2 = 0.0_dp
     coords2 = 0.0_dp
     effective_nuclear_charge2 = 0.0_dp
     READ(spin_x_FID,*)headerlines
     READ(spin_x_FID,*)headerlines
     DO i=1,4
       READ(spin_x_FID,*)(parameters2(i,j),j=1,4)
	   parameters2(i,2:4) = parameters2(i,2:4) * distance_scale
     END DO
     IF (sum(abs(parameters - parameters2)) > zero_tolerance) THEN
       spin_flag = 1
     END IF
     DO i=1,natoms
       READ(spin_x_FID,*)atomic_number2(i),effective_nuclear_charge2(i),(coords2(j,i),j=1,3)
     END DO
     coords2 = coords2 * distance_scale
     IF(sum(abs(coords - coords2)) + sum(abs(effective_nuclear_charge - effective_nuclear_charge2)) + sum(abs(atomic_number-&
     atomic_number2)) > zero_tolerance) THEN
       spin_flag = 1
       WRITE(output_FID,*)'The spin_density_x.cube file does not contain the same lattice vectors, grid spacing, or atom positions&
       & as the density files.'
       WRITE(output_FID,*)'Program will skip the spin moment analysis.'
     ELSE
       ALLOCATE(raw_spin_x(totnumC,totnumB,totnumA))
       DO i = 1, totnumA
         DO j = 1, totnumB
           READ(spin_x_FID,*) (raw_spin_x(k,j,i),k=1,totnumC)
         END DO
       END DO
       CLOSE(spin_x_FID)
       raw_spin_x = raw_spin_x*density_scale
     END IF
     IF (spin_flag == 1) THEN
       spin_available = .false.
     END IF
     IF (spin_flag == 0) THEN
       parameters2 = 0.0_dp
       atomic_number2= 0.0_dp
       coords2 = 0.0_dp
       effective_nuclear_charge2 = 0.0_dp
       !Load the spin_density_y.cube file
       WRITE(output_FID,*)'inputfile = "spin_density_y.cube"'
       FLUSH(output_FID)
       READ(spin_y_FID,*)headerlines
       READ(spin_y_FID,*)headerlines
       DO i=1,4
         READ(spin_y_FID,*)(parameters2(i,j),j=1,4)
		 parameters2(i,2:4) = parameters2(i,2:4) * distance_scale
       END DO
       IF (sum(abs(parameters - parameters2)) > zero_tolerance) THEN
         spin_flag = 1
       END IF
       DO i=1,natoms
         READ(spin_y_FID,*)atomic_number2(i),effective_nuclear_charge2(i),(coords2(j,i),j=1,3)
       END DO
       coords2 = coords2 * distance_scale
       IF(sum(abs(coords - coords2)) + sum(abs(effective_nuclear_charge - effective_nuclear_charge2)) + sum(abs(atomic_number-&
       atomic_number2)) > zero_tolerance) THEN
         spin_flag = 1
         WRITE(output_FID,*)'The spin_density_y.cube file does not contain the same lattice vectors, grid spacing, or atom posit&
         &ions as the density files.'
         WRITE(output_FID,*)'Program will skip the spin moment analysis.'
       ELSE
         ALLOCATE(raw_spin_y(totnumC,totnumB,totnumA))
         DO i = 1, totnumA
           DO j = 1, totnumB
             READ(spin_y_FID,*) (raw_spin_y(k,j,i),k=1,totnumC)
           END DO
         END DO
         CLOSE(spin_y_FID)
         raw_spin_y = raw_spin_y*density_scale
       END IF
     END IF
     IF (spin_flag == 1) THEN
       spin_available = .false.
     END IF
     IF (spin_flag == 0) THEN
       parameters2 = 0.0_dp
       atomic_number2= 0.0_dp
       coords2 = 0.0_dp
       effective_nuclear_charge2 = 0.0_dp
       !Load the spin_density_z.cube file
       WRITE(output_FID,*)'inputfile = spin_density_z.cube"'
       FLUSH(output_FID)
       READ(spin_z_FID,*)headerlines
       READ(spin_z_FID,*)headerlines
       DO i=1,4
         READ(spin_z_FID,*)(parameters2(i,j),j=1,4)
		 parameters2(i,2:4) = parameters2(i,2:4) * distance_scale
       END DO
       IF (sum(abs(parameters - parameters2)) > zero_tolerance) THEN
         spin_flag = 1
       END IF
       DO i=1,natoms
         READ(spin_z_FID,*)atomic_number2(i),effective_nuclear_charge2(i),(coords2(j,i),j=1,3)
       END DO
       coords2 = coords2 * distance_scale
       IF(sum(abs(coords - coords2)) + sum(abs(effective_nuclear_charge - effective_nuclear_charge2)) + sum(abs(atomic_number-&
       atomic_number2)) > zero_tolerance) THEN
         spin_flag = 1
         WRITE(output_FID,*)'The spin_density_z.cube file does not contain the same lattice vectors, grid spacing, or atom pos&
         &itions as the density files.'
         WRITE(output_FID,*)'Program will skip the spin moment analysis.'
       ELSE
         ALLOCATE(raw_spin_z(totnumC,totnumB,totnumA))
         DO i = 1, totnumA
           DO j = 1, totnumB
             READ(spin_z_FID,*) (raw_spin_z(k,j,i),k=1,totnumC)
           END DO
         END DO
         CLOSE(spin_z_FID)
         raw_spin_z = raw_spin_z*density_scale
       END IF
     END IF
	 IF (spin_flag == 1) THEN
       spin_available = .false.
     END IF
     WRITE(output_FID,*)'Noncollinear analysis has not been thoroughly debugged when reading cube files. Please perform your own &
     &check to make sure result is sensible. Program will proceed.'
     FLUSH(output_FID)
     !Construct the spin density vector arrays
     ALLOCATE(spin_density_vector(3,totnumA,totnumB,totnumC))
     spin_density_vector = 0.0_dp
!$omp parallel do default(none) &
!$omp private(ka,kb,kc,p,s) &
!$omp shared(M,totnumA,totnumB,totnumC,valence_density,core_density,raw_spin_x,raw_spin_y,&
!$omp raw_spin_z,spin_density_vector) &
!$omp schedule(dynamic,collapsed_chunk_size)
     DO ka = 1,totnumA
       DO kb = 1,totnumB
         DO kc = 1,totnumC
           p = valence_density(ka,kb,kc) + core_density(ka,kb,kc)
           s = sqrt(raw_spin_x(kc,kb,ka)**2 + raw_spin_y(kc,kb,ka)**2 + raw_spin_z(kc,kb,ka)**2)
           IF ((abs(p) < zero_tolerance) .or. (abs(s) < zero_tolerance)) THEN
             spin_density_vector(1,ka,kb,kc) = 0.0_dp
             spin_density_vector(2,ka,kb,kc) = 0.0_dp
             spin_density_vector(3,ka,kb,kc) = 0.0_dp
           ELSE    
             spin_density_vector(1,ka,kb,kc) = raw_spin_x(kc,kb,ka)*min(p/s,1.0_dp)
             spin_density_vector(2,ka,kb,kc) = raw_spin_y(kc,kb,ka)*min(p/s,1.0_dp)
             spin_density_vector(3,ka,kb,kc) = raw_spin_z(kc,kb,ka)*min(p/s,1.0_dp)
           END IF     
         END DO
       END DO
     END DO
!$omp end parallel do
     DEALLOCATE(raw_spin_x)
     DEALLOCATE(raw_spin_y)
     DEALLOCATE(raw_spin_z)
   END IF
 END IF

 END SUBROUTINE read_spin_density_cube_files

 END MODULE module_read_spin_density_cube_files