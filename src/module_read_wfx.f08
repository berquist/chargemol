 MODULE module_read_wfx
 !===================================================================================
 ! Module with subroutine to read wfx files
 ! Copyright (c) 2014, 2015, 2016, 2017 by Thomas A. Manz and Nidia Gabaldon Limas. Rights reserved.
 !===================================================================================

 USE module_precision
 USE module_global_parameters
 USE module_string_utilities
 USE module_common_variable_declarations
 
 IMPLICIT NONE

 CONTAINS
 
 SUBROUTINE read_wfx( )
 !===================================================================================

 CHARACTER(200) :: text
 CHARACTER(200) :: line_of_text
 CHARACTER(15) :: temporal_string
 INTEGER :: io,k,ios,min_center,temp_index,temp_count,fileSize
 REAL(kind=dp) :: temp_primitive_exponents_row(2,5),temp_exp
 REAL(kind=dp), ALLOCATABLE, DIMENSION(:) :: temp_array2,temp_array,temp_occupation
 REAL(kind=dp), ALLOCATABLE, DIMENSION(:,:) :: temp_orbital_coefficients_column

 !Open and read the wfx file
 OPEN(NEWUNIT = input_FID,FILE = long_input_filename,IOSTAT=iostat_value,STATUS='old')
 IF (iostat_value /= 0) THEN
   WRITE(output_FID,*)'Could not find input wfx file. Program will terminate.'
   FLUSH(output_FID)
 ELSE
   INQUIRE(unit=input_FID,size=fileSize)
   IF (fileSize < 10) THEN
     WRITE(output_FID,*)'Input wfx file appears to be empty. Program will terminate'
     FLUSH(output_FID)
     STOP
   END IF
 END IF
 
 num_periodic_directions = 0
 nprimitives = 0
 ncenters = 0
 n_edf_primitives = 0
 norbitals = 0
 netcharge = 0.001_dp
 edf_available = .false.
 core_available = .false.
 spin_available = .false.
 non_collinear= .false.
 DO
   READ(input_FID,'(a)',IOSTAT=io) line_of_text
   IF (io > 0) THEN
     WRITE(output_FID,*) "Could not read wfx file. Program will terminate."
     FLUSH(output_FID)
     STOP
   ELSE IF (io < 0) THEN
     EXIT
   ELSE
     !adjusting left, lowercase, trim data and read
     line_of_text = SimpleLine(line_of_text)
     line_of_text = StrLowCase(line_of_text)
     !assinging a value to each item
     IF(line_of_text == '<keywords>')THEN
       READ(input_FID,*)basis_set_type
       basis_set_type = SimpleLine(basis_set_type)
       basis_set_type = StrLowCase(basis_set_type) 
       IF(basis_set_type == 'gto') THEN
         WRITE(output_FID,*)'GTO basis set type is supported. Calculation&
         & will proceed.'
         FLUSH(output_FID)
       ELSE
         WRITE(output_FID,*)'Basis set type is not supported. Calculation&
         & will terminate.'
         FLUSH(output_FID)
         STOP
       END IF
     ELSE IF(line_of_text == '<net charge>') THEN
       READ(input_FID,*)netcharge
       WRITE(output_FID,'(a,f4.2)') " netcharge= ",netcharge
       FLUSH(output_FID)
     ELSE IF(line_of_text == '<number of translation vectors>') THEN
       READ(input_FID,*)num_periodic_directions 
       WRITE(output_FID,'(a,i5)')" num_periodic_directoins=",num_periodic_directions
       FLUSH(output_FID)
     ELSE IF(line_of_text =='<number of nuclei>') THEN
        READ(input_FID,*)ncenters
        WRITE(output_FID,'(a,i5)')" ncenters= ",ncenters
        FLUSH(output_FID)
     ELSE IF(line_of_text == '<number of primitives>')THEN
       READ(input_FID,*)nprimitives
       WRITE(output_FID,'(a,i5)')" nprimitives= ",nprimitives
       FLUSH(output_FID)
     ELSE IF(line_of_text == '<number of edf primitives>')THEN
       READ(input_FID,*)n_edf_primitives
       WRITE(output_FID,'(a,i7)')" n_edf_primitive= ",n_edf_primitives 
       FLUSH(output_FID)
       edf_available = .true.
     ELSE IF(line_of_text == '<number of occupied molecular orbitals>')THEN
       READ(input_FID,*)norbitals
       WRITE(output_FID,'(a,i7)')" norbitals= ",norbitals
       FLUSH(output_FID)
     END IF
   END IF
 END DO
 IF(len_trim(basis_set_type) == 0)THEN
   WRITE(output_FID,*)'Basis set type could not be read from wfx file.&
   & Program will terminate.'
   FLUSH(output_FID)
   STOP
 ELSE IF(netcharge == 0.001_dp)THEN
   WRITE(output_FID,*)'Net charge could not be read from wfx file. Pro&
   &gram will terminate.'
   FLUSH(output_FID)
   STOP
 ELSE IF(ncenters == 0)THEN
   WRITE(output_FID,*)'Number of basis set centers could not be read f&
   &rom wfx file. Program will terminate.'
   FLUSH(output_FID)
   STOP
 ELSE IF(nprimitives == 0)THEN
   WRITE(output_FID,*)'Number of basis set primitives could not be rea&
   &d from wfxm file. Program will terminate.'
   FLUSH(output_FID)
   STOP
 ELSE IF(norbitals == 0)THEN
   WRITE(output_FID,*)'Number of natural orbitals could not be read fr&
   &om wfx file. Program will terminate.'
   FLUSH(output_FID)
   STOP
 END IF
   
 ALLOCATE(periodic_vectors(num_periodic_directions,3))
 ALLOCATE(basis_set_centers(ncenters,6))
 ALLOCATE(primitive_exponents(nprimitives,5))
 ALLOCATE(edf_primitives(n_edf_primitives,15))
 ALLOCATE(alpha_beta_occupation(norbitals,2))
 ALLOCATE(orbital_coefficients(norbitals,nprimitives))
 ALLOCATE(temp_array2(n_edf_primitives))
 ALLOCATE(temp_array(nprimitives))
 ALLOCATE(temp_occupation(norbitals))
 periodic_vectors = 0.0_dp
 basis_set_centers = 0.0_dp
 primitive_exponents = 0.0_dp
 edf_primitives = 0.0_dp
 alpha_beta_occupation = 0.0_dp
 orbital_coefficients = 0.0_dp
 temp_array2 = 0.0_dp
 temp_array = 0.0_dp
 temp_occupation = 0.0_dp
 CLOSE(input_FID)
 OPEN(NEWUNIT = input_FID,FILE = long_input_filename)
 !Reading data
 DO
   READ(input_FID,'(a)',IOSTAT=ios) line_of_text
   IF (ios > 0) THEN
     WRITE(output_FID,*)"Could not open wfx file. Program will temrinate"
     FLUSH(output_FID)
     STOP
   ELSE IF (ios < 0) THEN
     EXIT
   ELSE
     line_of_text = SimpleLine(line_of_text)
     line_of_text = StrLowCase(line_of_text)
     IF (line_of_text == '<translation vectors>') THEN
       DO j=1,num_periodic_directions
         READ(input_FID,*)periodic_vectors(j,1:3)
       END DO
       WRITE(output_FID,*)'periodic_vectors'
       DO i=1,num_periodic_directions
         WRITE(output_FID,'(3f10.4)')(periodic_vectors(i,j),j=1,3)
       END DO
     ELSE IF (line_of_text == '<atomic numbers>') THEN
       DO i=1,ncenters
        READ(input_FID,*)basis_set_centers(i,1)
       END DO
     ELSE IF (line_of_text == '<nuclear charges>') THEN
       DO i=1,ncenters
         READ(input_FID,*)basis_set_centers(i,2)
       END DO
     ELSE IF (line_of_text == '<nuclear cartesian coordinates>') THEN
       DO i = 1,ncenters
         READ(input_FID,*)(basis_set_centers(i,j),j = 3,5)
       END DO
     ELSE IF (line_of_text == '<primitive centers>') THEN
       READ(input_FID,*)(primitive_exponents(i,1),i=1,nprimitives)
     ELSE IF (line_of_text == '<primitive exponents>') THEN
       READ(input_FID,*)(primitive_exponents(i,2),i=1,nprimitives)
     ELSE IF (line_of_text == '<primitive types>' ) THEN
       READ(input_FID,*)(temp_array(i),i=1,nprimitives)
       DO i=1,nprimitives
         IF (temp_array(i) == 1) THEN
           primitive_exponents(i,3) = 0
           primitive_exponents(i,4) = 0
           primitive_exponents(i,5) = 0
         ELSE IF (temp_array(i) == 2) THEN
           primitive_exponents(i,3) = 1
           primitive_exponents(i,4) = 0
           primitive_exponents(i,5) = 0
         ELSE IF (temp_array(i) == 3) THEN
           primitive_exponents(i,3) = 0
           primitive_exponents(i,4) = 1
           primitive_exponents(i,5) = 0
         ELSE IF (temp_array(i) == 4) THEN
           primitive_exponents(i,3) = 0
           primitive_exponents(i,4) = 0
           primitive_exponents(i,5) = 1
          ELSE IF (temp_array(i) == 5) THEN
           primitive_exponents(i,3) = 2
           primitive_exponents(i,4) = 0
           primitive_exponents(i,5) = 0
         ELSE IF (temp_array(i) == 6) THEN
           primitive_exponents(i,3) = 0
           primitive_exponents(i,4) = 2
           primitive_exponents(i,5) = 0
         ELSE IF (temp_array(i) == 7) THEN
           primitive_exponents(i,3) = 0
           primitive_exponents(i,4) = 0
           primitive_exponents(i,5) = 2
         ELSE IF (temp_array(i) == 8) THEN
           primitive_exponents(i,3) = 1
           primitive_exponents(i,4) = 1
           primitive_exponents(i,5) = 0
         ELSE IF (temp_array(i) == 9) THEN
           primitive_exponents(i,3) = 1
           primitive_exponents(i,4) = 0
           primitive_exponents(i,5) = 1
         ELSE IF (temp_array(i) == 10) THEN
           primitive_exponents(i,3) = 0
           primitive_exponents(i,4) = 1
           primitive_exponents(i,5) = 1
         ELSE IF (temp_array(i) == 11) THEN
           primitive_exponents(i,3) = 3
           primitive_exponents(i,4) = 0
           primitive_exponents(i,5) = 0
         ELSE IF (temp_array(i) == 12) THEN
           primitive_exponents(i,3) = 0
           primitive_exponents(i,4) = 3
           primitive_exponents(i,5) = 0
         ELSE IF (temp_array(i) == 13) THEN
           primitive_exponents(i,3) = 0
           primitive_exponents(i,4) = 0
           primitive_exponents(i,5) = 3
         ELSE IF (temp_array(i) == 14) THEN
           primitive_exponents(i,3) = 2
           primitive_exponents(i,4) = 1
           primitive_exponents(i,5) = 0
         ELSE IF (temp_array(i) == 15) THEN
           primitive_exponents(i,3) = 2
           primitive_exponents(i,4) = 0
           primitive_exponents(i,5) = 1
         ELSE IF (temp_array(i) == 16) THEN
           primitive_exponents(i,3) = 0
           primitive_exponents(i,4) = 2
           primitive_exponents(i,5) = 1
         ELSE IF (temp_array(i) == 17) THEN
           primitive_exponents(i,3) = 1
           primitive_exponents(i,4) = 2
           primitive_exponents(i,5) = 0
         ELSE IF (temp_array(i) == 18) THEN
           primitive_exponents(i,3) = 1
           primitive_exponents(i,4) = 0
           primitive_exponents(i,5) = 2
         ELSE IF (temp_array(i) == 19) THEN
           primitive_exponents(i,3) = 0
           primitive_exponents(i,4) = 1
           primitive_exponents(i,5) = 2
         ELSE IF (temp_array(i) == 20) THEN
           primitive_exponents(i,3) = 1
           primitive_exponents(i,4) = 1
           primitive_exponents(i,5) = 1
         ELSE IF (temp_array(i) == 21) THEN
           primitive_exponents(i,3) = 4
           primitive_exponents(i,4) = 0
           primitive_exponents(i,5) = 0
         ELSE IF (temp_array(i) == 22) THEN
           primitive_exponents(i,3) = 0
           primitive_exponents(i,4) = 4
           primitive_exponents(i,5) = 0
         ELSE IF (temp_array(i) == 23) THEN
           primitive_exponents(i,3) = 0
           primitive_exponents(i,4) = 0
           primitive_exponents(i,5) = 4
         ELSE IF (temp_array(i) == 24) THEN
           primitive_exponents(i,3) = 3
           primitive_exponents(i,4) = 1
           primitive_exponents(i,5) = 0
         ELSE IF (temp_array(i) == 25) THEN
           primitive_exponents(i,3) = 3
           primitive_exponents(i,4) = 0
           primitive_exponents(i,5) = 1
         ELSE IF (temp_array(i) == 26) THEN
           primitive_exponents(i,3) = 1
           primitive_exponents(i,4) = 3
           primitive_exponents(i,5) = 0
         ELSE IF (temp_array(i) == 27) THEN
           primitive_exponents(i,3) = 0
           primitive_exponents(i,4) = 3
           primitive_exponents(i,5) = 1
         ELSE IF (temp_array(i) == 28) THEN
           primitive_exponents(i,3) = 1
           primitive_exponents(i,4) = 0
           primitive_exponents(i,5) = 3
         ELSE IF (temp_array(i) == 29) THEN
           primitive_exponents(i,3) = 0
           primitive_exponents(i,4) = 1
           primitive_exponents(i,5) = 3
         ELSE IF (temp_array(i) == 30) THEN
           primitive_exponents(i,3) = 2
           primitive_exponents(i,4) = 2
           primitive_exponents(i,5) = 0
         ELSE IF (temp_array(i) == 31) THEN
           primitive_exponents(i,3) = 2
           primitive_exponents(i,4) = 0
           primitive_exponents(i,5) = 2
         ELSE IF (temp_array(i) == 32) THEN
           primitive_exponents(i,3) = 0
           primitive_exponents(i,4) = 2
           primitive_exponents(i,5) = 2
         ELSE IF (temp_array(i) == 33) THEN
           primitive_exponents(i,3) = 2
           primitive_exponents(i,4) = 1
           primitive_exponents(i,5) = 1
         ELSE IF (temp_array(i) == 34) THEN
           primitive_exponents(i,3) = 1
           primitive_exponents(i,4) = 2
           primitive_exponents(i,5) = 1
         ELSE IF (temp_array(i) == 35) THEN
           primitive_exponents(i,3) = 1
           primitive_exponents(i,4) = 1
           primitive_exponents(i,5) = 2
         ELSE IF (temp_array(i) == 36) THEN
           primitive_exponents(i,3) = 0
           primitive_exponents(i,4) = 0
           primitive_exponents(i,5) = 5
         ELSE IF (temp_array(i) == 37) THEN
           primitive_exponents(i,3) = 0
           primitive_exponents(i,4) = 1
           primitive_exponents(i,5) = 4
         ELSE IF (temp_array(i) == 38) THEN
           primitive_exponents(i,3) = 0
           primitive_exponents(i,4) = 2
           primitive_exponents(i,5) = 3
         ELSE IF (temp_array(i) == 39) THEN
           primitive_exponents(i,3) = 0
           primitive_exponents(i,4) = 3
           primitive_exponents(i,5) = 2
         ELSE IF (temp_array(i) == 40) THEN
           primitive_exponents(i,3) = 0
           primitive_exponents(i,4) = 4 
           primitive_exponents(i,5) = 1
         ELSE IF (temp_array(i) == 41) THEN
           primitive_exponents(i,3) = 0
           primitive_exponents(i,4) = 5
           primitive_exponents(i,5) = 0
         ELSE IF (temp_array(i) == 42) THEN
           primitive_exponents(i,3) = 1
           primitive_exponents(i,4) = 0
           primitive_exponents(i,5) = 4
         ELSE IF (temp_array(i) == 43) THEN
           primitive_exponents(i,3) = 1
           primitive_exponents(i,4) = 1
           primitive_exponents(i,5) = 3
         ELSE IF (temp_array(i) == 44) THEN
           primitive_exponents(i,3) = 1
           primitive_exponents(i,4) = 2
           primitive_exponents(i,5) = 2
         ELSE IF (temp_array(i) == 45) THEN
           primitive_exponents(i,3) = 1
           primitive_exponents(i,4) = 3
           primitive_exponents(i,5) = 1
         ELSE IF (temp_array(i) == 46) THEN
           primitive_exponents(i,3) = 1
           primitive_exponents(i,4) = 4
           primitive_exponents(i,5) = 0
         ELSE IF (temp_array(i) == 47) THEN
           primitive_exponents(i,3) = 2
           primitive_exponents(i,4) = 0
           primitive_exponents(i,5) = 3
         ELSE IF (temp_array(i) == 48) THEN
           primitive_exponents(i,3) = 2
           primitive_exponents(i,4) = 1
           primitive_exponents(i,5) = 2
         ELSE IF (temp_array(i) == 49) THEN
           primitive_exponents(i,3) = 2
           primitive_exponents(i,4) = 2
           primitive_exponents(i,5) = 1
         ELSE IF (temp_array(i) == 50) THEN
           primitive_exponents(i,3) = 2
           primitive_exponents(i,4) = 3
           primitive_exponents(i,5) = 0
         ELSE IF (temp_array(i) == 51) THEN
           primitive_exponents(i,3) = 3
           primitive_exponents(i,4) = 0
           primitive_exponents(i,5) = 2
         ELSE IF (temp_array(i) == 52) THEN
           primitive_exponents(i,3) = 3
           primitive_exponents(i,4) = 1
           primitive_exponents(i,5) = 1
         ELSE IF (temp_array(i) == 53) THEN
           primitive_exponents(i,3) = 3
           primitive_exponents(i,4) = 2
           primitive_exponents(i,5) = 0
         ELSE IF (temp_array(i) == 54) THEN
           primitive_exponents(i,3) = 4
           primitive_exponents(i,4) = 0
           primitive_exponents(i,5) = 1
         ELSE IF (temp_array(i) == 55) THEN
           primitive_exponents(i,3) = 4
           primitive_exponents(i,4) = 1
           primitive_exponents(i,5) = 0
         ELSE IF (temp_array(i) == 56) THEN
           primitive_exponents(i,3) = 5
           primitive_exponents(i,4) = 0
           primitive_exponents(i,5) = 0
         ELSE
           WRITE(output_FID,*)'Only basis set primitives up to H shells &
           &(L=5) are supported. Program will terminate.'
           FLUSH(output_FID)
           STOP 
         END IF
       END DO
     ELSE IF (line_of_text == '<edf primitive centers>') THEN
       READ(input_FID,*)(edf_primitives(i,1),i=1,n_edf_primitives)
     ELSE IF (line_of_text == '<edf primitive exponents>') THEN
       READ(input_FID,*)(edf_primitives(i,2),i=1,n_edf_primitives)
     ELSE IF (line_of_text == '<edf primitive coefficients>') THEN
       READ(input_FID,*)(edf_primitives(i,6),i=1,n_edf_primitives)
     ELSE IF (line_of_text == '<edf primitive types>') THEN
       READ(input_FID,*)(temp_array2(i),i=1,n_edf_primitives)
       DO i=1,n_edf_primitives
         IF (temp_array2(i) == 1) THEN
           edf_primitives(i,3) = 0
           edf_primitives(i,4) = 0
           edf_primitives(i,5) = 0 
         ELSE IF (temp_array2(i) == 2) THEN
           edf_primitives(i,3) = 1
           edf_primitives(i,4) = 0
           edf_primitives(i,5) = 0
         ELSE IF (temp_array2(i) == 3) THEN  
           edf_primitives(i,3) = 0
           edf_primitives(i,4) = 1
           edf_primitives(i,5) = 0 
         ELSE IF (temp_array2(i) == 4) THEN 
           edf_primitives(i,3) = 0
           edf_primitives(i,4) = 0
           edf_primitives(i,5) = 1
         ELSE IF (temp_array2(i) == 5) THEN
           edf_primitives(i,3) = 2
           edf_primitives(i,4) = 0
           edf_primitives(i,5) = 0
         ELSE IF (temp_array2(i) == 6) THEN
           edf_primitives(i,3) = 0
           edf_primitives(i,4) = 2
           edf_primitives(i,5) = 0  
         ELSE IF (temp_array2(i) == 7) THEN
           edf_primitives(i,3) = 0
           edf_primitives(i,4) = 0
           edf_primitives(i,5) = 2
         ELSE IF (temp_array2(i) == 8) THEN
           edf_primitives(i,3) = 1
           edf_primitives(i,4) = 1
           edf_primitives(i,5) = 0
         ELSE IF (temp_array2(i) == 9) THEN
          edf_primitives(i,3) = 1
          edf_primitives(i,4) = 0
          edf_primitives(i,5) = 1
         ELSE IF (temp_array2(i) == 10) THEN
           edf_primitives(i,3) = 0
           edf_primitives(i,4) = 1
           edf_primitives(i,5) = 1
         ELSE IF (temp_array2(i) == 11) THEN
           edf_primitives(i,3) = 3
           edf_primitives(i,4) = 0
           edf_primitives(i,5) = 0
         ELSE IF (temp_array2(i) == 12) THEN
           edf_primitives(i,3) = 0
           edf_primitives(i,4) = 3
           edf_primitives(i,5) = 0
         ELSE IF (temp_array2(i) == 13) THEN
           edf_primitives(i,3) = 0
           edf_primitives(i,4) = 0
           edf_primitives(i,5) = 3
         ELSE IF (temp_array2(i) == 14) THEN
           edf_primitives(i,3) = 2
           edf_primitives(i,4) = 1
           edf_primitives(i,5) = 0
         ELSE IF (temp_array2(i) == 15) THEN
           edf_primitives(i,3) = 2
           edf_primitives(i,4) = 0
           edf_primitives(i,5) = 1
         ELSE IF (temp_array2(i) == 16) THEN
           edf_primitives(i,3) = 0
           edf_primitives(i,4) = 2
           edf_primitives(i,5) = 1
         ELSE IF (temp_array2(i) == 17) THEN
           edf_primitives(i,3) = 1
           edf_primitives(i,4) = 2
           edf_primitives(i,5) = 0
         ELSE IF (temp_array2(i) == 18) THEN
           edf_primitives(i,3) = 1
           edf_primitives(i,4) = 0
           edf_primitives(i,5) = 2 
         ELSE IF (temp_array2(i) == 19) THEN
           edf_primitives(i,3) = 0
           edf_primitives(i,4) = 1
           edf_primitives(i,5) = 2
         ELSE IF (temp_array2(i) == 20) THEN
           edf_primitives(i,3) = 1
           edf_primitives(i,4) = 1
           edf_primitives(i,5) = 1
         ELSE IF (temp_array2(i) == 21) THEN
           edf_primitives(i,3) = 4
           edf_primitives(i,4) = 0
           edf_primitives(i,5) = 0
         ELSE IF (temp_array2(i) == 22) THEN
           edf_primitives(i,3) = 0
           edf_primitives(i,4) = 4
           edf_primitives(i,5) = 0
         ELSE IF (temp_array2(i) == 23) THEN
           edf_primitives(i,3) = 0
           edf_primitives(i,4) = 0
           edf_primitives(i,5) = 4
         ELSE IF (temp_array2(i) == 24) THEN
           edf_primitives(i,3) = 3
           edf_primitives(i,4) = 1
           edf_primitives(i,5) = 0
         ELSE IF (temp_array2(i) == 25) THEN
           edf_primitives(i,3) = 3
           edf_primitives(i,4) = 0
           edf_primitives(i,5) = 1
         ELSE IF (temp_array2(i) == 26) THEN
           edf_primitives(i,3) = 1
           edf_primitives(i,4) = 3
           edf_primitives(i,5) = 0
          ELSE IF (temp_array2(i) == 27) THEN
           edf_primitives(i,3) = 0
           edf_primitives(i,4) = 3
           edf_primitives(i,5) = 1
         ELSE IF (temp_array2(i) == 28) THEN
           edf_primitives(i,3) = 1
           edf_primitives(i,4) = 0
           edf_primitives(i,5) = 3
         ELSE IF (temp_array2(i) == 29) THEN
           edf_primitives(i,3) = 0
           edf_primitives(i,4) = 1
           edf_primitives(i,5) = 3
         ELSE IF (temp_array2(i) == 30) THEN
           edf_primitives(i,3) = 2
           edf_primitives(i,4) = 2
           edf_primitives(i,5) = 0
         ELSE IF (temp_array2(i) == 31) THEN
           edf_primitives(i,3) = 2 
           edf_primitives(i,4) = 0 
           edf_primitives(i,5) = 2 
         ELSE IF (temp_array2(i) == 32) THEN
           edf_primitives(i,3) = 0
           edf_primitives(i,4) = 2
           edf_primitives(i,5) = 2
         ELSE IF (temp_array2(i) == 33) THEN
           edf_primitives(i,3) = 2
           edf_primitives(i,4) = 1
           edf_primitives(i,5) = 1
         ELSE IF (temp_array2(i) == 34) THEN
           edf_primitives(i,3) = 1
           edf_primitives(i,4) = 2
           edf_primitives(i,5) = 1
         ELSE IF (temp_array2(i) == 35) THEN
           edf_primitives(i,3) = 1
           edf_primitives(i,4) = 1
           edf_primitives(i,5) = 2
         ELSE IF (temp_array2(i) == 36) THEN
           edf_primitives(i,3) = 0
           edf_primitives(i,4) = 0
           edf_primitives(i,5) = 5
         ELSE IF (temp_array2(i) == 37) THEN
           edf_primitives(i,3) = 0
           edf_primitives(i,4) = 1
           edf_primitives(i,5) = 4
         ELSE IF (temp_array2(i) == 38) THEN
           edf_primitives(i,3) = 0
           edf_primitives(i,4) = 2
           edf_primitives(i,5) = 3
         ELSE IF (temp_array2(i) == 39) THEN
           edf_primitives(i,3) = 0
           edf_primitives(i,4) = 3
           edf_primitives(i,5) = 2
         ELSE IF (temp_array2(i) == 40) THEN
           edf_primitives(i,3) = 0
           edf_primitives(i,4) = 4
           edf_primitives(i,5) = 1
         ELSE IF (temp_array2(i) == 41) THEN
           edf_primitives(i,3) = 0
           edf_primitives(i,4) = 5
           edf_primitives(i,5) = 0
         ELSE IF (temp_array2(i) == 42) THEN
           edf_primitives(i,3) = 1
           edf_primitives(i,4) = 0
           edf_primitives(i,5) = 4
         ELSE IF (temp_array2(i) == 43) THEN
           edf_primitives(i,3) = 1
           edf_primitives(i,4) = 1
           edf_primitives(i,5) = 3
         ELSE IF (temp_array2(i) == 44) THEN
           edf_primitives(i,3) = 1
           edf_primitives(i,4) = 2
           edf_primitives(i,5) = 2
         ELSE IF (temp_array2(i) == 45) THEN
           edf_primitives(i,3) = 1
           edf_primitives(i,4) = 3
           edf_primitives(i,5) = 1
         ELSE IF (temp_array2(i) == 46) THEN
           edf_primitives(i,3) = 1
           edf_primitives(i,4) = 4
           edf_primitives(i,5) = 0
         ELSE IF (temp_array2(i) == 47) THEN
           edf_primitives(i,3) = 2
           edf_primitives(i,4) = 0
           edf_primitives(i,5) = 3
         ELSE IF (temp_array2(i) == 48) THEN
           edf_primitives(i,3) = 2
           edf_primitives(i,4) = 1
           edf_primitives(i,5) = 2
         ELSE IF (temp_array2(i) == 49) THEN
           edf_primitives(i,3) = 2
           edf_primitives(i,4) = 2
           edf_primitives(i,5) = 1
         ELSE IF (temp_array2(i) == 50) THEN
           edf_primitives(i,3) = 2
           edf_primitives(i,4) = 3
           edf_primitives(i,5) = 0
         ELSE IF (temp_array2(i) == 51) THEN
           edf_primitives(i,3) = 3
           edf_primitives(i,4) = 0
           edf_primitives(i,5) = 2
         ELSE IF (temp_array2(i) == 52) THEN
           edf_primitives(i,3) = 3
           edf_primitives(i,4) = 1
           edf_primitives(i,5) = 1
         ELSE IF (temp_array2(i) == 53) THEN
           edf_primitives(i,3) = 3
           edf_primitives(i,4) = 2
           edf_primitives(i,5) = 0
         ELSE IF (temp_array2(i) == 54) THEN
           edf_primitives(i,3) = 4
           edf_primitives(i,4) = 0
           edf_primitives(i,5) = 1
         ELSE IF (temp_array2(i) == 55) THEN
           edf_primitives(i,3) = 4
           edf_primitives(i,4) = 1
           edf_primitives(i,5) = 0
         ELSE IF (temp_array2(i) == 56) THEN
           edf_primitives(i,3) = 5
           edf_primitives(i,4) = 0
           edf_primitives(i,5) = 0
         ELSE
           WRITE(output_FID,*)'Only basis set primitives up to H shells &
           &(L=5) are supported. Program will terminate.'
           FLUSH(output_FID)
           STOP
         END IF
       END DO
     ELSE IF (line_of_text == '<molecular orbital occupation numbers>') THEN
       READ(input_FID,*)(temp_occupation(i),i=1,norbitals)
     ELSE IF (line_of_text == '<molecular orbital spin types>') THEN
       DO i=1,norbitals
         READ(input_FID,'(a)')temporal_string
         temporal_string=SimpleLine(temporal_string)
         IF (temporal_string == 'Alpha') THEN  
           IF (temp_occupation(i) < 1.01_dp) THEN
             alpha_beta_occupation(i,1) = temp_occupation(i)
             alpha_beta_occupation(i,2) = 0.0_dp
             spin_available = .true.
             non_collinear = .false.
           ELSE
             WRITE(output_FID,*)"Alpha occupation numbers greater than one are &
	      &unphysical."
             WRITE(output_FID,*)"The program you used to generate the wfx file &
             &contains bug."
             WRITE(output_FID,*)"Please regenerate the wfx file using the lates&
             &t version of your quantum chemistry program."
             WRITE(output_FID,*)"If the error persists, please inform your quan&
             &tum chemistry program vendor of the bug."
             WRITE(output_FID,*)"Program will terminate."
             FLUSH(output_FID)
             STOP
           END IF
         ELSE IF (temporal_string == 'Beta') THEN 
           alpha_beta_occupation(i,1) = 0.0_dp
           alpha_beta_occupation(i,2) = temp_occupation(i)
         ELSE
           alpha_beta_occupation(i,1) = temp_occupation(i)/2.0_dp
           alpha_beta_occupation(i,2) = temp_occupation(i)/2.0_dp
         END IF
       END DO
       WRITE(output_FID,*)'Natural orbital occupation numbers for alpha (1st column) and beta (2nd column) electrons.'
       FLUSH(output_FID)
       DO i=1,norbitals
         WRITE(output_FID, '(2(F12.4))') alpha_beta_occupation(i,:)
         FLUSH(output_FID)
       END DO
       WRITE(output_FID,*) 'The total number of electrons in natural orbitals is'
       included_electrons=nint(sum(alpha_beta_occupation))
       WRITE(output_FID,'(a,i7)') ' included electrons= ',included_electrons
       FLUSH(output_FID)
       CYCLE
     ELSE IF (line_of_text == '<molecular orbital primitive coefficients>') THEN
       DO i=1,norbitals
         READ(input_FID,'(a)')text
         IF (text == "</molecular orbital primitive coefficients>") THEN
           EXIT
         ELSE
           READ(input_FID,*)text
           READ(input_FID,*)text
           READ(input_FID,*)(orbital_coefficients(i,j),j=1,nprimitives)
         END IF
       END DO
     END IF
   END IF
 END DO
 IF (sum(alpha_beta_occupation) == 0) THEN
   WRITE(output_FID,*)'Natural orbital occupation numbers could not be read from wfx file. Program will terminate.'
   FLUSH(output_FID)
   STOP
 ELSE IF (sum(orbital_coefficients) == 0)THEN
   WRITE(output_FID,*)'Natural orbital coefficients could not be read from wfx file. Program will terminate.'
   FLUSH(output_FID)
   STOP
 END IF
 !Sort the basis set primitives by center number and primitive exponent
 WRITE(output_FID,*)'Sorting the basis set primitives:'
 FLUSH(output_FID)
 !Sort by the center number
 temp_primitive_exponents_row = 0.0_dp
 ALLOCATE(temp_orbital_coefficients_column(norbitals,2))
 temp_orbital_coefficients_column = 0.0_dp
 min_center = 0
 DO i = 1,(nprimitives-1)
   min_center = nint(primitive_exponents(i,1))
   temp_index = i
   DO j = i,nprimitives
     IF (nint(primitive_exponents(j,1)) < min_center) THEN
       min_center = nint(primitive_exponents(j,1))
       temp_index = j
     END IF
   END DO
   IF (temp_index /= i) THEN
     temp_primitive_exponents_row(1,:) = primitive_exponents(i,:)
     temp_primitive_exponents_row(2,:) = primitive_exponents(temp_index,:)
     temp_orbital_coefficients_column(:,1) = orbital_coefficients(:,i)
     temp_orbital_coefficients_column(:,2) = orbital_coefficients(:,temp_index)
     primitive_exponents(temp_index,:) = temp_primitive_exponents_row(1,:)
     primitive_exponents(i,:) = temp_primitive_exponents_row(2,:)
     orbital_coefficients(:,temp_index) = temp_orbital_coefficients_column(:,1)
     orbital_coefficients(:,i) = temp_orbital_coefficients_column(:,2)
   END IF
 END DO
 !Count the number of primitives assigned to each basis set center
 ALLOCATE(nprimitives_per_center(ncenters))
 nprimitives_per_center = 0
 DO i = 1,nprimitives
   nprimitives_per_center(nint(primitive_exponents(i,1))) = nprimitives_per_center(nint(primitive_exponents(i,1))) + 1
 END DO
 !Subsort by the primitive exponent
 start_index = 1
 temp_exp = 0.0_dp
 DO k = 1,ncenters
   stop_index = start_index + nprimitives_per_center(k) - 1
   DO i = start_index,(stop_index - 1)
     temp_exp = primitive_exponents(i,2)
     temp_index = i
     DO j = i,stop_index
       IF (primitive_exponents(j,2) < temp_exp) THEN
         temp_exp = primitive_exponents(j,2)
         temp_index = j
       END IF
     END DO
     IF (temp_index /= i) THEN
       temp_primitive_exponents_row(1,:) = primitive_exponents(i,:)
       temp_primitive_exponents_row(2,:) = primitive_exponents(temp_index,:)
       temp_orbital_coefficients_column(:,1) = orbital_coefficients(:,i)
       temp_orbital_coefficients_column(:,2) = orbital_coefficients(:,temp_index)
       primitive_exponents(temp_index,:) = temp_primitive_exponents_row(1,:)
       primitive_exponents(i,:) = temp_primitive_exponents_row(2,:)
       orbital_coefficients(:,temp_index) = temp_orbital_coefficients_column(:,1)
       orbital_coefficients(:,i) = temp_orbital_coefficients_column(:,2)
     END IF
   END DO
   start_index = stop_index + 1
 END DO
 !Count the number of common exponent blocks
 nsingleblocks = 0
 temp_index = 0
 temp_exp = 0.0_dp
 DO i = 1,nprimitives
   IF ((nint(primitive_exponents(i,1)) /= temp_index) .or. (abs(temp_exp - primitive_exponents(i,2)) > zero_tolerance)) THEN
     nsingleblocks = nsingleblocks + 1
     temp_index = nint(primitive_exponents(i,1))
     temp_exp = primitive_exponents(i,2)
   END IF
 END DO
 !Store the range for each common exponent block
 ALLOCATE(single_block_ranges(nsingleblocks,2))
 single_block_ranges = 0
 nsingleblocks = 0
 temp_index = 0
 temp_exp = 0.0_dp
 DO i = 1,nprimitives
   IF ((nint(primitive_exponents(i,1)) /= temp_index) .or. (abs(temp_exp - primitive_exponents(i,2)) > zero_tolerance)) THEN
     nsingleblocks = nsingleblocks + 1
     temp_index = nint(primitive_exponents(i,1))
     temp_exp = primitive_exponents(i,2)
     single_block_ranges(nsingleblocks,1) = i
   END IF
 END DO
 single_block_ranges(nsingleblocks,2) = nprimitives
 DO i=1,(nsingleblocks - 1)
   single_block_ranges(i,2) = single_block_ranges((i+1),1) - 1
 END DO   
 !Determining number of atoms
 natoms = 0
 DO i=1,ncenters
   IF (basis_set_centers(i,2) > 0.5_dp) THEN
     natoms = natoms+1
   END IF
 END DO
 WRITE(output_FID,'(a,i5)')' natoms= ',natoms
 FLUSH(output_FID)
 !Set the atomic coordinates, core charges, etc.
 ALLOCATE(coords(3,natoms))
 ALLOCATE(atomic_number(natoms))
 ALLOCATE(core_electrons(natoms))
 ALLOCATE(missing_core(natoms))
 coords = 0.0_dp
 atomic_number = 0
 core_electrons = 0.0_dp
 missing_core = 0
 temp_count = 0
 DO i=1,ncenters
   IF (basis_set_centers(i,2) > 0.5_dp) THEN
     temp_count = temp_count + 1
     atomic_number(temp_count) = nint(basis_set_centers(i,1))
     core_electrons(temp_count) = basis_set_centers(i,1)-basis_set_centers(i,2)
     coords(1,temp_count) = basis_set_centers(i,3)
     coords(2,temp_count) = basis_set_centers(i,4)
     coords(3,temp_count) = basis_set_centers(i,5)
     basis_set_centers(i,6) = temp_count  !set the assigned atom in the last column of basis set centers
     IF (edf_available) THEN
       missing_core(temp_count) = 0
     ELSE
       missing_core(temp_count) = nint(core_electrons(temp_count))
     END IF
   END IF
 END DO
 !Run consistency checks for the number of valence electrons 
 IF (abs(sum(atomic_number) - netcharge - sum(core_electrons) - sum(temp_occupation))>0.001_dp) THEN
   WRITE(output_FID,*)'The number of valence electrons differs by more than 0.001e from the sum of orbital occupations. The &
   &quantum chemistry program you used to generate the wfx file contains a bug. Please contact the vendor of that program. &
   &Program will terminate.'
   FLUSH(output_FID)
   STOP
 ELSE
   WRITE(output_FID,*)'The number of valence electrons is consistent with the sum of orbital occupations.'
   FLUSH(output_FID)
 END IF
 CLOSE(input_FID)

END SUBROUTINE read_wfx

END MODULE module_read_wfx