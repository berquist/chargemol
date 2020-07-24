 MODULE module_format_vasp_densities
 ! Copyright (c) 2014, 2015, 2016, 2017 by Thomas A. Manz and Nidia Gabaldon Limas. Rights reserved.

 USE module_common_variable_declarations
 USE module_global_parameters
 USE module_matrix_operations
 USE module_atomic_symbol_to_number
 USE module_compute_dominant_atom_volumes
 USE module_initialize_atomic_densities
 USE module_add_missing_core_density

 IMPLICIT NONE
 
 CONTAINS

 SUBROUTINE format_vasp_densities( )
 !===================================================================================
 
 CHARACTER(20) :: CHG_file,totnumA_string,totnumB_string,totnumC_string
 CHARACTER(200) :: tempstring
 CHARACTER(50) :: headerlines,tline,atomic_symbol_string
 CHARACTER(2) :: atomic_symbol
 REAL(kind=dp) :: latticevectorfactor,vectors(3,3),latticevectorfactor2,vectors2(3,3),correction_size,v,&
 checkspin1_positive,checkspin2_positive,checkspin1_negative,checkspin2_negative,p,s,checkspin3_positive,&
 checkspin3_negative,checkspin3,totnum(3),checkspin1_positive_vector(3),&
 checkspin2_positive_vector(3),checkspin1_negative_vector(3),checkspin2_negative_vector(3),enhancement_factor,&
 checkspin3_positive_vector(3),checkspin3_negative_vector(3),s_trial(3),mag_s_trial,checkspin3_vector(3),checkspin3_mag,alpha,&
 beta,raw_valence_temp
 REAL(kind=dp), ALLOCATABLE, DIMENSION(:) :: raw_valence,nvalence_type,raw_spin_x,raw_spin_y,raw_spin_z,&
 raw_core,raw_valence_pseudodensity,raw_spin,raw_spin_1,raw_spin_2,raw_spin_3,temp_number
 REAL(kind=dp), ALLOCATABLE, DIMENSION(:,:) :: direct_coords,direct_coords2
 REAL(kind=dp), ALLOCATABLE, DIMENSION(:,:,:) :: valence_pseudodensity
 INTEGER :: AECCAR2_FID,iostat_AECCAR2,k,natomtypes,nheader,core_flag,natomtypes2,Q,M,numA,numB,numC,&
 temp,temp_count,totnumA2,totnumB2,totnumC2,valence_grid_correct_flag,z,AECCAR0_FID,CHG_FID,iostat_CHG,&
 col_position,eq_position,component,blank_position,start_index,end_index,first_blank_position,io&
 , totnumA_int, totnumB_int, totnumC_int,NaN_valence,NaN_core,NaN_spin,NaN_pseudodensity
 INTEGER(kind=8) :: fileSize
 INTEGER, ALLOCATABLE, DIMENSION(:) :: atomspertype,atomspertype2,atomic_number_type
 LOGICAL :: valence_grid_correct
 
 INTEGER :: clock_count,clock_count_rate,clock_count_max,clock_count2
 
 WRITE(output_FID,*)' '
 WRITE(output_FID,*)'******************** TIME CONTROL ********************'
 WRITE(output_FID,*)'Starting format_vasp_densities'
 FLUSH(output_FID)
 CALL system_clock(clock_count,clock_count_rate,clock_count_max)
 WRITE(output_FID,*)'******************************************************'
 WRITE(output_FID,*)' '
 FLUSH(output_FID)
 
 !Open and read the valence density file
 OPEN(newunit=AECCAR2_FID, file='AECCAR2', status='old',iostat=iostat_AECCAR2) 
 IF (iostat_AECCAR2 /= 0) THEN 
   WRITE(output_FID,*) 'Could not find AECCAR2 file. Program will terminate.'
   STOP
 ELSE
   WRITE(output_FID,*) 'inputfile = AECCAR2'
 END IF
 natomtypes=0
 READ(AECCAR2_FID,'(A)') headerlines
 READ(AECCAR2_FID,*) latticevectorfactor
 DO i=1,3
   READ(AECCAR2_FID,*)(vectors(i,j),j=1,3)
 END DO
 vectors=vectors*latticevectorfactor*bohrperangstrom
 WRITE(output_FID,*)'vectors'
 DO i=1,3
   WRITE(output_FID,'(3f10.6)')(vectors(i,j),j=1,3)
 END DO
 READ(AECCAR2_FID,'(A)')tempstring
 DO
   tempstring=adjustl(tempstring)
   blank_position=index(tempstring,' ')
   IF (len(trim(tempstring)) == 0) THEN
     EXIT
   ELSE
     natomtypes = natomtypes + 1
   END IF
   tempstring=tempstring(blank_position:)
 END DO
 CLOSE(AECCAR2_FID)
 OPEN(newunit=AECCAR2_FID, file='AECCAR2', status='old',iostat=iostat_AECCAR2) 
 ALLOCATE(atomspertype(natomtypes))
 DO i=1,5
  READ(AECCAR2_FID,*)headerlines
 END DO
 READ(AECCAR2_FID,*,IOSTAT=iostat_value)(atomspertype(i),i=1,natomtypes)
 IF (iostat_value /= 0) THEN
   nheader = 7
   READ(AECCAR2_FID,*)(atomspertype(i),i=1,natomtypes)
 ELSE 
   nheader = 6
 END IF
 natoms=sum(atomspertype)
 ALLOCATE(direct_coords(natoms,3))
 READ(AECCAR2_FID,*)headerlines
 DO i=1,natoms
   READ(AECCAR2_FID,'(3f10.6)')(direct_coords(i,j),j=1,3)
 END DO
 WRITE(output_FID,*)'direct_coords'
 DO i=1,natoms
   WRITE(output_FID,'(3f10.6)')(direct_coords(i,j),j=1,3)
 END DO 
 READ(AECCAR2_FID,*)totnumA,totnumB,totnumC
 WRITE(output_FID,'(A,I6)')' totnumA= ',totnumA
 WRITE(output_FID,'(A,I6)')' totnumB= ',totnumB
 WRITE(output_FID,'(A,I6)')' totnumC= ',totnumC
 FLUSH(output_FID)
 totnum=[totnumA,totnumB,totnumC] 
 M=totnumA*totnumB*totnumC 
 ALLOCATE(raw_valence(M))
 buffer = ceiling(real(M)/read_buffer_size)
 DO i=1,buffer
   start_index = (i-1)*read_buffer_size + 1
   end_index = min((i*read_buffer_size),M)
   READ(AECCAR2_FID,*)(raw_valence(j),j=start_index,end_index)
 END DO
 CLOSE(AECCAR2_FID)
 !Open and read the core density file
 OPEN(newunit=AECCAR0_FID,FILE='AECCAR0',IOSTAT=iostat_value,STATUS='old')
 IF (iostat_value == 0) THEN 
   core_available = .true.   
 ELSE
   core_available = .false.
 END IF
 core_flag=0
 IF (core_available) THEN
   WRITE(output_FID,*)'inputfile = AECCAR0'
   FLUSH(output_FID)
   READ(AECCAR0_FID,'(A)') headerlines
   READ(AECCAR0_FID,*) latticevectorfactor2
   IF (abs(latticevectorfactor-latticevectorfactor2) > zero_tolerance) THEN 
     core_flag=1
   END IF
   DO i=1,3
     READ(AECCAR0_FID,*)(vectors2(i,j),j=1,3)
   END DO
   IF (abs(sum(vectors)-sum(vectors2*latticevectorfactor*bohrperangstrom)) > zero_tolerance ) THEN
     core_flag=1
   END IF
   IF (nheader == 7) THEN
     READ(AECCAR0_FID,*)headerlines
   END IF
   ALLOCATE(atomspertype2(natomtypes))
   READ(AECCAR0_FID,*)(atomspertype2(i),i=1,natomtypes)
   IF (abs(sum(atomspertype)-sum(atomspertype2)) > zero_tolerance) THEN
     core_flag = 1
   END IF
   READ(AECCAR0_FID,*)headerlines
   ALLOCATE(direct_coords2(natoms,3))
   DO i=1,natoms
     READ(AECCAR0_FID,*)(direct_coords2(i,j),j=1,3)
   END DO
   IF (abs(sum(direct_coords - direct_coords2)) > zero_tolerance) THEN
     core_flag = 1
   END IF
   READ(AECCAR0_FID,*)totnumA2,totnumB2,totnumC2
   IF (abs(totnumA+totnumB+totnumC-(totnumA2+totnumB2+totnumC2)) > zero_tolerance) THEN
     core_flag = 1
   END IF
   ALLOCATE(raw_core(M))
 DO i=1,buffer
   start_index = (i-1)*read_buffer_size + 1
   end_index = min((i*read_buffer_size),M)
   READ(AECCAR0_FID,*)(raw_core(j),j=start_index,end_index)
 END DO
   CLOSE(AECCAR0_FID)
 END IF
 IF (core_flag == 1) THEN
   WRITE(output_FID,*)'The core density (AECCAR0) file does not contain the same lattice vectors, grid spacing, or atom &
   &positions as the valence density (AECCAR2) file.'
   WRITE(output_FID,*) 'This can occur if the AECCAR0 and AECCAR2 were produced from a geometry optimization rather than&
   &from a single-point calculation.'
   WRITE(output_FID,*) 'The program will generate its own core densities.'
   core_available = .false.
 END IF  
 WRITE(output_FID,*) 'core_available',core_available
 FLUSH(output_FID)
 !Open and read the valence pseudodensity
 OPEN(NEWUNIT=CHG_FID,FILE='CHGCAR',IOSTAT=iostat_CHG,STATUS='old')
 IF (iostat_CHG == 0) THEN
   CHG_file='CHGCAR'
   WRITE(output_FID,*)'inputfile = CHGCAR'
   FLUSH(output_FID)
 ELSE IF (iostat_CHG /= 0) THEN
   OPEN(NEWUNIT=CHG_FID,FILE='CHGCAR_noncollinear',IOSTAT=iostat_CHG,STATUS='old') 
   IF (iostat_CHG == 0) THEN
     CHG_file='CHGCAR_noncollinear'
     WRITE(output_FID,*)'inputfile = CHGCAR_noncollinear'
     FLUSH(output_FID)
   ELSE IF (iostat_CHG /= 0) THEN
     WRITE(output_FID,*)'Could not read CHGCAR file. Program will try to read CHG file.' 
     OPEN(NEWUNIT=CHG_FID,FILE='CHG',IOSTAT=iostat_CHG,STATUS='old')
     IF (iostat_CHG == 0) THEN
       CHG_file='CHG'
       WRITE(output_FID,*)'inputfile = CHG'
	   FLUSH(output_FID)
     ELSE IF (iostat_CHG /= 0) THEN
	   OPEN(NEWUNIT=CHG_FID,FILE='CHG_noncollinear',IOSTAT=iostat_CHG,STATUS='old')
       IF (iostat_CHG == 0) THEN
         CHG_file='CHG_noncollinear'
         WRITE(output_FID,*)'inputfile = CHG_noncollinear'
         FLUSH(output_FID)
       ELSE IF (iostat_CHG /= 0) THEN
         WRITE(output_FID,*)'Could not open CHGCAR or CHG file.'
       END IF
     END IF
   END IF 
 END IF
 fileSize = 0
 IF (iostat_CHG /= 0) THEN
   WRITE(output_FID,*)'Could not read CHGCAR or CHG file.'
   WRITE(output_FID,*)'This can occur if you forgot to specify LCHARG = .TRUE. in the & 
   &VASP INCAR file.'
   WRITE(output_FID,*)'Please regenerate the CHGCAR or CHG file. Program will terminate.'
   STOP
   valence_grid_correct = .false.
   spin_available = .false. 
 ELSE
   INQUIRE(unit=CHG_FID, size=fileSize)
   IF (fileSize > 10) THEN
     valence_grid_correct = .true.
   ELSE
     WRITE(output_FID,*)'CHGCAR or CHG file appears to be empty.'
	 WRITE(output_FID,*)'This can occur if you forgot to specify LCHARG = .TRUE. in the & 
	 &VASP INCAR file.'
     WRITE(output_FID,*)'Please regenerate the CHGCAR or CHG file. Program will terminate.'
     STOP
     valence_grid_correct = .false.
   END IF
 END IF
 latticevectorfactor2 = 0.0_dp
 atomspertype2 = 0
 direct_coords2 = 0.0_dp
 totnumA2 = 0
 totnumB2 = 0
 totnumC2 = 0
 valence_grid_correct_flag = 0

 IF (valence_grid_correct) THEN
  READ(CHG_FID,*)headerlines
  READ(CHG_FID,*)latticevectorfactor2
  IF (abs(latticevectorfactor-latticevectorfactor2) > zero_tolerance) THEN 
    valence_grid_correct_flag = 1
  END IF 
  DO i=1,3
     READ(CHG_FID,*)(vectors2(i,j),j=1,3)
   END DO
   IF (abs(sum(vectors)-sum(vectors2*latticevectorfactor*bohrperangstrom)) > zero_tolerance ) THEN
     valence_grid_correct_flag = 1
   END IF
   IF (nheader == 7) THEN
     READ(CHG_FID,*)headerlines
   END IF
   READ(CHG_FID,*)(atomspertype2(i),i=1,natomtypes)
   IF (abs(sum(atomspertype)-sum(atomspertype2)) > zero_tolerance) THEN
     valence_grid_correct_flag = 1
   END IF
   READ(CHG_FID,*)headerlines
   DO i=1,natoms
     READ(CHG_FID,*)(direct_coords2(i,j),j=1,3)
   END DO
   IF (abs(sum(direct_coords - direct_coords2)) > zero_tolerance) THEN
     valence_grid_correct_flag = 1
   END IF
   READ(CHG_FID,*)totnumA2,totnumB2,totnumC2
   IF (abs(totnumA+totnumB+totnumC-(totnumA2+totnumB2+totnumC2)) > zero_tolerance) THEN
     valence_grid_correct_flag = 1
   END IF
   ALLOCATE(raw_valence_pseudodensity(M))
   DO i=1,buffer
     start_index = (i-1)*read_buffer_size + 1
     end_index = min((i*read_buffer_size),M)
     READ(CHG_FID,*)(raw_valence_pseudodensity(j),j=start_index,end_index)
   END DO 
 END IF  
 IF (valence_grid_correct_flag == 1) THEN
   WRITE(output_FID,*)'The ',CHG_file,' does not contain the same lattice vectors, grid spacing, or atom positions as the &
   valence density (AECCAR2) file.'
   WRITE(output_FID,*) 'This can occur if the files were produced from a geometry optimization rather than from a fixed geometry &
   &calculation.'
   WRITE(output_FID,*) 'For best results, please regenerate the files using a fixed geometry calculation.'
   WRITE(output_FID,*)'Program will terminate.'
   STOP
   valence_grid_correct = .false.
 END IF 
 WRITE(output_FID,*)'valence_grid_correct =',valence_grid_correct
 FLUSH(output_FID)
 !Construct the density array:
 boundary(1,1)=vectors(1,1)/totnumA
 boundary(1,2)=vectors(1,2)/totnumA
 boundary(1,3)=vectors(1,3)/totnumA
 boundary(2,1)=vectors(2,1)/totnumB
 boundary(2,2)=vectors(2,2)/totnumB 
 boundary(2,3)=vectors(2,3)/totnumB 
 boundary(3,1)=vectors(3,1)/totnumC
 boundary(3,2)=vectors(3,2)/totnumC
 boundary(3,3)=vectors(3,3)/totnumC
 pixelvolume=determinant(boundary)
 vector1=[vectors(1,1),vectors(1,2),vectors(1,3)]
 vector2=[vectors(2,1),vectors(2,2),vectors(2,3)]
 vector3=[vectors(3,1),vectors(3,2),vectors(3,3)]
 ALLOCATE(coords(3,natoms))
 coords=0.0_dp
 DO k=1,natoms
   coords(1,k)=direct_coords(k,1)*vector1(1) + direct_coords(k,2)*vector2(1)+ direct_coords(k,3)*vector3(1)
   coords(2,k)=direct_coords(k,1)*vector1(2) + direct_coords(k,2)*vector2(2)+ direct_coords(k,3)*vector3(2)
   coords(3,k)=direct_coords(k,1)*vector1(3) + direct_coords(k,2)*vector2(3)+ direct_coords(k,3)*vector3(3)
 END DO
 origin = 0.0_dp
 sum_negative_density = 0.0_dp
 ALLOCATE(core_density(totnumA,totnumB,totnumC))
 core_density=0.0_dp
 ALLOCATE(valence_density(totnumA,totnumB,totnumC))
 valence_density=0.0_dp
 IF (valence_grid_correct) THEN
   ALLOCATE(valence_pseudodensity(totnumA,totnumB,totnumC))
   valence_pseudodensity=0.0_dp
 END IF
 max_correction_size = 0.0_dp
 num_corrected_pixels = 0

 WRITE(output_FID,*)' '
 WRITE(output_FID,*)'******************** TIME CONTROL ********************'
 CALL system_clock(clock_count2,clock_count_rate,clock_count_max)
 IF (clock_count_rate == 0) THEN
   WRITE(output_FID,*)'Processor clock not available, the time calculation wil not be performed.'
 ELSE
   seconds = modulo((clock_count2-clock_count),clock_count_max)/real(clock_count_rate)
   WRITE(output_FID,*)'Read input files in ',seconds, 'seconds'
   FLUSH(output_FID)
 END IF
 WRITE(output_FID,*)'******************************************************'
 WRITE(output_FID,*)' '
 FLUSH(output_FID)  
 
!$omp parallel do default(none) &
!$omp private(j,Q,numA,numB,numC,temp,v,correction_size) &
!$omp shared(M,totnumA,totnumB,totnumC,raw_valence,pixelvolume,valence_density,valence_grid_correct,valence_pseudodensity,&
!$omp raw_valence_pseudodensity,core_available,core_density,raw_core) &
!$omp reduction(+:sum_negative_density,num_corrected_pixels) &
!$omp reduction(max: max_correction_size) &
!$omp schedule(static,collapsed_chunk_size)
 DO j = 1,M
   Q = j - 1
   numA = modulo(Q,totnumA)
   temp = nint((Q - numA)/real(totnumA))
   numB = modulo(temp,totnumB)
   numC = nint((temp - numB)/real(totnumB))
   v = max(raw_valence(j),0.0_dp)/(pixelvolume*totnumA*totnumB*totnumC)
   IF (v*pixelvolume > pixel_integration_tolerance) THEN
     correction_size = v*pixelvolume - pixel_integration_tolerance
     max_correction_size = max(max_correction_size,correction_size)
     num_corrected_pixels = num_corrected_pixels +1
     v = pixel_integration_tolerance/pixelvolume
   END IF
   valence_density(numA+1,numB+1,numC+1) = v !valence density
   IF (valence_grid_correct) THEN
     valence_pseudodensity(numA+1,numB+1,numC+1) = raw_valence_pseudodensity(j)/(pixelvolume*totnumA*totnumB*totnumC)
   END IF   
   IF (core_available) THEN
     core_density(numA+1,numB+1,numC+1) = max(raw_core(j),0.0_dp)/(pixelvolume*totnumA*totnumB*totnumC)
   END IF
   IF (raw_valence(j) < 0.0_dp) THEN
     sum_negative_density = sum_negative_density - (pixelvolume*v - raw_valence(j)/(totnumA*totnumB*totnumC))
   END IF
 END DO
!$omp end parallel do
 DEALLOCATE(raw_core)
 DEALLOCATE(raw_valence)
 DEALLOCATE(raw_valence_pseudodensity)
 WRITE(output_FID,'(a,f10.6)')' The maximum pixel electron correction was: ', max_correction_size
 WRITE(output_FID,'(a,i6)')' The number of pixel corrections was: ', num_corrected_pixels
 WRITE(output_FID,'(a,f10.6)')' sum_negative_density:',sum_negative_density
 !Load the atomic numbers and numbers of valence and core electrons
 OPEN(NEWUNIT=potcar_FID,FILE='POTCAR',IOSTAT=iostat_value,STATUS="OLD")
 IF (iostat_value /= 0) THEN
   WRITE(output_FID,*)'Could not open POTCAR file. Program will terminate.'
   STOP
 END IF
 WRITE(output_FID,*)'POTCAR'
 ALLOCATE(nvalence_type(natomtypes))
 ALLOCATE(atomic_number_type(natomtypes))
 nvalence_type = 0.0_dp
 atomic_number_type = 0
 DO i = 1,natomtypes
   READ(potcar_FID,'(A)')tline
   READ(potcar_FID,*)nvalence_type(i)
   READ(potcar_FID,*)tline
   READ(potcar_FID,'(A)')atomic_symbol_string
   eq_position=index(atomic_symbol_string,'=')
   col_position=index(atomic_symbol_string,':')
   atomic_symbol=trim(adjustl(atomic_symbol_string(eq_position+1:col_position-1)))
   CALL atomic_symbol_to_number(atomic_symbol,output_FID,z)
   atomic_number_type(i) = z
   DO
     READ(potcar_FID,'(A)')tline
     IF (trim(adjustl(tline)) == 'End of Dataset') THEN
       EXIT
     END IF
   END DO
 END DO 
 CLOSE(potcar_FID)
 temp_count = 0
 ALLOCATE(atomic_number(natoms))
 ALLOCATE(core_electrons(natoms))
 atomic_number = 0
 core_electrons = 0
 DO i = 1,natomtypes
   DO j = 1,atomspertype(i)
     temp_count = temp_count + 1
     atomic_number(temp_count) = atomic_number_type(i)
     core_electrons(temp_count) = atomic_number(temp_count) - nint(nvalence_type(i))
   END DO
 END DO
 ALLOCATE(missing_core(natoms))
 IF (core_available) THEN
   missing_core = 0
 ELSE
   missing_core = core_electrons
 END IF
 !Add the missing core densities
 ALLOCATE(center_nabc(3,natoms))
 ALLOCATE(center_shift(3,natoms))
 center_nabc = 0
 center_shift = 0.0_dp
 CALL add_missing_core_density()
 !Read the charge and spin information from the CHGCAR or CHG file
 IF (iostat_CHG == 0) THEN
   INQUIRE(unit=CHG_FID, Size=fileSize)
   WRITE(output_FID,*)'The ',CHG_file, ' file size is: ',fileSize, 'bytes.'
   IF ((CHG_file == 'CHG') .or. (CHG_file == 'CHGCAR')) THEN
     IF (fileSize > 10) THEN
       spin_available = .true.
       non_collinear = .false.
     ELSE
       spin_available = .false.
     END IF
   ELSEIF ((CHG_file == 'CHG_noncollinear') .or. (CHG_file =='CHGCAR_noncollinear')) THEN
     IF (fileSize > 10) THEN
       spin_available = .true.
       non_collinear = .true.
     ELSE
       spin_available = .false.
     END IF
   END IF
 END IF
 WRITE(output_FID,*)'non_collinear =',non_collinear
 FLUSH(output_FID)
 !---------------------------------------------
 !Assign totnum to string
 write(totnumA_string,'(i5)')totnumA
 write(totnumB_string,'(i5)')totnumB
 write(totnumC_string,'(i5)')totnumC 
 IF (spin_available) THEN
   IF (non_collinear) THEN
     !Read the 'CHGCAR' file
     temp_count = 0
     DO
       READ(CHG_FID,'(a)',IOSTAT=iostat_value) tempstring
       IF (iostat_value < 0) THEN
         IF (temp_count == 0) THEN
           spin_available = .false.
         END IF
         EXIT
       ELSE IF (iostat_value > 0) THEN
         WRITE(output_FID,*)'Could not read CHGCAR or CHG file. Please re-generate the file and try again.'
		 WRITE(output_FID,*)'Program will terminate.'
         FLUSH(output_FID)
         STOP
       ELSE IF (iostat_value == 0) THEN
	     !adjusting left
         tempstring = ADJUSTL(tempstring)
	     !looking for totnumA
	     first_blank_position = index(tempstring,' ')
         read(tempstring(1:first_blank_position),'(i5)',iostat=io) totnumA_int
         IF (io == 0) THEN
           IF ( abs( totnumA_int-totnumA) < zero_tolerance) THEN
             tempstring = tempstring(first_blank_position:200)
             tempstring = ADJUSTL(tempstring)
	      !looking for totnumB
	      first_blank_position = index(tempstring,' ')
             read(tempstring(1:first_blank_position),'(i5)',iostat=io) totnumB_int
             IF (io==0) THEN
               IF ( abs( totnumB_int-totnumB) < zero_tolerance) THEN
                 tempstring = tempstring(first_blank_position:200)
                 tempstring = ADJUSTL(tempstring)
                 !looking for totnumC
                 first_blank_position = index(tempstring,' ')
                 read(tempstring(1:first_blank_position),'(i5)',iostat=io) totnumC_int
                 IF (io == 0) THEN
                   IF (abs( totnumC_int-totnumC) < zero_tolerance) THEN
                     temp_count = temp_count + 1
                     !Read CHGCAR 
                     IF (temp_count == 1) THEN
                       ALLOCATE(raw_spin_1(totnumA*totnumB*totnumC))
                       READ(CHG_FID,*,iostat=io)(raw_spin_1(i),i=1,totnumA*totnumB*totnumC)
                       IF (io /= 0) THEN
                           IF (CHG_file == 'CHG_noncollinear') THEN
                             WRITE(output_FID,*)'One of the values of the CHG_noncollinear file could not be read.'
                             WRITE(output_FID,*)'This could be due to a large number being printed in a small space in the file.'
                             WRITE(output_FID,*)'This error does not occur in the CHGCAR_noncollinear file.'
                             WRITE(output_FID,*)'Please use the CHGCAR_noncollinear instead of the CHG_&
                             &noncollinear file and try again.'
                             WRITE(output_FID,*)'Program will terminate.'  
                             FLUSH(output_FID)
                             STOP
                           END IF
                       END IF
			ELSE IF (temp_count == 2) THEN
                       ALLOCATE(raw_spin_2(totnumA*totnumB*totnumC))
			  READ(CHG_FID,*,iostat=io)(raw_spin_2(i),i=1,totnumA*totnumB*totnumC)
                       IF (io /= 0) THEN
                         IF (CHG_file == 'CHG_noncollinear') THEN
                             WRITE(output_FID,*)'One of the values of the CHG_noncollinear file could not be read.'
                             WRITE(output_FID,*)'This could be due to a large number being printed in a small space in the file.'
                             WRITE(output_FID,*)'This error does not occur in the CHGCAR_noncollinear file.'
                             WRITE(output_FID,*)'Please use the CHGCAR_noncollinear instead of the CHG_&
                             &noncollinear file and try again.'
                             WRITE(output_FID,*)'Program will terminate.' 
                             FLUSH(output_FID)
                             STOP
                          END IF
                       END IF
                     ELSE IF (temp_count == 3) THEN
                       ALLOCATE(raw_spin_3(totnumA*totnumB*totnumC))
                       READ(CHG_FID,*,iostat=io)(raw_spin_3(i),i=1,totnumA*totnumB*totnumC)
                       IF (io /= 0) THEN
                         IF (CHG_file == 'CHG_noncollinear') THEN
                             WRITE(output_FID,*)'One of the values of the CHG_noncollinear file could not be read.'
                             WRITE(output_FID,*)'This could be due to a large number being printed in a small space in the file.'
                             WRITE(output_FID,*)'This error does not occur in the CHGCAR_noncollinear file.'
                             WRITE(output_FID,*)'Please use the CHGCAR_noncollinear instead of the CHG_&
                             &noncollinear file and try again.'
                             WRITE(output_FID,*)'Program will terminate.' 
                             FLUSH(output_FID)
                             STOP
                         END IF
                       END IF
                     END IF
                   END IF
                 END IF
               END IF
             END IF
           END IF
         END IF
       END IF
     END DO
     IF (spin_available) THEN
       !Compute the x, y, and z spin components using the SAXIS value
       alpha = atan2(SAXIS(2),SAXIS(1))
       beta = atan2(sqrt(SAXIS(1)**2 + SAXIS(2)**2),SAXIS(3))
       IF ((abs(alpha) > 1.58_dp) .or. (abs(beta) > 1.58_dp)) THEN
         WRITE(output_FID,*)'Error reading SAXIS value. Program will terminate.'
         STOP
       END IF
       ALLOCATE(raw_spin_x(M))
       ALLOCATE(raw_spin_y(M))
       ALLOCATE(raw_spin_z(M))
       raw_spin_x = cos(beta)*cos(alpha)*raw_spin_1 - sin(alpha)*raw_spin_2 + sin(beta)*cos(alpha)*raw_spin_3
       raw_spin_y = cos(beta)*sin(alpha)*raw_spin_1 + cos(alpha)*raw_spin_2 + sin(beta)*sin(alpha)*raw_spin_3
       raw_spin_z = -sin(beta)*raw_spin_1 + cos(beta)*raw_spin_3
       DEALLOCATE(raw_spin_1)  
       DEALLOCATE(raw_spin_2)
       DEALLOCATE(raw_spin_3)
       !Construct the spin density matrices
       ALLOCATE(spin_density_vector(3,totnumA,totnumB,totnumC))
       spin_density_vector=0.0_dp
       checkspin1_positive_vector = 0.0_dp
       checkspin2_positive_vector = 0.0_dp
       checkspin1_negative_vector = 0.0_dp
       checkspin2_negative_vector = 0.0_dp
!$omp parallel do default(none) &
!$omp private(j,Q,numA,numB,numC,temp,p,s,enhancement_factor) &
!$omp shared(M,totnumA,totnumB,totnumC,core_density,valence_density,raw_spin_x,raw_spin_y,raw_spin_z,valence_pseudodensity,&
!$omp spin_density,spin_density_vector,pixelvolume) &
!$omp schedule(dynamic,collapsed_chunk_size) &
!$omp reduction(+:checkspin1_positive_vector,checkspin2_positive_vector,checkspin1_negative_vector,checkspin2_negative_vector)
       DO j = 1,M
         Q = j - 1
         numA = modulo(Q,totnumA)
         temp = nint((Q - numA)/real(totnumA))
         numB = modulo(temp,totnumB)
         numC = nint((temp - numB)/real(totnumB))
         p = core_density(numA+1,numB+1,numC+1) + valence_density(numA+1,numB+1,numC+1)
         s = sqrt(raw_spin_x(j)**2 + raw_spin_y(j)**2 + raw_spin_z(j)**2)/(pixelvolume*totnumA*totnumB*totnumC)
         IF ((p <= 0.0_dp) .or. (s == 0.0_dp) .or. (valence_pseudodensity(numA+1,numB+1,numC+1) <= 0.0_dp)) THEN
           spin_density_vector(1,numA+1,numB+1,numC+1) = 0.0_dp
           spin_density_vector(2,numA+1,numB+1,numC+1) = 0.0_dp
           spin_density_vector(3,numA+1,numB+1,numC+1) = 0.0_dp
         ELSE
           enhancement_factor = min(abs(p/valence_pseudodensity(numA+1,numB+1,numC+1)),10.0_dp)
           spin_density_vector(1,numA+1,numB+1,numC+1) = enhancement_factor*raw_spin_x(j)*min(p/(enhancement_factor*s),1.0_dp)/&
           (pixelvolume*totnumA*totnumB*totnumC)
           spin_density_vector(2,numA+1,numB+1,numC+1) = enhancement_factor*raw_spin_y(j)*min(p/(enhancement_factor*s),1.0_dp)/&
           (pixelvolume*totnumA*totnumB*totnumC)
           spin_density_vector(3,numA+1,numB+1,numC+1) = enhancement_factor*raw_spin_z(j)*min(p/(enhancement_factor*s),1.0_dp)/&
           (pixelvolume*totnumA*totnumB*totnumC)
         END IF
         IF (raw_spin_x(j) > 0.0_dp) THEN
           checkspin1_positive_vector(1) = checkspin1_positive_vector(1) + raw_spin_x(j)/(totnumA*totnumB*totnumC)
           checkspin2_positive_vector(1) = checkspin2_positive_vector(1) + spin_density_vector(1,numA+1,numB+1,numC+1)*pixelvolume
         ELSE
           checkspin1_negative_vector(1) = checkspin1_negative_vector(1) - raw_spin_x(j)/(totnumA*totnumB*totnumC)
           checkspin2_negative_vector(1) = checkspin2_negative_vector(1) - spin_density_vector(1,numA+1,numB+1,numC+1)*pixelvolume 
         END IF
         IF (raw_spin_y(j) > 0.0_dp) THEN
           checkspin1_positive_vector(2) = checkspin1_positive_vector(2) + raw_spin_y(j)/(totnumA*totnumB*totnumC)
           checkspin2_positive_vector(2) = checkspin2_positive_vector(2) + spin_density_vector(2,numA+1,numB+1,numC+1)*pixelvolume
         ELSE
           checkspin1_negative_vector(2) = checkspin1_negative_vector(2) - raw_spin_y(j)/(totnumA*totnumB*totnumC)
           checkspin2_negative_vector(2) = checkspin2_negative_vector(2) - spin_density_vector(2,numA+1,numB+1,numC+1)*pixelvolume
         END IF
         IF (raw_spin_z(j) > 0.0_dp) THEN
           checkspin1_positive_vector(3) = checkspin1_positive_vector(3) + raw_spin_z(j)/(totnumA*totnumB*totnumC)
           checkspin2_positive_vector(3) = checkspin2_positive_vector(3) + spin_density_vector(3,numA+1,numB+1,numC+1)*pixelvolume
         ELSE
           checkspin1_negative_vector(3) = checkspin1_negative_vector(3) - raw_spin_z(j)/(totnumA*totnumB*totnumC)
           checkspin2_negative_vector(3) = checkspin2_negative_vector(3) - spin_density_vector(3,numA+1,numB+1,numC+1)*pixelvolume
         END IF 
       END DO
!$omp end parallel do 
       DO i=1,3
         checkspin3_positive_vector=0.0_dp
         checkspin3_negative_vector=0.0_dp
!$omp parallel do default(none) &
!$omp private(ka,kb,kc,s_trial,p,mag_s_trial,component,s) &
!$omp shared(totnumA,totnumB,totnumC,core_density,valence_density,spin_density_vector,checkspin1_positive_vector,&
!$omp checkspin1_negative_vector,checkspin2_positive_vector,checkspin2_negative_vector,pixelvolume) &
!$omp schedule(dynamic,chunk_size) &
!$omp reduction(+:checkspin3_negative_vector,checkspin3_positive_vector)
         DO kc=1,totnumC
           DO kb=1,totnumB
             DO ka=1,totnumA
               p = core_density(ka,kb,kc) + valence_density(ka,kb,kc)
               s = sqrt(spin_density_vector(1,ka,kb,kc)**2 + spin_density_vector(2,ka,kb,kc)**2 + spin_density_vector(3,ka,kb,kc)&
               **2)
               s_trial(1)=spin_density_vector(1,ka,kb,kc)
               s_trial(2)=spin_density_vector(2,ka,kb,kc)
               s_trial(3)=spin_density_vector(3,ka,kb,kc)
               DO component = 1,3
                 IF ((spin_density_vector(component,ka,kb,kc) > 0.0_dp) .and. (checkspin2_positive_vector(component) > &
                 zero_tolerance)) THEN
                   s_trial(component) = spin_density_vector(component,ka,kb,kc)*checkspin1_positive_vector(component)/&
                   checkspin2_positive_vector(component)
                 ELSE IF ((spin_density_vector(component,ka,kb,kc) < 0.0_dp) .and. (checkspin2_negative_vector(component) > &
                 zero_tolerance)) THEN
                   s_trial(component) = spin_density_vector(component,ka,kb,kc)*checkspin1_negative_vector(component)/&
                   checkspin2_negative_vector(component)
                 END IF
               END DO
               mag_s_trial = sqrt(s_trial(1)**2 + s_trial(2)**2 + s_trial(3)**2)
               IF (mag_s_trial > 0.0_dp) THEN
                 spin_density_vector(1,ka,kb,kc) = s_trial(1)*min(p/mag_s_trial,1.0_dp)
                 spin_density_vector(2,ka,kb,kc) = s_trial(2)*min(p/mag_s_trial,1.0_dp)
                 spin_density_vector(3,ka,kb,kc) = s_trial(3)*min(p/mag_s_trial,1.0_dp)
               END IF
               DO component=1,3
                 IF (spin_density_vector(component,ka,kb,kc) > 0.0_dp ) THEN
                   checkspin3_positive_vector(component) = checkspin3_positive_vector(component) + &
                   spin_density_vector(component,ka,kb,kc)*pixelvolume
                 ELSE
                   checkspin3_negative_vector(component) = checkspin3_negative_vector(component) - &
                   spin_density_vector(component,ka,kb,kc)*pixelvolume
                 END IF
               END DO
             END DO
           END DO
         END DO
!$omp end parallel do
         checkspin2_positive_vector=checkspin3_positive_vector
         checkspin2_negative_vector=checkspin3_negative_vector
       END DO
       checkspin3_vector=checkspin3_positive_vector - checkspin3_negative_vector
       WRITE(output_FID,*)'The renormalized spin magnetization vector integrated over the unit cell is: '
       WRITE(output_FID,*)checkspin3_vector
       FLUSH(output_FID)
       checkspin3_mag = sqrt(checkspin3_vector(1)**2 + checkspin3_vector(2)**2 + checkspin3_vector(3)**2)
       WRITE(output_FID,*)'The magnitude of the spin magnetization vector integrated over the unit cell is: '
       WRITE(output_FID,*)checkspin3_mag
       FLUSH(output_FID)
       DEALLOCATE(raw_spin_x)
       DEALLOCATE(raw_spin_y)
       DEALLOCATE(raw_spin_z)
     END IF  
   ELSE
     temp_count = 0
     DO
       READ(CHG_FID,'(a)',IOSTAT=iostat_value) tempstring
       IF (iostat_value < 0) THEN
         IF (temp_count == 0) THEN
           spin_available = .false.
         END IF
         EXIT
       ELSE IF (iostat_value > 0) THEN
         WRITE(output_FID,*)'Could not read CHGCAR or CHG file. Please re-generate the file and try again.'
		 WRITE(output_FID,*)'Program will terminate.'
         FLUSH(output_FID)
         STOP
       ELSE IF (iostat_value == 0) THEN
	     !adjusting left
         tempstring = ADJUSTL(tempstring)
	     !looking for totnumA
	     first_blank_position = index(tempstring,' ')
         read(tempstring(1:first_blank_position),'(i5)',iostat=io) totnumA_int
         IF (io == 0) THEN
           IF ( abs( totnumA_int-totnumA) < zero_tolerance) THEN
             tempstring = tempstring(first_blank_position:200)
             tempstring = ADJUSTL(tempstring)
	         !looking for totnumB
	         first_blank_position = index(tempstring,' ')
             read(tempstring(1:first_blank_position),'(i5)',iostat=io) totnumB_int
             IF (io==0) THEN
               IF ( abs( totnumB_int-totnumB) < zero_tolerance) THEN
                 tempstring = tempstring(first_blank_position:200)
                 tempstring = ADJUSTL(tempstring)
                 !looking for totnumC
                 first_blank_position = index(tempstring,' ')
                 read(tempstring(1:first_blank_position),'(I5)',iostat=io) totnumC_int
                 IF (io == 0) THEN
                   IF (abs( totnumC_int-totnumC) < zero_tolerance) THEN
                     temp_count = temp_count + 1
                     !Read CHGCAR raw valence pseudodensity
                     IF (temp_count == 1) THEN
                       ALLOCATE(raw_spin(M))
                       DO i=1,buffer
                         start_index = (i-1)*read_buffer_size + 1
                         end_index = min((i*read_buffer_size),M)
                         READ(CHG_FID,*,iostat=io)(raw_spin(j),j=start_index,end_index)
                         IF (io /= 0) THEN
                           IF (CHG_file == 'CHG') THEN
                             WRITE(output_FID,*)'One of the values of the CHG file could not be read.' 
                             WRITE(output_FID,*)'This could be due to a large number being printed in a small space in the file.' 
                             WRITE(output_FID,*)'This type of error does not occur in CHGCAR files.'
                             WRITE(output_FID,*)'Please use the CHGCAR instead of the CHG file and &
                             &try again.'
                             WRITE(output_FID,*)'Program will terminate.' 
                             FLUSH(output_FID)
                             STOP
                           END IF
                         END IF
                       END DO
                     ELSE
                       IF (CHG_file == 'CHG') THEN
                         WRITE(output_FID,*)'The file read has noncollinear magnetism. Please set the file name to &
                         &CHG_noncollinear and try again. Program will terminate.'
                       ELSE IF (CHG_file == 'CHGCAR') THEN
                         WRITE(output_FID,*)'The file read has noncollinear magnetism. Please set the file name to &
                         &CHGCAR_noncollinear and try again. Program will terminate.'
                       END IF
                       STOP
                     END IF
                   END IF
                 END IF
               END IF
             END IF
           END IF
         END IF
       END IF
     END DO
     IF (spin_available) THEN
       ALLOCATE(spin_density(totnumA,totnumB,totnumC))
       spin_density = 0.0_dp
       checkspin1_positive = 0.0_dp
       checkspin2_positive = 0.0_dp
       checkspin1_negative = 0.0_dp
       checkspin2_negative = 0.0_dp
!$omp parallel do default(none) &
!$omp private(j,Q,numA,numB,numC,temp,p,s) &
!$omp shared(M,totnumA,totnumB,totnumC,core_density,valence_density,raw_spin,valence_pseudodensity,spin_density,pixelvolume) &
!$omp schedule(dynamic,collapsed_chunk_size) &
!$omp reduction(+:checkspin1_positive,checkspin2_positive,checkspin1_negative,checkspin2_negative)     
       DO j = 1,M
         Q = j - 1
         numA = modulo(Q,totnumA)
         temp = nint((Q - numA)/real(totnumA))
         numB = modulo(temp,totnumB)
         numC = nint((temp - numB)/real(totnumB))
         p = core_density(numA+1,numB+1,numC+1) + valence_density(numA+1,numB+1,numC+1)
         IF ((p <= 0.0_dp) .or. (raw_spin(j) == 0.0_dp) .or. (valence_pseudodensity(numA+1,numB+1,numC+1) <= 0.0_dp)) THEN
           spin_density(numA+1,numB+1,numC+1) = 0.0_dp
         ELSE
           s = raw_spin(j)/(pixelvolume*totnumA*totnumB*totnumC)
           s = s*min(abs(p/valence_pseudodensity(numA+1,numB+1,numC+1)),10.0_dp)
           spin_density(numA+1,numB+1,numC+1) = s*min(abs(p)/abs(s),1.0_dp)
         END IF
         IF (raw_spin(j) > 0.0_dp) THEN
           checkspin1_positive = checkspin1_positive + raw_spin(j)/(totnumA*totnumB*totnumC)
           checkspin2_positive = checkspin2_positive + spin_density(numA+1,numB+1,numC+1)*pixelvolume
         ELSE
           checkspin1_negative = checkspin1_negative - raw_spin(j)/(totnumA*totnumB*totnumC)
           checkspin2_negative = checkspin2_negative - spin_density(numA+1,numB+1,numC+1)*pixelvolume   
         END IF
       END DO
!$omp end parallel do
       DO i=1,3
         checkspin3_positive = 0.0_dp
         checkspin3_negative = 0.0_dp
!$omp parallel do default(none) &
!$omp private(ka,kb,kc,s,p) &
!$omp shared(totnumA,totnumB,totnumC,core_density,valence_density,spin_density,checkspin1_positive,checkspin1_negative,&
!$omp checkspin2_positive,checkspin2_negative,pixelvolume) &
!$omp schedule(dynamic,chunk_size) &
!$omp reduction(+:checkspin3_negative,checkspin3_positive)
         DO kc=1,totnumC
           DO kb=1,totnumB
             DO ka=1,totnumA
               p = core_density(ka,kb,kc) + valence_density(ka,kb,kc)
               IF ((spin_density(ka,kb,kc) > 0.0_dp) .and. (checkspin2_positive > zero_tolerance)) THEN
                 s = spin_density(ka,kb,kc)*checkspin1_positive/checkspin2_positive
                 spin_density(ka,kb,kc) = s*min(abs(p)/abs(s),1.0_dp)
                 checkspin3_positive = checkspin3_positive + spin_density(ka,kb,kc)*pixelvolume
               ELSE IF ((spin_density(ka,kb,kc) < 0.0_dp) .and. (checkspin2_negative > zero_tolerance)) THEN
                 s = spin_density(ka,kb,kc)*checkspin1_negative/checkspin2_negative
                 spin_density(ka,kb,kc) = s*min(abs(p)/abs(s),1.0_dp)
                 checkspin3_negative = checkspin3_negative - spin_density(ka,kb,kc)*pixelvolume
               END IF
             END DO
           END DO
         END DO
!$omp end parallel do
         checkspin2_positive=checkspin3_positive
         checkspin2_negative=checkspin3_negative
       END DO
       checkspin3=checkspin3_positive - checkspin3_negative
       WRITE(output_FID,*)'The renormalized spin magnetization scalar integrated over the unit cell is: '
       WRITE(output_FID,*)checkspin3
       DEALLOCATE(raw_spin)
     END IF
   END IF
 END IF
 WRITE(output_FID,*)'spin_available =',spin_available
 FLUSH(output_FID)
 CLOSE(CHG_FID)
 !Check if any of the input VASP files contains NaN 
 NaN_valence = 0
 NaN_core = 0
 NaN_spin = 0
 NaN_pseudodensity = 0
 DO kc = 1,totnumC
   DO kb = 1,totnumB
     DO ka = 1,totnumA
       IF (isnan(valence_density(ka,kb,kc)) .or. (abs(valence_density(ka,kb,kc)) > 1.0D12) ) THEN
         NaN_valence = NaN_valence + 1
       END IF
       IF (isnan(core_density(ka,kb,kc)) .or. (abs(core_density(ka,kb,kc)) > 1.0D12) ) THEN
         NaN_core = NaN_core + 1
       END IF
	   IF (spin_available) THEN
         IF (isnan(valence_pseudodensity(ka,kb,kc)) .or. (abs(valence_pseudodensity(ka,kb,kc)) > 1.0D12) ) THEN
           NaN_pseudodensity = NaN_pseudodensity + 1
         END IF
         IF (non_collinear) THEN
           DO i=1,3
             IF (isnan(spin_density_vector(i,ka,kb,kc)) .or. (abs(spin_density_vector(i,ka,kb,kc)) > 1.0D12) ) THEN
               NaN_spin = NaN_spin + 1
             END IF
           END DO
         ELSE
           IF (isnan(spin_density(ka,kb,kc)) .or. (abs(spin_density(ka,kb,kc)) > 1.0D12) ) THEN
             NaN_spin = NaN_spin + 1
           END IF
         END IF
       END IF
     END DO
   END DO
 END DO
 IF (NaN_valence > 0) THEN
   WRITE(output_FID,*)'Your AECCAR2 file contains ',NaN_valence,' NaN entries or numbers with absolute values greater than 1.0D+12.'
   WRITE(output_FID,*)'Please generate a new AECCAR2 file.'
   WRITE(output_FID,*)'Program will terminate.'
   STOP
 END IF
 IF (NaN_core > 0) THEN
   WRITE(output_FID,*)'Your AECCAR0 file contains ',NaN_core,' NaN entries or numbers with absolute values greater than 1.0D+12.'
   WRITE(output_FID,*)'Please generate a new AECCAR0 file.'
   WRITE(output_FID,*)'Program will terminate.'
   STOP
 END IF 
 IF (NaN_pseudodensity > 0) THEN
   WRITE(output_FID,*)'Your ',CHG_file,' file contains ',NaN_pseudodensity,' NaN entries &
   &(or numbers with absolute values greater than 1.0D+12) in the pseudodensity grids.'
   WRITE(output_FID,*)'Please generate a new ',CHG_file,' file.'
   WRITE(output_FID,*)'Program will terminate.'
   STOP
 END IF 
 IF (NaN_spin > 0) THEN
   WRITE(output_FID,*)'Your ',CHG_file,' file contains ',NaN_pseudodensity,' NaN entries &
   &(or numbers with absolute values greater than 1.0D+12) in the spin density grids.'
   WRITE(output_FID,*)'Please generate a new ',CHG_file,' file.'
   WRITE(output_FID,*)'Program will terminate.'
   STOP
 END IF 
 CALL initialize_atomic_densities()
 !Compute the valence occupancy corrections
 ALLOCATE(occupancy_correction(11,natoms))
 occupancy_correction=0.0_dp
 IF (valence_grid_correct) THEN
   ALLOCATE(dominant_atom_weight(nshells,natoms))
   dominant_atom_weight = neutral_density
   CALL compute_dominant_atom_volumes() 
   !Free some space
   DEALLOCATE(dominant_atom_weight)
!$omp parallel do default(none) &
!$omp private(kc,kb,ka) &
!$omp shared(totnumA,totnumB,totnumC,dominant_atom_points,pixelvolume,valence_pseudodensity,&
!$omp valence_density,chunk_size,output_FID) &
!$omp reduction(+:occupancy_correction) &
!$omp schedule(static,chunk_size)
   DO kc = 1,totnumC
     DO kb = 1,totnumB
       DO ka = 1,totnumA
         IF (dominant_atom_points(ka,kb,kc) > 0) THEN
           occupancy_correction(1,dominant_atom_points(ka,kb,kc)) = occupancy_correction(1,dominant_atom_points&
         (ka,kb,kc)) + pixelvolume*(valence_pseudodensity(ka,kb,kc) - valence_density(ka,kb,kc))
         END IF   
       END DO
     END DO
   END DO
!$omp end parallel do
   WRITE(output_FID,*)'The largest occupancy correction was'
   WRITE(output_FID,*) maxval(abs(occupancy_correction))
 END IF
 !Free up some memory by clearing unneeded arrays.    
 DEALLOCATE(valence_pseudodensity)
 DEALLOCATE(dominant_atom_points)
 DEALLOCATE(direct_coords)

 WRITE(output_FID,*)' '
 WRITE(output_FID,*)'******************** TIME CONTROL ********************'
 CALL system_clock(clock_count2,clock_count_rate,clock_count_max)
 IF (clock_count_rate == 0) THEN
   WRITE(output_FID,*)'Processor clock not available, the time calculation wil not be performed.'
 ELSE
   seconds = modulo((clock_count2-clock_count),clock_count_max)/real(clock_count_rate)
   WRITE(output_FID,*)'Finished format_vasp_densities in ',seconds, 'seconds'
   FLUSH(output_FID)
 END IF
 WRITE(output_FID,*)'******************************************************'
 WRITE(output_FID,*)' '
 FLUSH(output_FID)

 END SUBROUTINE format_vasp_densities
 
 END MODULE module_format_vasp_densities