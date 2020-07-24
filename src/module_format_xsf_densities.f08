 MODULE module_format_xsf_densities
 ! Copyright (c) 2014, 2015, 2016, 2017 by Thomas A. Manz and Nidia Gabaldon Limas. Rights reserved.
 
 USE module_common_variable_declarations
 USE module_global_parameters
 USE module_matrix_operations
 USE module_add_missing_core_density
 USE module_initialize_atomic_densities

 IMPLICIT NONE
 
 CONTAINS
 
 SUBROUTINE format_xsf_densities( )
 !=================================================================================== 
 
 CHARACTER(200) :: line_of_text
 REAL(kind=dp) :: vectors(3,3),v,correction_size,p,s,checkspin
 REAL(kind=dp),ALLOCATABLE,DIMENSION(:) :: raw_valence,raw_spin,raw_spin_1,raw_spin_2,raw_spin_3,raw_spin_4,raw_spin_x,raw_spin_y&
 ,raw_spin_z,temp
 INTEGER :: xsf_FID,fileSize,spins_loaded,found_primvec,natomtypes,iostat_value1,npoints,counter,Q,numA,numB,&
 numC,temp_var,start_index,end_index
 
 INTEGER :: clock_count,clock_count_rate,clock_count_max,clock_count2
 REAL(kind=dp) :: seconds
 
 WRITE(output_FID,*)' '
 WRITE(output_FID,*)'******************** TIME CONTROL ********************'
 WRITE(output_FID,*)'Starting format_xsf_densities'
 FLUSH(output_FID)
 CALL system_clock(clock_count,clock_count_rate,clock_count_max)
 WRITE(output_FID,*)'******************************************************'
 WRITE(output_FID,*)' '
 FLUSH(output_FID)

 !Initialize variables
 !Open and read the xsf file containing the valence density
 OPEN(NEWUNIT = xsf_FID,FILE = long_input_filename,IOSTAT=iostat_value,STATUS='old')
 IF (iostat_value /= 0) THEN
   WRITE(output_FID,*)'Could not find input xsf file. Program will terminate.'
   FLUSH(output_FID)
   STOP
 ELSE
   INQUIRE(unit=xsf_FID, size=fileSize)
   IF (fileSize < 10) THEN
     WRITE(output_FID,*)'Input xsf file appears to be empty. Program will terminate.'
     FLUSH(output_FID)
     STOP
   END IF
 END IF
 DO
   READ(xsf_FID,'(a)',IOSTAT=iostat_value) line_of_text
   IF (iostat_value > 0) THEN
     WRITE(output_FID,*) "Could not read xsf file. Program will terminate."
     FLUSH(output_FID)
     STOP
   ELSE IF (iostat_value < 0) THEN
    EXIT
   ELSE
     !adjusting left, lowercase, trim data and read
     line_of_text = SimpleLine(line_of_text)
     line_of_text = StrLowCase(line_of_text)
     !assinging a value to each item
     IF(line_of_text == 'atoms') THEN
       num_periodic_directions = 0
       periodicA = .false.
       periodicB = .false.
       periodicC = .false.
     ELSE IF (line_of_text == 'polymer') THEN
       num_periodic_directions = 1
       periodicA = .true.
       periodicB = .false.
       periodicC = .false.
     ELSE IF (line_of_text == 'slab') THEN
       num_periodic_directions = 2
       periodicA = .true.
       periodicB = .true.
       periodicC = .false.
     ELSE IF(line_of_text == 'crystal') THEN 
       periodicA = .true.
       periodicB = .true.
       periodicC = .true.
     END IF
   END IF
 END DO
 CLOSE(xsf_FID)
 OPEN(NEWUNIT = xsf_FID,FILE = long_input_filename,IOSTAT=iostat_value,STATUS='old')
 ALLOCATE(periodic_vectors(num_periodic_directions,3))
 spins_loaded = 0  
 DO
   READ(xsf_FID,'(a)',IOSTAT=iostat_value) line_of_text
   IF (iostat_value > 0) THEN
     WRITE(output_FID,*) "Could not read xsf file. Program will terminate."
     FLUSH(output_FID)
     STOP
   ELSE IF (iostat_value < 0) THEN
    EXIT
   ELSE
     !adjusting left, lowercase, trim data and read
     line_of_text = SimpleLine(line_of_text)
     line_of_text = StrLowCase(line_of_text)
     !assinging a value to each item
     IF(line_of_text == 'primvec') THEN
       DO i=1,num_periodic_directions
         READ(xsf_FID,*)(periodic_vectors(i,j),j=1,3)
       END DO
       periodic_vectors=periodic_vectors*bohrperangstrom
       WRITE(output_FID,*)'periodic_vectors: '
       DO i=1,num_periodic_directions
         WRITE(output_FID,*)(periodic_vectors(i,j),j=1,3)
       END DO
       found_primvec = 1
     ELSE IF (line_of_text == 'primcoord') THEN
       READ(xsf_FID,*)natoms
       ALLOCATE(atomic_number(natoms))
       ALLOCATE(coords(3,natoms))
       DO i=1,natoms
         READ(xsf_FID,*)atomic_number(i),(coords(j,i),j=1,3)
       END DO
       coords = coords * bohrperangstrom
       CYCLE
     ELSE IF (line_of_text == 'atoms') THEN
       ALLOCATE(atomic_number(10000))
       ALLOCATE(coords(3,10000))
       natoms = 0
       DO i=1,10000
         READ(xsf_FID,*,IOSTAT=iostat_value1)atomic_number(i),(coords(j,i),j=1,3)
         IF (iostat_value1 /= 0) THEN
           EXIT
         END IF
         natoms = natoms + 1
       END DO
     ELSE IF(line_of_text == 'begin_block_datagrid_3d') THEN 
       READ(xsf_FID,*)line_of_text
     ELSE IF (line_of_text == 'begin_datagrid_3d_rho:spin_1') THEN
         READ(xsf_FID,*,IOSTAT=iostat_value1)totnumA,totnumB,totnumC
         IF (iostat_value1 /= 0) THEN
         WRITE(output_FID,*)'Error reading xsf file. Program will terminate.'
           FLUSH(output_FID)
           STOP
         ELSE
           totnumA = totnumA - 1
           totnumB = totnumB - 1
           totnumC = totnumC - 1
         END IF
         READ(xsf_FID,*)(origin(i),i=1,3)
         origin = origin *bohrperangstrom
         DO i=1,3
           READ(xsf_FID,*)(vectors(i,j),j=1,3)
         END DO
         vectors = vectors * bohrperangstrom
         IF ((num_periodic_directions == 3) .and. (found_primvec == 1)) THEN
           IF (sum(abs(periodic_vectors - vectors)) > zero_tolerance) THEN
             WRITE(output_FID,*)'The periodic vectors and the grid vectors must match for a periodic system.'
             WRITE(output_FID,*)'Program will terminate. Please prepare a new input file and try again.'
             FLUSH(output_FID)
             STOP
           END IF
         END IF
         npoints = (totnumA+1)*(totnumB+1)*(totnumC+1)
         ALLOCATE(raw_spin_1(npoints))
         buffer = ceiling(real(npoints)/read_buffer_size)
         DO i=1,buffer
           start_index = (i-1)*read_buffer_size + 1
           end_index = min((i*read_buffer_size),npoints)
           READ(xsf_FID,*)(raw_spin_1(j),j=start_index,end_index)
         END DO
         spins_loaded = spins_loaded + 1
         WRITE(output_FID,*)'spins_loaded= ',spins_loaded
         FLUSH(output_FID)
     ELSE IF(line_of_text == 'begin_datagrid_3d_rho:spin_2') THEN
         READ(xsf_FID,*)totnumA,totnumB,totnumC
         totnumA = totnumA - 1
         totnumB = totnumB - 1
         totnumC = totnumC - 1
         READ(xsf_FID,*)(origin(i),i=1,3)
         origin = origin * bohrperangstrom
         DO i=1,3
           READ(xsf_FID,*)(vectors(i,j),j=1,3)
         END DO
         vectors = vectors * bohrperangstrom
         IF ((num_periodic_directions ==3) .and. (found_primvec == 1)) THEN
           IF (sum(abs(periodic_vectors - vectors)) > zero_tolerance) THEN
             WRITE(output_FID,*)'The periodic vectors and the grid vectors must match for a periodic system.'
             FLUSH(output_FID)
             WRITE(output_FID,*)'Program will terminate. Please prepare a new input file and try again.'
             FLUSH(output_FID)
             STOP
           END IF
         END IF
         npoints=(totnumA+1)*(totnumB+1)*(totnumC+1)
         ALLOCATE(raw_spin_2(npoints))
         DO i=1,buffer
           start_index = (i-1)*read_buffer_size + 1
           end_index = min((i*read_buffer_size),npoints)
           READ(xsf_FID,*)(raw_spin_2(j),j=start_index,end_index)
         END DO
         spins_loaded = spins_loaded +1
         WRITE(output_FID,*)'spins_loaded= ',spins_loaded
         FLUSH(output_FID)
     ELSE IF (line_of_text == 'begin_datagrid_3d_rho:spin_3') THEN
         READ(xsf_FID,*)totnumA,totnumB,totnumC
         totnumA = totnumA - 1
         totnumB = totnumB - 1
         totnumC = totnumC - 1
         READ(xsf_FID,*)(origin(i),i=1,3)
         origin = origin * bohrperangstrom
         DO i=1,3
           READ(xsf_FID,*)(vectors(i,j),j=1,3)
         END DO
         vectors = vectors * bohrperangstrom
         IF ((num_periodic_directions ==3) .and. (found_primvec == 1)) THEN
           IF (sum(abs(periodic_vectors - vectors)) > zero_tolerance) THEN
             WRITE(output_FID,*)'The periodic vectors and the grid vectors must match for a periodic system.'
             FLUSH(output_FID)
             WRITE(output_FID,*)'Program will terminate. Please prepare a new input file and try again.'
             FLUSH(output_FID)
             STOP
           END IF
         END IF
         npoints=(totnumA+1)*(totnumB+1)*(totnumC+1)
         ALLOCATE(raw_spin_3(npoints))
         DO i=1,buffer
           start_index = (i-1)*read_buffer_size + 1
           end_index = min((i*read_buffer_size),npoints)
           READ(xsf_FID,*)(raw_spin_3(j),j=start_index,end_index)
         END DO
         spins_loaded = spins_loaded +1
         WRITE(output_FID,*)'spins_loaded= ',spins_loaded
         FLUSH(output_FID)
     ELSE IF (line_of_text == 'begin_datagrid_3d_rho:spin_4') THEN
         READ(xsf_FID,*)totnumA,totnumB,totnumC
         totnumA = totnumA - 1
         totnumB = totnumB - 1
         totnumC = totnumC - 1
         READ(xsf_FID,*)(origin(i),i=1,3)
         origin = origin * bohrperangstrom
         DO i=1,3
           READ(xsf_FID,*)(vectors(i,j),j=1,3)
         END DO
         vectors = vectors * bohrperangstrom
         IF ((num_periodic_directions ==3) .and. (found_primvec == 1)) THEN
           IF (sum(abs(periodic_vectors - vectors)) > zero_tolerance) THEN
             WRITE(output_FID,*)'The periodic vectors and the grid vectors must match for a periodic system.'
             FLUSH(output_FID)
             WRITE(output_FID,*)'Program will terminate. Please prepare a new input file and try again.'
             FLUSH(output_FID)
             STOP
           END IF
         END IF
         npoints=(totnumA+1)*(totnumB+1)*(totnumC+1)
         ALLOCATE(raw_spin_4(npoints))
         DO i=1,buffer
           start_index = (i-1)*read_buffer_size + 1
           end_index = min((i*read_buffer_size),npoints)
           READ(xsf_FID,*)(raw_spin_4(j),j=start_index,end_index)
         END DO
         spins_loaded = spins_loaded +1
         WRITE(output_FID,*)'spins_loaded= ',spins_loaded
         FLUSH(output_FID)
     END IF
   END IF
 END DO
 CLOSE(xsf_FID)
 WRITE(output_FID,*)'totnumA= ',totnumA
 FLUSH(output_FID)
 WRITE(output_FID,*)'totnumB= ',totnumB
 FLUSH(output_FID)
 WRITE(output_FID,*)'totnumC= ',totnumC
 FLUSH(output_FID)
 WRITE(output_FID,*)'origin= ',origin
 FLUSH(output_FID)
 WRITE(output_FID,*)'vectors= '
 FLUSH(output_FID)
 DO i=1,3
   WRITE(output_FID,*)(vectors(i,j),j=1,3)
   FLUSH(output_FID)
 END DO
 WRITE(output_FID,*)'natoms= ',natoms
 WRITE(output_FID,*)'coords: '
 DO i=1,natoms
   WRITE(output_FID,*)(coords(j,i),j=1,3)
 END DO
 FLUSH(output_FID)
 temp = [natoms,3]   
 natoms = temp(1)
 !Construct the spin density matrices from the individual spin components
 IF (spins_loaded == 1) THEN
   spin_available = .false.
   WRITE(output_FID,*)'No spin moment analysis will be performed.'
   FLUSH(output_FID)
   ALLOCATE(raw_valence(npoints))
   raw_valence = raw_spin_1
   DEALLOCATE(raw_spin_1)
 ELSE IF(spins_loaded == 2) THEN
   spin_available = .true.
   non_collinear = .false.
   WRITE(output_FID,*)'Collinear spin moment analysis will be performed.'
   FLUSH(output_FID)
   ALLOCATE(raw_valence(npoints))
   raw_valence = raw_spin_1 + raw_spin_2
   ALLOCATE(raw_spin(npoints))
   raw_spin = raw_spin_1 - raw_spin_2
   DEALLOCATE(raw_spin_1)
   DEALLOCATE(raw_spin_2)
 ELSE IF (spins_loaded == 4) THEN
   WRITE(output_FID,*)'Warning: Non-collinear spin feature is still being tested. Results may not be accurate!'
   FLUSH(output_FID)
   spin_available = .true.
   non_collinear = .true.
   WRITE(output_FID,*)'Noncollinear spin moment analysis will be performed.'
   FLUSH(output_FID)
   ALLOCATE(raw_valence(npoints))
   ALLOCATE(raw_spin_x(npoints))
   ALLOCATE(raw_spin_y(npoints))
   ALLOCATE(raw_spin_z(npoints))
   raw_valence = raw_spin_1 + raw_spin_2
   raw_spin_x = 2.0_dp*raw_spin_3
   raw_spin_y = -2.0_dp*raw_spin_4
   raw_spin_z = raw_spin_1 - raw_spin_2
   DEALLOCATE(raw_spin_1)
   DEALLOCATE(raw_spin_2)
   DEALLOCATE(raw_spin_3)
   DEALLOCATE(raw_spin_4)
 ELSE
   WRITE(output_FID,*)'Error reading the density and spin information in the xsf file.'
   FLUSH(output_FID)
   WRITE(output_FID,*)'Program will terminate.'
   FLUSH(output_FID)
   STOP
 END IF    
 !Initialize the valence density arrays
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
 sum_negative_density = 0.0_dp
 ALLOCATE(valence_density(totnumA,totnumB,totnumC))
 valence_density = 0.0_dp
 max_correction_size = 0.0_dp
 num_corrected_pixels = 0
!$omp parallel do default(none) &
!$omp private(j,Q,numA,numB,numC,temp_var,v,correction_size) &
!$omp shared(totnumA,totnumB,totnumC,raw_valence,pixelvolume,npoints,valence_density) &
!$omp reduction(+:sum_negative_density,num_corrected_pixels) &
!$omp reduction(max:max_correction_size) &
!$omp schedule(static,collapsed_chunk_size)
 DO j = 1,npoints
    Q = j - 1
    numA = modulo(Q,(totnumA+1))
    temp_var = nint((Q - numA)/real(totnumA + 1))
    numB = modulo(temp_var,(totnumB+1))
    numC = nint((temp_var - numB)/real(totnumB + 1))
    v=max(raw_valence(j),0.0_dp)
    IF (v*pixelvolume > pixel_integration_tolerance) THEN
      correction_size = v*pixelvolume - pixel_integration_tolerance
      max_correction_size = max(max_correction_size,correction_size)
      num_corrected_pixels = num_corrected_pixels + 1
    END IF
    IF (raw_valence(j) < 0.0_dp) THEN
      sum_negative_density = sum_negative_density - (pixelvolume*v - raw_valence(j)/(totnumA*totnumB*totnumC))
    END IF
    numA = modulo(numA,totnumA)
    numB = modulo(numB,totnumB)
    numC = modulo(numC,totnumC)
    valence_density(numA+1,numB+1,numC+1) = v !valence density   
 END DO
!$omp end parallel do
 DEALLOCATE(raw_valence)
 WRITE(output_FID,*)'The maximum pixel electron correction was: ', max_correction_size
 WRITE(output_FID,*)'The number of pixel corrections was: ',num_corrected_pixels
 WRITE(output_FID,*)'sum_negative_density: ',sum_negative_density
 FLUSH(output_FID)
 ALLOCATE(core_electrons(natoms))
 DO j = 1,natoms
  core_electrons(j) = num_core(atomic_number(j))
 END DO
 core_available = .false.
 ALLOCATE(core_density(totnumA,totnumB,totnumC))
 core_density = 0.0_dp
 !Assign the number of core electrons and effective nuclear charge if core fitting is not performed
 IF (.not. core_available) THEN
   DO j = 1,natoms
     core_electrons(j) = num_core(atomic_number(j))
   END DO
   ALLOCATE(missing_core(natoms))
   missing_core=core_electrons
 END IF
 ALLOCATE(center_nabc(3,natoms))
 ALLOCATE(center_shift(3,natoms))
 center_nabc = 0.0_dp
 center_shift = 0.0_dp
 CALL add_missing_core_density( )
 !Process the spin densities
 IF (spin_available) THEN
   IF (non_collinear) THEN
     ALLOCATE(spin_density_vector(3,totnumA,totnumB,totnumC))
     spin_density_vector = 0.0_dp
!$omp parallel do default(none) &
!$omp private(j,Q,numA,numB,numC,temp_var,p,s) &
!$omp shared(npoints,totnumA,totnumB,totnumC,core_density,valence_density,raw_spin_x,raw_spin_y,raw_spin_z,spin_density_vector) &
!$omp schedule(static,collapsed_chunk_size)
     DO j = 1,npoints
       Q = j - 1
       numA = modulo(Q,(totnumA+1))
       temp_var = nint((Q - numA)/real((totnumA + 1)))
       numB = modulo(temp_var,(totnumB+1))
       numC = nint((temp_var - numB)/real((totnumB + 1)))
       numA = modulo(numA,totnumA)
       numB = modulo(numB,totnumB)
       numC = modulo(numC,totnumC)
       p = core_density(numA+1,numB+1,numC+1) + valence_density(numA+1,numB+1,numC+1)
       s = sqrt(raw_spin_x(j)**2 + raw_spin_y(j)**2 + raw_spin_z(j)**2)
       IF ((abs(p) < zero_tolerance) .or. (abs(s) < zero_tolerance)) THEN
         spin_density_vector(1,numA+1,numB+1,numC+1) = 0.0_dp
         spin_density_vector(2,numA+1,numB+1,numC+1) = 0.0_dp
         spin_density_vector(3,numA+1,numB+1,numC+1) = 0.0_dp
       ELSE           
         spin_density_vector(1,numA+1,numB+1,numC+1) = raw_spin_x(j)*min(p/s,1.0_dp)
         spin_density_vector(2,numA+1,numB+1,numC+1) = raw_spin_y(j)*min(p/s,1.0_dp)
         spin_density_vector(3,numA+1,numB+1,numC+1) = raw_spin_z(j)*min(p/s,1.0_dp)
       END IF   
     END DO
!$omp end parallel do
     DEALLOCATE(raw_spin_x)
     DEALLOCATE(raw_spin_y)
     DEALLOCATE(raw_spin_z) 
   ELSE
     ALLOCATE(spin_density(totnumA,totnumB,totnumC))
     spin_density = 0.0_dp
!$omp parallel do default(none) &
!$omp private(j,Q,numA,numB,numC,temp_var,p,s) &
!$omp shared(npoints,totnumA,totnumB,totnumC,core_density,valence_density,raw_spin,spin_density) &
!$omp schedule(static,collapsed_chunk_size)
    DO j = 1,npoints
        Q = j - 1
        numA = modulo(Q,(totnumA+1))
        temp_var = nint((Q - numA)/real((totnumA + 1)))
        numB = modulo(temp_var,(totnumB+1))
        numC = nint((temp_var - numB)/real((totnumB + 1)))
        numA = modulo(numA,totnumA)
        numB = modulo(numB,totnumB)
        numC = modulo(numC,totnumC)
        p = core_density(numA+1,numB+1,numC+1) + valence_density(numA+1,numB+1,numC+1)
        s = raw_spin(j)
        IF ((abs(p) < zero_tolerance) .or. (abs(s) < zero_tolerance)) THEN
            spin_density(numA+1,numB+1,numC+1) = 0.0_dp
        ELSE
            spin_density(numA+1,numB+1,numC+1) = s*min(abs(p)/abs(s),1.0_dp)
        END IF    
    END DO
!$omp end parallel do
    checkspin=sum(spin_density)*pixelvolume
    WRITE(output_FID,*)'checkspin: ',checkspin
    FLUSH(output_FID)
    DEALLOCATE(raw_spin)
   END IF
 END IF
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
   WRITE(output_FID,*)'Finished format_xsf_densities in ',seconds, 'seconds'
   FLUSH(output_FID)
 END IF
 WRITE(output_FID,*)'******************************************************'
 WRITE(output_FID,*)' '
 FLUSH(output_FID) 

 END SUBROUTINE format_xsf_densities

 END MODULE module_format_xsf_densities

