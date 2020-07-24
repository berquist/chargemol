 MODULE module_read_job_control
 ! Copyright (c) 2014, 2015, 2016, 2017 by Thomas A. Manz and Nidia Gabaldon Limas. Rights reserved.
 
 USE module_precision
 USE module_common_variable_declarations
 USE module_global_parameters
 USE module_string_utilities

 IMPLICIT NONE
 
 CONTAINS
 
 SUBROUTINE read_job_control()
 !--------------------------------------------------------------
 ! Reads job_control.txt to find the atomic densities' directory
 ! and re-write default variables
 !--------------------------------------------------------------
 
 CHARACTER(200) :: line_of_text,headerlines
 INTEGER :: iostat_value1,i,read_core_electrons
 output_filename = 'output.'
 long_input_filename = ' '
 netcharge = 0.0_dp
 distance_scale = 1.0_dp
 density_scale = 1.0_dp
 OPEN(newunit=job_control_FID,file='job_control.txt',status='old')
 DO
   READ(job_control_FID,'(a)',IOSTAT=iostat_value) line_of_text
   IF (iostat_value > 0) THEN
     OPEN(newunit = output_FID,file='output.txt',status='replace')
     WRITE(output_FID,*) 'Could not read job_control.txt file. Program will terminate.'
     FLUSH(output_FID)
     CLOSE(output_FID)
     STOP
   ELSE IF (iostat_value < 0) THEN
     EXIT
   ELSE
     !adjusting left, lowercase, trim data and read
     line_of_text = SimpleLine(line_of_text)
     line_of_text = StrLowCase(line_of_text)
     !assinging a value to each item
     IF(line_of_text == '<input filename>')THEN
       READ(job_control_FID,'(A)')long_input_filename 
       long_input_filename = TRIM(ADJUSTL(long_input_filename)) 
       base_length = index(long_input_filename,'.',back = .true.) - 1
       output_filename = long_input_filename(1:base_length)//'.output'
     ELSE IF (line_of_text == '<charge type>')THEN
       READ(job_control_FID,'(A)')charge_type
     ELSE IF (line_of_text == '<net charge>')THEN
       READ(job_control_FID,*)netcharge
     ELSE IF (line_of_text == '<density format>')THEN
       READ(job_control_FID,*)density_format
       density_format = StrLowCase(density_format)
       IF (density_format == 'e_per_angs3') THEN
         distance_scale = bohrperangstrom
         density_scale = 1.0_dp/bohrperangstrom**3
       ELSE IF (density_format == 'onetep') THEN
         distance_scale = 1.0_dp
         density_scale = 1.0_dp/bohrperangstrom**3
       ELSE 
         distance_scale = 1.0_dp
         density_scale = 1.0_dp
       END IF
     ELSE IF (line_of_text == '<compute bos>')THEN
       READ(job_control_FID,*)compute_BOs
     ELSE IF (line_of_text == '<print atomic densities>')THEN
       READ(job_control_FID,*)print_atomic_densities
     ELSE IF (line_of_text == '<number of core electrons>')THEN
       DO
         READ(job_control_FID,'(a)',IOSTAT=iostat_value) line_of_text
         IF ((adjustl(line_of_text(1:1)) == "<") .or. (iostat_value/=0) ) THEN    
           EXIT
         ELSE
           backspace(job_control_FID)        
           READ(job_control_FID,*) i,read_core_electrons
           num_core(i) =  read_core_electrons
         END IF
       END DO
     END IF
   END IF
 END DO
 CLOSE(job_control_FID)
  !Check if it is wfx or some other kind of file
 input_type = 'nan'
 IF (StrLowCase(long_input_filename((base_length + 2) : LEN(long_input_filename)))=='wfx') THEN
   input_type = 'wfx'
   periodicA = .false.
   periodicB = .false.
   periodicC = .false.
 ELSE IF (StrLowCase(long_input_filename((base_length + 2) : LEN(long_input_filename)))=='xsf') THEN
   input_type = 'xsf'
   periodicA = .true.
   periodicB = .true.
   periodicC = .true.
 END IF
 IF (input_type == 'nan') THEN
   OPEN(NEWUNIT=potcar_FID, file='POTCAR', status='OLD', iostat=iostat_value)
   IF (iostat_value == 0) THEN
      input_type = 'vsp'
      periodicA = .true.
      periodicB = .true.
      periodicC = .true.
      output_filename='VASP_DDEC_analysis.output'
      CLOSE(potcar_FID)
   END IF
 END IF
 IF (input_type == 'nan') THEN
   OPEN(NEWUNIT=total_density_FID,file='total_density.cube', status='OLD', iostat=iostat_value)
   IF (iostat_value == 0) THEN
     OPEN(NEWUNIT=valence_FID,file='valence_density.cube',IOSTAT=iostat_value1,STATUS='old')
     IF (iostat_value1 == 0) THEN
       input_type = 'vtc'
       periodicA = .false.
       periodicB = .false.
       periodicC = .false.
       output_filename='valence_total_cube_DDEC_analysis.output'
       CLOSE(valence_FID)
     ELSE
       input_type = 'cub'
       periodicA = .false.
       periodicB = .false.
       periodicC = .false.
       output_filename='total_cube_DDEC_analysis.output'
     END IF
     CLOSE(total_density_FID)
   ELSE
     OPEN(NEWUNIT=valence_FID,file='valence_density.cube',IOSTAT=iostat_value1,STATUS='old')
     IF (iostat_value1 == 0) THEN
       input_type = 'val'
       periodicA = .false.
       periodicB = .false.
       periodicC = .false.
       output_filename='valence_cube_DDEC_analysis.output'
     END IF
     CLOSE(valence_FID)
   END IF
 END IF
 OPEN(newunit=job_control_FID,file='job_control.txt',status='old',iostat=iostat_value)
 IF (iostat_value /= 0 ) THEN
   OPEN(newunit=output_FID,file='output.txt',status='replace')
   WRITE(output_FID,*)'job_control.txt could not be found. '
   WRITE(output_FID,*)'Program will terminate.'
   FLUSH(output_FID)
   CLOSE(output_FID)
   STOP
 END IF
 DO
   READ(job_control_FID,'(a)',IOSTAT=iostat_value) line_of_text
   IF (iostat_value > 0) THEN
     OPEN(newunit=output_FID,file='output.txt',status='replace')
     WRITE(output_FID,*) 'Could not read job_control.txt file. Program will terminate.'
     FLUSH(output_FID)
     CLOSE(output_FID)
     STOP
   ELSE IF (iostat_value < 0) THEN
     EXIT
   ELSE
     !adjusting left, lowercase, trim data and read
     line_of_text = SimpleLine(line_of_text)
     line_of_text = StrLowCase(line_of_text)
     !assinging a value to each item
     IF(line_of_text == '<atomic densities directory complete path>')THEN
       READ(job_control_FID,'(a)')atomic_densities_directory
       atomic_densities_directory=trim(atomic_densities_directory)
     ELSE IF(line_of_text == '<periodicity along a, b, and c vectors>') THEN
       READ(job_control_FID,*)periodicA
       READ(job_control_FID,*)periodicB
       READ(job_control_FID,*)periodicC
     END IF
   END IF
 END DO
 CLOSE(job_control_FID)
 OPEN(newunit = output_FID,file=trim(output_filename),status='replace')

 END SUBROUTINE read_job_control
 
 END MODULE module_read_job_control
   