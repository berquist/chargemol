 MODULE module_gen_dens_grids_from_gaussian_basis_set_coeff
 ! Copyright (c) 2014, 2015, 2016, 2017 by Thomas A. Manz and Nidia Gabaldon Limas. Rights reserved.

 USE module_precision
 USE module_global_parameters
 USE module_matrix_operations
 USE module_common_variable_declarations
 USE module_charge_center_positions_and_parallelpiped
 USE module_gaussian_functions
 USE module_add_missing_core_density
 USE module_initialize_atomic_densities
 USE module_align_periodic_atoms

 IMPLICIT NONE

 CONTAINS
 
 SUBROUTINE generate_density_grids_from_gaussian_basis_set_coefficients()
 !===================================================================================
 ! Generates valence, spin, and core densities from natural orbital coefficients with
 ! gaussian primitive basis set
 !===================================================================================
 
 REAL(kind=dp) :: temp_vector1(3),temp_vector2(3),temp_vector3(3),check_included_electrons,temp1,temp2,temp3,temp4,temp5,&
 X1,Y1,Z1,X2,Y2,Z2,alpha1,alpha2,max_coeff,test,partial_pair_weight_i,partial_pair_weight_j,full_1PDM,&
 radius,r12,temp_center_shift(3),XD,YD,ZD,xyz(3),temp_vector(3),temp,analytic_pair_overlap,overlap_sum,&
 abs_overlap_sum,X,Y,Z,temp_scalar,spin_contribution,distance,alpha,Ax,Ay,Az,temp_nabc(3),alpha_sum,basis_center_separation,&
 dsquared,Rsquared_cutoff,common_gaussian_factor,center_1_nuclear_charge,center_2_nuclear_charge,r12_squared,&
 temp_radius_squared,temp_included_primitive_data_row(40),temp6,temp7,temp8,max_pair_correction,dX1,dY1,dZ1,dX2,dY2,dZ2,&
 fa,fb,fc,amplitude,a_dot_a,a_dot_b,a_dot_c,b_dot_b,b_dot_c,c_dot_c,cos_a_b,cos_a_c,cos_b_c,sin_a_b,sin_a_c,sin_b_c,&
 radial_moment(6)
 INTEGER :: min_na,min_nb,min_nc,max_na,max_nb,max_nc, periodicsumlimitA,periodicsumlimitB,periodicsumlimitC,&
 nimages,k,analytic_pair_count,numeric_pair_count,total_count,repeata,repeatb,repeatc,Lx1,Ly1,Lz1,Lx2,Ly2,Lz2,&
 num_analytic_primitive_pairs,num_numeric_primitive_pairs,num_included_primitive_pairs,atom_index_for_i,&
 atom_index_for_j,renormalization_step,delta_na_pair,delta_nb_pair,delta_nc_pair,na,nb,nc,ka,kb,kc,basis_index_for_i,&
 basis_index_for_j,center_index_for_i,center_index_for_j,temp_index,Lx,Ly,Lz,k2a,k2b,k2c,k4a,k4b,k4c,k6a,k6b,k6c,k3a,k3b,k3c,&
 include_pair_block,block_count,center_1,center_2,pair_count,block_i,block_j,num_pair_blocks,pairs_in_block,max_delta_na,&
 max_delta_nb,max_delta_nc,max_grid_spacing,delta_na_block,delta_nb_block,delta_nc_block,Na_block_center,&
 Nb_block_center,Nc_block_center,gridspacing,k8a,k8b,k8c,k12a,k12b,k12c,k12a_lower,k12a_upper,k12b_lower,k12b_upper,k12c_lower,&
 k12c_upper,k8a_lower,k8a_upper,k8b_lower,k8b_upper,k8c_lower,k8c_upper,k6a_lower,k6a_upper,k6b_lower,k6b_upper,k6c_lower,&
 k6c_upper,k4a_lower,k4a_upper,k4b_lower,k4b_upper,k4c_lower,k4c_upper,k3a_lower,k3a_upper,k3b_lower,k3b_upper,k3c_lower,k3c_upper&
 ,k2a_lower,k2a_upper,k2b_lower,k2b_upper,k2c_lower,k2c_upper
 INTEGER :: clock_count,clock_count_rate,clock_count_max,clock_count2 
 INTEGER, ALLOCATABLE, DIMENSION (:) :: included_primitive_pairs
 REAL(kind=dp), ALLOCATABLE, DIMENSION(:,:) :: alpha_1PDM,beta_1PDM,included_primitive_data,pair_block_data
 REAL(kind=dp), ALLOCATABLE, DIMENSION (:) :: vector_length, pair_overlap,MO_normalization_check
 REAL(kind=dp), ALLOCATABLE, DIMENSION (:,:,:) :: valence_density_grid2,valence_density_grid3,valence_density_grid4,&
 valence_density_grid6,spin_density_grid2,spin_density_grid3,spin_density_grid4,spin_density_grid6,valence_density_grid8,&
 valence_density_grid12,spin_density_grid8,spin_density_grid12
 LOGICAL :: print_single_atom_radial_moments
  
 WRITE(output_FID,*)' '
 WRITE(output_FID,*)'******************** TIME CONTROL ********************'
 WRITE(output_FID,*)'Starting generate_density_grids_from_gaussian_basis_set_coefficients'
 FLUSH(output_FID)
 CALL system_clock(clock_count,clock_count_rate,clock_count_max)
 WRITE(output_FID,*)'******************************************************'
 WRITE(output_FID,*)' '
 FLUSH(output_FID)  
  
 !----------------------------------------------------------------------------------
 ! Step 1: Choose the Grid Points
 !----------------------------------------------------------------------------------
 ! Choose the pixel boundary used to generate the grid
 ALLOCATE(vector_length(num_periodic_directions))
 vector_length=0.0_dp
 IF (num_periodic_directions == 0) THEN
  boundary(1,:) = [preferred_grid_spacing,0.0_dp,0.0_dp]
  boundary(2,:) = [0.0_dp,preferred_grid_spacing,0.0_dp]
  boundary(3,:) = [0.0_dp,0.0_dp,preferred_grid_spacing]
 ELSE
   DO i=1,num_periodic_directions
     vector_length(i)= sqrt(periodic_vectors(i,1)**2 + periodic_vectors(i,2)**2 + periodic_vectors(i,3)**2)
     boundary(i,:) = periodic_vectors(i,:)/(24*ceiling(vector_length(i)/(24*preferred_grid_spacing)))
   END DO
   IF (num_periodic_directions == 1) THEN
     IF (abs(periodic_vectors(1,1)) > abs(periodic_vectors(1,3))) THEN
        temp_vector2=cross([0.0_dp, 0.0_dp, 1.0_dp],periodic_vectors)
       ELSE
           temp_vector2=cross(periodic_vectors,[1.0_dp, 0.0_dp, 0.0_dp])
     END IF
     boundary(2,:) = temp_vector2(:)*preferred_grid_spacing/norm2(temp_vector2)
   END IF       
   IF (num_periodic_directions < 3) THEN
     temp_vector1(:)=boundary(1,:)
     temp_vector2(:)=boundary(2,:)
     temp_vector3(:)=cross(temp_vector1,temp_vector2)
     boundary(3,:)= temp_vector3(:)*preferred_grid_spacing/norm2(temp_vector3)
   END IF
 END IF
 DO i=1,3
   WRITE(output_FID,'(3f10.4)')(boundary(i,j),j=1,3)
 END DO
 periodicA = .false.
 periodicB = .false.
 periodicC = .false.
 IF (num_periodic_directions > 0) THEN
   periodicA=.true.
 END IF
 IF (num_periodic_directions > 1) THEN
   periodicB=.true.
 END IF
 IF (num_periodic_directions > 2) THEN
   periodicC=.true.
 END IF
 temp_vector1(:)=boundary(1,:)
 temp_vector2(:)=boundary(2,:)
 temp_vector3(:)=boundary(3,:)
 pixelvolume=determinant(boundary)
 WRITE(output_FID,*)'Initial setup of atoms using temporary origin (0,0,0)'
 FLUSH(output_FID)
 origin = [0.0_dp, 0.0_dp, 0.0_dp]
 ALLOCATE(center_nabc(3,natoms))
 ALLOCATE(center_shift(3,natoms))
 center_nabc = 0.0_dp
 center_shift = 0.0_dp
 CALL charge_center_positions( )
 min_na=minval(center_nabc(1,1:natoms))
 min_nb=minval(center_nabc(2,1:natoms))
 min_nc=minval(center_nabc(3,1:natoms))
 max_na=maxval(center_nabc(1,1:natoms))
 max_nb=maxval(center_nabc(2,1:natoms))
 max_nc=maxval(center_nabc(3,1:natoms))
 CALL parallelpiped( )
 WRITE(output_FID,*)'Compute a better origin and number of grid points along each grid direction'
 FLUSH(output_FID)
 IF (.not. periodicA) THEN
   totnumA = max_na - min_na + 2*delta_na
   !Eleven buffer pixels ensure no wrap around errors
   totnumA = totnumA + 11
   !Make sure totnumA is divisible by twentyfour
   totnumA = totnumA + 23.99_dp - modulo(totnumA + 23.99_dp,24.0_dp)
   WRITE(output_FID,'(a,i7)')' totnumA= ',totnumA
   FLUSH(output_FID)
   origin = origin + (min_na - delta_na - 6)*temp_vector1
 ELSE
   totnumA = 24*ceiling(vector_length(1)/(24*preferred_grid_spacing))
   WRITE(output_FID,'(a,i7)')' totnumA',totnumA
   FLUSH(output_FID)
 END IF
 IF (.not. periodicB) THEN
   totnumB = max_nb - min_nb + 2*delta_nb
   !Eleven buffer pixels ensure no wrap around errors
   totnumB = totnumB + 11
   !Make sure totnumB is divisible by twentyfour
   totnumB = totnumB + 23.99_dp - modulo(totnumB+23.99_dp,24.0_dp)
   WRITE(output_FID,'(a,i7)')' totnumB= ',totnumB
   FLUSH(output_FID)
   origin = origin + (min_nb - delta_nb - 6)*temp_vector2
 ELSE
   totnumB = 24*ceiling(vector_length(2)/(24*preferred_grid_spacing))
   WRITE(output_FID,'(a,i7)')' totnumB= ',totnumB
   FLUSH(output_FID)  
 END IF
 IF (.not. periodicC) THEN
   totnumC = max_nc - min_nc + 2*delta_nc
   !Eleven buffer pixels ensure no wrap around errors
   totnumC = totnumC + 11
   !Make sure totnumC is divisible by twentyfour
   totnumC = totnumC + 23.99_dp - modulo(totnumC+23.99_dp,24.0_dp)
   WRITE(output_FID,'(a,i7)')' totnumC= ',totnumC
   FLUSH(output_FID)
   origin = origin + (min_nc - delta_nc - 6)*temp_vector3
 ELSE
   totnumC = 24*ceiling(vector_length(3)/(24*preferred_grid_spacing))
   WRITE(output_FID,'(a,i7)')' totnumC= ',totnumC
   FLUSH(output_FID)  
 END IF
 WRITE(output_FID,*)'The new origin is'
 FLUSH(output_FID)
 WRITE(output_FID,'(3f10.4)')origin
 FLUSH(output_FID)
 DO i=1,3
  vector1(i)=temp_vector1(i)*totnumA
 END DO
 DO i=1,3
  vector2(i)=temp_vector2(i)*totnumB
 END DO
 DO i=1,3
  vector3(i)=temp_vector3(i)*totnumC
 END DO
 WRITE(output_FID,*)'Compute atomic positions using the better origin'
 FLUSH(output_FID)
 CALL charge_center_positions( )
 CALL align_periodic_atoms( )
 !Compute the number of basis center images to include along each periodic direction
 a_dot_a = vector1(1)**2 + vector1(2)**2 + vector1(3)**2
 a_dot_b = vector1(1)*vector2(1) + vector1(2)*vector2(2) + vector1(3)*vector2(3)
 a_dot_c = vector1(1)*vector3(1) + vector1(2)*vector3(2) + vector1(3)*vector3(3)
 b_dot_b = vector2(1)**2 + vector2(2)**2 + vector2(3)**2
 b_dot_c = vector2(1)*vector3(1) + vector2(2)*vector3(2) + vector2(3)*vector3(3)
 c_dot_c = vector3(1)**2 + vector3(2)**2 + vector3(3)**2
 cos_a_b = (a_dot_b)/sqrt(a_dot_a*b_dot_b)
 cos_a_c = (a_dot_c)/sqrt(a_dot_a*c_dot_c)
 cos_b_c = (b_dot_c)/sqrt(b_dot_b*c_dot_c)
 sin_a_b = sqrt(1.0_dp - cos_a_b*cos_a_b)
 sin_a_c = sqrt(1.0_dp - cos_a_c*cos_a_c)
 sin_b_c = sqrt(1.0_dp - cos_b_c*cos_b_c)
 IF (num_periodic_directions > 0.5_dp) THEN
    periodicsumlimitA = ceiling(max((periodic_cutoff_length/sin_a_b),(periodic_cutoff_length/sin_a_c))/sqrt(a_dot_a)) + 1 
 ELSE
    periodicsumlimitA = 0
 END IF
 IF (num_periodic_directions > 1.5_dp) THEN
    periodicsumlimitB = ceiling(max((periodic_cutoff_length/sin_a_b),(periodic_cutoff_length/sin_b_c))/sqrt(b_dot_b)) + 1
 ELSE
     periodicsumlimitB = 0
 END IF
 IF (num_periodic_directions > 2.5_dp) THEN
    periodicsumlimitC = ceiling(max((periodic_cutoff_length/sin_a_c),(periodic_cutoff_length/sin_b_c))/sqrt(c_dot_c)) + 1
 ELSE
    periodicsumlimitC = 0
 END IF
 nimages=(2*periodicsumlimitA + 1)*(2*periodicsumlimitB + 1)*(2*periodicsumlimitC + 1)
 WRITE(output_FID,'(a,i7)')' nimages= ',nimages
 FLUSH(output_FID)
 
 WRITE(output_FID,*)' '
 WRITE(output_FID,*)'******************** TIME CONTROL ********************'
 CALL system_clock(clock_count2,clock_count_rate,clock_count_max)
 IF (clock_count_rate == 0) THEN
   WRITE(output_FID,*)'Processor clock not available, the time calculation wil not be performed.'
 ELSE
   seconds = modulo((clock_count2-clock_count),clock_count_max)/real(clock_count_rate)
   WRITE(output_FID,*)'Finished Step 1 in ',seconds, 'seconds.'
   FLUSH(output_FID)
 END IF
 WRITE(output_FID,*)'******************************************************'
 WRITE(output_FID,*)' '
 FLUSH(output_FID)  
 
 !----------------------------------------------------------------------------------
 ! Step 2: Compute the density matrices for alpha and beta electrons
 !----------------------------------------------------------------------------------
 ALLOCATE(alpha_1PDM(nprimitives,nprimitives))
 ALLOCATE(beta_1PDM(nprimitives,nprimitives))
 alpha_1PDM = 0.0_dp
 beta_1PDM = 0.0_dp
 DO j=1,nprimitives
   DO k = j,nprimitives
     DO i=1,norbitals
       alpha_1PDM(j,k) = alpha_1PDM(j,k) + alpha_beta_occupation(i,1)*orbital_coefficients(i,j)*orbital_coefficients(i,k)
       beta_1PDM(j,k) = beta_1PDM(j,k) + alpha_beta_occupation(i,2)*orbital_coefficients(i,j)*orbital_coefficients(i,k)
     END DO
   END DO
 END DO 
 WRITE(output_FID,*)'Setup of grid points was successful.'
 FLUSH(output_FID)
 
 WRITE(output_FID,*)' '
 WRITE(output_FID,*)'******************** TIME CONTROL ********************'
 CALL system_clock(clock_count2,clock_count_rate,clock_count_max)
 IF (clock_count_rate == 0) THEN
   WRITE(output_FID,*)'Processor clock not available, the time calculation wil not be performed.'
 ELSE
   seconds = modulo((clock_count2-clock_count),clock_count_max)/real(clock_count_rate)
   WRITE(output_FID,*)'Finished Step 2 in ',seconds, 'seconds'
   FLUSH(output_FID)
 END IF
 WRITE(output_FID,*)'******************************************************'
 WRITE(output_FID,*)' '
 FLUSH(output_FID)  

 !If a single atom and non-periodic, then compute and print the radial moments
 IF ((ncenters == 1) .AND. (num_periodic_directions == 0)) THEN
   print_single_atom_radial_moments = .true.
   WRITE(output_FID,*)' '
   WRITE(output_FID,*)'Since the system is a single isolated atom, the radial moments will be computed and printed.'
   WRITE(output_FID,*)'These radial moments include only the explicit electrons.'
   WRITE(output_FID,*)'If your calculation used an effective core potential, the core electrons replaced &
   &by that potential are not included in these computed radial moments.'
   WRITE(output_FID,*)' '
   DO block_i = 1,nsingleblocks
      alpha1 = primitive_exponents(single_block_ranges(block_i,1),2)
      DO block_j = block_i,nsingleblocks
         alpha2 = primitive_exponents(single_block_ranges(block_j,1),2)
         alpha_sum = alpha1 + alpha2
         DO i = single_block_ranges(block_i,1),single_block_ranges(block_i,2)
            DO j = single_block_ranges(block_j,1),single_block_ranges(block_j,2)
               IF (j < i) THEN
                 CYCLE
               END IF 
               Lx1=primitive_exponents(i,3) + primitive_exponents(j,3)
               Ly1=primitive_exponents(i,4) + primitive_exponents(j,4)
               Lz1=primitive_exponents(i,5) + primitive_exponents(j,5)
               IF (j > i) THEN
                  temp1 = 2.0_dp
               ELSE
                  temp1 = 1.0_dp
               END IF
               IF ((modulo(Lx1,2) == 1) .OR. (modulo(Ly1,2) == 1) .OR. (modulo(Lz1,2) == 1)) THEN
                  CYCLE
               END IF
               IF (Lx1+Ly1+Lz1 > 8) THEN
                  print_single_atom_radial_moments = .false.
               END IF
               temp2 = alpha_1PDM(i,j) + beta_1PDM(i,j)
               temp3 = spherical_multiplier(Lx1,Ly1,Lz1)
               temp4 = temp1*temp2*temp3
               k = (Lx1+Ly1+Lz1)/2 + 1
               radial_moment(1) = radial_moment(1) + temp4 * (factorial(k-1)/2.0_dp) /(alpha_sum**k)
               radial_moment(2) = radial_moment(2) + temp4 * gaussian_integration(alpha_sum,2*k)
               radial_moment(3) = radial_moment(3) + temp4 * (factorial(k)/2.0_dp) /(alpha_sum**(k+1))
               radial_moment(4) = radial_moment(4) + temp4 * gaussian_integration(alpha_sum,2*k+2)
               radial_moment(5) = radial_moment(5) + temp4 * (factorial(k+1)/2.0_dp) /(alpha_sum**(k+2))  
               radial_moment(6) = radial_moment(6) + temp4 * gaussian_integration(alpha_sum,2*k+4)
            END DO
         END DO
      END DO
   END DO
   radial_moment(1) = 4.0_dp*pi*radial_moment(1)
   radial_moment(2) = 2.0_dp*pi*radial_moment(2)
   radial_moment(3) = 4.0_dp*pi*radial_moment(3)  
   radial_moment(4) = 2.0_dp*pi*radial_moment(4)
   radial_moment(5) = 4.0_dp*pi*radial_moment(5)
   radial_moment(6) = 2.0_dp*pi*radial_moment(6)
   IF (print_single_atom_radial_moments) THEN
      DO k = 1,6
         WRITE(output_FID,*)'The ',k-2,' radial moment is ',radial_moment(k)
      END DO
   ELSE
      WRITE(output_FID,*)'High angular momentum basis function detected. Isolated atom radial moments printing will be skipped.'
   END IF   
 END IF
 
 !----------------------------------------------------------------------------------
 ! Step 3: Use cutoff condition to determine which primitive pairs (i,j) to include
 !----------------------------------------------------------------------------------
 ALLOCATE(included_primitive_pairs((nprimitives*(nprimitives+1)*nimages/2)))
 included_primitive_pairs = 0.0_dp
 total_count = 0
 include_pair_block = 0
 block_count = 0
 center_1 = 0
 center_2 = 0
 pair_count = 0
 analytic_pair_count=0
 numeric_pair_count=0
 alpha_sum = 0.0_dp
 DO repeata=-periodicsumlimitA,periodicsumlimitA
   DO repeatb=-periodicsumlimitB,periodicsumlimitB
     DO repeatc=-periodicsumlimitC,periodicsumlimitC
       DO block_i = 1,nsingleblocks
         alpha1 = primitive_exponents(single_block_ranges(block_i,1),2)
         center_1 = nint(primitive_exponents(single_block_ranges(block_i,1),1))
         X1=basis_set_centers(center_1,3)
         Y1=basis_set_centers(center_1,4)
         Z1=basis_set_centers(center_1,5)
         DO block_j = block_i,nsingleblocks
           alpha2 = primitive_exponents(single_block_ranges(block_j,1),2)
           center_2 = nint(primitive_exponents(single_block_ranges(block_j,1),1))
           X2=basis_set_centers(center_2,3) + repeata*vector1(1) + repeatb*vector2(1) + repeatc*vector3(1)
           Y2=basis_set_centers(center_2,4) + repeata*vector1(2) + repeatb*vector2(2) + repeatc*vector3(2)
           Z2=basis_set_centers(center_2,5) + repeata*vector1(3) + repeatb*vector2(3) + repeatc*vector3(3)
           include_pair_block = 0
           alpha_sum = alpha1 + alpha2
           DO i = single_block_ranges(block_i,1),single_block_ranges(block_i,2)
             ! Round powers to even number so overlap test is positive definite
             DO j = single_block_ranges(block_j,1),single_block_ranges(block_j,2)
               IF (j < i) THEN
                 CYCLE
               END IF 
               total_count = total_count + 1
               !Round powers to even number so overlap test is positive definite
               IF (center_1 == center_2) THEN
                 IF (alpha_sum > 1.0_dp) THEN
                   !Round sum down to even number
                   Lx1=nint((primitive_exponents(i,3) + primitive_exponents(j,3)) - modulo((nint(primitive_exponents(i,3) + &
                   primitive_exponents(j,3))),2))
                   Ly1=nint((primitive_exponents(i,4) + primitive_exponents(j,4)) - modulo((nint(primitive_exponents(i,4) + &
                   primitive_exponents(j,4))),2))
                   Lz1=nint((primitive_exponents(i,5) + primitive_exponents(j,5)) - modulo((nint(primitive_exponents(i,5) + &
                   primitive_exponents(j,5))),2))
                   Lx2 = 0
                   Ly2 = 0
                   Lz2 = 0
                 ELSE
                   !Round sum up to even number
                   Lx1=2*ceiling((primitive_exponents(i,3) + primitive_exponents(j,3))/2.0_dp)
                   Ly1=2*ceiling((primitive_exponents(i,4) + primitive_exponents(j,4))/2.0_dp)
                   Lz1=2*ceiling((primitive_exponents(i,5) + primitive_exponents(j,5))/2.0_dp)
                   Lx2 = 0
                   Ly2 = 0
                   Lz2 = 0
                 END IF                                
               ELSE    
                 IF (alpha_sum > 1.0_dp) THEN
                   !Round each down to even number
                   Lx1=nint(primitive_exponents(i,3) - modulo(nint(primitive_exponents(i,3)),2))
                   Ly1=nint(primitive_exponents(i,4) - modulo(nint(primitive_exponents(i,4)),2))
                   Lz1=nint(primitive_exponents(i,5) - modulo(nint(primitive_exponents(i,5)),2))
                   Lx2=nint(primitive_exponents(j,3) - modulo(nint(primitive_exponents(j,3)),2))
                   Ly2=nint(primitive_exponents(j,4) - modulo(nint(primitive_exponents(j,4)),2))
                   Lz2=nint(primitive_exponents(j,5) - modulo(nint(primitive_exponents(j,5)),2))
                 ELSE
                   !Round each up to even number
                   Lx1=2*ceiling(primitive_exponents(i,3)/2.0_dp)
                   Ly1=2*ceiling(primitive_exponents(i,4)/2.0_dp)
                   Lz1=2*ceiling(primitive_exponents(i,5)/2.0_dp)
                   Lx2=2*ceiling(primitive_exponents(j,3)/2.0_dp)
                   Ly2=2*ceiling(primitive_exponents(j,4)/2.0_dp)
                   Lz2=2*ceiling(primitive_exponents(j,5)/2.0_dp)
                 END IF
               END IF
               max_coeff = abs(alpha_1PDM(i,j)) + abs(beta_1PDM(i,j))
               test = max_coeff*gaussian_overlap_integral(alpha1,Lx1,Ly1,Lz1,X1,Y1,Z1,alpha2,Lx2,Ly2,Lz2,X2,Y2,Z2)
               IF (test < gaussian_overlap_tolerance) THEN
                 CYCLE 
               ELSE 
                 include_pair_block = 1
                 pair_count = pair_count + 1
                 IF ((alpha1*min(basis_set_centers(nint(primitive_exponents(i,1)),2),1.0_dp) + alpha2*min(basis_set_centers(nint&
                 (primitive_exponents(j,1)),2),1.0_dp)) > analytic_alpha_cutoff) THEN
                   analytic_pair_count = analytic_pair_count + 1
                   included_primitive_pairs(total_count)= 1.0_dp
                 ELSE
                   numeric_pair_count = numeric_pair_count + 1
                   included_primitive_pairs(total_count)= 2.0_dp
                 END IF       
               END IF
             END DO
           END DO
           IF (include_pair_block == 1) THEN
             block_count = block_count + 1
           END IF
         END DO
       END DO
     END DO
   END DO
 END DO
 num_pair_blocks = block_count
 num_included_primitive_pairs = pair_count
 num_analytic_primitive_pairs = analytic_pair_count
 num_numeric_primitive_pairs = numeric_pair_count
 WRITE(output_FID,'(a,i7)')' The number of primitive pair blocks are ',num_pair_blocks
 FLUSH(output_FID)
 WRITE(output_FID,'(a,i7)')' The total number of included primitives is: ',num_included_primitive_pairs
 FLUSH(output_FID)
 WRITE(output_FID,'(a,i7)')' The number of primitive products to treat analytically is: ',num_analytic_primitive_pairs
 FLUSH(output_FID)
 WRITE(output_FID,'(a,i7)')' The number of primitive products to treat numerically is: ',num_numeric_primitive_pairs
 WRITE(output_FID,*)' '
 FLUSH(output_FID)
 
 WRITE(output_FID,*)' '
 WRITE(output_FID,*)'******************** TIME CONTROL ********************'
 CALL system_clock(clock_count2,clock_count_rate,clock_count_max)
 IF (clock_count_rate == 0) THEN
   WRITE(output_FID,*)'Processor clock not available, the time calculation wil not be performed.'
 ELSE
   seconds = modulo((clock_count2-clock_count),clock_count_max)/real(clock_count_rate)
   WRITE(output_FID,*)'Finished Step 3 in ',seconds, 'seconds'
   FLUSH(output_FID)
 END IF
 WRITE(output_FID,*)'******************************************************'
 WRITE(output_FID,*)' '
 FLUSH(output_FID)
 
 
 !----------------------------------------------------------------------------------
 ! Step 4: Check whether the primitive pairs give the correct total number of electons
 !---------------------------------------------------------------------------------- 
 check_included_electrons = 0.0_dp
 ALLOCATE(pair_overlap(nprimitives*(nprimitives+1)*nimages/2))
 pair_overlap = 0.0_dp
 total_count = 0
 center_1 = 0
 center_2 = 0
 DO repeata=-periodicsumlimitA,periodicsumlimitA
   DO repeatb=-periodicsumlimitB,periodicsumlimitB
     DO repeatc=-periodicsumlimitC,periodicsumlimitC
       DO block_i = 1,nsingleblocks
         alpha1 = primitive_exponents(single_block_ranges(block_i,1),2)
         center_1 = nint(primitive_exponents(single_block_ranges(block_i,1),1))
         X1=basis_set_centers(center_1,3)
         Y1=basis_set_centers(center_1,4)
         Z1=basis_set_centers(center_1,5)
         DO block_j = block_i,nsingleblocks
           alpha2 = primitive_exponents(single_block_ranges(block_j,1),2)
           center_2 = nint(primitive_exponents(single_block_ranges(block_j,1),1))
           X2=basis_set_centers(center_2,3) + repeata*vector1(1) + repeatb*vector2(1) + repeatc*vector3(1)
           Y2=basis_set_centers(center_2,4) + repeata*vector1(2) + repeatb*vector2(2) + repeatc*vector3(2)
           Z2=basis_set_centers(center_2,5) + repeata*vector1(3) + repeatb*vector2(3) + repeatc*vector3(3)
           DO i = single_block_ranges(block_i,1),single_block_ranges(block_i,2)
             Lx1 = nint(primitive_exponents(i,3))
             Ly1 = nint(primitive_exponents(i,4))
             Lz1 = nint(primitive_exponents(i,5))
             DO j = single_block_ranges(block_j,1),single_block_ranges(block_j,2)
               IF (j < i) THEN
                 CYCLE
               END IF
               total_count = total_count + 1
               Lx2 = nint(primitive_exponents(j,3))
               Ly2 = nint(primitive_exponents(j,4))
               Lz2 = nint(primitive_exponents(j,5))
               pair_overlap(total_count) = gaussian_overlap_integral(alpha1,Lx1,Ly1,Lz1,X1,Y1,Z1,alpha2,Lx2,Ly2,Lz2,X2,Y2,Z2)
               IF (j > i) THEN
                 check_included_electrons = check_included_electrons + 2.0_dp*(alpha_1PDM(i,j) + beta_1PDM(i,j))&
                 *pair_overlap(total_count)
               ELSE
                 check_included_electrons = check_included_electrons + (alpha_1PDM(i,j) + beta_1PDM(i,j))*&
                 pair_overlap(total_count)
               END IF
             END DO
           END DO
         END DO
      END DO
     END DO
   END DO
 END DO
 total_count = 0
 IF (abs(check_included_electrons - included_electrons) > 0.001) THEN
   WRITE(output_FID,*) 'The one-particle density matrix is not properly normalized.'
   FLUSH(output_FID)
   WRITE(output_FID,*) check_included_electrons - included_electrons
   FLUSH(output_FID)
   WRITE(output_FID,*) 'Computing the self-overlap (should equal one) for each natural orbital.'
   FLUSH(output_FID)
   ALLOCATE(MO_normalization_check(norbitals))
   MO_normalization_check = 0.0_dp
   total_count = 0
   DO repeata=-periodicsumlimitA,periodicsumlimitA
     DO repeatb=-periodicsumlimitB,periodicsumlimitB
       DO repeatc=-periodicsumlimitC,periodicsumlimitC
         DO block_i = 1,nsingleblocks
           DO block_j = block_i,nsingleblocks
             DO i = single_block_ranges(block_i,1),single_block_ranges(block_i,2)
               DO j = single_block_ranges(block_j,1),single_block_ranges(block_j,2)
                 IF (j < i) THEN
                   CYCLE
                 END IF    
                 total_count = total_count + 1
                 IF (included_primitive_pairs(total_count)== 0) THEN
                   CYCLE
                 END IF 
                 DO k=1,norbitals
                   IF (j > i ) THEN
                     MO_normalization_check(k) = MO_normalization_check(k) + 2*orbital_coefficients(k,i)*orbital_coefficients&
                     (k,j)*pair_overlap(total_count)
                   ELSE
                     MO_normalization_check(k) = MO_normalization_check(k) + orbital_coefficients(k,i)*orbital_coefficients(k,j)&
                     *pair_overlap(total_count)
                   END IF                                
                 END DO
               END DO
             END DO
           END DO
         END DO
       END DO
     END DO
   END DO
   WRITE(output_FID,*)'MO_normalization_check= ',MO_normalization_check
   FLUSH(output_FID)
   WRITE(output_FID,*)'Program will terminate.'  
   FLUSH(output_FID)
   DEALLOCATE(MO_normalization_check)
   STOP
 ELSE
   WRITE(output_FID,*)'The one-particle density matrix is properly normalized.'
   FLUSH(output_FID)
 END IF
 
 WRITE(output_FID,*)' '
 WRITE(output_FID,*)'******************** TIME CONTROL ********************'
 CALL system_clock(clock_count2,clock_count_rate,clock_count_max)
 IF (clock_count_rate == 0) THEN
   WRITE(output_FID,*)'Processor clock not available, the time calculation wil not be performed.'
 ELSE
   seconds = modulo((clock_count2-clock_count),clock_count_max)/real(clock_count_rate)
   WRITE(output_FID,*)'Finished Step 4 in ',seconds, 'seconds'
   FLUSH(output_FID)
 END IF
 WRITE(output_FID,*)'******************************************************'
 WRITE(output_FID,*)' '
 FLUSH(output_FID) 
 
 !----------------------------------------------------------------------------------
 ! Step 5: Collect key data for each primitive pair
 !---------------------------------------------------------------------------------- 
 ALLOCATE(included_primitive_data(num_included_primitive_pairs,40))
 included_primitive_data = 0.0_dp
 pair_count = 0
 analytic_pair_count = 0
 total_count = 0
 temp_nabc = 0.0_dp
 temp_center_shift = 0.0_dp
 ALLOCATE(occupancy_correction(11,natoms))
 occupancy_correction = 0.0_dp
 pairs_in_block = 0
 center_1_nuclear_charge = 0.0_dp
 center_2_nuclear_charge = 0.0_dp
 r12_squared = 0.0_dp
 block_count = 0
 full_1PDM = 0.0_dp
 atom_index_for_i = 0 
 atom_index_for_j = 0 
 radius = 0.0_dp
 max_delta_na = 0
 max_delta_nb = 0
 max_delta_nc = 0
 ALLOCATE(pair_block_data(num_pair_blocks,2))
 pair_block_data = 0.0_dp
 max_grid_spacing = 0
 DO repeata=-periodicsumlimitA,periodicsumlimitA
   DO repeatb=-periodicsumlimitB,periodicsumlimitB
     DO repeatc=-periodicsumlimitC,periodicsumlimitC
       DO block_i = 1,nsingleblocks
         alpha1 = primitive_exponents(single_block_ranges(block_i,1),2)
         center_1 = nint(primitive_exponents(single_block_ranges(block_i,1),1))
         X1=basis_set_centers(center_1,3)
         Y1=basis_set_centers(center_1,4)
         Z1=basis_set_centers(center_1,5)
         DO block_j = block_i,nsingleblocks
           alpha2 = primitive_exponents(single_block_ranges(block_j,1),2)
           center_2 = nint(primitive_exponents(single_block_ranges(block_j,1),1))
           X2=basis_set_centers(center_2,3) + repeata*vector1(1) + repeatb*vector2(1) + repeatc*vector3(1)
           Y2=basis_set_centers(center_2,4) + repeata*vector1(2) + repeatb*vector2(2) + repeatc*vector3(2)
           Z2=basis_set_centers(center_2,5) + repeata*vector1(3) + repeatb*vector2(3) + repeatc*vector3(3)
           alpha_sum = alpha1 + alpha2
           pairs_in_block = 0
           max_delta_na = 0
           max_delta_nb = 0
           max_delta_nc = 0
           DO i = single_block_ranges(block_i,1),single_block_ranges(block_i,2)
             Lx1 = nint(primitive_exponents(i,3))
             Ly1 = nint(primitive_exponents(i,4))
             Lz1 = nint(primitive_exponents(i,5))
             DO j = single_block_ranges(block_j,1),single_block_ranges(block_j,2)
               IF (j < i) THEN
                 CYCLE
               END IF    
               total_count = total_count + 1
               IF (included_primitive_pairs(total_count) == 0) THEN
                 CYCLE
               END IF
               pair_count = pair_count + 1
               pairs_in_block = pairs_in_block + 1
               IF (pairs_in_block == 1) THEN
                 block_count = block_count + 1
                 pair_block_data(block_count,1) = pair_count
               END IF
               Lx2 = nint(primitive_exponents(j,3))
               Ly2 = nint(primitive_exponents(j,4))
               Lz2 = nint(primitive_exponents(j,5))
               IF (i == j) THEN
                 included_primitive_data(pair_count,3) = (alpha_1PDM(i,j) + beta_1PDM(i,j))
                 included_primitive_data(pair_count,4) = (alpha_1PDM(i,j) - beta_1PDM(i,j))
               ELSE
                 included_primitive_data(pair_count,3) = 2.0_dp*(alpha_1PDM(i,j) + beta_1PDM(i,j))
                 included_primitive_data(pair_count,4) = 2.0_dp*(alpha_1PDM(i,j) - beta_1PDM(i,j))
               END IF
               IF (included_primitive_pairs(total_count) == 1) THEN !Analytic pairs
                 analytic_pair_count = analytic_pair_count + 1
                 IF ((basis_set_centers(nint(primitive_exponents(i,1)),2) > 0.0_dp) .and. (basis_set_centers(nint&
                 (primitive_exponents(j,1)),2) > 0.0_dp)) THEN
                   partial_pair_weight_i = alpha1**10/(alpha1**10 + alpha2**10)
                   partial_pair_weight_j = alpha2**10/(alpha1**10 + alpha2**10)
                 ELSE IF ((basis_set_centers(nint(primitive_exponents(i,1)),2) > 0.0_dp)) THEN
                   partial_pair_weight_i = 1.0_dp
                   partial_pair_weight_j = 0.0_dp
                 ELSE IF ((basis_set_centers(nint(primitive_exponents(j,1)),2) > 0.0_dp)) THEN
                   partial_pair_weight_i = 0.0_dp
                   partial_pair_weight_j = 1.0_dp
                 ELSE
                   partial_pair_weight_i = 0.5_dp
                   partial_pair_weight_j = 0.5_dp
                 END IF
                 included_primitive_data(pair_count,37) = partial_pair_weight_i
                 included_primitive_data(pair_count,38) = partial_pair_weight_j
                 full_1PDM = included_primitive_data(pair_count,3)
                 atom_index_for_i = basis_set_centers(nint(primitive_exponents(i,1)),6)
                 atom_index_for_j = basis_set_centers(nint(primitive_exponents(j,1)),6)
                 IF (atom_index_for_i > 0.5_dp) THEN
                   !Number of valence electron correction
                   core_electrons(atom_index_for_i) = core_electrons(atom_index_for_i) + full_1PDM*pair_overlap(total_count)&
                   *partial_pair_weight_i
                   !Spin magnetic moment correction
                   occupancy_correction(2,atom_index_for_i) = occupancy_correction(2,atom_index_for_i) + included_primitive_data&
                   (pair_count,4)*pair_overlap(total_count)*partial_pair_weight_i
                   !Dipole corrections
                   occupancy_correction(3,atom_index_for_i) = occupancy_correction(3,atom_index_for_i) - &
                   full_1PDM*partial_pair_weight_i*gaussian_overlap_integral(alpha1,Lx1+1,Ly1,Lz1,X1,Y1,Z1,alpha2,Lx2,Ly2,Lz2,&
                   X2,Y2,Z2)
                   occupancy_correction(4,atom_index_for_i) = occupancy_correction(4,atom_index_for_i) - &
                   full_1PDM*partial_pair_weight_i*gaussian_overlap_integral(alpha1,Lx1,Ly1+1,Lz1,X1,Y1,Z1,alpha2,Lx2,Ly2,Lz2,&
                   X2,Y2,Z2)
                   occupancy_correction(5,atom_index_for_i) = occupancy_correction(5,atom_index_for_i) - &
                   full_1PDM*partial_pair_weight_i*gaussian_overlap_integral(alpha1,Lx1,Ly1,Lz1+1,X1,Y1,Z1,alpha2,Lx2,Ly2,Lz2,&
                   X2,Y2,Z2)
                   !Quadrupole corrections
                   occupancy_correction(6,atom_index_for_i) = occupancy_correction(6,atom_index_for_i) - &
                   full_1PDM*partial_pair_weight_i*gaussian_overlap_integral(alpha1,Lx1+2,Ly1,Lz1,X1,Y1,Z1,alpha2,Lx2,Ly2,Lz2,&
                   X2,Y2,Z2) ! x**2
                   occupancy_correction(7,atom_index_for_i) = occupancy_correction(7,atom_index_for_i) - &
                   full_1PDM*partial_pair_weight_i*gaussian_overlap_integral(alpha1,Lx1,Ly1+2,Lz1,X1,Y1,Z1,alpha2,Lx2,Ly2,Lz2,&
                   X2,Y2,Z2) ! y**2
                   occupancy_correction(8,atom_index_for_i) = occupancy_correction(8,atom_index_for_i) - &
                   full_1PDM*partial_pair_weight_i*gaussian_overlap_integral(alpha1,Lx1,Ly1,Lz1+2,X1,Y1,Z1,alpha2,Lx2,Ly2,Lz2,&
                   X2,Y2,Z2) ! z**2
                   occupancy_correction(9,atom_index_for_i) = occupancy_correction(9,atom_index_for_i) - &
                   full_1PDM*partial_pair_weight_i*gaussian_overlap_integral(alpha1,Lx1+1,Ly1+1,Lz1,X1,Y1,Z1,alpha2,Lx2,Ly2,Lz2,&
                   X2,Y2,Z2) ! xy
                   occupancy_correction(10,atom_index_for_i) = occupancy_correction(10,atom_index_for_i) - &
                   full_1PDM*partial_pair_weight_i*gaussian_overlap_integral(alpha1,Lx1+1,Ly1,Lz1+1,X1,Y1,Z1,alpha2,Lx2,Ly2,Lz2,&
                   X2,Y2,Z2) ! xz
                   occupancy_correction(11,atom_index_for_i) = occupancy_correction(11,atom_index_for_i) - &
                   full_1PDM*partial_pair_weight_i*gaussian_overlap_integral(alpha1,Lx1,Ly1+1,Lz1+1,X1,Y1,Z1,alpha2,Lx2,Ly2,Lz2,&
                   X2,Y2,Z2) ! yz
                 END IF
                 IF (atom_index_for_j > 0.5_dp) THEN
                   !Number of valence electron correction
                   core_electrons(atom_index_for_j) = core_electrons(atom_index_for_j) + full_1PDM*pair_overlap(total_count)&
                   *partial_pair_weight_j
                   !Spin magnetic moment correction
                   occupancy_correction(2,atom_index_for_j) = occupancy_correction(2,atom_index_for_j) + included_primitive_data&
                   (pair_count,4)*pair_overlap(total_count)*partial_pair_weight_j
                   !Dipole corrections
                   occupancy_correction(3,atom_index_for_j) = occupancy_correction(3,atom_index_for_j) - &
                   full_1PDM*partial_pair_weight_j*gaussian_overlap_integral(alpha1,Lx1,Ly1,Lz1,X1,Y1,Z1,alpha2,Lx2+1,Ly2,Lz2,X2&
                   ,Y2,Z2)
                   occupancy_correction(4,atom_index_for_j) = occupancy_correction(4,atom_index_for_j) - &
                   full_1PDM*partial_pair_weight_j*gaussian_overlap_integral(alpha1,Lx1,Ly1,Lz1,X1,Y1,Z1,alpha2,Lx2,Ly2+1,Lz2,X2&
                   ,Y2,Z2)
                   occupancy_correction(5,atom_index_for_j) = occupancy_correction(5,atom_index_for_j) - &
                   full_1PDM*partial_pair_weight_j*gaussian_overlap_integral(alpha1,Lx1,Ly1,Lz1,X1,Y1,Z1,alpha2,Lx2,Ly2,Lz2+1,X2&
                   ,Y2,Z2)
                   !Quadrupole corrections
                   occupancy_correction(6,atom_index_for_j) = occupancy_correction(6,atom_index_for_j) - &
                   full_1PDM*partial_pair_weight_j*gaussian_overlap_integral(alpha1,Lx1,Ly1,Lz1,X1,Y1,Z1,alpha2,Lx2+2,Ly2,Lz2,X2&
                   ,Y2,Z2) ! x**2
                   occupancy_correction(7,atom_index_for_j) = occupancy_correction(7,atom_index_for_j) - full_1PDM*&
                   partial_pair_weight_j*gaussian_overlap_integral(alpha1,Lx1,Ly1,Lz1,X1,Y1,Z1,alpha2,Lx2,Ly2+2,Lz2,X2,Y2,Z2) !y**2
                   occupancy_correction(8,atom_index_for_j) = occupancy_correction(8,atom_index_for_j) - full_1PDM*&
                   partial_pair_weight_j*gaussian_overlap_integral(alpha1,Lx1,Ly1,Lz1,X1,Y1,Z1,alpha2,Lx2,Ly2,Lz2+2,X2,Y2,Z2) ! z**2
                   occupancy_correction(9,atom_index_for_j) = occupancy_correction(9,atom_index_for_j) - full_1PDM*&
                   partial_pair_weight_j*gaussian_overlap_integral(alpha1,Lx1,Ly1,Lz1,X1,Y1,Z1,alpha2,Lx2+1,Ly2+1,Lz2,X2,Y2,Z2) ! xy
                   occupancy_correction(10,atom_index_for_j) = occupancy_correction(10,atom_index_for_j) - full_1PDM*&
                   partial_pair_weight_j*gaussian_overlap_integral(alpha1,Lx1,Ly1,Lz1,X1,Y1,Z1,alpha2,Lx2+1,Ly2,Lz2+1,X2,Y2,Z2) ! xz
                   occupancy_correction(11,atom_index_for_j) = occupancy_correction(11,atom_index_for_j) - full_1PDM*&
                   partial_pair_weight_j*gaussian_overlap_integral(alpha1,Lx1,Ly1,Lz1,X1,Y1,Z1,alpha2,Lx2,Ly2+1,Lz2+1,X2,Y2,Z2) ! yz 
                 END IF                               
               END IF
               included_primitive_data(pair_count,1) = i ! original basis function number
               included_primitive_data(pair_count,2) = j ! original basis function number
               included_primitive_data(pair_count,5) = pair_overlap(total_count)
               included_primitive_data(pair_count,6) = alpha1
               included_primitive_data(pair_count,7) = Lx1
               included_primitive_data(pair_count,8) = Ly1
               included_primitive_data(pair_count,9) = Lz1
               included_primitive_data(pair_count,10) = X1
               included_primitive_data(pair_count,11) = Y1
               included_primitive_data(pair_count,12) = Z1
               included_primitive_data(pair_count,13) = alpha2
               included_primitive_data(pair_count,14) = Lx2
               included_primitive_data(pair_count,15) = Ly2
               included_primitive_data(pair_count,16) = Lz2
               included_primitive_data(pair_count,17) = X2
               included_primitive_data(pair_count,18) = Y2
               included_primitive_data(pair_count,19) = Z2
               r12 = sqrt((X1-X2)**2 + (Y1-Y2)**2 + (Z1-Z2)**2)
               radius = sqrt(max((-log(gaussian_overlap_tolerance)/alpha_sum - alpha1*alpha2*(r12**2)/(alpha_sum**2)),0.25_dp))
               included_primitive_data(pair_count,20) = ceiling(radius*100*delta_na/(cutoff_radius*bohrperangstrom)) + 1.0_dp ! delta_na value for pair
               included_primitive_data(pair_count,21) = ceiling(radius*100*delta_nb/(cutoff_radius*bohrperangstrom)) + 1.0_dp ! delta_nb value for pair
               included_primitive_data(pair_count,22) = ceiling(radius*100*delta_nc/(cutoff_radius*bohrperangstrom)) + 1.0_dp ! delta_nc value for pair
               !Compute the center (na,nb,nc) position and the center shifts
               XD = (alpha1*X1 + alpha2*X2)/alpha_sum !Coordinates of pair product center
               YD = (alpha1*Y1 + alpha2*Y2)/alpha_sum
               ZD = (alpha1*Z1 + alpha2*Z2)/alpha_sum
               included_primitive_data(pair_count,23) = XD
               included_primitive_data(pair_count,24) = YD
               included_primitive_data(pair_count,25) = ZD        
               xyz=[XD,YD,ZD]
               temp_vector = matmul(inv_boundary,(xyz - origin))
               included_primitive_data(pair_count,26) = nint(temp_vector(1)) !Pair product center Na, Nb, Nc values
               included_primitive_data(pair_count,27) = nint(temp_vector(2))
               included_primitive_data(pair_count,28) = nint(temp_vector(3))
               temp_center_shift(:) = matmul(TRANSPOSE(boundary),(temp_vector(:)) - nint(temp_vector(:))) !Pair product xyz center shifts
               included_primitive_data(pair_count,29) = temp_center_shift(1) 
               included_primitive_data(pair_count,30) = temp_center_shift(2)
               included_primitive_data(pair_count,31) = temp_center_shift(3)
               !Pairwise renormalization factor goes in column 32. Set initial estimate to 1.
               included_primitive_data(pair_count,32) = 1.0_dp
               !Whether the pair is treated numerically or analytically
               IF (included_primitive_pairs(total_count) == 1) THEN
                 included_primitive_data(pair_count,33) = 1.0_dp
               END IF
               !The integration cutoff radius squared (in bohr) for the primitive pair
               included_primitive_data(pair_count,34) = radius*radius
               !The common factor for the gaussian
               included_primitive_data(pair_count,35) = exp(-alpha1*alpha2*(r12**2)/alpha_sum)
               IF ((included_primitive_data(pair_count,33) == 0.0_dp)) THEN
                 !Grid spacing criteria based on the alpha_sum
                 IF (( alpha_sum > coarser_grid_alpha_cutoff )) THEN
                   included_primitive_data(pair_count,36) = 1.0_dp
                 ELSE IF (alpha_sum > coarser_grid_alpha_cutoff*4.0_dp/9.0_dp) THEN
                   included_primitive_data(pair_count,36) = 2.0_dp
                 ELSE IF (alpha_sum > coarser_grid_alpha_cutoff/4.0_dp) THEN
                   included_primitive_data(pair_count,36) = 3.0_dp
                 ELSE IF (alpha_sum > coarser_grid_alpha_cutoff/9.0_dp) THEN
                   included_primitive_data(pair_count,36) = 4.0_dp
                 ELSE IF (alpha_sum > coarser_grid_alpha_cutoff/16.0_dp) THEN
                   included_primitive_data(pair_count,36) = 6.0_dp
                 ELSE IF (alpha_sum > coarser_grid_alpha_cutoff/36.0_dp) THEN 
                   included_primitive_data(pair_count,36) = 8.0_dp
                 ELSE 
                   included_primitive_data(pair_count,36) = 12.0_dp
                 END IF
                 max_grid_spacing = max(max_grid_spacing,nint(included_primitive_data(pair_count,36)))
               END IF
             END DO
           END DO
           IF(pairs_in_block > 0) THEN
             pair_block_data(block_count,2) = pair_block_data(block_count,1) + pairs_in_block - 1
           END IF
         END DO      
       END DO
     END DO
   END DO
 END DO
 WRITE(output_FID,*)'The transformation of gaussian primitives was succesful.'
 FLUSH(output_FID)
 
 WRITE(output_FID,*)' '
 WRITE(output_FID,*)'******************** TIME CONTROL ********************'
 CALL system_clock(clock_count2,clock_count_rate,clock_count_max)
 IF (clock_count_rate == 0) THEN
   WRITE(output_FID,*)'Processor clock not available, the time calculation wil not be performed.'
 ELSE
   seconds = modulo((clock_count2-clock_count),clock_count_max)/real(clock_count_rate)
   WRITE(output_FID,*)'Finished Step 5 in ',seconds, 'seconds'
   FLUSH(output_FID)
 END IF
 WRITE(output_FID,*)'******************************************************'
 WRITE(output_FID,*)' '
 FLUSH(output_FID) 
 
 !----------------------------------------------------------------------------------
 ! Step 6: Sort pairs within a block according to cutoff radius
 !---------------------------------------------------------------------------------- 
 temp_radius_squared = 0.0_dp
 temp_included_primitive_data_row = 0.0_dp
 DO block_count = 1,num_pair_blocks
   start_index = pair_block_data(block_count,1)
   stop_index = pair_block_data(block_count,2)
   DO i = start_index,(stop_index - 1)
     temp_radius_squared = included_primitive_data(i,34)
     temp_index = i
     DO j = i,stop_index
       IF (included_primitive_data(j,34) > temp_radius_squared) THEN
         temp_radius_squared = included_primitive_data(j,34)
         temp_index = j
       END IF
     END DO
     IF (temp_index /= i) THEN
       temp_included_primitive_data_row(:) = included_primitive_data(i,:)
       included_primitive_data(i,:) = included_primitive_data(temp_index,:)
       included_primitive_data(temp_index,:) = temp_included_primitive_data_row(:)
     END IF
   END DO
 END DO
 
 WRITE(output_FID,*)' '
 WRITE(output_FID,*)'******************** TIME CONTROL ********************'
 CALL system_clock(clock_count2,clock_count_rate,clock_count_max)
 IF (clock_count_rate == 0) THEN
   WRITE(output_FID,*)'Processor clock not available, the time calculation wil not be performed.'
 ELSE
   seconds = modulo((clock_count2-clock_count),clock_count_max)/real(clock_count_rate)
   WRITE(output_FID,*)'Finished Step 6 in ',seconds, 'seconds'
   FLUSH(output_FID)
 END IF
 WRITE(output_FID,*)'******************************************************'
 WRITE(output_FID,*)' '
 FLUSH(output_FID) 
 
 !----------------------------------------------------------------------------------
 ! Step 7:  Compute the valence_density and spin_density values at each grid point
 !The first iteration computes the pairwise normalization factor, which is
 !used to renormalize the pair contributions in the second pass, but only
 !if the adjustment would lead to a change greater than the gaussian overlap tolerance.
 !---------------------------------------------------------------------------------- 
 temp1 = 0.0_dp
 temp2 = 0.0_dp
 temp3 = 0.0_dp
 temp4 = 0.0_dp
 temp5 = 0.0_dp
 temp6 = 0.0_dp
 temp7 = 0.0_dp
 temp8 = 0.0_dp
 dsquared = 0.0_dp
 Rsquared_cutoff = 0.0_dp
 common_gaussian_factor = 1.0_dp
 delta_na_block = 0
 delta_nb_block = 0
 delta_nc_block = 0
 Na_block_center = 0
 Nb_block_center = 0
 Nc_block_center = 0
 gridspacing = 0
 max_pair_correction = 0.0_dp
 IF (spin_available) THEN
   ALLOCATE(valence_density(totnumA,totnumB,totnumC))
   valence_density = 0.0_dp
   ALLOCATE(spin_density(totnumA,totnumB,totnumC))
   spin_density = 0.0_dp
   IF (max_grid_spacing > 1) THEN
     ALLOCATE(valence_density_grid2((totnumA/2),(totnumB/2),(totnumC/2)))
     valence_density_grid2 = 0.0_dp
     ALLOCATE(spin_density_grid2((totnumA/2),(totnumB/2),(totnumC/2)))
     spin_density_grid2 = 0.0_dp
   END IF
   IF (max_grid_spacing > 2 ) THEN
     ALLOCATE(valence_density_grid3((totnumA/3),(totnumB/3),(totnumC/3)))
     valence_density_grid3 = 0.0_dp
     ALLOCATE(spin_density_grid3((totnumA/3),(totnumB/3),(totnumC/3)))
     spin_density_grid3 = 0.0_dp   
   END IF
   IF (max_grid_spacing > 3) THEN
     ALLOCATE(valence_density_grid4((totnumA/4),(totnumB/4),(totnumC/4)))
     valence_density_grid4 = 0.0_dp
     ALLOCATE(spin_density_grid4((totnumA/4),(totnumB/4),(totnumC/4)))
     spin_density_grid4 = 0.0_dp
   END IF
   IF (max_grid_spacing > 4) THEN
     ALLOCATE(valence_density_grid6((totnumA/6),(totnumB/6),(totnumC/6)))
     valence_density_grid6 = 0.0_dp
     ALLOCATE(spin_density_grid6((totnumA/6),(totnumB/6),(totnumC/6)))
     spin_density_grid6 = 0.0_dp   
   END IF
   IF (max_grid_spacing > 6) THEN
     ALLOCATE(valence_density_grid8((totnumA/8),(totnumB/8),(totnumC/8)))
     valence_density_grid8 = 0.0_dp
     ALLOCATE(spin_density_grid8((totnumA/8),(totnumB/8),(totnumC/8)))
     spin_density_grid8 = 0.0_dp
   END IF
   IF (max_grid_spacing > 8) THEN
     ALLOCATE(valence_density_grid12((totnumA/12),(totnumB/12),(totnumC/12)))
     valence_density_grid12 = 0.0_dp
     ALLOCATE(spin_density_grid12((totnumA/12),(totnumB/12),(totnumC/12)))
     spin_density_grid12 = 0.0_dp
   END IF
   ALLOCATE(core_density(totnumA,totnumB,totnumC))
   core_density=0.0_dp
   temp = 0.0_dp
   alpha_sum = 0.0_dp
   basis_center_separation = 0.0_dp
   dsquared = 0.0_dp
   Rsquared_cutoff = 0.0_dp
   common_gaussian_factor = 1.0_dp
   
 WRITE(output_FID,*)' '
 WRITE(output_FID,*)'******************** TIME CONTROL ********************'
 CALL system_clock(clock_count2,clock_count_rate,clock_count_max)
 IF (clock_count_rate == 0) THEN
   WRITE(output_FID,*)'Processor clock not available, the time calculation wil not be performed.'
 ELSE
   seconds = modulo((clock_count2-clock_count),clock_count_max)/real(clock_count_rate)
   WRITE(output_FID,*)'Before loop ',seconds, 'seconds'
   FLUSH(output_FID)
 END IF
 WRITE(output_FID,*)'******************************************************'
 WRITE(output_FID,*)' '
 FLUSH(output_FID)    
   
!$omp parallel do default(none) &
!$omp private(block_count,start_index,stop_index,alpha1,X1,Y1,Z1,alpha2,alpha_sum,X2,Y2,Z2,delta_na_block,delta_nb_block,&
!$omp delta_nc_block,Apoints,Bpoints,Cpoints,XD,YD,ZD,Na_block_center,Nb_block_center,Nc_block_center,temp1,temp2,temp3,temp4,&
!$omp temp5,temp6,temp7,temp8,Rsquared_cutoff,common_gaussian_factor,gridspacing,partial_pair_weight_i,partial_pair_weight_j,&
!$omp basis_index_for_i,basis_index_for_j,center_index_for_i,center_index_for_j,atom_index_for_i,atom_index_for_j,min_na,max_na,&
!$omp na,ka,min_nb,max_nb,nb,min_nc,max_nc,nc,renormalization_step,X,Y,Z,dsquared,dX1,dY1,dZ1,dX2,dY2,dZ2,pair_count,k2a,k2b,k2c,&
!$omp k3a,k3b,k3c,k4a,k4b,k4c,k6a,k6b,k6c,k8a,k8b,k8c,k12a,k12b,k12c,spin_contribution,max_coeff,analytic_pair_overlap,&
!$omp overlap_sum,abs_overlap_sum,kb,kc)&
!$omp shared(num_pair_blocks,pair_block_data,included_primitive_data,periodicA,totnumA,periodicB,totnumB,periodicC,totnumC,&
!$omp boundary,pixelvolume,valence_density,spin_density,valence_density_grid2,spin_density_grid2,valence_density_grid3,&
!$omp spin_density_grid3,valence_density_grid4,spin_density_grid4,valence_density_grid6,spin_density_grid6,valence_density_grid8,&
!$omp spin_density_grid8,valence_density_grid12,spin_density_grid12,core_density,occupancy_correction,primitive_exponents,&
!$omp basis_set_centers) &
!$omp reduction(max:max_pair_correction) &
!$omp schedule(dynamic,1)
   DO block_count = 1,num_pair_blocks
     !Data common to all primitives in a block
     start_index = pair_block_data(block_count,1)
     stop_index = pair_block_data(block_count,2)
     alpha1=included_primitive_data(start_index,6)
     X1 = included_primitive_data(start_index,10)
     Y1 = included_primitive_data(start_index,11)
     Z1 = included_primitive_data(start_index,12)
     alpha2=included_primitive_data(start_index,13)
     alpha_sum = alpha1 + alpha2
     X2 = included_primitive_data(start_index,17)
     Y2 = included_primitive_data(start_index,18)
     Z2 = included_primitive_data(start_index,19)
     delta_na_block = nint(included_primitive_data(start_index,20))
     delta_nb_block = nint(included_primitive_data(start_index,21))
     delta_nc_block = nint(included_primitive_data(start_index,22))
     ALLOCATE(Apoints(2*delta_na_block + 1))
     ALLOCATE(Bpoints(2*delta_nb_block + 1))
     ALLOCATE(Cpoints(2*delta_nc_block + 1))   
     XD = included_primitive_data(start_index,23)
     YD = included_primitive_data(start_index,24)
     ZD = included_primitive_data(start_index,25)
     Na_block_center = nint(included_primitive_data(start_index,26))
     Nb_block_center = nint(included_primitive_data(start_index,27))
     Nc_block_center = nint(included_primitive_data(start_index,28))
     temp1 = included_primitive_data(start_index,23) - included_primitive_data(start_index,29)
     temp2 = included_primitive_data(start_index,24) - included_primitive_data(start_index,30)
     temp3 = included_primitive_data(start_index,25) - included_primitive_data(start_index,31)
     Rsquared_cutoff = included_primitive_data(start_index,34)
     common_gaussian_factor = included_primitive_data(start_index,35)
     gridspacing = nint(included_primitive_data(start_index,36))
     partial_pair_weight_i = included_primitive_data(start_index,37)
     partial_pair_weight_j = included_primitive_data(start_index,38)
     basis_index_for_i = nint(included_primitive_data(start_index,1))
     basis_index_for_j = nint(included_primitive_data(start_index,2))
     center_index_for_i = nint(primitive_exponents(basis_index_for_i,1))
     center_index_for_j = nint(primitive_exponents(basis_index_for_j,1))
     atom_index_for_i = nint(basis_set_centers(center_index_for_i,6))
     atom_index_for_j = nint(basis_set_centers(center_index_for_j,6))
     min_na = delta_na_block
     max_na = -delta_na_block
     Apoints = 0
     DO na = -delta_na_block,delta_na_block
       IF (periodicA) THEN
         ka = modulo((na + Na_block_center),totnumA) + 1
       ELSE
         ka = na + Na_block_center + 1
       END IF
       Apoints(delta_na_block + na + 1) = ka
       IF ((ka < 1) .or. (ka > totnumA)) THEN
         CYCLE
       ELSE IF ((gridspacing > 1) .and. (modulo(ka,gridspacing) /= (gridspacing - 1))) THEN
         CYCLE
       ELSE
         min_na = min(min_na,na)
         max_na = max(max_na,na)
       END IF  
     END DO
     min_nb = delta_nb_block
     max_nb = -delta_nb_block
     Bpoints = 0
     DO nb = -delta_nb_block,delta_nb_block
       IF (periodicB) THEN
         kb = modulo((nb + Nb_block_center),totnumB) + 1
       ELSE
         kb = nb + Nb_block_center + 1
       END IF
       Bpoints(delta_nb_block + nb + 1) = kb
       IF ((kb < 1) .or. (kb > totnumB)) THEN
         CYCLE
       ELSE IF ((gridspacing > 1) .and. (modulo(kb,gridspacing) /= (gridspacing - 1))) THEN
         CYCLE
       ELSE
         min_nb = min(min_nb,nb)
         max_nb = max(max_nb,nb)
       END IF     
     END DO
     min_nc = delta_nc_block
     max_nc = -delta_nc_block
     Cpoints = 0
     DO nc = -delta_nc_block,delta_nc_block
       IF (periodicC) THEN
         kc = modulo((nc + Nc_block_center),totnumC) + 1
       ELSE
         kc = nc + Nc_block_center + 1
       END IF
       Cpoints(delta_nc_block + nc + 1) = kc
       IF ((kc < 1) .or. (kc > totnumC)) THEN
         CYCLE
       ELSE IF ((gridspacing > 1) .and. (modulo(kc,gridspacing) /= (gridspacing - 1))) THEN
         CYCLE
       ELSE
         min_nc = min(min_nc,nc)
         max_nc = max(max_nc,nc)
       END IF 
     END DO
     IF ((min_na > max_na) .or. (min_nb > max_nb) .or. (min_nc > max_nc)) THEN
       !No valid point found
       CYCLE
     END IF
     DO renormalization_step = 1,2
       IF ((renormalization_step == 2) .and. (gridspacing == 0)) THEN
         CYCLE
       END IF  
       IF (gridspacing == 1) THEN
         DO nc = min_nc,max_nc
           kc = Cpoints(delta_nc_block + nc + 1)
           DO nb = min_nb,max_nb
             kb = Bpoints(delta_nb_block + nb + 1)
             DO na = min_na,max_na
               ka = Apoints(delta_na_block + na + 1)
               X = na*boundary(1,1) + nb*boundary(2,1) + nc*boundary(3,1) + temp1
               Y = na*boundary(1,2) + nb*boundary(2,2) + nc*boundary(3,2) + temp2
               Z = na*boundary(1,3) + nb*boundary(2,3) + nc*boundary(3,3) + temp3
               dsquared = (X - XD)**2 + (Y - YD)**2 + (Z - ZD)**2
               IF (dsquared > Rsquared_cutoff) THEN
                 CYCLE
               END IF   
               dX1 = X-X1
               dY1 = Y-Y1
               dZ1 = Z-Z1
               dX2 = X-X2
               dY2 = Y-Y2
               dZ2 = Z-Z2
               temp4 = exp(-alpha_sum*dsquared)*common_gaussian_factor
               temp5 = temp4*pixelvolume
               temp7 = 0.0_dp
               temp8 = 0.0_dp
               DO pair_count = start_index,stop_index 
                 IF ((dsquared > included_primitive_data(pair_count,34))) THEN
                   EXIT
                 ELSE IF (included_primitive_data(pair_count,32) == 0.0_dp) THEN
                   CYCLE 
                 END IF
                 temp6 = dX1**nint(included_primitive_data(pair_count,7))*dY1**nint(included_primitive_data(pair_count,8))*dZ1**&
                 nint(included_primitive_data(pair_count,9))*dX2**nint(included_primitive_data(pair_count,14))*dY2**nint&
                 (included_primitive_data(pair_count,15))*dZ2**nint(included_primitive_data(pair_count,16))
                 IF (renormalization_step == 1) THEN
                   included_primitive_data(pair_count,39) = included_primitive_data(pair_count,39) + temp5*temp6
                   included_primitive_data(pair_count,40) = included_primitive_data(pair_count,40) + temp5*abs(temp6)
                 ELSE    
                   temp6 = abs(temp6)
                 END IF
                 temp7 = temp7 + temp6*included_primitive_data(pair_count,32)*included_primitive_data(pair_count,3)
                 temp8 = temp8 + temp6*included_primitive_data(pair_count,32)*included_primitive_data(pair_count,4)
               END DO 
!$omp atomic
               valence_density(ka,kb,kc) = valence_density(ka,kb,kc) + temp7*temp4 
!$omp atomic
               spin_density(ka,kb,kc) = spin_density(ka,kb,kc) + temp8*temp4
             END DO
           END DO
         END DO      
       ELSE IF (gridspacing == 2) THEN
         !Code for computing the grid2 densities
         DO nc = min_nc,max_nc,2
           kc = Cpoints(delta_nc_block + nc + 1)
           DO nb = min_nb,max_nb,2
             kb = Bpoints(delta_nb_block + nb + 1)
             DO na = min_na,max_na,2
               ka = Apoints(delta_na_block + na + 1)
               X = na*boundary(1,1) + nb*boundary(2,1) + nc*boundary(3,1) + temp1
               Y = na*boundary(1,2) + nb*boundary(2,2) + nc*boundary(3,2) + temp2
               Z = na*boundary(1,3) + nb*boundary(2,3) + nc*boundary(3,3) + temp3
               dsquared = (X - XD)**2 + (Y - YD)**2 + (Z - ZD)**2
               IF (dsquared > Rsquared_cutoff) THEN
                 CYCLE
               END IF
               dX1 = X-X1
               dY1 = Y-Y1
               dZ1 = Z-Z1
               dX2 = X-X2
               dY2 = Y-Y2
               dZ2 = Z-Z2
               temp4 = exp(-alpha_sum*dsquared)*common_gaussian_factor
               temp5 = temp4*pixelvolume*8.0_dp
               temp7 = 0.0_dp
               temp8 = 0.0_dp
               DO pair_count = start_index,stop_index  
                 IF (dsquared > included_primitive_data(pair_count,34)) THEN
                   EXIT
                 ELSE IF (included_primitive_data(pair_count,32) == 0.0_dp) THEN
                   CYCLE 
                 END IF
                 temp6 = dX1**nint(included_primitive_data(pair_count,7))*dY1**nint(included_primitive_data(pair_count,8))*dZ1**&
                 nint(included_primitive_data(pair_count,9))*dX2**nint(included_primitive_data(pair_count,14))*dY2**nint&
                 (included_primitive_data(pair_count,15))*dZ2**nint(included_primitive_data(pair_count,16))
                 IF (renormalization_step == 1) THEN
                   included_primitive_data(pair_count,39) = included_primitive_data(pair_count,39) + temp5*temp6
                   included_primitive_data(pair_count,40) = included_primitive_data(pair_count,40) + temp5*abs(temp6)
                 ELSE    
                   temp6 = abs(temp6)
                 END IF
                 temp7 = temp7 + temp6*included_primitive_data(pair_count,32)*included_primitive_data(pair_count,3)
                 temp8 = temp8 + temp6*included_primitive_data(pair_count,32)*included_primitive_data(pair_count,4)
               END DO            
               k2a = nint((ka+1)/2.0_dp)
               k2b = nint((kb+1)/2.0_dp)
               k2c = nint((kc+1)/2.0_dp)
!$omp atomic
               valence_density_grid2(k2a,k2b,k2c) = valence_density_grid2(k2a,k2b,k2c) + temp7*temp4
!$omp atomic
               spin_density_grid2(k2a,k2b,k2c) = spin_density_grid2(k2a,k2b,k2c) + temp8*temp4
             END DO
           END DO
         END DO
       ELSE IF (gridspacing == 3) THEN
         !Code for computing the grid3 densities
         DO nc = min_nc,max_nc,3
           kc = Cpoints(delta_nc_block + nc + 1)
           DO nb = min_nb,max_nb,3
             kb = Bpoints(delta_nb_block + nb + 1)
             DO na = min_na,max_na,3
               ka = Apoints(delta_na_block + na + 1)
               X = na*boundary(1,1) + nb*boundary(2,1) + nc*boundary(3,1) + temp1
               Y = na*boundary(1,2) + nb*boundary(2,2) + nc*boundary(3,2) + temp2
               Z = na*boundary(1,3) + nb*boundary(2,3) + nc*boundary(3,3) + temp3
               dsquared = (X - XD)**2 + (Y - YD)**2 + (Z - ZD)**2
               IF (dsquared > Rsquared_cutoff) THEN
                 CYCLE
               END IF
               dX1 = X-X1
               dY1 = Y-Y1
               dZ1 = Z-Z1
               dX2 = X-X2
               dY2 = Y-Y2
               dZ2 = Z-Z2
               temp4 = exp(-alpha_sum*dsquared)*common_gaussian_factor
               temp5 = temp4*pixelvolume*27.0_dp
               temp7 = 0.0_dp
               temp8 = 0.0_dp
               DO pair_count = start_index,stop_index  
                 IF (dsquared > included_primitive_data(pair_count,34)) THEN
                   EXIT
                 ELSE IF (included_primitive_data(pair_count,32) == 0.0_dp) THEN
                   CYCLE 
                 END IF
                 temp6 = dX1**nint(included_primitive_data(pair_count,7))*dY1**nint(included_primitive_data(pair_count,8))*dZ1**&
                 nint(included_primitive_data(pair_count,9))*dX2**nint(included_primitive_data(pair_count,14))*dY2**nint&
                 (included_primitive_data(pair_count,15))*dZ2**nint(included_primitive_data(pair_count,16))
                 IF (renormalization_step == 1) THEN
                   included_primitive_data(pair_count,39) = included_primitive_data(pair_count,39) + temp5*temp6
                   included_primitive_data(pair_count,40) = included_primitive_data(pair_count,40) + temp5*abs(temp6)
                 ELSE    
                   temp6 = abs(temp6)
                 END IF
                 temp7 = temp7 + temp6*included_primitive_data(pair_count,32)*included_primitive_data(pair_count,3)
                 temp8 = temp8 + temp6*included_primitive_data(pair_count,32)*included_primitive_data(pair_count,4)
               END DO             
               k3a = nint((ka+1)/3.0_dp)
               k3b = nint((kb+1)/3.0_dp)
               k3c = nint((kc+1)/3.0_dp)
!$omp atomic
               valence_density_grid3(k3a,k3b,k3c) = valence_density_grid3(k3a,k3b,k3c) + temp7*temp4
!$omp atomic
               spin_density_grid3(k3a,k3b,k3c) = spin_density_grid3(k3a,k3b,k3c) + temp8*temp4
             END DO
           END DO
         END DO   
       ELSE IF (gridspacing == 4) THEN
         !code for computing the grid4 densities
         DO nc = min_nc,max_nc,4
           kc = Cpoints(delta_nc_block + nc + 1)
           DO nb = min_nb,max_nb,4
             kb = Bpoints(delta_nb_block + nb + 1)
             DO na = min_na,max_na,4
               ka = Apoints(delta_na_block + na + 1)
               X = na*boundary(1,1) + nb*boundary(2,1) + nc*boundary(3,1) + temp1
               Y = na*boundary(1,2) + nb*boundary(2,2) + nc*boundary(3,2) + temp2
               Z = na*boundary(1,3) + nb*boundary(2,3) + nc*boundary(3,3) + temp3
               dsquared = (X - XD)**2 + (Y - YD)**2 + (Z - ZD)**2
               IF (dsquared > Rsquared_cutoff) THEN
                 CYCLE
               END IF
               dX1 = X-X1
               dY1 = Y-Y1
               dZ1 = Z-Z1
               dX2 = X-X2
               dY2 = Y-Y2
               dZ2 = Z-Z2
               temp4 = exp(-alpha_sum*dsquared)*common_gaussian_factor
               temp5 = temp4*pixelvolume*64.0_dp
               temp7 = 0.0_dp
               temp8 = 0.0_dp
               DO pair_count = start_index,stop_index  
                 IF ((dsquared > included_primitive_data(pair_count,34))) THEN
                   EXIT
                 ELSE IF (included_primitive_data(pair_count,32) == 0.0_dp) THEN
                   CYCLE 
                 END IF
                 temp6 = dX1**nint(included_primitive_data(pair_count,7))*dY1**nint(included_primitive_data(pair_count,8))*dZ1**&
                 nint(included_primitive_data(pair_count,9))*dX2**nint(included_primitive_data(pair_count,14))*dY2**nint&
                 (included_primitive_data(pair_count,15))*dZ2**nint(included_primitive_data(pair_count,16))
                 IF (renormalization_step == 1) THEN
                   included_primitive_data(pair_count,39) = included_primitive_data(pair_count,39) + temp5*temp6
                   included_primitive_data(pair_count,40) = included_primitive_data(pair_count,40) + temp5*abs(temp6)
                 ELSE    
                   temp6 = abs(temp6)
                 END IF
                 temp7 = temp7 + temp6*included_primitive_data(pair_count,32)*included_primitive_data(pair_count,3)
                 temp8 = temp8 + temp6*included_primitive_data(pair_count,32)*included_primitive_data(pair_count,4)
               END DO            
               k4a = nint((ka+1)/4.0_dp)
               k4b = nint((kb+1)/4.0_dp)
               k4c = nint((kc+1)/4.0_dp)
!$omp atomic
               valence_density_grid4(k4a,k4b,k4c) = valence_density_grid4(k4a,k4b,k4c) + temp7*temp4
!$omp atomic
               spin_density_grid4(k4a,k4b,k4c) = spin_density_grid4(k4a,k4b,k4c) + temp8*temp4
             END DO
           END DO
         END DO  
       ELSE IF (gridspacing == 6) THEN
         !Code for computing the grid6 densities
         DO nc = min_nc,max_nc,6
           kc = Cpoints(delta_nc_block + nc + 1)
           DO nb = min_nb,max_nb,6
             kb = Bpoints(delta_nb_block + nb + 1)
             DO na = min_na,max_na,6
               ka = Apoints(delta_na_block + na + 1)
               X = na*boundary(1,1) + nb*boundary(2,1) + nc*boundary(3,1) + temp1
               Y = na*boundary(1,2) + nb*boundary(2,2) + nc*boundary(3,2) + temp2
               Z = na*boundary(1,3) + nb*boundary(2,3) + nc*boundary(3,3) + temp3
               dsquared = (X - XD)**2 + (Y - YD)**2 + (Z - ZD)**2
               IF (dsquared > Rsquared_cutoff) THEN
                 CYCLE
               END IF
               dX1 = X-X1
               dY1 = Y-Y1
               dZ1 = Z-Z1
               dX2 = X-X2
               dY2 = Y-Y2
               dZ2 = Z-Z2
               temp4 = exp(-alpha_sum*dsquared)*common_gaussian_factor
               temp5 = temp4*pixelvolume*216.0_dp
               temp7 = 0.0_dp
               temp8 = 0.0_dp
               DO pair_count = start_index,stop_index  
                 IF (dsquared > included_primitive_data(pair_count,34)) THEN
                   EXIT
                 ELSE IF (included_primitive_data(pair_count,32) == 0.0_dp) THEN
                    CYCLE 
                 END IF
                 temp6 = dX1**nint(included_primitive_data(pair_count,7))*dY1**nint(included_primitive_data(pair_count,8))*dZ1**&
                 nint(included_primitive_data(pair_count,9))*dX2**nint(included_primitive_data(pair_count,14))*dY2**nint&
                 (included_primitive_data(pair_count,15))*dZ2**nint(included_primitive_data(pair_count,16))
                 IF (renormalization_step == 1) THEN
                   included_primitive_data(pair_count,39) = included_primitive_data(pair_count,39) + temp5*temp6
                   included_primitive_data(pair_count,40) = included_primitive_data(pair_count,40) + temp5*abs(temp6)
                 ELSE    
                   temp6 = abs(temp6)
                 END IF
                 temp7 = temp7 + temp6*included_primitive_data(pair_count,32)*included_primitive_data(pair_count,3)
                 temp8 = temp8 + temp6*included_primitive_data(pair_count,32)*included_primitive_data(pair_count,4)
               END DO             
               k6a = nint((ka+1)/6.0_dp)
               k6b = nint((kb+1)/6.0_dp)
               k6c = nint((kc+1)/6.0_dp)
!$omp atomic
               valence_density_grid6(k6a,k6b,k6c) = valence_density_grid6(k6a,k6b,k6c) + temp7*temp4
!$omp atomic
               spin_density_grid6(k6a,k6b,k6c) = spin_density_grid6(k6a,k6b,k6c) + temp8*temp4
             END DO
           END DO
         END DO
       ELSE IF (gridspacing == 8) THEN
         !Code for computing the grid8 densities
         DO nc = min_nc,max_nc,8
           kc = Cpoints(delta_nc_block + nc + 1)
           DO nb = min_nb,max_nb,8
             kb = Bpoints(delta_nb_block + nb + 1)
             DO na = min_na,max_na,8
               ka = Apoints(delta_na_block + na + 1)
               X = na*boundary(1,1) + nb*boundary(2,1) + nc*boundary(3,1) + temp1
               Y = na*boundary(1,2) + nb*boundary(2,2) + nc*boundary(3,2) + temp2
               Z = na*boundary(1,3) + nb*boundary(2,3) + nc*boundary(3,3) + temp3
               dsquared = (X - XD)**2 + (Y - YD)**2 + (Z - ZD)**2
               IF (dsquared > Rsquared_cutoff) THEN
                 CYCLE
               END IF
               dX1 = X-X1
               dY1 = Y-Y1
               dZ1 = Z-Z1
               dX2 = X-X2
               dY2 = Y-Y2
               dZ2 = Z-Z2
               temp4 = exp(-alpha_sum*dsquared)*common_gaussian_factor
               temp5 = temp4*pixelvolume*512.0_dp
               temp7 = 0.0_dp
               temp8 = 0.0_dp
               DO pair_count = start_index,stop_index  
                 IF (dsquared > included_primitive_data(pair_count,34)) THEN
                   EXIT
                 ELSE IF (included_primitive_data(pair_count,32) == 0.0_dp) THEN
                   CYCLE 
                 END IF
                 temp6 = dX1**nint(included_primitive_data(pair_count,7))*dY1**nint(included_primitive_data(pair_count,8))*dZ1&
                 **nint(included_primitive_data(pair_count,9))*dX2**nint(included_primitive_data(pair_count,14))*dY2**nint&
                 (included_primitive_data(pair_count,15))*dZ2**nint(included_primitive_data(pair_count,16))
                 IF (renormalization_step == 1) THEN
                   included_primitive_data(pair_count,39) = included_primitive_data(pair_count,39) + temp5*temp6
                   included_primitive_data(pair_count,40) = included_primitive_data(pair_count,40) + temp5*abs(temp6)
                 ELSE    
                   temp6 = abs(temp6)
                 END IF
                 temp7 = temp7 + temp6*included_primitive_data(pair_count,32)*included_primitive_data(pair_count,3)
                 temp8 = temp8 + temp6*included_primitive_data(pair_count,32)*included_primitive_data(pair_count,4)
               END DO             
               k8a = nint((ka+1)/8.0_dp)
               k8b = nint((kb+1)/8.0_dp)
               k8c = nint((kc+1)/8.0_dp)
!$omp atomic
               valence_density_grid8(k8a,k8b,k8c) = valence_density_grid8(k8a,k8b,k8c) + temp7*temp4
!$omp atomic
               spin_density_grid8(k8a,k8b,k8c) = spin_density_grid8(k8a,k8b,k8c) + temp8*temp4
             END DO
           END DO
         END DO
       ELSE IF (gridspacing == 12)  THEN
         !Code for computing the grid12 densities
         DO nc = min_nc,max_nc,12
           kc = Cpoints(delta_nc_block + nc + 1)
           DO nb = min_nb,max_nb,12
             kb = Bpoints(delta_nb_block + nb + 1)
             DO na = min_na,max_na,12
               ka = Apoints(delta_na_block + na + 1)
               X = na*boundary(1,1) + nb*boundary(2,1) + nc*boundary(3,1) + temp1
               Y = na*boundary(1,2) + nb*boundary(2,2) + nc*boundary(3,2) + temp2
               Z = na*boundary(1,3) + nb*boundary(2,3) + nc*boundary(3,3) + temp3
               dsquared = (X - XD)**2 + (Y - YD)**2 + (Z - ZD)**2
               IF (dsquared > Rsquared_cutoff) THEN
                 CYCLE
               END IF
               dX1 = X-X1
               dY1 = Y-Y1
               dZ1 = Z-Z1
               dX2 = X-X2
               dY2 = Y-Y2
               dZ2 = Z-Z2
               temp4 = exp(-alpha_sum*dsquared)*common_gaussian_factor
               temp5 = temp4*pixelvolume*1728.0_dp
               temp7 = 0.0_dp
               temp8 = 0.0_dp
               DO pair_count = start_index,stop_index  
                 IF (dsquared > included_primitive_data(pair_count,34)) THEN
                   EXIT
                 ELSE IF (included_primitive_data(pair_count,32) == 0.0_dp) THEN
                   CYCLE 
                 END IF
                 temp6 = dX1**nint(included_primitive_data(pair_count,7))*dY1**nint(included_primitive_data(pair_count,8))*dZ1&
                 **nint(included_primitive_data(pair_count,9))*dX2**nint(included_primitive_data(pair_count,14))*dY2**nint&
                 (included_primitive_data(pair_count,15))*dZ2**nint(included_primitive_data(pair_count,16))
                 IF (renormalization_step == 1) THEN
                   included_primitive_data(pair_count,39) = included_primitive_data(pair_count,39) + temp5*temp6
                   included_primitive_data(pair_count,40) = included_primitive_data(pair_count,40) + temp5*abs(temp6)
                 ELSE    
                   temp6 = abs(temp6)
                 END IF
                 temp7 = temp7 + temp6*included_primitive_data(pair_count,32)*included_primitive_data(pair_count,3)
                 temp8 = temp8 + temp6*included_primitive_data(pair_count,32)*included_primitive_data(pair_count,4)
               END DO             
               k12a = nint((ka+1)/12.0_dp)
               k12b = nint((kb+1)/12.0_dp)
               k12c = nint((kc+1)/12.0_dp) 
!$omp atomic
               valence_density_grid12(k12a,k12b,k12c) = valence_density_grid12(k12a,k12b,k12c) + temp7*temp4
!$omp atomic
               spin_density_grid12(k12a,k12b,k12c) = spin_density_grid12(k12a,k12b,k12c) + temp8*temp4
             END DO
           END DO
         END DO    
       ELSE IF (gridspacing == 0) THEN
         DO nc = min_nc,max_nc
           kc = Cpoints(delta_nc_block + nc + 1)
           DO nb = min_nb,max_nb
             kb = Bpoints(delta_nb_block + nb + 1)
             DO na = min_na,max_na
               ka = Apoints(delta_na_block + na + 1)
               X = na*boundary(1,1) + nb*boundary(2,1) + nc*boundary(3,1) + temp1
               Y = na*boundary(1,2) + nb*boundary(2,2) + nc*boundary(3,2) + temp2
               Z = na*boundary(1,3) + nb*boundary(2,3) + nc*boundary(3,3) + temp3
               dsquared = (X - XD)**2 + (Y - YD)**2 + (Z - ZD)**2
               IF (dsquared > Rsquared_cutoff) THEN
                 CYCLE
               END IF
               dX1 = X-X1
               dY1 = Y-Y1
               dZ1 = Z-Z1
               dX2 = X-X2
               dY2 = Y-Y2
               dZ2 = Z-Z2
               temp4 = exp(-alpha_sum*dsquared)*common_gaussian_factor
               temp7 = 0.0_dp
               temp8 = 0.0_dp
               DO pair_count = start_index,stop_index  
                 IF (dsquared > included_primitive_data(pair_count,34)) THEN
                   EXIT
                 END IF
                 temp6 = dX1**nint(included_primitive_data(pair_count,7))*dY1**nint(included_primitive_data(pair_count,8))*dZ1**&
                 nint(included_primitive_data(pair_count,9))*dX2**nint(included_primitive_data(pair_count,14))*dY2**nint&
                 (included_primitive_data(pair_count,15))*dZ2**nint(included_primitive_data(pair_count,16))
                 temp7 = temp7 + temp6*included_primitive_data(pair_count,32)*included_primitive_data(pair_count,3)
                 temp8 = temp8 + temp6*included_primitive_data(pair_count,32)*included_primitive_data(pair_count,4)
               END DO   
!$omp atomic
               core_density(ka,kb,kc) = core_density(ka,kb,kc) + temp7*temp4
               spin_contribution = temp8*temp4
!$omp atomic
               spin_density(ka,kb,kc) = spin_density(ka,kb,kc) + spin_contribution
               IF (atom_index_for_i > 0) THEN
                 !Spin magnetic moment correction
!$omp atomic
                 occupancy_correction(2,atom_index_for_i) = occupancy_correction(2,atom_index_for_i) - spin_contribution*&
                 pixelvolume*partial_pair_weight_i
               END IF
               IF (atom_index_for_j > 0) THEN
                 !Spin magnetic moment correction
!$omp atomic
                 occupancy_correction(2,atom_index_for_j) = occupancy_correction(2,atom_index_for_j) - spin_contribution&
                 *pixelvolume*partial_pair_weight_j
               END IF                                       
             END DO                           
           END DO    
         END DO
       END IF
       !Update the renormalization factor
       IF ((gridspacing > 0) .and. (renormalization_step == 1)) THEN
         DO pair_count = start_index,stop_index
           !The following line makes sure the correction is not more than 5 percent of the overlap sum
           max_coeff=abs((included_primitive_data(pair_count,3) + included_primitive_data(pair_count,4))/2.0_dp) + abs((&
           included_primitive_data(pair_count,3) - included_primitive_data(pair_count,4))/2.0_dp)
           analytic_pair_overlap = included_primitive_data(pair_count,5)
           overlap_sum = included_primitive_data(pair_count,39)
           abs_overlap_sum = included_primitive_data(pair_count,40)
           IF (max_coeff*abs_overlap_sum > gaussian_overlap_tolerance) THEN
             included_primitive_data(pair_count,32) = max(min(((analytic_pair_overlap - overlap_sum)/abs_overlap_sum),0.05_dp),&
             -0.05_dp)
             max_pair_correction = max(abs((analytic_pair_overlap - overlap_sum)/abs_overlap_sum), max_pair_correction)
           ELSE
             included_primitive_data(pair_count,32) = 0.0_dp
           END IF
         END DO     
       END IF
     END DO
     DEALLOCATE(Apoints)
     DEALLOCATE(Bpoints)
     DEALLOCATE(Cpoints)
   END DO 
!$omp end parallel do
   !Interpolate the coarser grids back onto the regular grid

 WRITE(output_FID,*)' '
 WRITE(output_FID,*)'******************** TIME CONTROL ********************'
 CALL system_clock(clock_count2,clock_count_rate,clock_count_max)
 IF (clock_count_rate == 0) THEN
   WRITE(output_FID,*)'Processor clock not available, the time calculation wil not be performed.'
 ELSE
   seconds = modulo((clock_count2-clock_count),clock_count_max)/real(clock_count_rate)
   WRITE(output_FID,*)'Before interpolation ',seconds, 'seconds'
   FLUSH(output_FID)
 END IF
 WRITE(output_FID,*)'******************************************************'
 WRITE(output_FID,*)' '
 FLUSH(output_FID)   
   IF (max_grid_spacing > 8) THEN
!$omp parallel do default(none) &
!$omp private(k6c,k6b,k6a,k12a_lower,k12a_upper,k12b_lower,k12b_upper,k12c_lower,k12c_upper) &
!$omp shared(totnumC,totnumB,totnumA,valence_density_grid12,spin_density_grid12,valence_density_grid6,spin_density_grid6) &
!$omp schedule(static,1)
     DO k6c = 1,totnumC/6
       DO k6b = 1,totnumB/6
         DO k6a = 1,totnumA/6
           k12a_lower = modulo((nint((k6a - 0.5_dp)/2.0_dp)-1),totnumA/12)+1
           k12a_upper = modulo((nint((k6a + 0.5_dp)/2.0_dp)-1),totnumA/12)+1
           k12b_lower = modulo((nint((k6b - 0.5_dp)/2.0_dp)-1),totnumB/12)+1
           k12b_upper = modulo((nint((k6b + 0.5_dp)/2.0_dp)-1),totnumB/12)+1
           k12c_lower = modulo((nint((k6c - 0.5_dp)/2.0_dp)-1),totnumC/12)+1
           k12c_upper = modulo((nint((k6c + 0.5_dp)/2.0_dp)-1),totnumC/12)+1
           valence_density_grid6(k6a,k6b,k6c) = valence_density_grid6(k6a,k6b,k6c)+(1.0_dp/8.0_dp)*(valence_density_grid12(&
           k12a_lower,k12b_lower,k12c_lower)+valence_density_grid12(k12a_upper,k12b_lower,k12c_lower)+valence_density_grid12(&
           k12a_lower,k12b_upper,k12c_lower)+valence_density_grid12(k12a_lower,k12b_lower,k12c_upper)+valence_density_grid12(&
           k12a_upper,k12b_upper,k12c_lower)+valence_density_grid12(k12a_upper,k12b_lower,k12c_upper)+valence_density_grid12(&
           k12a_lower,k12b_upper,k12c_upper)+valence_density_grid12(k12a_upper,k12b_upper,k12c_upper))
           spin_density_grid6(k6a,k6b,k6c) = spin_density_grid6(k6a,k6b,k6c)+(1.0_dp/8.0_dp)*(spin_density_grid12(k12a_lower,&
           k12b_lower,k12c_lower)+spin_density_grid12(k12a_upper,k12b_lower,k12c_lower)+spin_density_grid12(k12a_lower,k12b_upper&
           ,k12c_lower)+spin_density_grid12(k12a_lower,k12b_lower,k12c_upper)+spin_density_grid12(k12a_upper,k12b_upper,k12c_lower&
           )+spin_density_grid12(k12a_upper,k12b_lower,k12c_upper)+spin_density_grid12(k12a_lower,k12b_upper,k12c_upper)+&
           spin_density_grid12(k12a_upper,k12b_upper,k12c_upper))
         END DO
       END DO
     END DO
!$omp end parallel do
     DEALLOCATE(valence_density_grid12)
     DEALLOCATE(spin_density_grid12)
   END IF
   IF (max_grid_spacing > 6) THEN
!$omp parallel do default(none) &
!$omp private(k4c,k4b,k4a,k8a_lower,k8a_upper,k8b_lower,k8b_upper,k8c_lower,k8c_upper) &
!$omp shared(totnumC,totnumB,totnumA,valence_density_grid8,spin_density_grid8,valence_density_grid4,spin_density_grid4) &
!$omp schedule(static,1)
     DO k4c = 1,totnumC/4
       DO k4b = 1,totnumB/4
         DO k4a = 1,totnumA/4
           k8a_lower = modulo((nint((k4a - 0.5_dp)/2.0_dp)-1),totnumA/8)+1
           k8a_upper = modulo((nint((k4a + 0.5_dp)/2.0_dp)-1),totnumA/8)+1
           k8b_lower = modulo((nint((k4b - 0.5_dp)/2.0_dp)-1),totnumB/8)+1
           k8b_upper = modulo((nint((k4b + 0.5_dp)/2.0_dp)-1),totnumB/8)+1
           k8c_lower = modulo((nint((k4c - 0.5_dp)/2.0_dp)-1),totnumC/8)+1
           k8c_upper = modulo((nint((k4c + 0.5_dp)/2.0_dp)-1),totnumC/8)+1
           valence_density_grid4(k4a,k4b,k4c) = valence_density_grid4(k4a,k4b,k4c)+(1.0_dp/8.0_dp)*(valence_density_grid8(&
           k8a_lower,k8b_lower,k8c_lower)+valence_density_grid8(k8a_upper,k8b_lower,k8c_lower)+valence_density_grid8(k8a_lower,&
           k8b_upper,k8c_lower)+valence_density_grid8(k8a_lower,k8b_lower,k8c_upper)+valence_density_grid8(k8a_upper,k8b_upper,&
           k8c_lower)+valence_density_grid8(k8a_upper,k8b_lower,k8c_upper)+valence_density_grid8(k8a_lower,k8b_upper,k8c_upper)+&
           valence_density_grid8(k8a_upper,k8b_upper,k8c_upper))
           spin_density_grid4(k4a,k4b,k4c) = spin_density_grid4(k4a,k4b,k4c)+(1.0_dp/8.0_dp)*(spin_density_grid8(k8a_lower,&
           k8b_lower,k8c_lower)+spin_density_grid8(k8a_upper,k8b_lower,k8c_lower)+spin_density_grid8(k8a_lower,k8b_upper,&
           k8c_lower)+spin_density_grid8(k8a_lower,k8b_lower,k8c_upper)+spin_density_grid8(k8a_upper,k8b_upper,k8c_lower)+&
           spin_density_grid8(k8a_upper,k8b_lower,k8c_upper)+spin_density_grid8(k8a_lower,k8b_upper,k8c_upper)+&
           spin_density_grid8(k8a_upper,k8b_upper,k8c_upper))
         END DO
       END DO
     END DO
!$omp end parallel do
     DEALLOCATE(valence_density_grid8)
     DEALLOCATE(spin_density_grid8)
   END IF
   IF (max_grid_spacing > 4) THEN
!$omp parallel do default(none) &
!$omp private(k3c,k3b,k3a,k6a_lower,k6a_upper,k6b_lower,k6b_upper,k6c_lower,k6c_upper) &
!$omp shared(totnumC,totnumB,totnumA,valence_density_grid6,spin_density_grid6,valence_density_grid3,spin_density_grid3) &
!$omp schedule(static,1)
     DO k3c = 1,totnumC/3
       DO k3b = 1,totnumB/3
         DO k3a = 1,totnumA/3
           k6a_lower = modulo((nint((k3a - 0.5_dp)/2.0_dp)-1),totnumA/6)+1
           k6a_upper = modulo((nint((k3a + 0.5_dp)/2.0_dp)-1),totnumA/6)+1
           k6b_lower = modulo((nint((k3b - 0.5_dp)/2.0_dp)-1),totnumB/6)+1
           k6b_upper = modulo((nint((k3b + 0.5_dp)/2.0_dp)-1),totnumB/6)+1
           k6c_lower = modulo((nint((k3c - 0.5_dp)/2.0_dp)-1),totnumC/6)+1
           k6c_upper = modulo((nint((k3c + 0.5_dp)/2.0_dp)-1),totnumC/6)+1
           valence_density_grid3(k3a,k3b,k3c) = valence_density_grid3(k3a,k3b,k3c)+(1.0_dp/8.0_dp)*(valence_density_grid6(&
           k6a_lower,k6b_lower,k6c_lower)+valence_density_grid6(k6a_upper,k6b_lower,k6c_lower)+valence_density_grid6(&
           k6a_lower,k6b_upper,k6c_lower)+valence_density_grid6(k6a_lower,k6b_lower,k6c_upper)+valence_density_grid6(k6a_upper,&
           k6b_upper,k6c_lower)+valence_density_grid6(k6a_upper,k6b_lower,k6c_upper)+valence_density_grid6(k6a_lower,k6b_upper,&
           k6c_upper)+valence_density_grid6(k6a_upper,k6b_upper,k6c_upper))
           spin_density_grid3(k3a,k3b,k3c) = spin_density_grid3(k3a,k3b,k3c)+(1.0_dp/8.0_dp)*(spin_density_grid6(k6a_lower,&
           k6b_lower,k6c_lower)+spin_density_grid6(k6a_upper,k6b_lower,k6c_lower)+spin_density_grid6(k6a_lower,k6b_upper,&
           k6c_lower)+spin_density_grid6(k6a_lower,k6b_lower,k6c_upper)+spin_density_grid6(k6a_upper,k6b_upper,k6c_lower)+&
           spin_density_grid6(k6a_upper,k6b_lower,k6c_upper)+spin_density_grid6(k6a_lower,k6b_upper,k6c_upper)+&
           spin_density_grid6(k6a_upper,k6b_upper,k6c_upper))
         END DO
       END DO
     END DO
!$omp end parallel do
     DEALLOCATE(valence_density_grid6)
     DEALLOCATE(spin_density_grid6)
   END IF
   IF (max_grid_spacing > 3) THEN
!$omp parallel do default(none) &
!$omp private(k2c,k2b,k2a,k4a_lower,k4a_upper,k4b_lower,k4b_upper,k4c_lower,k4c_upper) &
!$omp shared(totnumC,totnumB,totnumA,valence_density_grid4,spin_density_grid4,valence_density_grid2,spin_density_grid2) &
!$omp schedule(static,1)
     DO k2c = 1,totnumC/2
       DO k2b = 1,totnumB/2
         DO k2a = 1,totnumA/2
           k4a_lower = modulo((nint((k2a - 0.5_dp)/2.0_dp)-1),totnumA/4)+1
           k4a_upper = modulo((nint((k2a + 0.5_dp)/2.0_dp)-1),totnumA/4)+1
           k4b_lower = modulo((nint((k2b - 0.5_dp)/2.0_dp)-1),totnumB/4)+1
           k4b_upper = modulo((nint((k2b + 0.5_dp)/2.0_dp)-1),totnumB/4)+1
           k4c_lower = modulo((nint((k2c - 0.5_dp)/2.0_dp)-1),totnumC/4)+1
           k4c_upper = modulo((nint((k2c + 0.5_dp)/2.0_dp)-1),totnumC/4)+1
           valence_density_grid2(k2a,k2b,k2c) = valence_density_grid2(k2a,k2b,k2c)+(1.0_dp/8.0_dp)*(valence_density_grid4(&
           k4a_lower,k4b_lower,k4c_lower)+valence_density_grid4(k4a_upper,k4b_lower,k4c_lower)+valence_density_grid4(k4a_lower,&
           k4b_upper,k4c_lower)+valence_density_grid4(k4a_lower,k4b_lower,k4c_upper)+valence_density_grid4(k4a_upper,k4b_upper,&
           k4c_lower)+valence_density_grid4(k4a_upper,k4b_lower,k4c_upper)+valence_density_grid4(k4a_lower,k4b_upper,k4c_upper)+&
           valence_density_grid4(k4a_upper,k4b_upper,k4c_upper))
           spin_density_grid2(k2a,k2b,k2c) = spin_density_grid2(k2a,k2b,k2c)+(1.0_dp/8.0_dp)*(spin_density_grid4(k4a_lower,&
           k4b_lower,k4c_lower)+spin_density_grid4(k4a_upper,k4b_lower,k4c_lower)+spin_density_grid4(k4a_lower,k4b_upper,&
           k4c_lower)+spin_density_grid4(k4a_lower,k4b_lower,k4c_upper)+spin_density_grid4(k4a_upper,k4b_upper,k4c_lower)+&
           spin_density_grid4(k4a_upper,k4b_lower,k4c_upper)+spin_density_grid4(k4a_lower,k4b_upper,k4c_upper)+&
           spin_density_grid4(k4a_upper,k4b_upper,k4c_upper))
         END DO
       END DO
     END DO
!$omp end parallel do
     DEALLOCATE(valence_density_grid4)
     DEALLOCATE(spin_density_grid4)
   END IF
   IF (max_grid_spacing > 2) THEN
!$omp parallel do default(none) &
!$omp private(kc,kb,ka,k3a_lower,k3a_upper,k3b_lower,k3b_upper,k3c_lower,k3c_upper,fa,fb,fc) &
!$omp shared(totnumC,totnumB,totnumA,valence_density_grid3,spin_density_grid3,valence_density,spin_density) &
!$omp schedule(static,1)
     DO kc = 1,totnumC
       DO kb = 1,totnumB
         DO ka = 1,totnumA
           k3a_lower = modulo((nint(ka/3.0_dp)-1),totnumA/3)+1
           k3a_upper = modulo((nint((ka+2.0_dp)/3.0_dp)-1),totnumA/3)+1
           k3b_lower = modulo((nint(kb/3.0_dp)-1),totnumB/3)+1
           k3b_upper = modulo((nint((kb+2.0_dp)/3.0_dp)-1),totnumB/3)+1
           k3c_lower = modulo((nint(kc/3.0_dp)-1),totnumC/3)+1
           k3c_upper = modulo((nint((kc+2.0_dp)/3.0_dp)-1),totnumC/3)+1
           fa = modulo((ka+1),3)/3.0_dp
           fb = modulo((kb+1),3)/3.0_dp
           fc = modulo((kc+1),3)/3.0_dp
           valence_density(ka,kb,kc)=valence_density(ka,kb,kc)+((1.0_dp-fa)*(1.0_dp-fb)*(1.0_dp-fc))*valence_density_grid3(&
           k3a_lower,k3b_lower,k3c_lower)+(fa*(1.0_dp-fb)*(1.0_dp-fc))*valence_density_grid3(k3a_upper,k3b_lower,k3c_lower)+((&
           1.0_dp-fa)*fb*(1.0_dp-fc))*valence_density_grid3(k3a_lower,k3b_upper,k3c_lower)+((1.0_dp-fa)*(1.0_dp-fb)*fc)*&
           valence_density_grid3(k3a_lower,k3b_lower,k3c_upper)+(fa*fb*(1.0_dp-fc))*valence_density_grid3(k3a_upper,k3b_upper,&
           k3c_lower)+(fa*(1-fb)*fc)*valence_density_grid3(k3a_upper,k3b_lower,k3c_upper)+((1.0_dp-fa)*fb*fc)*&
           valence_density_grid3(k3a_lower,k3b_upper,k3c_upper)+(fa*fb*fc)*valence_density_grid3(k3a_upper,k3b_upper,k3c_upper)
           spin_density(ka,kb,kc)=spin_density(ka,kb,kc)+((1.0_dp-fa)*(1.0_dp-fb)*(1.0_dp-fc))*spin_density_grid3(k3a_lower,&
           k3b_lower,k3c_lower)+(fa*(1.0_dp-fb)*(1.0_dp-fc))*spin_density_grid3(k3a_upper,k3b_lower,k3c_lower)+((1.0_dp-fa)*fb*(&
           1.0_dp-fc))*spin_density_grid3(k3a_lower,k3b_upper,k3c_lower)+((1.0_dp-fa)*(1.0_dp-fb)*fc)*spin_density_grid3(&
           k3a_lower,k3b_lower,k3c_upper)+(fa*fb*(1.0_dp-fc))*spin_density_grid3(k3a_upper,k3b_upper,k3c_lower)+(fa*(1.0_dp-fb)*&
           fc)*spin_density_grid3(k3a_upper,k3b_lower,k3c_upper)+((1.0_dp-fa)*fb*fc)*spin_density_grid3(k3a_lower,k3b_upper,&
           k3c_upper)+(fa*fb*fc)*spin_density_grid3(k3a_upper,k3b_upper,k3c_upper)
         END DO
       END DO
     END DO
!$omp end parallel do
     DEALLOCATE(valence_density_grid3)
     DEALLOCATE(spin_density_grid3)
   END IF
   IF (max_grid_spacing > 1) THEN
!$omp parallel do default(none) &
!$omp private(kc,kb,ka,k2a_lower,k2a_upper,k2b_lower,k2b_upper,k2c_lower,k2c_upper) &
!$omp shared(totnumC,totnumB,totnumA,valence_density_grid2,spin_density_grid2,valence_density,spin_density) &
!$omp schedule(static,1)
     DO kc = 1,totnumC
       DO kb = 1,totnumB
         DO ka = 1,totnumA
           k2a_lower = modulo((nint((ka+0.5_dp)/2.0_dp)-1),totnumA/2)+1
           k2a_upper = modulo((nint((ka+1.5_dp)/2.0_dp)-1),totnumA/2)+1
           k2b_lower = modulo((nint((kb+0.5_dp)/2.0_dp)-1),totnumB/2)+1
           k2b_upper = modulo((nint((kb+1.5_dp)/2.0_dp)-1),totnumB/2)+1
           k2c_lower = modulo((nint((kc+0.5_dp)/2.0_dp)-1),totnumC/2)+1
           k2c_upper = modulo((nint((kc+1.5_dp)/2.0_dp)-1),totnumC/2)+1   
           valence_density(ka,kb,kc) = valence_density(ka,kb,kc)+(1.0_dp/8.0_dp)*(valence_density_grid2(k2a_lower,k2b_lower,&
           k2c_lower)+valence_density_grid2(k2a_upper,k2b_lower,k2c_lower)+valence_density_grid2(k2a_lower,k2b_upper,k2c_lower)+&
           valence_density_grid2(k2a_lower,k2b_lower,k2c_upper)+valence_density_grid2(k2a_upper,k2b_upper,k2c_lower)+&
           valence_density_grid2(k2a_upper,k2b_lower,k2c_upper)+valence_density_grid2(k2a_lower,k2b_upper,k2c_upper)+&
           valence_density_grid2(k2a_upper,k2b_upper,k2c_upper))
           spin_density(ka,kb,kc) = spin_density(ka,kb,kc)+1.0_dp/8.0_dp*(spin_density_grid2(k2a_lower,k2b_lower,&
           k2c_lower)+spin_density_grid2(k2a_upper,k2b_lower,k2c_lower)+spin_density_grid2(k2a_lower,k2b_upper,k2c_lower)+&
           spin_density_grid2(k2a_lower,k2b_lower,k2c_upper)+spin_density_grid2(k2a_upper,k2b_upper,k2c_lower)+&
           spin_density_grid2(k2a_upper,k2b_lower,k2c_upper)+spin_density_grid2(k2a_lower,k2b_upper,k2c_upper)+&
           spin_density_grid2(k2a_upper,k2b_upper,k2c_upper))
         END DO
       END DO
     END DO
!$omp end parallel do
     DEALLOCATE(valence_density_grid2)
     DEALLOCATE(spin_density_grid2)
   END IF    
   WRITE(output_FID,*)'Generation of the valence and spin density for each grid point was successful.'
   FLUSH(output_FID)
 ELSE
   ALLOCATE(valence_density(totnumA,totnumB,totnumC))
   valence_density = 0.0_dp
   IF (max_grid_spacing > 1) THEN
     ALLOCATE(valence_density_grid2((totnumA/2),(totnumB/2),(totnumC/2)))
     valence_density_grid2 = 0.0_dp
   END IF
   IF (max_grid_spacing > 2 ) THEN
     ALLOCATE(valence_density_grid3((totnumA/3),(totnumB/3),(totnumC/3)))
     valence_density_grid3 = 0.0_dp
   END IF
   IF (max_grid_spacing > 3) THEN
     ALLOCATE(valence_density_grid4((totnumA/4),(totnumB/4),(totnumC/4)))
     valence_density_grid4 = 0.0_dp
   END IF
   IF (max_grid_spacing > 4) THEN
     ALLOCATE(valence_density_grid6((totnumA/6),(totnumB/6),(totnumC/6)))
     valence_density_grid6 = 0.0_dp
   END IF
   IF (max_grid_spacing > 6) THEN
     ALLOCATE(valence_density_grid8((totnumA/8),(totnumB/8),(totnumC/8)))
     valence_density_grid8 = 0.0_dp
   END IF
   IF (max_grid_spacing > 8) THEN
     ALLOCATE(valence_density_grid12((totnumA/12),(totnumB/12),(totnumC/12)))
     valence_density_grid12 = 0.0_dp
   END IF
   ALLOCATE(core_density(totnumA,totnumB,totnumC))
   core_density=0.0_dp
   temp = 0.0_dp
   alpha_sum = 0.0_dp
   basis_center_separation = 0.0_dp
   dsquared = 0.0_dp
   Rsquared_cutoff = 0.0_dp
   common_gaussian_factor = 1.0_dp
   
 WRITE(output_FID,*)' '
 WRITE(output_FID,*)'******************** TIME CONTROL ********************'
 CALL system_clock(clock_count2,clock_count_rate,clock_count_max)
 IF (clock_count_rate == 0) THEN
   WRITE(output_FID,*)'Processor clock not available, the time calculation wil not be performed.'
 ELSE
   seconds = modulo((clock_count2-clock_count),clock_count_max)/real(clock_count_rate)
   WRITE(output_FID,*)'Before loop ',seconds, 'seconds'
   FLUSH(output_FID)
 END IF
 WRITE(output_FID,*)'******************************************************'
 WRITE(output_FID,*)' '
 FLUSH(output_FID)    
   
!$omp parallel do default(none) &
!$omp private(block_count,start_index,stop_index,alpha1,X1,Y1,Z1,alpha2,alpha_sum,X2,Y2,Z2,delta_na_block,delta_nb_block,&
!$omp delta_nc_block,Apoints,Bpoints,Cpoints,XD,YD,ZD,Na_block_center,Nb_block_center,Nc_block_center,temp1,temp2,temp3,temp4,&
!$omp temp5,temp6,temp7,temp8,Rsquared_cutoff,common_gaussian_factor,gridspacing,partial_pair_weight_i,partial_pair_weight_j,&
!$omp basis_index_for_i,basis_index_for_j,center_index_for_i,center_index_for_j,atom_index_for_i,atom_index_for_j,min_na,max_na,&
!$omp na,ka,min_nb,max_nb,nb,min_nc,max_nc,nc,renormalization_step,X,Y,Z,dsquared,dX1,dY1,dZ1,dX2,dY2,dZ2,pair_count,k2a,k2b,k2c,&
!$omp k3a,k3b,k3c,k4a,k4b,k4c,k6a,k6b,k6c,k8a,k8b,k8c,k12a,k12b,k12c,max_coeff,analytic_pair_overlap,overlap_sum,abs_overlap_sum,&
!$omp kb,kc)&
!$omp shared(num_pair_blocks,pair_block_data,included_primitive_data,periodicA,totnumA,periodicB,totnumB,periodicC,totnumC,&
!$omp boundary,pixelvolume,valence_density,valence_density_grid2,valence_density_grid3,valence_density_grid4,valence_density_grid6&
!$omp ,valence_density_grid8,valence_density_grid12,core_density,occupancy_correction,primitive_exponents,&
!$omp basis_set_centers) &
!$omp reduction(max:max_pair_correction) &
!$omp schedule(dynamic,1) 
   DO block_count = 1,num_pair_blocks
     !Data common to all primitives in a block
     start_index = pair_block_data(block_count,1)
     stop_index = pair_block_data(block_count,2)
     alpha1=included_primitive_data(start_index,6)
     X1 = included_primitive_data(start_index,10)
     Y1 = included_primitive_data(start_index,11)
     Z1 = included_primitive_data(start_index,12)
     alpha2=included_primitive_data(start_index,13)
     alpha_sum = alpha1 + alpha2
     X2 = included_primitive_data(start_index,17)
     Y2 = included_primitive_data(start_index,18)
     Z2 = included_primitive_data(start_index,19)
     delta_na_block = nint(included_primitive_data(start_index,20))
     delta_nb_block = nint(included_primitive_data(start_index,21))
     delta_nc_block = nint(included_primitive_data(start_index,22))
     ALLOCATE(Apoints(2*delta_na_block + 1))
     ALLOCATE(Bpoints(2*delta_nb_block + 1))
     ALLOCATE(Cpoints(2*delta_nc_block + 1))   
     XD = included_primitive_data(start_index,23)
     YD = included_primitive_data(start_index,24)
     ZD = included_primitive_data(start_index,25)
     Na_block_center = nint(included_primitive_data(start_index,26))
     Nb_block_center = nint(included_primitive_data(start_index,27))
     Nc_block_center = nint(included_primitive_data(start_index,28))
     temp1 = included_primitive_data(start_index,23) - included_primitive_data(start_index,29)
     temp2 = included_primitive_data(start_index,24) - included_primitive_data(start_index,30)
     temp3 = included_primitive_data(start_index,25) - included_primitive_data(start_index,31)
     Rsquared_cutoff = included_primitive_data(start_index,34)
     common_gaussian_factor = included_primitive_data(start_index,35)
     gridspacing = nint(included_primitive_data(start_index,36))
     partial_pair_weight_i = included_primitive_data(start_index,37)
     partial_pair_weight_j = included_primitive_data(start_index,38)
     basis_index_for_i = nint(included_primitive_data(start_index,1))
     basis_index_for_j = nint(included_primitive_data(start_index,2))
     center_index_for_i = nint(primitive_exponents(basis_index_for_i,1))
     center_index_for_j = nint(primitive_exponents(basis_index_for_j,1))
     atom_index_for_i = nint(basis_set_centers(center_index_for_i,6))
     atom_index_for_j = nint(basis_set_centers(center_index_for_j,6))
     min_na = delta_na_block
     max_na = -delta_na_block
     Apoints =  0
     DO na = -delta_na_block,delta_na_block
       IF (periodicA) THEN
         ka = modulo((na + Na_block_center),totnumA) + 1
       ELSE
         ka = na + Na_block_center + 1
       END IF
       Apoints(delta_na_block + na + 1) = ka
       IF ((ka < 1) .or. (ka > totnumA)) THEN
         CYCLE
       ELSE IF ((gridspacing > 1) .and. (modulo(ka,gridspacing) /= (gridspacing - 1))) THEN
         CYCLE
       ELSE
         min_na = min(min_na,na)
         max_na = max(max_na,na)
       END IF  
     END DO
     min_nb = delta_nb_block
     max_nb = -delta_nb_block
     Bpoints = 0
     DO nb = -delta_nb_block,delta_nb_block
       IF (periodicB) THEN
         kb = modulo((nb + Nb_block_center),totnumB) + 1
       ELSE
         kb = nb + Nb_block_center + 1
       END IF
       Bpoints(delta_nb_block + nb + 1) = kb
       IF ((kb < 1) .or. (kb > totnumB)) THEN
         CYCLE
       ELSE IF ((gridspacing > 1) .and. (modulo(kb,gridspacing) /= (gridspacing - 1))) THEN
         CYCLE
       ELSE
         min_nb = min(min_nb,nb)
         max_nb = max(max_nb,nb)
       END IF     
     END DO
     min_nc = delta_nc_block
     max_nc = -delta_nc_block
     Cpoints = 0
     DO nc = -delta_nc_block,delta_nc_block
       IF (periodicC) THEN
         kc = modulo((nc + Nc_block_center),totnumC) + 1
       ELSE
         kc = nc + Nc_block_center + 1
       END IF
       Cpoints(delta_nc_block + nc + 1) = kc
       IF ((kc < 1) .or. (kc > totnumC)) THEN
         CYCLE
       ELSE IF ((gridspacing > 1) .and. (modulo(kc,gridspacing) /= (gridspacing - 1))) THEN
         CYCLE
       ELSE
         min_nc = min(min_nc,nc)
         max_nc = max(max_nc,nc)
       END IF 
     END DO
     IF ((min_na > max_na) .or. (min_nb > max_nb) .or. (min_nc > max_nc)) THEN
       !No valid point found
       CYCLE
     END IF
     DO renormalization_step = 1,2
       IF ((renormalization_step == 2) .and. (gridspacing == 0)) THEN
         CYCLE
       END IF  
       IF (gridspacing == 1) THEN
         DO nc = min_nc,max_nc
           kc = Cpoints(delta_nc_block + nc + 1)
           DO nb = min_nb,max_nb
             kb = Bpoints(delta_nb_block + nb + 1)
             DO na = min_na,max_na
               ka = Apoints(delta_na_block + na + 1)
               X = na*boundary(1,1) + nb*boundary(2,1) + nc*boundary(3,1) + temp1
               Y = na*boundary(1,2) + nb*boundary(2,2) + nc*boundary(3,2) + temp2
               Z = na*boundary(1,3) + nb*boundary(2,3) + nc*boundary(3,3) + temp3
               dsquared = (X - XD)**2 + (Y - YD)**2 + (Z - ZD)**2
               IF (dsquared > Rsquared_cutoff) THEN
                 CYCLE
               END IF   
               dX1 = X-X1
               dY1 = Y-Y1
               dZ1 = Z-Z1
               dX2 = X-X2
               dY2 = Y-Y2
               dZ2 = Z-Z2
               temp4 = exp(-alpha_sum*dsquared)*common_gaussian_factor
               temp5 = temp4*pixelvolume
               temp7 = 0.0_dp
               DO pair_count = start_index,stop_index 
                 IF ((dsquared > included_primitive_data(pair_count,34))) THEN
                   EXIT
                 ELSE IF (included_primitive_data(pair_count,32) == 0.0_dp) THEN
                   CYCLE 
                 END IF
                 temp6 = dX1**nint(included_primitive_data(pair_count,7))*dY1**nint(included_primitive_data(pair_count,8))*dZ1**&
                 nint(included_primitive_data(pair_count,9))*dX2**nint(included_primitive_data(pair_count,14))*dY2**nint&
                 (included_primitive_data(pair_count,15))*dZ2**nint(included_primitive_data(pair_count,16))
                 IF (renormalization_step == 1) THEN
                   included_primitive_data(pair_count,39) = included_primitive_data(pair_count,39) + temp5*temp6
                   included_primitive_data(pair_count,40) = included_primitive_data(pair_count,40) + temp5*abs(temp6)
                 ELSE    
                   temp6 = abs(temp6)
                 END IF
                 temp7 = temp7 + temp6*included_primitive_data(pair_count,32)*included_primitive_data(pair_count,3)
               END DO 
!$omp atomic
               valence_density(ka,kb,kc) = valence_density(ka,kb,kc) + temp7*temp4 
             END DO
           END DO
         END DO      
       ELSE IF (gridspacing == 2) THEN
         !Code for computing the grid2 densities
         DO nc = min_nc,max_nc,2
           kc = Cpoints(delta_nc_block + nc + 1)
           DO nb = min_nb,max_nb,2
             kb = Bpoints(delta_nb_block + nb + 1)
             DO na = min_na,max_na,2
               ka = Apoints(delta_na_block + na + 1)
               X = na*boundary(1,1) + nb*boundary(2,1) + nc*boundary(3,1) + temp1
               Y = na*boundary(1,2) + nb*boundary(2,2) + nc*boundary(3,2) + temp2
               Z = na*boundary(1,3) + nb*boundary(2,3) + nc*boundary(3,3) + temp3
               dsquared = (X - XD)**2 + (Y - YD)**2 + (Z - ZD)**2
               IF (dsquared > Rsquared_cutoff) THEN
                 CYCLE
               END IF
               dX1 = X-X1
               dY1 = Y-Y1
               dZ1 = Z-Z1
               dX2 = X-X2
               dY2 = Y-Y2
               dZ2 = Z-Z2
               temp4 = exp(-alpha_sum*dsquared)*common_gaussian_factor
               temp5 = temp4*pixelvolume*8.0_dp
               temp7 = 0.0_dp
               DO pair_count = start_index,stop_index  
                 IF (dsquared > included_primitive_data(pair_count,34)) THEN
                   EXIT
                 ELSE IF (included_primitive_data(pair_count,32) == 0.0_dp) THEN
                   CYCLE 
                 END IF
                 temp6 = dX1**nint(included_primitive_data(pair_count,7))*dY1**nint(included_primitive_data(pair_count,8))*dZ1&
                 **nint(included_primitive_data(pair_count,9))*dX2**nint(included_primitive_data(pair_count,14))*dY2**nint&
                 (included_primitive_data(pair_count,15))*dZ2**nint(included_primitive_data(pair_count,16))
                 IF (renormalization_step == 1) THEN
                   included_primitive_data(pair_count,39) = included_primitive_data(pair_count,39) + temp5*temp6
                   included_primitive_data(pair_count,40) = included_primitive_data(pair_count,40) + temp5*abs(temp6)
                 ELSE    
                   temp6 = abs(temp6)
                 END IF
                 temp7 = temp7 + temp6*included_primitive_data(pair_count,32)*included_primitive_data(pair_count,3)
               END DO            
               k2a = nint((ka+1)/2.0_dp)
               k2b = nint((kb+1)/2.0_dp)
               k2c = nint((kc+1)/2.0_dp)
!$omp atomic
               valence_density_grid2(k2a,k2b,k2c) = valence_density_grid2(k2a,k2b,k2c) + temp7*temp4
             END DO
           END DO
         END DO
       ELSE IF (gridspacing == 3) THEN
         !Code for computing the grid3 densities
         DO nc = min_nc,max_nc,3
           kc = Cpoints(delta_nc_block + nc + 1)
           DO nb = min_nb,max_nb,3
             kb = Bpoints(delta_nb_block + nb + 1)
             DO na = min_na,max_na,3
               ka = Apoints(delta_na_block + na + 1)
               X = na*boundary(1,1) + nb*boundary(2,1) + nc*boundary(3,1) + temp1
               Y = na*boundary(1,2) + nb*boundary(2,2) + nc*boundary(3,2) + temp2
               Z = na*boundary(1,3) + nb*boundary(2,3) + nc*boundary(3,3) + temp3
               dsquared = (X - XD)**2 + (Y - YD)**2 + (Z - ZD)**2
               IF (dsquared > Rsquared_cutoff) THEN
                 CYCLE
               END IF
               dX1 = X-X1
               dY1 = Y-Y1
               dZ1 = Z-Z1
               dX2 = X-X2
               dY2 = Y-Y2
               dZ2 = Z-Z2
               temp4 = exp(-alpha_sum*dsquared)*common_gaussian_factor
               temp5 = temp4*pixelvolume*27.0_dp
               temp7 = 0.0_dp
               DO pair_count = start_index,stop_index  
                 IF (dsquared > included_primitive_data(pair_count,34)) THEN
                   EXIT
                 ELSE IF (included_primitive_data(pair_count,32) == 0.0_dp) THEN
                   CYCLE 
                 END IF
                 temp6 = dX1**nint(included_primitive_data(pair_count,7))*dY1**nint(included_primitive_data(pair_count,8))*dZ1**&
                 nint(included_primitive_data(pair_count,9))*dX2**nint(included_primitive_data(pair_count,14))*dY2**nint&
                 (included_primitive_data(pair_count,15))*dZ2**nint(included_primitive_data(pair_count,16))
                 IF (renormalization_step == 1) THEN
                   included_primitive_data(pair_count,39) = included_primitive_data(pair_count,39) + temp5*temp6
                   included_primitive_data(pair_count,40) = included_primitive_data(pair_count,40) + temp5*abs(temp6)
                 ELSE    
                   temp6 = abs(temp6)
                 END IF
                 temp7 = temp7 + temp6*included_primitive_data(pair_count,32)*included_primitive_data(pair_count,3)
               END DO             
               k3a = nint((ka+1)/3.0_dp)
               k3b = nint((kb+1)/3.0_dp)
               k3c = nint((kc+1)/3.0_dp)
!$omp atomic
               valence_density_grid3(k3a,k3b,k3c) = valence_density_grid3(k3a,k3b,k3c) + temp7*temp4
             END DO
           END DO
         END DO   
       ELSE IF (gridspacing == 4) THEN
         !Code for computing the grid4 densities
         DO nc = min_nc,max_nc,4
           kc = Cpoints(delta_nc_block + nc + 1) 
           DO nb = min_nb,max_nb,4
             kb = Bpoints(delta_nb_block + nb + 1)
             DO na = min_na,max_na,4
               ka = Apoints(delta_na_block + na + 1)
               X = na*boundary(1,1) + nb*boundary(2,1) + nc*boundary(3,1) + temp1
               Y = na*boundary(1,2) + nb*boundary(2,2) + nc*boundary(3,2) + temp2
               Z = na*boundary(1,3) + nb*boundary(2,3) + nc*boundary(3,3) + temp3
               dsquared = (X - XD)**2 + (Y - YD)**2 + (Z - ZD)**2
               IF (dsquared > Rsquared_cutoff) THEN
                 CYCLE
               END IF
               dX1 = X-X1
               dY1 = Y-Y1
               dZ1 = Z-Z1
               dX2 = X-X2
               dY2 = Y-Y2
               dZ2 = Z-Z2
               temp4 = exp(-alpha_sum*dsquared)*common_gaussian_factor
               temp5 = temp4*pixelvolume*64.0_dp
               temp7 = 0.0_dp
               DO pair_count = start_index,stop_index  
                 IF ((dsquared > included_primitive_data(pair_count,34))) THEN
                   EXIT
                 ELSE IF (included_primitive_data(pair_count,32) == 0.0_dp) THEN
                   CYCLE 
                 END IF
                 temp6 = dX1**nint(included_primitive_data(pair_count,7))*dY1**nint(included_primitive_data(pair_count,8))*dZ1**&
                 nint(included_primitive_data(pair_count,9))*dX2**nint(included_primitive_data(pair_count,14))*dY2**nint&
                 (included_primitive_data(pair_count,15))*dZ2**nint(included_primitive_data(pair_count,16))
                 IF (renormalization_step == 1) THEN
                   included_primitive_data(pair_count,39) = included_primitive_data(pair_count,39) + temp5*temp6
                   included_primitive_data(pair_count,40) = included_primitive_data(pair_count,40) + temp5*abs(temp6)
                 ELSE    
                   temp6 = abs(temp6)
                 END IF
                 temp7 = temp7 + temp6*included_primitive_data(pair_count,32)*included_primitive_data(pair_count,3)
               END DO            
               k4a = nint((ka+1)/4.0_dp)
               k4b = nint((kb+1)/4.0_dp)
               k4c = nint((kc+1)/4.0_dp)
!$omp atomic
               valence_density_grid4(k4a,k4b,k4c) = valence_density_grid4(k4a,k4b,k4c) + temp7*temp4
             END DO
           END DO
         END DO  
       ELSE IF (gridspacing == 6) THEN
         !Code for computing the grid6 densities
         DO nc = min_nc,max_nc,6
           kc = Cpoints(delta_nc_block + nc + 1)
           DO nb = min_nb,max_nb,6
             kb = Bpoints(delta_nb_block + nb + 1)
             DO na = min_na,max_na,6
               ka = Apoints(delta_na_block + na + 1)
               X = na*boundary(1,1) + nb*boundary(2,1) + nc*boundary(3,1) + temp1
               Y = na*boundary(1,2) + nb*boundary(2,2) + nc*boundary(3,2) + temp2
               Z = na*boundary(1,3) + nb*boundary(2,3) + nc*boundary(3,3) + temp3
               dsquared = (X - XD)**2 + (Y - YD)**2 + (Z - ZD)**2
               IF (dsquared > Rsquared_cutoff) THEN
                 CYCLE
               END IF
               dX1 = X-X1
               dY1 = Y-Y1
               dZ1 = Z-Z1
               dX2 = X-X2
               dY2 = Y-Y2
               dZ2 = Z-Z2
               temp4 = exp(-alpha_sum*dsquared)*common_gaussian_factor
               temp5 = temp4*pixelvolume*216.0_dp
               temp7 = 0.0_dp
               DO pair_count = start_index,stop_index  
                 IF (dsquared > included_primitive_data(pair_count,34)) THEN
                   EXIT
                 ELSE IF (included_primitive_data(pair_count,32) == 0.0_dp) THEN
                    CYCLE 
                 END IF
                 temp6 = dX1**nint(included_primitive_data(pair_count,7))*dY1**nint(included_primitive_data(pair_count,8))*dZ1**&
                 nint(included_primitive_data(pair_count,9))*dX2**nint(included_primitive_data(pair_count,14))*dY2**nint&
                 (included_primitive_data(pair_count,15))*dZ2**nint(included_primitive_data(pair_count,16))
                 IF (renormalization_step == 1) THEN
                   included_primitive_data(pair_count,39) = included_primitive_data(pair_count,39) + temp5*temp6
                   included_primitive_data(pair_count,40) = included_primitive_data(pair_count,40) + temp5*abs(temp6)
                 ELSE    
                   temp6 = abs(temp6)
                 END IF
                 temp7 = temp7 + temp6*included_primitive_data(pair_count,32)*included_primitive_data(pair_count,3)
               END DO             
               k6a = nint((ka+1)/6.0_dp)
               k6b = nint((kb+1)/6.0_dp)
               k6c = nint((kc+1)/6.0_dp)
!$omp atomic
               valence_density_grid6(k6a,k6b,k6c) = valence_density_grid6(k6a,k6b,k6c) + temp7*temp4
             END DO
           END DO
         END DO
       ELSE IF (gridspacing == 8) THEN
         !Code for computing the grid8 densities
         DO nc = min_nc,max_nc,8
           kc = Cpoints(delta_nc_block + nc + 1)
           DO nb = min_nb,max_nb,8
             kb = Bpoints(delta_nb_block + nb + 1)
             DO na = min_na,max_na,8
               ka = Apoints(delta_na_block + na + 1)
               X = na*boundary(1,1) + nb*boundary(2,1) + nc*boundary(3,1) + temp1
               Y = na*boundary(1,2) + nb*boundary(2,2) + nc*boundary(3,2) + temp2
               Z = na*boundary(1,3) + nb*boundary(2,3) + nc*boundary(3,3) + temp3
               dsquared = (X - XD)**2 + (Y - YD)**2 + (Z - ZD)**2
               IF (dsquared > Rsquared_cutoff) THEN
                 CYCLE
               END IF
               dX1 = X-X1
               dY1 = Y-Y1
               dZ1 = Z-Z1
               dX2 = X-X2
               dY2 = Y-Y2
               dZ2 = Z-Z2
               temp4 = exp(-alpha_sum*dsquared)*common_gaussian_factor
               temp5 = temp4*pixelvolume*512.0_dp
               temp7 = 0.0_dp
               DO pair_count = start_index,stop_index  
                 IF (dsquared > included_primitive_data(pair_count,34)) THEN
                   EXIT
                 ELSE IF (included_primitive_data(pair_count,32) == 0.0_dp) THEN
                   CYCLE 
                 END IF
                 temp6 = dX1**nint(included_primitive_data(pair_count,7))*dY1**nint(included_primitive_data(pair_count,8))*dZ1&
                 **nint(included_primitive_data(pair_count,9))*dX2**nint(included_primitive_data(pair_count,14))*dY2**nint&
                 (included_primitive_data(pair_count,15))*dZ2**nint(included_primitive_data(pair_count,16))
                 IF (renormalization_step == 1) THEN
                   included_primitive_data(pair_count,39) = included_primitive_data(pair_count,39) + temp5*temp6
                   included_primitive_data(pair_count,40) = included_primitive_data(pair_count,40) + temp5*abs(temp6)
                 ELSE    
                   temp6 = abs(temp6)
                 END IF
                 temp7 = temp7 + temp6*included_primitive_data(pair_count,32)*included_primitive_data(pair_count,3)
               END DO             
               k8a = nint((ka+1)/8.0_dp)
               k8b = nint((kb+1)/8.0_dp)
               k8c = nint((kc+1)/8.0_dp)
!$omp atomic
               valence_density_grid8(k8a,k8b,k8c) = valence_density_grid8(k8a,k8b,k8c) + temp7*temp4
             END DO
           END DO
         END DO
       ELSE IF (gridspacing == 12)  THEN
         !Code for computing the grid12 densities
         DO nc = min_nc,max_nc,12
           kc = Cpoints(delta_nc_block + nc + 1)
           DO nb = min_nb,max_nb,12
             kb = Bpoints(delta_nb_block + nb + 1)
             DO na = min_na,max_na,12
               ka = Apoints(delta_na_block + na + 1)
               X = na*boundary(1,1) + nb*boundary(2,1) + nc*boundary(3,1) + temp1
               Y = na*boundary(1,2) + nb*boundary(2,2) + nc*boundary(3,2) + temp2
               Z = na*boundary(1,3) + nb*boundary(2,3) + nc*boundary(3,3) + temp3
               dsquared = (X - XD)**2 + (Y - YD)**2 + (Z - ZD)**2
               IF (dsquared > Rsquared_cutoff) THEN
                 CYCLE
               END IF
               dX1 = X-X1
               dY1 = Y-Y1
               dZ1 = Z-Z1
               dX2 = X-X2
               dY2 = Y-Y2
               dZ2 = Z-Z2
               temp4 = exp(-alpha_sum*dsquared)*common_gaussian_factor
               temp5 = temp4*pixelvolume*1728.0_dp
               temp7 = 0.0_dp
               DO pair_count = start_index,stop_index  
                 IF (dsquared > included_primitive_data(pair_count,34)) THEN
                   EXIT
                 ELSE IF (included_primitive_data(pair_count,32) == 0.0_dp) THEN
                   CYCLE 
                 END IF
                 temp6 = dX1**nint(included_primitive_data(pair_count,7))*dY1**nint(included_primitive_data(pair_count,8))*dZ1&
                 **nint(included_primitive_data(pair_count,9))*dX2**nint(included_primitive_data(pair_count,14))*dY2**nint&
                 (included_primitive_data(pair_count,15))*dZ2**nint(included_primitive_data(pair_count,16))
                 IF (renormalization_step == 1) THEN
                   included_primitive_data(pair_count,39) = included_primitive_data(pair_count,39) + temp5*temp6
                   included_primitive_data(pair_count,40) = included_primitive_data(pair_count,40) + temp5*abs(temp6)
                 ELSE    
                   temp6 = abs(temp6)
                 END IF
                 temp7 = temp7 + temp6*included_primitive_data(pair_count,32)*included_primitive_data(pair_count,3)
               END DO             
               k12a = nint((ka+1)/12.0_dp)
               k12b = nint((kb+1)/12.0_dp)
               k12c = nint((kc+1)/12.0_dp) 
!$omp atomic
               valence_density_grid12(k12a,k12b,k12c) = valence_density_grid12(k12a,k12b,k12c) + temp7*temp4
             END DO
           END DO
         END DO    
       ELSE IF (gridspacing == 0) THEN
         DO nc = min_nc,max_nc
           kc = Cpoints(delta_nc_block + nc + 1)
           DO nb = min_nb,max_nb
             kb = Bpoints(delta_nb_block + nb + 1)
             DO na = min_na,max_na
               ka = Apoints(delta_na_block + na + 1)
               X = na*boundary(1,1) + nb*boundary(2,1) + nc*boundary(3,1) + temp1
               Y = na*boundary(1,2) + nb*boundary(2,2) + nc*boundary(3,2) + temp2
               Z = na*boundary(1,3) + nb*boundary(2,3) + nc*boundary(3,3) + temp3
               dsquared = (X - XD)**2 + (Y - YD)**2 + (Z - ZD)**2
               IF (dsquared > Rsquared_cutoff) THEN
                 CYCLE
               END IF
               dX1 = X-X1
               dY1 = Y-Y1
               dZ1 = Z-Z1
               dX2 = X-X2
               dY2 = Y-Y2
               dZ2 = Z-Z2
               temp4 = exp(-alpha_sum*dsquared)*common_gaussian_factor
               temp7 = 0.0_dp
               DO pair_count = start_index,stop_index  
                 IF (dsquared > included_primitive_data(pair_count,34)) THEN
                   EXIT
                 END IF
                 temp6 = dX1**nint(included_primitive_data(pair_count,7))*dY1**nint(included_primitive_data(pair_count,8))*dZ1**&
                 nint(included_primitive_data(pair_count,9))*dX2**nint(included_primitive_data(pair_count,14))*dY2**nint&
                 (included_primitive_data(pair_count,15))*dZ2**nint(included_primitive_data(pair_count,16))
                 temp7 = temp7 + temp6*included_primitive_data(pair_count,32)*included_primitive_data(pair_count,3)
               END DO
!$omp atomic
               core_density(ka,kb,kc) = core_density(ka,kb,kc) + temp7*temp4
             END DO                           
           END DO    
         END DO
       END IF
       !Update the renormalization factor
       IF ((gridspacing > 0) .and. (renormalization_step == 1)) THEN
         DO pair_count = start_index,stop_index
           !The following line makes sure the correction is not more than 5 percent of the overlap sum
           max_coeff=abs((included_primitive_data(pair_count,3) + included_primitive_data(pair_count,4))/2.0_dp) + abs((&
           included_primitive_data(pair_count,3) - included_primitive_data(pair_count,4))/2.0_dp)
           analytic_pair_overlap = included_primitive_data(pair_count,5)
           overlap_sum = included_primitive_data(pair_count,39)
           abs_overlap_sum = included_primitive_data(pair_count,40)
           IF (max_coeff*abs_overlap_sum > gaussian_overlap_tolerance) THEN
             included_primitive_data(pair_count,32) = max(min(((analytic_pair_overlap - overlap_sum)/abs_overlap_sum),0.05_dp),&
             -0.05_dp)
             max_pair_correction = max(abs((analytic_pair_overlap - overlap_sum)/abs_overlap_sum), max_pair_correction)
           ELSE
             included_primitive_data(pair_count,32) = 0.0_dp
           END IF
         END DO     
       END IF
     END DO
     DEALLOCATE(Apoints)
     DEALLOCATE(Bpoints)
     DEALLOCATE(Cpoints)
   END DO   
!$omp end parallel do
   !Interpolate the coarser grids back onto the regular grid
   
 WRITE(output_FID,*)' '
 WRITE(output_FID,*)'******************** TIME CONTROL ********************'
 CALL system_clock(clock_count2,clock_count_rate,clock_count_max)
 IF (clock_count_rate == 0) THEN
   WRITE(output_FID,*)'Processor clock not available, the time calculation wil not be performed.'
 ELSE
   seconds = modulo((clock_count2-clock_count),clock_count_max)/real(clock_count_rate)
   WRITE(output_FID,*)'Before interpolation ',seconds, 'seconds'
   FLUSH(output_FID)
 END IF
 WRITE(output_FID,*)'******************************************************'
 WRITE(output_FID,*)' '
 FLUSH(output_FID)    

   IF (max_grid_spacing > 8) THEN
!$omp parallel do default(none) &
!$omp private(k6c,k6b,k6a,k12a_lower,k12a_upper,k12b_lower,k12b_upper,k12c_lower,k12c_upper) &
!$omp shared(totnumC,totnumB,totnumA,valence_density_grid12,valence_density_grid6) &
!$omp schedule(static,1)
     DO k6c = 1,totnumC/6
       DO k6b = 1,totnumB/6
         DO k6a = 1,totnumA/6
           k12a_lower = modulo((nint((k6a - 0.5_dp)/2.0_dp)-1),totnumA/12)+1
           k12a_upper = modulo((nint((k6a + 0.5_dp)/2.0_dp)-1),totnumA/12)+1
           k12b_lower = modulo((nint((k6b - 0.5_dp)/2.0_dp)-1),totnumB/12)+1
           k12b_upper = modulo((nint((k6b + 0.5_dp)/2.0_dp)-1),totnumB/12)+1
           k12c_lower = modulo((nint((k6c - 0.5_dp)/2.0_dp)-1),totnumC/12)+1
           k12c_upper = modulo((nint((k6c + 0.5_dp)/2.0_dp)-1),totnumC/12)+1
           valence_density_grid6(k6a,k6b,k6c) = valence_density_grid6(k6a,k6b,k6c)+1.0_dp/8.0_dp*(valence_density_grid12(&
           k12a_lower,k12b_lower,k12c_lower)+valence_density_grid12(k12a_upper,k12b_lower,k12c_lower)+valence_density_grid12(&
           k12a_lower,k12b_upper,k12c_lower)+valence_density_grid12(k12a_lower,k12b_lower,k12c_upper)+valence_density_grid12(&
           k12a_upper,k12b_upper,k12c_lower)+valence_density_grid12(k12a_upper,k12b_lower,k12c_upper)+valence_density_grid12(&
           k12a_lower,k12b_upper,k12c_upper)+valence_density_grid12(k12a_upper,k12b_upper,k12c_upper))
         END DO
       END DO
     END DO
!$omp end parallel do
     DEALLOCATE(valence_density_grid12)
   END IF
   IF (max_grid_spacing > 6) THEN
!$omp parallel do default(none) &
!$omp private(k4c,k4b,k4a,k8a_lower,k8a_upper,k8b_lower,k8b_upper,k8c_lower,k8c_upper) &
!$omp shared(totnumC,totnumB,totnumA,valence_density_grid8,valence_density_grid4) &
!$omp schedule(static,1)
     DO k4c = 1,totnumC/4
       DO k4b = 1,totnumB/4
         DO k4a = 1,totnumA/4
           k8a_lower = modulo((nint((k4a - 0.5_dp)/2.0_dp)-1),totnumA/8)+1
           k8a_upper = modulo((nint((k4a + 0.5_dp)/2.0_dp)-1),totnumA/8)+1
           k8b_lower = modulo((nint((k4b - 0.5_dp)/2.0_dp)-1),totnumB/8)+1
           k8b_upper = modulo((nint((k4b + 0.5_dp)/2.0_dp)-1),totnumB/8)+1
           k8c_lower = modulo((nint((k4c - 0.5_dp)/2.0_dp)-1),totnumC/8)+1
           k8c_upper = modulo((nint((k4c + 0.5_dp)/2.0_dp)-1),totnumC/8)+1
           valence_density_grid4(k4a,k4b,k4c) = valence_density_grid4(k4a,k4b,k4c)+1.0_dp/8.0_dp*(valence_density_grid8(k8a_lower&
           ,k8b_lower,k8c_lower)+valence_density_grid8(k8a_upper,k8b_lower,k8c_lower)+valence_density_grid8(k8a_lower,k8b_upper,&
           k8c_lower)+valence_density_grid8(k8a_lower,k8b_lower,k8c_upper)+valence_density_grid8(k8a_upper,k8b_upper,k8c_lower)+&
           valence_density_grid8(k8a_upper,k8b_lower,k8c_upper)+valence_density_grid8(k8a_lower,k8b_upper,k8c_upper)+&
           valence_density_grid8(k8a_upper,k8b_upper,k8c_upper))
         END DO
       END DO
     END DO
!$omp end parallel do
     DEALLOCATE(valence_density_grid8)
   END IF
   IF (max_grid_spacing > 4) THEN
!$omp parallel do default(none) &
!$omp private(k3c,k3b,k3a,k6a_lower,k6a_upper,k6b_lower,k6b_upper,k6c_lower,k6c_upper) &
!$omp shared(totnumC,totnumB,totnumA,valence_density_grid6,valence_density_grid3) &
!$omp schedule(static,1)
     DO k3c = 1,totnumC/3
       DO k3b = 1,totnumB/3
         DO k3a = 1,totnumA/3
           k6a_lower = modulo((nint((k3a - 0.5_dp)/2.0_dp)-1),totnumA/6)+1
           k6a_upper = modulo((nint((k3a + 0.5_dp)/2.0_dp)-1),totnumA/6)+1
           k6b_lower = modulo((nint((k3b - 0.5_dp)/2.0_dp)-1),totnumB/6)+1
           k6b_upper = modulo((nint((k3b + 0.5_dp)/2.0_dp)-1),totnumB/6)+1
           k6c_lower = modulo((nint((k3c - 0.5_dp)/2.0_dp)-1),totnumC/6)+1
           k6c_upper = modulo((nint((k3c + 0.5_dp)/2.0_dp)-1),totnumC/6)+1
           valence_density_grid3(k3a,k3b,k3c) = valence_density_grid3(k3a,k3b,k3c)+1.0_dp/8.0_dp*(valence_density_grid6(k6a_lower&
           ,k6b_lower,k6c_lower)+valence_density_grid6(k6a_upper,k6b_lower,k6c_lower)+valence_density_grid6(k6a_lower,k6b_upper,&
           k6c_lower)+valence_density_grid6(k6a_lower,k6b_lower,k6c_upper)+valence_density_grid6(k6a_upper,k6b_upper,k6c_lower)+&
           valence_density_grid6(k6a_upper,k6b_lower,k6c_upper)+valence_density_grid6(k6a_lower,k6b_upper,k6c_upper)+&
           valence_density_grid6(k6a_upper,k6b_upper,k6c_upper))
         END DO
       END DO
     END DO
!$omp end parallel do
     DEALLOCATE(valence_density_grid6)
   END IF
   IF (max_grid_spacing > 3) THEN
!$omp parallel do default(none) &
!$omp private(k2c,k2b,k2a,k4a_lower,k4a_upper,k4b_lower,k4b_upper,k4c_lower,k4c_upper) &
!$omp shared(totnumC,totnumB,totnumA,valence_density_grid4,valence_density_grid2) &
!$omp schedule(static,1)
     DO k2c = 1,totnumC/2
       DO k2b = 1,totnumB/2
         DO k2a = 1,totnumA/2
           k4a_lower = modulo((nint((k2a - 0.5_dp)/2.0_dp)-1),totnumA/4)+1
           k4a_upper = modulo((nint((k2a + 0.5_dp)/2.0_dp)-1),totnumA/4)+1
           k4b_lower = modulo((nint((k2b - 0.5_dp)/2.0_dp)-1),totnumB/4)+1
           k4b_upper = modulo((nint((k2b + 0.5_dp)/2.0_dp)-1),totnumB/4)+1
           k4c_lower = modulo((nint((k2c - 0.5_dp)/2.0_dp)-1),totnumC/4)+1
           k4c_upper = modulo((nint((k2c + 0.5_dp)/2.0_dp)-1),totnumC/4)+1
           valence_density_grid2(k2a,k2b,k2c) = valence_density_grid2(k2a,k2b,k2c)+1.0_dp/8.0_dp*(valence_density_grid4(k4a_lower&
           ,k4b_lower,k4c_lower)+valence_density_grid4(k4a_upper,k4b_lower,k4c_lower)+valence_density_grid4(k4a_lower,k4b_upper,&
           k4c_lower)+valence_density_grid4(k4a_lower,k4b_lower,k4c_upper)+valence_density_grid4(k4a_upper,k4b_upper,k4c_lower)+&
           valence_density_grid4(k4a_upper,k4b_lower,k4c_upper)+valence_density_grid4(k4a_lower,k4b_upper,k4c_upper)+&
           valence_density_grid4(k4a_upper,k4b_upper,k4c_upper))
         END DO
       END DO
     END DO
!$omp end parallel do
     DEALLOCATE(valence_density_grid4)
   END IF
   IF (max_grid_spacing > 2) THEN
!$omp parallel do default(none) &
!$omp private(kc,kb,ka,k3a_lower,k3a_upper,k3b_lower,k3b_upper,k3c_lower,k3c_upper,fa,fb,fc) &
!$omp shared(totnumC,totnumB,totnumA,valence_density_grid3,valence_density) &
!$omp schedule(static,1)
     DO kc = 1,totnumC
       DO kb = 1,totnumB
         DO ka = 1,totnumA
           k3a_lower = modulo((nint(ka/3.0_dp)-1),totnumA/3)+1
           k3a_upper = modulo((nint((ka+2.0_dp)/3.0_dp)-1),totnumA/3)+1
           k3b_lower = modulo((nint(kb/3.0_dp)-1),totnumB/3)+1
           k3b_upper = modulo((nint((kb+2.0_dp)/3.0_dp)-1),totnumB/3)+1
           k3c_lower = modulo((nint(kc/3.0_dp)-1),totnumC/3)+1
           k3c_upper = modulo((nint((kc+2.0_dp)/3.0_dp)-1),totnumC/3)+1
           fa = modulo((ka+1),3)/3.0_dp
           fb = modulo((kb+1),3)/3.0_dp
           fc = modulo((kc+1),3)/3.0_dp
           valence_density(ka,kb,kc)=valence_density(ka,kb,kc)+((1.0_dp-fa)*(1.0_dp-fb)*(1.0_dp-fc))*valence_density_grid3(&
           k3a_lower,k3b_lower,k3c_lower)+(fa*(1.0_dp-fb)*(1.0_dp-fc))*valence_density_grid3(k3a_upper,k3b_lower,k3c_lower)+((&
           1.0_dp-fa)*fb*(1.0_dp-fc))*valence_density_grid3(k3a_lower,k3b_upper,k3c_lower)+((1.0_dp-fa)*(1.0_dp-fb)*fc)*&
           valence_density_grid3(k3a_lower,k3b_lower,k3c_upper)+(fa*fb*(1.0_dp-fc))*valence_density_grid3(k3a_upper,k3b_upper,&
           k3c_lower)+(fa*(1.0_dp-fb)*fc)*valence_density_grid3(k3a_upper,k3b_lower,k3c_upper)+((1.0_dp-fa)*fb*fc)*&
           valence_density_grid3(k3a_lower,k3b_upper,k3c_upper)+(fa*fb*fc)*valence_density_grid3(k3a_upper,k3b_upper,k3c_upper)
         END DO
       END DO
     END DO
!$omp end parallel do
     DEALLOCATE(valence_density_grid3)
   END IF
   IF (max_grid_spacing > 1) THEN
!$omp parallel do default(none) &
!$omp private(kc,kb,ka,k2a_lower,k2a_upper,k2b_lower,k2b_upper,k2c_lower,k2c_upper) &
!$omp shared(totnumC,totnumB,totnumA,valence_density_grid2,valence_density) &
!$omp schedule(static,1)
     DO kc = 1,totnumC
       DO kb = 1,totnumB
         DO ka = 1,totnumA
           k2a_lower = modulo((nint((ka+0.5_dp)/2.0_dp)-1),totnumA/2)+1
           k2a_upper = modulo((nint((ka+1.5_dp)/2.0_dp)-1),totnumA/2)+1
           k2b_lower = modulo((nint((kb+0.5_dp)/2.0_dp)-1),totnumB/2)+1
           k2b_upper = modulo((nint((kb+1.5_dp)/2.0_dp)-1),totnumB/2)+1
           k2c_lower = modulo((nint((kc+0.5_dp)/2.0_dp)-1),totnumC/2)+1
           k2c_upper = modulo((nint((kc+1.5_dp)/2.0_dp)-1),totnumC/2)+1   
           valence_density(ka,kb,kc) = valence_density(ka,kb,kc)+(1.0_dp/8.0_dp)*(valence_density_grid2(k2a_lower,k2b_lower,&
           k2c_lower)+valence_density_grid2(k2a_upper,k2b_lower,k2c_lower)+valence_density_grid2(k2a_lower,k2b_upper,k2c_lower)+&
           valence_density_grid2(k2a_lower,k2b_lower,k2c_upper)+valence_density_grid2(k2a_upper,k2b_upper,k2c_lower)+&
           valence_density_grid2(k2a_upper,k2b_lower,k2c_upper)+valence_density_grid2(k2a_lower,k2b_upper,k2c_upper)+&
           valence_density_grid2(k2a_upper,k2b_upper,k2c_upper))
         END DO
       END DO
     END DO
!$omp end parallel do
     DEALLOCATE(valence_density_grid2)
   END IF    
   WRITE(output_FID,*)'Generation of the valence and spin density for each grid point was successful.'
   FLUSH(output_FID)
 END IF
 WRITE(output_FID,'(a,es13.4)')' The maximum pair correction was: ',max_pair_correction
 FLUSH(output_FID)

 WRITE(output_FID,*)' '
 WRITE(output_FID,*)'******************** TIME CONTROL ********************'
 CALL system_clock(clock_count2,clock_count_rate,clock_count_max)
 IF (clock_count_rate == 0) THEN
   WRITE(output_FID,*)'Processor clock not available, the time calculation wil not be performed.'
 ELSE
   seconds = modulo((clock_count2-clock_count),clock_count_max)/real(clock_count_rate)
   WRITE(output_FID,*)'Finished Step 7 (loop + interpolation) in ',seconds, 'seconds'
   FLUSH(output_FID)
 END IF
 WRITE(output_FID,*)'******************************************************'
 WRITE(output_FID,*)' '
 FLUSH(output_FID)

 !----------------------------------------------------------------------------------
 ! Step 8:  Generate the core density grid
 !---------------------------------------------------------------------------------- 
 IF (edf_available) THEN
   ALLOCATE(Apoints(2*delta_na + 1))
   ALLOCATE(Bpoints(2*delta_nb + 1))
   ALLOCATE(Cpoints(2*delta_nc + 1))
   DO j=1,n_edf_primitives
     edf_primitives(j,7)= basis_set_centers(nint(edf_primitives(j,1)),3) !edf center X
     edf_primitives(j,8)= basis_set_centers(nint(edf_primitives(j,1)),4) !edf center Y
     edf_primitives(j,9)= basis_set_centers(nint(edf_primitives(j,1)),5) !edf center Z
     xyz=[edf_primitives(j,7),edf_primitives(j,8),edf_primitives(j,9)]
     temp_vector = matmul(inv_boundary,(xyz - origin))
     edf_primitives(j,10)=nint(temp_vector(1)) !edf center na
     edf_primitives(j,11)=nint(temp_vector(2)) !edf center nb
     edf_primitives(j,12)=nint(temp_vector(3)) !edf center nc
     temp_center_shift = matmul(TRANSPOSE(boundary),(temp_vector) - nint(temp_vector))
     edf_primitives(j,13)= temp_center_shift(1) !edf center shift x
     edf_primitives(j,14)= temp_center_shift(2) !edf center shift y
     edf_primitives(j,15)= temp_center_shift(3) !edf center shift z
   END DO
   DO j=1,n_edf_primitives
     Apoints = 0
     DO na = -delta_na,delta_na
       IF (periodicA) THEN
         Apoints(delta_na + na + 1) = modulo((na + nint(edf_primitives(j,10))),totnumA) + 1
       ELSE
         Apoints(delta_na + na + 1) = nint(na + edf_primitives(j,10) + 1)
       END IF
     END DO
     Bpoints = 0
     DO nb = -delta_nb,delta_nb
       IF (periodicB) THEN
         Bpoints(delta_nb + nb + 1) = modulo((nb + nint(edf_primitives(j,11))),totnumB) + 1
       ELSE
         Bpoints(delta_nb + nb + 1) = nint(nb + edf_primitives(j,11) + 1)
       END IF
     END DO
     Cpoints = 0
     DO nc = -delta_nc,delta_nc
       IF (periodicC) THEN
         Cpoints(delta_nc + nc + 1) = modulo((nc + nint(edf_primitives(j,12))),totnumC) + 1
       ELSE
         Cpoints(delta_nc + nc + 1) = nint(nc + edf_primitives(j,12) + 1)
       END IF
     END DO
     DO nc = -delta_nc,delta_nc
       kc = Cpoints(delta_nc + nc + 1)
       IF (kc < 1 .or. kc > totnumC) THEN
         CYClE
       END IF
       DO nb = -delta_nb,delta_nb
         kb = Bpoints(delta_nb + nb + 1)
         IF (kb < 1 .or. kb > totnumB) THEN
           CYCLE
         END IF
         DO na = -delta_na,delta_na
           ka = Apoints(delta_na + na + 1)
           IF (ka < 1 .or. ka > totnumA) THEN 
             CYCLE
           END IF
           temp_vector(1) = na*boundary(1,1) + nb*boundary(2,1) + nc*boundary(3,1) - edf_primitives(j,13)
           temp_vector(2) = na*boundary(1,2) + nb*boundary(2,2) + nc*boundary(3,2) - edf_primitives(j,14)
           temp_vector(3) = na*boundary(1,3) + nb*boundary(2,3) + nc*boundary(3,3) - edf_primitives(j,15)
           distance = sqrt(temp_vector(1)*temp_vector(1) + temp_vector(2)*temp_vector(2) + temp_vector(3)*temp_vector(3))
           temp_index = ceiling(scalefactor*distance + zero_tolerance)
           IF (temp_index <= nshells) THEN
             alpha=edf_primitives(j,2)
             Lx=nint(edf_primitives(j,3))
             Ly=nint(edf_primitives(j,4))
             Lz=nint(edf_primitives(j,5))
             Ax=edf_primitives(j,7)
             Ay=edf_primitives(j,8)
             Az=edf_primitives(j,9)
             X = temp_vector(1) + edf_primitives(j,7)
             Y = temp_vector(2) + edf_primitives(j,8)
             Z = temp_vector(3) + edf_primitives(j,9)
             temp_scalar = gaussian_value(edf_primitives(j,6),alpha,Lx,Ly,Lz,Ax,Ay,Az,X,Y,Z)
             core_density(ka,kb,kc) = core_density(ka,kb,kc) + temp_scalar
           END IF
         END DO    
       END DO
     END DO
   END DO
   DEALLOCATE(Apoints)
   DEALLOCATE(Bpoints)
   DEALLOCATE(Cpoints)     
 END IF
 !Make sure the core and total electron densities are positive definite
 DO kc=1,totnumC
   DO kb=1,totnumB
      DO ka=1,totnumA
       core_density(ka,kb,kc) = max(core_density(ka,kb,kc),0.0_dp)
       valence_density(ka,kb,kc) = max(valence_density(ka,kb,kc),-core_density(ka,kb,kc))
     END DO
   END DO
 END DO
 DO j=1,natoms
   IF (core_electrons(j) < zero_tolerance) THEN
     occupancy_correction(1,j) = core_electrons(j)
     core_electrons(j) = 0.0_dp
   END IF
 END DO
 WRITE(output_FID,'(a,es13.4)')' The raw number of core electrons is: ', sum(core_density)*pixelvolume
 FLUSH(output_FID)
 CALL add_missing_core_density()
 CALL initialize_atomic_densities()
 
 sum_negative_density=0.0_dp
 
 DEALLOCATE(included_primitive_pairs)
 DEALLOCATE(alpha_1PDM)
 DEALLOCATE(beta_1PDM)
 DEALLOCATE(included_primitive_data)
 DEALLOCATE(pair_block_data)
 DEALLOCATE(pair_overlap)

 WRITE(output_FID,*)' '
 WRITE(output_FID,*)'******************** TIME CONTROL ********************'
 CALL system_clock(clock_count2,clock_count_rate,clock_count_max)
 IF (clock_count_rate == 0) THEN
   WRITE(output_FID,*)'Processor clock not available, the time calculation wil not be performed.'
 ELSE
   seconds = modulo((clock_count2-clock_count),clock_count_max)/real(clock_count_rate)
   WRITE(output_FID,*)'Finished Step 8 in ',seconds, 'seconds'
   FLUSH(output_FID)
 END IF
 WRITE(output_FID,*)'******************************************************'
 WRITE(output_FID,*)' '
 FLUSH(output_FID)

 WRITE(output_FID,*)' '
 WRITE(output_FID,*)'******************** TIME CONTROL ********************'
 CALL system_clock(clock_count2,clock_count_rate,clock_count_max)
 IF (clock_count_rate == 0) THEN
   WRITE(output_FID,*)'Processor clock not available, the time calculation wil not be performed.'
 ELSE
   seconds = modulo((clock_count2-clock_count),clock_count_max)/real(clock_count_rate)
   WRITE(output_FID,*)'Finished generate_density_grids_from_gaussian_basis_set_coefficients in ',seconds, 'seconds'
   FLUSH(output_FID)
 END IF
 WRITE(output_FID,*)'******************************************************'
 WRITE(output_FID,*)' '
 FLUSH(output_FID) 
 
 END SUBROUTINE generate_density_grids_from_gaussian_basis_set_coefficients

 END MODULE module_gen_dens_grids_from_gaussian_basis_set_coeff