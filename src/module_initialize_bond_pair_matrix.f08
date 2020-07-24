 MODULE module_initialize_bond_pair_matrix
 ! Copyright (c) 2014, 2015, 2016, 2017 by Thomas A. Manz and Nidia Gabaldon Limas. Rights reserved.
 
 USE module_common_variable_declarations
 USE module_global_parameters
 USE module_determine_whether_pair_overlap_is_significant

 IMPLICIT NONE
 
 CONTAINS
 
 SUBROUTINE initialize_bond_pair_matrix()
 !===================================================================================
 
 REAL(kind=dp) :: bond_cutoff_radius,X1,Y1,Z1,X2,Y2,Z2,distance,a_dot_a,a_dot_b,a_dot_c,b_dot_b,b_dot_c,c_dot_c,cos_a_b,&
 cos_a_c,cos_b_c,sin_a_b,sin_a_c,sin_b_c,multiplier,mid_na,mid_nb,mid_nc
 INTEGER :: k,periodicsumlimitA,periodicsumlimitB,periodicsumlimitC,temp_count,repeata,repeatb,repeatc
 
 INTEGER :: clock_count,clock_count_rate,clock_count_max,clock_count2
 
 WRITE(output_FID,*)' '
 WRITE(output_FID,*)'******************** TIME CONTROL ********************'
 WRITE(output_FID,*)'Starting initialize_bond_pair_matrix'
 FLUSH(output_FID)
 CALL system_clock(clock_count,clock_count_rate,clock_count_max)
 WRITE(output_FID,*)'******************************************************'
 WRITE(output_FID,*)' '
 FLUSH(output_FID) 
 
 WRITE(output_FID,*)'Initializing the bond pair matrix ...'
 !Determine the number of unit cell replications for initializing the bond pair matrix
 bond_cutoff_radius = 2.0_dp*cutoff_radius*bohrperangstrom/100.0_dp !in bohr
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
 IF (periodicA) THEN
   periodicsumlimitA = ceiling(max((bond_cutoff_radius/sin_a_b),(bond_cutoff_radius/sin_a_c))/sqrt(a_dot_a)) + 1
 ELSE
   periodicsumlimitA = 0
 END IF
 IF (periodicB) THEN
   periodicsumlimitB = ceiling(max((bond_cutoff_radius/sin_a_b),(bond_cutoff_radius/sin_b_c))/sqrt(b_dot_b)) + 1
 ELSE
   periodicsumlimitB = 0
 END IF
 IF (periodicC) THEN
   periodicsumlimitC = ceiling(max((bond_cutoff_radius/sin_a_c),(bond_cutoff_radius/sin_b_c))/sqrt(c_dot_c)) + 1
 ELSE
   periodicsumlimitC = 0
 END IF
 !nA_outside is used in the cutoffs below to determine which atom
 !Pairs are important to include in the bond pair matrix
 ALLOCATE(nA_outside(nshells,natoms))
 nA_outside = 0.0_dp
 DO j = 1,natoms
   nA_outside(nshells,j) = corrected_spherical_avg_density(nshells,j)*sum_points(nshells,j)*pixelvolume
   DO k = 1,(nshells - 1)
     nA_outside((nshells - k),j) = nA_outside((nshells - k + 1),j) + corrected_spherical_avg_density((nshells - k),j)*sum_points&
     ((nshells - k),j)*pixelvolume
   END DO
 END DO 
 !Determine the number of symmetry unique atom pairs to include in the bond pair matrix
 temp_count = 0
 DO i=1,natoms
   DO j=i,natoms
     DO repeatc = -periodicsumlimitC,periodicsumlimitC
       DO repeatb = -periodicsumlimitB,periodicsumlimitB
         DO repeata = -periodicsumlimitA,periodicsumlimitA
           active_pair = 0
           X1=coords(1,i)
           Y1=coords(2,i)
           Z1=coords(3,i)
           X2=coords(1,j) + repeata*vector1(1) + repeatb*vector2(1) + repeatc*vector3(1)
           Y2=coords(2,j) + repeata*vector1(2) + repeatb*vector2(2) + repeatc*vector3(2)
           Z2=coords(3,j) + repeata*vector1(3) + repeatb*vector2(3) + repeatc*vector3(3)
           distance = sqrt((X1-X2)**2 + (Y1-Y2)**2 + (Z1-Z2)**2)
           IF (distance > bond_cutoff_radius) THEN
             CYCLE
             !Determine if the pair is symmetrically unique
           ELSE IF (j > i) THEN
             active_pair = 1
           ELSE IF (repeata > 0) THEN
             active_pair = 1
           ELSE IF ((repeata == 0) .and. (repeatb > 0)) THEN
             active_pair = 1
           ELSE IF (((repeata == 0) .and. (repeatb == 0)) .and. (repeatc > 0)) THEN
             active_pair = 1
           ELSE
             CYCLE
           END IF
           CALL determine_whether_pair_overlap_is_significant(distance)
           IF (active_pair == 1)  THEN
             temp_count = temp_count + 1
           ELSE
             CYCLE
           END IF
         END DO
       END DO
     END DO
   END DO
 END DO
 WRITE(output_FID,*)'The total number of symmetry unique atom pairs to include in the bond pair matrix is'
 num_sym_unique_bond_pairs = temp_count
 WRITE(output_FID,'(A,i6)')' num_sym_unique_bond_pairs: ',num_sym_unique_bond_pairs
 FLUSH(output_FID)
 IF (temp_count == 0) THEN
   WRITE(output_FID,*)'The are no symmetry unique atom pairs to include in the bond pair matrix.'
   WRITE(output_FID,*)'Because there are no significant BOs to compute, BO analysis will be skipped.'
   FLUSH(output_FID)
   flag = 1
   RETURN  !Retruns to module_perform_BO_analysis
 END IF
 !Initialize the bond pair matrix
 ALLOCATE(bond_pair_matrix(20,temp_count))
 bond_pair_matrix = 0.0_dp
 temp_count = 0
 DO i=1,natoms
   DO j=i,natoms
     DO repeata = -periodicsumlimitA,periodicsumlimitA
       DO repeatb = -periodicsumlimitB,periodicsumlimitB
         DO repeatc = -periodicsumlimitC,periodicsumlimitC
           X1=coords(1,i)
           Y1=coords(2,i)
           Z1=coords(3,i)
           X2=coords(1,j) + repeata*vector1(1) + repeatb*vector2(1) + repeatc*vector3(1)
           Y2=coords(2,j) + repeata*vector1(2) + repeatb*vector2(2) + repeatc*vector3(2)
           Z2=coords(3,j) + repeata*vector1(3) + repeatb*vector2(3) + repeatc*vector3(3)
           distance = sqrt((X1-X2)**2 + (Y1-Y2)**2 + (Z1-Z2)**2)
           IF (distance > bond_cutoff_radius) THEN
             CYCLE
             !Determine if the pair is symmetrically unique
           ELSE IF (j > i) THEN
             active_pair = 1
           ELSE IF (repeata > 0) THEN
             active_pair = 1
           ELSE IF (((repeata == 0) .and. (repeatb > 0))) THEN
             active_pair = 1
           ELSE IF (((repeata == 0) .and. (repeatb == 0)) .and. (repeatc > 0)) THEN
             active_pair = 1
           ELSE
             CYCLE
           END IF          
           CALL determine_whether_pair_overlap_is_significant(distance)
           IF (active_pair == 1) THEN
             temp_count = temp_count + 1
           ELSE
             CYCLE
           END IF
           !Initialize the bond pair matrix
           bond_pair_matrix(1,temp_count) = i
           bond_pair_matrix(2,temp_count) = j
           bond_pair_matrix(3,temp_count) = repeata
           bond_pair_matrix(4,temp_count) = repeatb
           bond_pair_matrix(5,temp_count) = repeatc
           !Compute the relevant shared parallelpiped for the atom pair
           !Compute the minimum na value for shared parallelpiped
           multiplier = sqrt(1.0_dp - ((0.5_dp*distance*100.0_dp)/(bohrperangstrom*cutoff_radius))**2)
           mid_na = (center_nabc(1,i) + center_nabc(1,j) + repeata*totnumA)*0.5_dp
           mid_nb = (center_nabc(2,i) + center_nabc(2,j) + repeatb*totnumB)*0.5_dp
           mid_nc = (center_nabc(3,i) + center_nabc(3,j) + repeatc*totnumC)*0.5_dp
           bond_pair_matrix(6,temp_count) = -delta_na + max(center_nabc(1,i),(center_nabc(1,j) + repeata*totnumA))
           bond_pair_matrix(6,temp_count) = max(bond_pair_matrix(6,temp_count), nint(-delta_na*multiplier + mid_na) - 2.0_dp)
           !Compute the maximum na value for shared parallelpiped
           bond_pair_matrix(7,temp_count) = delta_na + min(center_nabc(1,i),(center_nabc(1,j) + repeata*totnumA))
           bond_pair_matrix(7,temp_count) = min(bond_pair_matrix(7,temp_count), nint(delta_na*multiplier + mid_na) + 2.0_dp)
           !Compute the minimum nb value for shared parallelpiped
           bond_pair_matrix(8,temp_count) = -delta_nb + max(center_nabc(2,i),(center_nabc(2,j) + repeatb*totnumB))
           bond_pair_matrix(8,temp_count) = max(bond_pair_matrix(8,temp_count), nint(-delta_nb*multiplier + mid_nb) - 2.0_dp)
           !Compute the maximum nb value for shared parallelpiped
           bond_pair_matrix(9,temp_count) = delta_nb + min(center_nabc(2,i),(center_nabc(2,j) + repeatb*totnumB))
           bond_pair_matrix(9,temp_count) = min(bond_pair_matrix(9,temp_count), nint(delta_nb*multiplier + mid_nb) + 2.0_dp)
           !Compute the minimum nc value for shared parallelpiped
           bond_pair_matrix(10,temp_count) = -delta_nc + max(center_nabc(3,i),(center_nabc(3,j) + repeatc*totnumC))
           bond_pair_matrix(10,temp_count) = max(bond_pair_matrix(10,temp_count), nint(-delta_nc*multiplier + mid_nc) - 2.0_dp) 
           !Compute the maximum nc value for shared parallelpiped
           bond_pair_matrix(11,temp_count) = delta_nc + min(center_nabc(3,i),(center_nabc(3,j) + repeatc*totnumC))
           bond_pair_matrix(11,temp_count) = min(bond_pair_matrix(11,temp_count), nint(delta_nc*multiplier + mid_nc) + 2.0_dp)
         END DO
       END DO
     END DO
   END DO
 END DO
 WRITE(output_FID,*)'The bond pair matrix has been initialized.' 
 
  WRITE(output_FID,*)' '
 WRITE(output_FID,*)'******************** TIME CONTROL ********************'
 CALL system_clock(clock_count2,clock_count_rate,clock_count_max)
 IF (clock_count_rate == 0) THEN
   WRITE(output_FID,*)'Processor clock not available, the time calculation wil not be performed.'
 ELSE
   seconds = modulo((clock_count2-clock_count),clock_count_max)/real(clock_count_rate)
   WRITE(output_FID,*)'Finished initialize_bond_pair_matrix in ',seconds, 'seconds'
   FLUSH(output_FID)
 END IF
 WRITE(output_FID,*)'******************************************************'
 WRITE(output_FID,*)' '
 FLUSH(output_FID) 
 
 END SUBROUTINE initialize_bond_pair_matrix
 
 END MODULE module_initialize_bond_pair_matrix