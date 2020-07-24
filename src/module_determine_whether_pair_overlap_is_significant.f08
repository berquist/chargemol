 MODULE module_determine_whether_pair_overlap_is_significant
 ! Copyright (c) 2014, 2015, 2016, 2017 by Thomas A. Manz and Nidia Gabaldon Limas. Rights reserved.
 
 USE module_common_variable_declarations
 USE module_global_parameters

 IMPLICIT NONE
 
 CONTAINS
 
 SUBROUTINE determine_whether_pair_overlap_is_significant(distance)
 !=================================================================================== 
 
 REAL(kind=dp) :: temp1,temp2,overlap_upper_bound,estimated_overlap,distance
 INTEGER :: shells_of_separation,h,k,x,index_i,index_j
 
 !Compute upper bound of overlap between the two atoms
 temp1 = 0.0_dp
 overlap_upper_bound = 1.0_dp !initialize to some large value
 shells_of_separation = ceiling(scalefactor*distance + zero_tolerance) !number of radial shells separating the two atoms
 DO k = 1,(shells_of_separation - 1)
   IF ((k > nshells) .or. ((shells_of_separation - k) > nshells)) THEN
     CYCLE
   END IF    
   temp1 = 2.0_dp*(nA_outside(k,i) + nA_outside((shells_of_separation - k),j))
   overlap_upper_bound = min(overlap_upper_bound,temp1)
 END DO
 IF (overlap_upper_bound >= BO_print_cutoff) THEN
   !Approximate the overlap term
   estimated_overlap = 0.0_dp
   IF (shells_of_separation <= 2*nshells) THEN
     DO x = (shells_of_separation - nshells),nshells
       DO h = 0,ceiling(sqrt(real(nshells**2 - max(x**2,(shells_of_separation - x)**2))))
         index_i = nint(sqrt(real(x**2 + h**2))) + 1
         index_j = nint(sqrt(real((shells_of_separation - x)**2 + h**2))) + 1
         IF ((index_i <= nshells) .and. (index_j <= nshells)) THEN
            temp1 = max(min_density(index_i,i),min_density(index_j,j))
            temp1 = max(temp1,(corrected_spherical_avg_density(index_i,i)+corrected_spherical_avg_density(index_j,j)))
            IF (temp1 > zero_tolerance) THEN
               temp2 = 4.0_dp*pi*h*corrected_spherical_avg_density(index_i,i)*corrected_spherical_avg_density(index_j,j)&
	       /(temp1*scalefactor**3)
               estimated_overlap = estimated_overlap + temp2
            END IF 
         END IF    
       END DO
     END DO
   END IF
   IF (estimated_overlap > 0.5_dp*BO_print_cutoff) THEN
     active_pair = 1
   ELSE
     active_pair = 0
   END IF
 ELSE
   active_pair = 0
 END IF

 END SUBROUTINE determine_whether_pair_overlap_is_significant
 
 END MODULE module_determine_whether_pair_overlap_is_significant
 