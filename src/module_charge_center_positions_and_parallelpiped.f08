 MODULE module_charge_center_positions_and_parallelpiped
 !=====================================================================================
 ! Module used in run_valence_core to calculate the denstity grids from the gaussian 
 ! basis set coefficients
 ! Copyright (c) 2014, 2015, 2016, 2017 by Thomas A. Manz and Nidia Gabaldon Limas. Rights reserved.
 !=====================================================================================
  
 USE module_precision
 USE module_global_parameters
 USE module_common_variable_declarations
 USE module_matrix_operations
!$ USE omp_lib
  
 IMPLICIT NONE
 
 CONTAINS
 
 SUBROUTINE charge_center_positions()
 !====================================================================================
 ! Compute the positions of the charge centers in terms of na, nb, nc
 !====================================================================================

 REAL(kind=dp) :: xyz(3),T(3,3),temp_vector(3)
 
 i = 0
 j = 0
 T=transpose(boundary)
 inv_boundary=inverse(T)
 DO i = 1,natoms
   DO j = 1,3
       xyz(j)=coords(j,i)
   END DO  
   temp_vector = matmul(inv_boundary,(xyz - origin))
   DO j=1,3
    center_nabc(j,i) = nint(temp_vector(j))
   END DO
   center_shift(:,i) = matmul(T,(temp_vector(:) - nint(temp_vector(:))))
 END DO

 WRITE(output_FID,*)'center_nabc ='
 FLUSH(output_FID)
 DO i=1,natoms
   WRITE(output_FID, *) (center_nabc(j,i),j=1,3)
 END DO
 FLUSH(output_FID)

 END SUBROUTINE charge_center_positions
 

 
 SUBROUTINE parallelpiped()
 !====================================================================================
 ! Compute the three side lengths of a parallelepiped that enclosed each atom, where 
 ! each side of the parallelpiped lies along the chosen axes.
 !====================================================================================
 
 REAL(kind=dp) :: temp_vector1(3), temp_vector2(3), temp_vector3(3),a_dot_a,a_dot_b,&
 a_dot_c,b_dot_b, b_dot_c,c_dot_c,R,cos_a_b,cos_a_c,cos_b_c,sin_a_b,sin_a_c,sin_b_c

 temp_vector1(:) = boundary(1,:)
 temp_vector2(:) = boundary(2,:)
 temp_vector3(:) = boundary(3,:)
 
 a_dot_a = dot_product(temp_vector1,temp_vector1)
 a_dot_b = dot_product(temp_vector1,temp_vector2)
 a_dot_c = dot_product(temp_vector1,temp_vector3)
 b_dot_b = dot_product(temp_vector2,temp_vector2)
 b_dot_c = dot_product(temp_vector2,temp_vector3)
 c_dot_c = dot_product(temp_vector3,temp_vector3)
 R = nshells/scalefactor
 cos_a_b = (a_dot_b)/sqrt(a_dot_a*b_dot_b)
 cos_a_c = (a_dot_c)/sqrt(a_dot_a*c_dot_c)
 cos_b_c = (b_dot_c)/sqrt(b_dot_b*c_dot_c)
 sin_a_b = sqrt(1.0_dp - cos_a_b*cos_a_b)
 sin_a_c = sqrt(1.0_dp - cos_a_c*cos_a_c)
 sin_b_c = sqrt(1.0_dp - cos_b_c*cos_b_c)
 delta_na = ceiling(max((R/sin_a_b), (R/sin_a_c))/sqrt(a_dot_a)) + 1.0_dp
 delta_nb = ceiling(max((R/sin_a_b), (R/sin_b_c))/sqrt(b_dot_b)) + 1.0_dp
 delta_nc = ceiling(max((R/sin_a_c), (R/sin_b_c))/sqrt(c_dot_c)) + 1.0_dp
!$ number_of_threads = omp_get_max_threads()
!$ chunk_size = ceiling((2*delta_nc+1)/real(number_of_threads*6)) 
 WRITE(output_FID,'(a,i6)')' delta_na= ',delta_na
 FLUSH(output_FID)
 WRITE(output_FID,'(a,i6)')' delta_nb= ',delta_nb
 FLUSH(output_FID)
 WRITE(output_FID,'(a,i6)')' delta_nc= ', delta_nc
 FLUSH(output_FID)
 
 END SUBROUTINE parallelpiped
 
 END MODULE module_charge_center_positions_and_parallelpiped