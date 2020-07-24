 MODULE module_cloud_penetration
 ! Copyright (c) 2014, 2015, 2016, 2017 by Thomas A. Manz and Nidia Gabaldon Limas. Rights reserved.
 
 USE module_precision
 USE module_global_parameters
 USE module_common_variable_declarations

 IMPLICIT NONE

 CONTAINS

 SUBROUTINE cloud_penetration()
 !===================================================================================
 ! Calculation of electron cloud penetration terms.
 !===================================================================================
 
 INTEGER :: index_min_penetration,k,n
 REAL(kind=dp) :: Sx,Sxx,Sxy,Sy,Syy,x,y
 
 index_min_penetration = ceiling(rmin_cloud_penetration*nshells/cutoff_radius)
 WRITE(output_FID,*) 'Calculation of the electron cloud penetration terms'
 ALLOCATE(fitted_tail_slope(natoms))
 ALLOCATE(fitted_tail_intercept(natoms))
 ALLOCATE(fitted_tail_Rsquared(natoms))
 DO j = 1,natoms
   Sx = 0.0_dp
   Sxx = 0.0_dp
   Sxy = 0.0_dp
   Sy = 0.0_dp
   Syy = 0.0_dp
   n = 0
   DO k = index_min_penetration,nshells
     x = (k-0.5_dp)*(cutoff_radius*bohrperangstrom/(100.0_dp*nshells))
     IF ((spherical_average_density(k,j) < zero_tolerance**1.5_dp)) THEN
       CYCLE
     END IF
     y = log(spherical_average_density(k,j))
     Sx = Sx + x
     Sxx = Sxx + x*x
     Sxy = Sxy + x*y
     Sy = Sy + y
     Syy = Syy + y*y
     n = n + 1
   END DO
    fitted_tail_slope(j) = (n*Sxy - Sx*Sy)/(n*Sxx - Sx*Sx)
    fitted_tail_intercept(j) = (Sy - fitted_tail_slope(j)*Sx)/n
    fitted_tail_Rsquared(j) = (n*Sxy - Sx*Sy)**2/((n*Sxx - Sx*Sx)*(n*Syy - Sy*Sy))
 END DO   
 
 END SUBROUTINE cloud_penetration

 END MODULE module_cloud_penetration