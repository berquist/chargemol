 MODULE module_generate_kApoints_kBpoints_kCpoints
 ! Copyright (c) 2014, 2015, 2016, 2017 by Thomas A. Manz and Nidia Gabaldon Limas. Rights reserved.
 
 USE module_common_variable_declarations

 IMPLICIT NONE

 CONTAINS
 
 SUBROUTINE generate_kApoints_kBpoints_kCpoints()
 !===================================================================================  

 ! Set up the kApoints, kBpoints, and kCpoints arrays
 ALLOCATE(lower_na(natoms))
 ALLOCATE(upper_na(natoms))
 ALLOCATE(lower_nb(natoms))
 ALLOCATE(upper_nb(natoms))
 ALLOCATE(lower_nc(natoms))
 ALLOCATE(upper_nc(natoms))
 ALLOCATE(kApoints(natoms,(2*delta_na + 1)))
 ALLOCATE(kBpoints(natoms,(2*delta_nb + 1)))
 ALLOCATE(kCpoints(natoms,(2*delta_nc + 1)))
 lower_na = 0
 upper_na = 0
 lower_nb = 0
 upper_nb = 0
 lower_nc = 0
 upper_nc = 0
 kApoints = 0
 kBpoints = 0
 kCpoints = 0
 DO j=1,natoms
   lower_na(j) = delta_na
   upper_na(j) = -delta_na
   DO na = -delta_na,delta_na
     IF (periodicA) THEN
       kApoints(j,delta_na + na + 1) = modulo((na + center_nabc(1,j)),totnumA) + 1
     ELSE
       kApoints(j,delta_na + na + 1) = na + center_nabc(1,j) + 1
     END IF
     ka = kApoints(j,delta_na + na + 1)
     IF ((ka < 1) .or. (ka > totnumA)) THEN
       CYCLE
     ELSE
       lower_na(j) = min(lower_na(j),na)
       upper_na(j) = max(upper_na(j),na)
     END IF  
   END DO
   lower_nb(j) = delta_nb
   upper_nb(j) = -delta_nb
   DO nb = -delta_nb,delta_nb
     IF (periodicB) THEN
       kBpoints(j,delta_nb + nb + 1) = modulo((nb + center_nabc(2,j)),totnumB) + 1
     ELSE
       kBpoints(j,delta_nb + nb + 1) = nb + center_nabc(2,j) + 1
     END IF
     kb = kBpoints(j,delta_nb + nb + 1)
     IF ((kb < 1) .or. (kb > totnumB)) THEN
       CYCLE
     ELSE
       lower_nb(j) = min(lower_nb(j),nb)
       upper_nb(j) = max(upper_nb(j),nb)
     END IF 
   END DO
   lower_nc(j) = delta_nc
   upper_nc(j) = -delta_nc
   DO nc = -delta_nc,delta_nc
     IF (periodicC) THEN
       kCpoints(j,delta_nc + nc + 1) = modulo((nc + center_nabc(3,j)),totnumC) + 1
     ELSE
       kCpoints(j,delta_nc + nc + 1) = nc + center_nabc(3,j) + 1
     END IF
     kc = kCpoints(j,delta_nc + nc + 1)
     IF ((kc < 1) .or. (kc > totnumC)) THEN
       CYCLE
     ELSE
       lower_nc(j) = min(lower_nc(j),nc)
       upper_nc(j) = max(upper_nc(j),nc)
     END IF 
   END DO
   !Check to see whether the atom has any valid grid points
   IF ((lower_na(j) > upper_na(j)) .or. (lower_nb(j) > upper_nb(j)) .or. (lower_nc(j) > upper_nc(j))) THEN
     WRITE(output_FID,*)'One of the atoms does not have any integration grid points.'
     WRITE(output_FID,*)'Program will terminate.'
     FLUSH(output_FID)
     STOP
   END IF
 END DO
 
 END SUBROUTINE generate_kApoints_kBpoints_kCpoints
 
 END MODULE module_generate_kApoints_kBpoints_kCpoints