 MODULE module_compute_CM5
 ! Copyright (c) 2014, 2015, 2016, 2017 by Thomas A. Manz and Nidia Gabaldon Limas. Rights reserved.

 USE module_precision
 USE module_global_parameters
 USE module_common_variable_declarations
 USE module_CM5_parameters

 IMPLICIT NONE
 
 CONTAINS
 
 SUBROUTINE compute_CM5_NACs()
 !===================================================================================
 ! Computes the CM5 net atomic charges
 !===================================================================================
 
 REAL(kind=dp) :: Bij,Tij,bond_cutoff_radius,a_dot_a,a_dot_b,a_dot_c,b_dot_b,b_dot_c,c_dot_c,cos_a_b,cos_a_c,cos_b_c,sin_a_b,&
 sin_a_c,sin_b_c,X1,Y1,Z1,X2,Y2,Z2,distance
 REAL(kind=dp),ALLOCATABLE,DIMENSION(:) :: CM5_net_atomic_charge
 INTEGER :: periodicsumlimitA,periodicsumlimitB,periodicsumlimitC,atomic_number_j,repeata,repeatb,repeatc,active_pair,&
 atomic_number_i
  
 !CM5_parameters
 Bij=0.0_dp
 Tij=0.0_dp
 !Determine the number of unit cell replications for compute CM5 NACs
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
   periodicsumlimitA = ceiling(max((bond_cutoff_radius/sin_a_b), (bond_cutoff_radius/sin_a_c))/sqrt(a_dot_a)) + 1.0_dp
 ELSE
   periodicsumlimitA = 0
 END IF
 IF (periodicB) THEN
   periodicsumlimitB = ceiling(max((bond_cutoff_radius/sin_a_b), (bond_cutoff_radius/sin_b_c))/sqrt(b_dot_b)) + 1.0_dp
 ELSE
   periodicsumlimitB = 0
 END IF
 IF (periodicC) THEN
   periodicsumlimitC = ceiling(max((bond_cutoff_radius/sin_a_c), (bond_cutoff_radius/sin_b_c))/sqrt(c_dot_c)) + 1.0_dp
 ELSE
   periodicsumlimitC = 0
 END IF
 ALLOCATE(CM5_net_atomic_charge(natoms))
 CM5_net_atomic_charge = Hirshfeld_net_atomic_charge
 DO i=1,natoms
   atomic_number_i = atomic_number(i)
   DO j=1,natoms
     atomic_number_j = atomic_number(j)
     DO repeata = -periodicsumlimitA,periodicsumlimitA
       DO repeatb = -periodicsumlimitB,periodicsumlimitB
         DO repeatc = -periodicsumlimitC,periodicsumlimitC
           active_pair = 0
           X1=coords(1,i)
           Y1=coords(2,i)
           Z1=coords(3,i)
           X2=coords(1,j) + repeata*vector1(1) + repeatb*vector2(1) + repeatc*vector3(1)
           Y2=coords(2,j) + repeata*vector1(2) + repeatb*vector2(2) + repeatc*vector3(2)
           Z2=coords(3,j) + repeata*vector1(3) + repeatb*vector2(3) + repeatc*vector3(3)
           distance = sqrt((X1-X2)**2 + (Y1-Y2)**2 + (Z1-Z2)**2)/bohrperangstrom
           Bij = exp(-CM5_alpha*(distance - covalent_radii(atomic_number_i) - covalent_radii(atomic_number_j)))
           IF (atomic_number_i == atomic_number_j) THEN
             Tij = 0.0_dp
           ELSEIF ((atomic_number_i == 1) .and. (atomic_number_j == 6)) THEN
             Tij = CM5_D_HC
           ELSEIF ((atomic_number_i == 6) .and. (atomic_number_j == 1)) THEN
             Tij = -CM5_D_HC
           ELSEIF ((atomic_number_i == 1) .and. (atomic_number_j == 7)) THEN
             Tij = CM5_D_HN
           ELSEIF ((atomic_number_i == 7) .and. (atomic_number_j == 1)) THEN
             Tij = -CM5_D_HN
           ELSEIF ((atomic_number_i == 1) .and. (atomic_number_j == 8)) THEN
             Tij = CM5_D_HO
           ELSEIF ((atomic_number_i == 8) .and. (atomic_number_j == 1)) THEN
             Tij = -CM5_D_HO
           ELSEIF ((atomic_number_i == 6) .and. (atomic_number_j == 7)) THEN
             Tij = CM5_D_CN
           ELSEIF ((atomic_number_i == 7) .and. (atomic_number_j == 6)) THEN
             Tij = -CM5_D_CN
           ELSE IF ((atomic_number_i == 6) .and. (atomic_number_j == 8)) THEN
             Tij = CM5_D_CO
           ELSEIF ((atomic_number_i == 8) .and. (atomic_number_j == 6)) THEN
             Tij = -CM5_D_CO
           ELSEIF ((atomic_number_i == 7) .and. (atomic_number_j == 8)) THEN
             Tij = CM5_D_NO
           ELSEIF ((atomic_number_i == 8) .and. (atomic_number_j == 7)) THEN
             Tij = -CM5_D_NO
           ELSE
             Tij = CM5_DZ(atomic_number_i) - CM5_DZ(atomic_number_j)
           END IF
           CM5_net_atomic_charge(i) = CM5_net_atomic_charge(i) + Tij*Bij
         END DO
       END DO
     END DO
   END DO
 END DO 
 WRITE(output_FID,*) 'The computed CM5 net atomic charges are:'
 WRITE(output_FID,'(10f13.6)') CM5_net_atomic_charge
 
 END SUBROUTINE compute_CM5_NACs
 END MODULE module_compute_CM5