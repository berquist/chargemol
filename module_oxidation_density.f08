 MODULE module_oxidation_density
 ! Copyright (c) 2014, 2015, 2016, 2017 by Thomas A. Manz and Nidia Gabaldon Limas. Rights reserved.
 
 USE module_precision
 USE module_string_utilities

 IMPLICIT NONE

 CONTAINS

 FUNCTION oxidation_density(output_FID,atomic_densities_directory,density_set_prefix,Atomic_Num,atom_charge,cutoff_radius,nshells&
 ,radial_shell_volume)
 !=================================================================================== 

 CHARACTER(200) :: atomic_densities_directory
 CHARACTER(2) :: density_set_prefix
 INTEGER :: atom_charge,sign_charge,q,nelectron,Atomic_Num,atomic_fid,io,i,refine_step,k,openstat,nshells,output_FID
 REAL(kind=dp),DIMENSION(nshells) :: radial_shell_volume
 REAL(kind=dp), ALLOCATABLE, DIMENSION (:) :: temp_oxidation_density,temp,oxidation_density
 CHARACTER(200) :: combinedstring,line_of_text
 REAL(kind=dp) :: summed_oxidation_density,summed_delta_oxidation_density,temp1,delta_rho_min,density_shift
 REAL(kind=dp) :: cutoff_radius
 ALLOCATE(temp_oxidation_density(nshells))
 ALLOCATE(temp(nshells))
 ALLOCATE(oxidation_density(nshells))
 temp_oxidation_density=0.0_dp
 temp=0.0_dp
 IF (atom_charge < 0) THEN
    sign_charge = -1
 ELSE
    sign_charge = 1
 END IF    
 DO q = 0,atom_charge,sign_charge
   !Construct the file name to read
   nelectron=Atomic_Num - q
   WRITE(combinedstring, '(A,A,I3.3,A,I3.3,A,I3.3,A,I3.3,A,I3.3,A)') density_set_prefix,'_',Atomic_Num,'_',Atomic_Num,'_',&
   &nelectron,'_',nint(cutoff_radius),'_',nshells,'.txt'
   OPEN(NEWUNIT=atomic_fid, FILE=(TRIM(ADJUSTL(atomic_densities_directory)))//combinedstring, STATUS='OLD',IOSTAT=openstat)
   IF (openstat > 0) THEN  
     WRITE(output_FID,*) combinedstring
     FLUSH(output_FID)
     WRITE(output_FID,*) 'Could not find a suitable reference density. Program will terminate.'
     FLUSH(output_FID)
     STOP
   END IF
   !Read in the density
   DO
     READ(atomic_fid,'(a)',IOSTAT=io) line_of_text
     IF (io > 0) THEN
       WRITE(output_FID,*) 'Could not find a suitable reference density for ',combinedstring,'. Program will terminate.'
       FLUSH(output_FID)
       STOP
     ELSE IF (io < 0) THEN
       EXIT
     ELSE
       !adjusting left, lowercase, trim data and read
       line_of_text = SimpleLine(line_of_text)
       line_of_text = StrLowCase(line_of_text)
       !assinging a value to each item
       IF(line_of_text == 'density in atomic units:')THEN
         READ(atomic_fid,*)(temp(i),i=1,nshells)  
       END IF
     END IF
   END DO
   CLOSE(atomic_fid)
   !Compute the normalization factor for the reference density
   summed_oxidation_density = 0.0_dp
   summed_delta_oxidation_density = 0.0_dp
   DO k = 2,nshells
     temp(k) = max(temp(k),1.0e-16_dp)
     summed_oxidation_density = summed_oxidation_density + radial_shell_volume(k)*temp(k)
   END DO 
   !Correct the reference density cusp at the nucleus
   temp(1) = max( (nelectron - summed_oxidation_density)/radial_shell_volume(1) , 1.0e-16_dp)
   summed_oxidation_density = summed_oxidation_density + radial_shell_volume(1)*temp(1)
   !Make the oxidation density monotonically decreasing
   temp1 = 1.0e-16_dp
   summed_oxidation_density = 0.0_dp
   DO k=nshells,1,-1
     IF (temp(k) > temp1) THEN
       temp1 = temp(k)
     ELSE
       temp(k) = temp1
     END IF
     summed_oxidation_density = summed_oxidation_density + radial_shell_volume(k)*temp(k)
   END DO
   !Make sure the normalization is exact
   DO k = 1,nshells
     temp(k) = max(temp(k)*nelectron/summed_oxidation_density,1.0e-16_dp)
   END DO
   summed_oxidation_density = 0.0_dp
   DO k=nshells,1,-1
     summed_oxidation_density = summed_oxidation_density + radial_shell_volume(k)*temp(k)
   END DO
   !Make the density difference monotonically increasing
   IF (q /= 0) THEN
     DO refine_step = 1,100
       summed_delta_oxidation_density = 0.0_dp
       delta_rho_min = 0.0_dp
       DO k = nshells,1,-1
         delta_rho_min = max(delta_rho_min,sign_charge*(temp_oxidation_density(k) - temp(k)))
         density_shift = temp_oxidation_density(k) - sign_charge*delta_rho_min - temp(k)
         temp(k) = temp(k) + density_shift
         summed_delta_oxidation_density = summed_delta_oxidation_density + density_shift*radial_shell_volume(k)
       END DO
       DO k = 1,nshells
         temp(k) = temp(k)*summed_oxidation_density/(summed_oxidation_density + summed_delta_oxidation_density)
       END DO
       IF (abs(summed_delta_oxidation_density) < 1.0e-16_dp) THEN
         EXIT
       END IF   
     END DO    
   END IF
   temp_oxidation_density = temp
 END DO
 oxidation_density = temp_oxidation_density

 END FUNCTION oxidation_density
 
 END MODULE module_oxidation_density