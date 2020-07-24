 MODULE module_gaussian_functions 
 ! Copyright (c) 2014, 2015, 2016, 2017 by Thomas A. Manz and Nidia Gabaldon Limas. Rights reserved.
 
 USE module_precision
 USE module_matrix_operations
 USE module_global_parameters

 IMPLICIT NONE
 
 CONTAINS
 
 PURE FUNCTION gaussian_overlap_integral(alpha1,Lx1,Ly1,Lz1,X1,Y1,Z1,alpha2,Lx2,Ly2,Lz2,X2,Y2,Z2)
 !==========================================================================================
 ! Computes the spatial overlap integral between two gaussians f1 and f2
 ! f1 =(X-X1)^Lx1*(Y-Y1)^Ly1*(Z-Z1)^Lz1*exp(-alpha1*((X-X1)^2 +(Y-Y1)^2 +(Z-Z1)^2))
 ! f2 =(X-X2)^Lx2*(Y-Y2)^Ly2*(Z-Z2)^Lz2*exp(-alpha2*((X-X2)^2 +(Y-Y2)^2 +(Z-Z2)^2)) 
 !==========================================================================================
 
 REAL(kind=dp) :: factor1,factor2,XD,YD,ZD,temp,SART,factor3,factor4,gaussian_overlap_integral
 REAL(kind=dp),INTENT(IN) :: alpha1,alpha2,X2,X1,Y2,Y1,Z2,Z1
 INTEGER :: i,j
 INTEGER,INTENT(IN) :: Lx1,Lx2,Ly1,Ly2,Lz1,Lz2

 factor1=exp(-alpha1*alpha2/(alpha1 + alpha2)*((X2-X1)**2 + (Y2-Y1)**2 + (Z2-Z1)**2))
 XD=(alpha1*X1 + alpha2*X2)/(alpha1 + alpha2)
 YD=(alpha1*Y1 + alpha2*Y2)/(alpha1 + alpha2)
 ZD=(alpha1*Z1 + alpha2*Z2)/(alpha1 + alpha2)
 factor2=0.0_dp
 DO i=0,Lx1
   DO j=0,Lx2
     IF (modulo((i+j),2) == 1) THEN 
       CYCLE
     ELSE
       temp=(factorial(Lx1)/(factorial(Lx1-i)*factorial(i)))
       temp=temp*(factorial(Lx2)/(factorial(Lx2-j)*factorial(j)))
       temp=temp*((XD - X1)**(Lx1-i))*((XD - X2)**(Lx2-j))
       SART = gaussian_integration((alpha1+alpha2),(i+j))
     END IF
       factor2 = factor2 + temp*SART
   END DO
 END DO
 factor3=0.0_dp
 DO i=0,Ly1
   DO j=0,Ly2
     IF (modulo((i+j),2) == 1) THEN
       CYCLE
     ELSE
       temp=(factorial(Ly1)/(factorial(Ly1-i)*factorial(i)))
       temp=temp*(factorial(Ly2)/(factorial(Ly2-j)*factorial(j)))
       temp=temp*((YD - Y1)**(Ly1-i))*((YD - Y2)**(Ly2-j))
       SART = gaussian_integration((alpha1+alpha2),(i+j))
     END IF
     factor3 = factor3 + temp*SART
   END DO
 END DO
 factor4=0.0_dp
 DO i=0,Lz1
    DO j=0,Lz2
        IF (modulo((i+j),2) == 1) THEN
           CYCLE
        ELSE
           temp=(factorial(Lz1)/(factorial(Lz1-i)*factorial(i)))
           temp=temp*(factorial(Lz2)/(factorial(Lz2-j)*factorial(j)))
           temp=temp*((ZD - Z1)**(Lz1-i))*((ZD - Z2)**(Lz2-j))
           SART = gaussian_integration((alpha1+alpha2),(i+j))
        END IF
        factor4 = factor4 + temp*SART
    END DO
 END DO
 gaussian_overlap_integral=factor1*factor2*factor3*factor4

 END FUNCTION gaussian_overlap_integral


 PURE FUNCTION gaussian_integration(alpha,Lx_fun)
 !==========================================================================================

 INTEGER :: i_fun,temp_fun
 INTEGER,INTENT(IN) :: Lx_fun
 REAL(kind=dp) :: gaussian_integration
 REAL(kind=dp),INTENT(IN) :: alpha

 temp_fun=1
 IF (modulo(Lx_fun,2) == 1) THEN
   gaussian_integration = 0.0_dp
 ELSE
   DO i_fun=1,Lx_fun,2
     temp_fun=temp_fun*i_fun
   END DO
   gaussian_integration =sqrt(pi/alpha)*temp_fun/((2*alpha)**((Lx_fun)/2.0_dp))
 END IF
 
 END FUNCTION gaussian_integration

 
 PURE FUNCTION gaussian_value(constant,alpha,Lx,Ly,Lz,Ax,Ay,Az,x,y,z)
 !==========================================================================================
 
 REAL(kind=dp) :: gaussian_value
 REAL(kind=dp),INTENT(IN) :: constant,alpha,Ax,Ay,Az,x,y,z
 INTEGER,INTENT(IN) :: Lx,Ly,Lz

 gaussian_value = constant*((x-Ax)**Lx)*((y-Ay)**Ly)*((z-Az)**Lz)*exp(-alpha*((x-Ax)**2 + (y-Ay)**2 + (z-Az)**2))
  
 END FUNCTION gaussian_value


 !---------------------------------------------------------------------------------------------------------------------
 PURE FUNCTION spherical_multiplier(Lx_sum,Ly_sum,Lz_sum) RESULT(spherical_avg_multiplier)
 IMPLICIT NONE
 
 INTEGER,INTENT(IN) :: Lx_sum,Ly_sum,Lz_sum
 INTEGER :: L_sum,L_plus_product
 REAL(kind=dp) :: spherical_avg_multiplier
 
 L_sum = Lx_sum + Ly_sum + Lz_sum
 L_plus_product = (Lx_sum+1)*(Ly_sum+1)*(Lz_sum+1)
 IF ((Lx_sum < 0) .or. (Ly_sum < 0) .or. (Lz_sum < 0)) THEN
   spherical_avg_multiplier = 0.0_dp  
 ELSE IF((modulo(Lx_sum,2) == 1) .or. (modulo(Ly_sum,2)==1) .or. (modulo(Lz_sum,2)==1)) THEN
   spherical_avg_multiplier = 0.0_dp
 ELSE IF (L_sum == 0) THEN
   spherical_avg_multiplier = 1.0_dp
 ELSE IF (L_sum == 2) THEN
   spherical_avg_multiplier = 1.0_dp/3.0_dp
 ELSE IF (L_sum > 8) THEN
   spherical_avg_multiplier = 0.0_dp  
 ELSE IF (L_plus_product == 5) THEN
   spherical_avg_multiplier = 1.0_dp/5.0_dp
 ELSE IF ((L_sum == 4) .AND. (L_plus_product == 9)) THEN
   spherical_avg_multiplier = 1.0_dp/15.0_dp
 ELSE IF (L_plus_product == 7) THEN
   spherical_avg_multiplier = 1.0_dp/7.0_dp
 ELSE IF (L_plus_product == 15) THEN
   spherical_avg_multiplier = 1.0_dp/35.0_dp
 ELSE IF (L_plus_product == 27) THEN
   spherical_avg_multiplier = 1.0_dp/105.0_dp
 ELSE IF ((L_sum == 8) .AND. (L_plus_product == 9)) THEN
   spherical_avg_multiplier = 1.0_dp/9.0_dp
 ELSE IF (L_plus_product == 21) THEN
   spherical_avg_multiplier = 1.0_dp/63.0_dp
 ELSE IF (L_plus_product == 25) THEN
   spherical_avg_multiplier = 1.0_dp/105.0_dp
 ELSE IF (L_plus_product == 45) THEN
   spherical_avg_multiplier = 1.0_dp/315.0_dp 
 END IF
 
 END FUNCTION spherical_multiplier

  
 END MODULE module_gaussian_functions
 