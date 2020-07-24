 MODULE module_CM5_parameters
 !===================================================================================
 ! These are the parameters to compute CM5 net atomic charges.
 !===================================================================================
 USE module_precision
 
 IMPLICIT NONE
 
 REAL(kind=dp),PARAMETER,DIMENSION(109) :: covalent_radii=[0.32_dp,0.37_dp,1.30_dp,0.99_dp,0.84_dp,0.75_dp,0.71_dp,0.64_dp,& 
 0.60_dp,0.62_dp,1.60_dp,1.40_dp,1.24_dp,1.14_dp,1.09_dp,1.04_dp,1.00_dp,1.01_dp,2.00_dp,1.74_dp,1.59_dp,1.48_dp,1.44_dp,& 
 1.30_dp,1.29_dp,1.24_dp,1.18_dp,1.17_dp,1.22_dp,1.20_dp,1.23_dp,1.20_dp,1.20_dp,1.18_dp,1.17_dp,1.16_dp,2.15_dp,1.90_dp,&
 1.76_dp,1.64_dp,1.56_dp,1.46_dp,1.38_dp,1.36_dp,1.34_dp,1.30_dp,1.36_dp,1.40_dp,1.42_dp,1.40_dp,1.40_dp,1.37_dp,1.36_dp,&
 1.36_dp,2.38_dp,2.06_dp,1.94_dp,1.84_dp,1.90_dp,1.88_dp,1.86_dp,1.85_dp,1.83_dp,1.82_dp,1.81_dp,1.80_dp,1.79_dp,1.77_dp,&
 1.77_dp,1.78_dp,1.74_dp,1.64_dp,1.58_dp,1.50_dp,1.41_dp,1.36_dp,1.32_dp,1.30_dp,1.30_dp,1.32_dp,1.44_dp,1.45_dp,1.50_dp,&
 1.42_dp,1.48_dp,1.46_dp,2.42_dp,2.11_dp,2.01_dp,1.90_dp,1.84_dp,1.83_dp,1.80_dp,1.80_dp,1.73_dp,1.68_dp,1.68_dp,1.68_dp,&
 1.65_dp,1.67_dp,1.73_dp,1.76_dp,1.61_dp,1.57_dp,1.49_dp,1.43_dp,1.41_dp,1.34_dp,1.29_dp]

 REAL(kind=dp),PARAMETER,DIMENSION(109) :: CM5_DZ=[0.0056_dp,-0.1543_dp,0.0000_dp,0.0333_dp,-0.1030_dp,-0.0446_dp,-0.1072_dp,&
 -0.0802_dp,-0.0629_dp,-0.1088_dp,0.0184_dp,0.0000_dp,-0.0726_dp,-0.0790_dp,-0.0756_dp,-0.0565_dp,-0.0444_dp,-0.0767_dp,&
 0.0130_dp,0.0000_dp,0.0000_dp,0.0000_dp,0.0000_dp,0.0000_dp,0.0000_dp,0.0000_dp,0.0000_dp,0.0000_dp,0.0000_dp,0.0000_dp,&
 -0.0512_dp,-0.0557_dp,-0.0533_dp,-0.0399_dp,-0.0313_dp,-0.0541_dp,0.0092_dp,0.0000_dp,0.0000_dp,0.0000_dp,0.0000_dp,0.0000_dp,&
 0.0000_dp,0.0000_dp,0.0000_dp,0.0000_dp,0.0000_dp,0.0000_dp,-0.0361_dp,-0.0393_dp,-0.0376_dp,-0.0281_dp,-0.0220_dp,-0.0381_dp,&
 0.0065_dp,0.0000_dp,0.0000_dp,0.0000_dp,0.0000_dp,0.0000_dp,0.0000_dp,0.0000_dp,0.0000_dp,0.0000_dp,0.0000_dp,0.0000_dp,&
 0.0000_dp,0.0000_dp,0.0000_dp,0.0000_dp,0.0000_dp,0.0000_dp,0.0000_dp,0.0000_dp,0.0000_dp,0.0000_dp,0.0000_dp,0.0000_dp,&
 0.0000_dp,0.0000_dp,-0.0255_dp,-0.0277_dp,-0.0265_dp,-0.0198_dp,-0.0155_dp,-0.0269_dp,0.0046_dp,0.0000_dp,0.0000_dp,0.0000_dp,&
 0.0000_dp,0.0000_dp,0.0000_dp,0.0000_dp,0.0000_dp,0.0000_dp,0.0000_dp,0.0000_dp,0.0000_dp,0.0000_dp,0.0000_dp,0.0000_dp,&
 0.0000_dp,0.0000_dp,0.0000_dp,0.0000_dp,0.0000_dp,0.0000_dp,0.0000_dp] 

 REAL(kind=dp),PARAMETER :: CM5_alpha = 2.474_dp
 REAL(kind=dp),PARAMETER :: CM5_D_HC = 0.0502_dp
 REAL(kind=dp),PARAMETER :: CM5_D_HN = 0.1747_dp
 REAL(kind=dp),PARAMETER :: CM5_D_HO = 0.1671_dp
 REAL(kind=dp),PARAMETER :: CM5_D_CN = 0.0556_dp
 REAL(kind=dp),PARAMETER :: CM5_D_CO = 0.0234_dp
 REAL(kind=dp),PARAMETER :: CM5_D_NO = -0.0346_dp
 
 END MODULE module_CM5_parameters