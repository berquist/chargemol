 MODULE module_matrix_operations

 USE module_precision

 IMPLICIT NONE
  
 CONTAINS
  
 PURE FUNCTION adj(a)
 !============================================================
 ! adj function to solve inverse of a 3x3 matrix.
 !------------------------------------------------------------
 ! http://scicomp.stackexchange.com/questions/5306/computing-
 ! the-pseudoinverse-of-a-3x3-matrix
 ! May 2nd 2013 
 !=========================================================== 

 !  Adjugate of a 3x3 matrix (the transpose of the cofactor matrix).
   REAL(kind=dp), dimension(3,3) :: adj
   REAL(kind=dp), dimension(3,3), intent(in) :: a

   adj(1,1) = a(2,2)*a(3,3) - a(2,3)*a(3,2)
   adj(1,2) = a(1,3)*a(3,2) - a(1,2)*a(3,3)
   adj(1,3) = a(1,2)*a(2,3) - a(1,3)*a(2,2)
   adj(2,1) = a(2,3)*a(3,1) - a(2,1)*a(3,3)
   adj(2,2) = a(1,1)*a(3,3) - a(1,3)*a(3,1)
   adj(2,3) = a(1,3)*a(2,1) - a(1,1)*a(2,3)
   adj(3,1) = a(2,1)*a(3,2) - a(2,2)*a(3,1)
   adj(3,2) = a(1,2)*a(3,1) - a(1,1)*a(3,2)
   adj(3,3) = a(1,1)*a(2,2) - a(1,2)*a(2,1)

 END FUNCTION adj
   
 
 PURE FUNCTION inverse(a)
 !============================================================
 ! Solves the inverse of a 3x3 matrix.
 !------------------------------------------------------------
 ! http://scicomp.stackexchange.com/questions/5306/computing-
 ! the-pseudoinverse-of-a-3x3-matrix
 ! May 2nd 2013 
 !=========================================================== 

   REAL(kind=dp), dimension(3,3) :: inverse
   REAL(kind=dp), dimension(3,3), intent(in) :: a
   REAL(kind=dp) :: detr
   detr = 1./determinant(a)
   inverse = adj(a)*detr
   
 END FUNCTION inverse
   

 PURE FUNCTION cross(a,b)
 !============================================================
 ! Cross Product of a 3x3 Matrix for REAL numbers
 !------------------------------------------------------------
 ! http://stackoverflow.com/questions/6511711/computing-the-cr
 ! oss-product-of-two-vectors-in-fortran-90
 ! February 11th 2013 
 !===========================================================
 REAL(kind=dp), DIMENSION(3) :: cross
 REAL(kind=dp), DIMENSION(3), INTENT(IN) :: a, b

 cross(1) = a(2) * b(3) - a(3) * b(2)
 cross(2) = a(3) * b(1) - a(1) * b(3)
 cross(3) = a(1) * b(2) - a(2) * b(1)

 END FUNCTION cross

 PURE FUNCTION determinant(a)
!============================================================
! Matrix Determinant
! Method: Computes the determinant of a 3x3 matrix.
!============================================================
 REAL(kind=dp),INTENT(IN) :: a(3,3)
 REAL(kind=dp) :: determinant
 
 determinant=a(1,1)*a(2,2)*a(3,3)+a(1,2)*a(2,3)*a(3,1)+a(1,3)*a(2,1)*&
 a(3,2)-a(1,3)*a(2,2)*a(3,1)-a(1,2)*a(2,1)*a(3,3)-a(1,1)*&
 a(2,3)*a(3,2)
 
 END FUNCTION determinant
 
 
 PURE FUNCTION factorial(a)
!============================================================
! Integer factorial.
! Method: Computes the factorial of an integer.
!============================================================
 INTEGER :: i,factorial
 INTEGER,INTENT(IN) :: a
 
 factorial=1
 
 IF (a == 0) THEN
   factorial=1
 ELSE IF(a > 0) THEN
   DO i=1,a
    factorial=factorial*i
   END DO
 END IF
 
 END FUNCTION factorial
 
 PURE FUNCTION eig(A)
!============================================================
!This function finds and sorts the eigenvalues of a 3x3 
!matrix.
!http://en.wikipedia.org/wiki/Eigenvalue_algorithm
!============================================================
 
 REAL(kind=dp),INTENT(IN) :: A(3,3)
 REAL(kind=dp) :: I(3,3),p1,q,p2,phi,eig1,eig2,eig3,p,B(3,3),r,eig(3)
 REAL(kind=dp), PARAMETER :: pi = 3.14159265358979323846_dp
 
 !Given a real symmetric 3x3 matrix A, compute the eigenvalues
 p1 = A(1,2)**2 + A(1,3)**2 + A(2,3)**2
 I(1,:) = [1.0_dp,0.0_dp,0.0_dp]
 I(2,:) = [0.0_dp,1.0_dp,0.0_dp]
 I(3,:) = [0.0_dp,0.0_dp,1.0_dp]
 q = (A(1,1)+A(2,2)+A(3,3))/3.0_dp
 IF (p1 < 0.0000000001) THEN
   !A is diagonal.
   eig1 = min(A(1,1),A(2,2),A(3,3))
   eig3 = max(A(1,1),A(2,2),A(3,3))
   eig2 = (q*3.0_dp)-eig1-eig3
 ELSE
   p2 = (A(1,1)-q)**2 + (A(2,2)-q)**2 + (A(3,3)-q)**2+2.0_dp*p1
   p = sqrt(p2/6.0_dp)
   B = (1.0_dp/p) * (A-q*I) 
   r = determinant(B)/2.0_dp
   !In exact arithmetic for a symmetric matrix  -1 <= r <= 1
   !but computation error can leave it slightly outside this range.
   IF (r <= -1.0_dp) THEN
      phi = pi/3.0_dp
   ELSE IF (r >= 1.0_dp) THEN
      phi = 0.0_dp
   ELSE
      phi = acos(r)/3.0_dp
   END IF
   !The eigenvalues satisfy eig1 <= eig2 <= eig3
   eig3 = q+2.0_dp*p*cos(phi)
   eig1 = q+2.0_dp*p*cos(phi+(2.0_dp*pi/3.0_dp))
   eig2 = 3.0_dp*q-eig1-eig3     !since trace(A) = eig1 + eig2 + eig3
 END IF
 eig = [eig1,eig2,eig3] 
 
 END FUNCTION eig
 
 END MODULE module_matrix_operations