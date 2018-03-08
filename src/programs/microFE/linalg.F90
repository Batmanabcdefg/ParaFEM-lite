MODULE linalg

  USE precision

CONTAINS
  
  FUNCTION DOUBLE_CONTRACTION_22(vector_input_1,vector_input_2)
    ! This function performs a double contraction of two second order tensors, 
    ! in vector notation, as following, alpha=AijCij
    IMPLICIT NONE
    
    real(iwp), INTENT(IN) :: vector_input_1(:), vector_input_2(:)
    real(iwp) :: double_contraction_22
    integer :: i
    
    double_contraction_22=0._iwp
    
    DO i=1,3
      double_contraction_22=double_contraction_22+vector_input_1(i)* &
       vector_input_2(i)
    END DO
    
    DO i=4,6
      double_contraction_22=double_contraction_22+2._iwp*vector_input_1(i)* &
       vector_input_2(i)
    END DO
    
  RETURN
  END FUNCTION DOUBLE_CONTRACTION_22
  
  FUNCTION TENSOR_PRODUCT_22(vector_input_1,vector_input_2)
    ! This function calculates a fourth order tensor through the tensorial 
    ! product of two second order tensor, in matrix notation, as following
    ! Aijkl=BijCkl
    IMPLICIT NONE
    
    real(iwp), INTENT(IN) :: vector_input_1(:), vector_input_2(:)
    real(iwp) :: tensor_product_22(6,6)
    integer :: i, j
    
    DO i=1,6
      DO j=1,6
        tensor_product_22(i,j)=vector_input_1(i)*vector_input_2(j)
      END DO
    END DO
    
  RETURN
  END FUNCTION TENSOR_PRODUCT_22
  
  FUNCTION DOT_PROD(vector_input_1,vector_input_2,dimen)
    ! This function performs a double contraction of two second order tensors, 
    ! in vector notation, as following, dot_prod=uivi, with dimen being the 
    ! dimension of both vectors
    IMPLICIT NONE
    
    real(iwp), INTENT(IN) :: vector_input_1(:), vector_input_2(:)
    integer, INTENT(IN) :: dimen
    real(iwp) :: dot_prod
    integer :: i
    
    dot_prod=0._iwp
    
    DO i=1,dimen
      dot_prod=dot_prod+vector_input_1(i)*vector_input_2(i)
    END DO
       
  RETURN
  END FUNCTION DOT_PROD
  
  SUBROUTINE INVERSE(a,c,n)
    !This subroutine calculates the inverse of a nxn matrix
    ! Input is a(n,n)
    ! n is the dimension
    ! c is the inverse of a    
    IMPLICIT NONE

    INTEGER :: n, i, j, k  
    REAL(iwp) :: a(n,n), c(n,n), l(n,n), u(n,n), b(n), d(n), x(n)
    REAL(iwp) :: coeff
    
    ! step 0: initialization for matrices L and U and b
    ! Fortran 90/95 allows such operations on matrices
    l=0._iwp
    u=0._iwp
    b=0._iwp

    ! step 1: forward elimination
    DO k=1,(n-1)
      DO i=(k+1),n
        coeff=a(i,k)/a(k,k)
        l(i,k)=coeff
        DO j=(k+1),n
          a(i,j)=a(i,j)-coeff*a(k,j)
        END DO
      END DO
    END DO

    ! Step 2: prepare L and U matrices 
    ! L matrix is a matrix of the elimination coefficient
    ! + the diagonal elements are 1.0
    DO i=1,n
      l(i,i)=1._iwp
    END DO
    
    ! U matrix is the upper triangular part of A
    DO j=1,n
      DO i=1,j
        u(i,j)=a(i,j)
      END DO
    END DO

    ! Step 3: compute columns of the inverse matrix C
    DO k=1,n
      b(k)=1._iwp
      d(1)=b(1)
      
      ! Step 3a: Solve Ld=b using the forward substitution
      DO i=2,n
        d(i)=b(i)
        DO j=1,(i-1)
          d(i)=d(i)-l(i,j)*d(j)
        END DO
      END DO
      
      ! Step 3b: Solve Ux=d using the back substitution
      x(n)=d(n)/u(n,n)
      DO i=(n-1),1,-1
        x(i)=d(i)
        DO j=n,(i+1),-1
        x(i)=x(i)-u(i,j)*x(j)
        END DO
        x(i)=x(i)/u(i,i)
      END DO
      
      ! Step 3c: fill the solutions x(n) into column k of C
      DO i=1,n
        c(i,k)=x(i)
      END DO
      b(k)=0._iwp
    END DO
    
  RETURN
  END SUBROUTINE INVERSE  

END MODULE linalg
