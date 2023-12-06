program number_conditionaly
   use Environment
   
   implicit none
   character(*), parameter    :: output_file = "output.txt"
   integer                    :: Out = 0, i = 0, j = 0, N = 8, iP = 0
   integer, allocatable       :: IPVT(:)
   real(R_), allocatable      :: Matrix(:,:), UnitMatrix(:,:), InverseMatrix(:,:), WORK(:), ResidualMatrix(:,:)
   real(R_)                   :: P(5), COND, DET, NORM, tmp

   P = (/1.0, 0.1, 0.01, 0.0001, 0.000001/)

   allocate(Matrix(N, N), UnitMatrix(N, N), InverseMatrix(N,N), ResidualMatrix(N,N))
   allocate(WORK(N))
   allocate(IPVT(N)) 

   DO i = 1, N
      DO j = 1, N
         IF (i == j) THEN
            UnitMatrix(i,j) = 1
         ELSE
            UnitMatrix(i,j) = 0
         END IF
      END DO
   END DO

   open (file=output_file, encoding=E_, newunit=Out)
      write (Out, *) "2 практическое задание"
   close (Out)

   DO iP = 1, 5
      
      call Init_matrix(Matrix, P(iP))
      
      open (file=output_file, encoding=E_, newunit=Out, position='append')
         write(Out, *)
         write(Out, '(a, T20, "= ", e10.4)') "P: ", P(iP)
         write(Out, "(a)") "Исходная матрица:"

         DO i = 1, N
            write(Out, '('//N//'(f12.6))') (Matrix(i, j), j = 1, N)
         END DO
      close (Out)

      call DECOMP(N, N, Matrix, COND, IPVT, WORK)

      !Вычисление определителя матрицы Matrix
      DET = IPVT(N)
      DO i = 1, N
         DET = DET*(Matrix(i, i))
      END DO 

      call Inverse_matrix(Matrix, InverseMatrix, UnitMatrix, N, IPVT)
      call Init_matrix(Matrix, P(iP))
      call Residual_matrix(Matrix, InverseMatrix, UnitMatrix, ResidualMatrix, N)
      
      open (file=output_file, encoding=E_, newunit=Out, position='append')
         write(Out, *)
         write(Out, '(a, T60, "= ", f18.6)') "Определитель исходной матрицы", DET
         write(Out, "(a)") "Обратная матрица:"

         DO i = 1, N
            write(Out, '('//N//'(f18.6))') (InverseMatrix(i, j), j = 1, N)
         END DO
      
         write(Out, *)
         write(Out, "(a)") "Матрица невязки:"

         DO i = 1, N
            write(Out, '('//N//'(e12.4))') (ResidualMatrix(i, j), j = 1, N)
         END DO
      close (Out)

      open (file=output_file, encoding=E_, newunit=Out, position='append')
         write(Out, '(a, T46, "= ", f18.6)') "COND", COND
      close(Out)

      !Вычисление нормы матрицы невязки
      NORM=0.0
      DO j=1,N
         tmp=0.0
         DO i=1,N
            tmp=tmp+ABS(ResidualMatrix(i,j))
         END DO
         IF(tmp > NORM) NORM=tmp
      END DO

      open (file=output_file, encoding=E_, newunit=Out, position='append')
         write(Out, '(a, T60, "= ", e13.6)') "NORM матрицы невязки", NORM
      close(Out)
   END DO
contains
   !Заполнение матрицы началыными значениями
   !На входе матрица - Matrix и p для первого элемента матрицы
   !На выходе Matrix заполненная начальными значениями. Первый элемент - p + 6
   SUBROUTINE Init_matrix(Matrix, p)
      real(R_), intent(in)                     :: p 
      real(R_), dimension(8,8), intent(out)    :: Matrix           

      Matrix(1, 1:8) = (/0, 2, 6, 8, -2, 1, 8, -5/)
      Matrix(2, 1:8) = (/6, -22, -2, -1, 0, 5, -6, 4/)
      Matrix(3, 1:8) = (/-2, -3, -16, 0, 0, -4, 2, -5/)
      Matrix(4, 1:8) = (/1, 1, 4, 9, 1,  0, 0, -6/)
      Matrix(5, 1:8) = (/0, 2, 0, 2, -3, -5, 7, 5/)
      Matrix(6, 1:8) = (/6, -2, -4, 2, -8, -12, 3, -3/)
      Matrix(7, 1:8) = (/-6, -6, 0, -8, 0, 5, -15, 0/)
      Matrix(8, 1:8) = (/0, 7, 6, 0, -5, -8, -5, -3/)

      Matrix(1, 1) = p + 6
   END SUBROUTINE Init_matrix

   !Получение обратной матрицы
   SUBROUTINE Inverse_matrix(Matrix, InverseMatrix, UnitMatrix, N, IPVT)
      integer, intent(in)                    :: N
      real(R_), dimension(N, N), intent(in)  :: Matrix, UnitMatrix
      real(R_), dimension(N, N), intent(out) :: InverseMatrix
      integer, dimension(N), intent(in)      :: IPVT

      real(R_), dimension(N)                 :: WORK_tmp
      integer                                :: i

      DO i = 1, N
         WORK_tmp = UnitMatrix(:, i)
         call SOLVE(N, N, Matrix, WORK_tmp, IPVT)
         InverseMatrix(:, i) = WORK_tmp
      END DO
   END SUBROUTINE Inverse_matrix

   !Получение матрицы невязки
   SUBROUTINE Residual_matrix(Matrix, InverseMatrix, UnitMatrix, ResidualMatrix, N)
      integer, intent(in)                    :: N
      real(R_), dimension(N, N), intent(in)  :: Matrix, InverseMatrix, UnitMatrix
      real(R_), dimension(N, N), intent(out) :: ResidualMatrix

      integer                                :: i, j

      ResidualMatrix = matmul(Matrix, InverseMatrix)
      DO i = 1, N
         DO j = 1, N
            ResidualMatrix(i, j) = UnitMatrix(i, j) - ResidualMatrix(i, j)       
         END DO
      END DO
   END SUBROUTINE Residual_matrix

   SUBROUTINE DECOMP(NDIM,N,A,COND,IPVT,WORK)
         INTEGER NDIM,N
         REAL(R_) A(NDIM,N),COND,WORK(N)
         INTEGER IPVT(N)
         REAL(R_) EK,T,ANORM,YNORM,ZNORM
         INTEGER NM1,I,J,K,KP1,KB,KM1,M
   
         IPVT(N)=1
         IF(N.EQ.1)GO TO 80
         NM1=N-1
   
         ANORM=0.0
         DO 10 J=1,N
            T=0.0
            DO 5 I=1,N
               T=T+ABS(A(I,J))
         5   CONTINUE
            IF(T.GT.ANORM) ANORM=T
      10 CONTINUE
   
         DO 35 K=1,NM1
            KP1=K+1
   
            M=K
            DO 15 I=KP1,N
               IF(ABS(A(I,K)).GT.ABS(A(M,K))) M=I
      15   CONTINUE
            IPVT(K)=M
            IF(M.NE.K)IPVT(N)=-IPVT(N)
            T=A(M,K)
            A(M,K)=A(K,K)
            A(K,K)=T
   
            IF(T.EQ.0.0)GO TO 35
   
            DO 20 I=KP1,N
               A(I,K)=-A(I,K)/T
      20   CONTINUE
   
            DO 30 J=KP1,N
               T=A(M,J)
               A(M,J)=A(K,J)
               A(K,J)=T
               IF(T.EQ.0.0)GO TO 30
               DO 25 I=KP1,N
               A(I,J)=A(I,J)+A(I,K)*T
      25     CONTINUE
      30   CONTINUE
      35 CONTINUE
   
         DO 50 K=1,N
            T=0.0
            IF(K.EQ.1)GO TO 45
            KM1=K-1
            DO 40 I=1,KM1
               T=T+A(I,K)*WORK(I)
      40   CONTINUE
      45   EK=1.0
            IF(T.LT.0.0)EK=-1.0
            IF(A(K,K).EQ.0.0)GO TO 90
            WORK(K)=-(EK+T)/A(K,K)
      50 CONTINUE
         DO 60 KB=1,NM1
            K=N-KB
            T=WORK(K)
            KP1=K+1
            DO 55 I=KP1,N
               T=T+A(I,K)*WORK(I)
      55   CONTINUE
            WORK(K)=T
            M=IPVT(K)
            IF(M.EQ.K)GO TO 60
            T=WORK(M)
            WORK(M)=WORK(K)
            WORK(K)=T
      60 CONTINUE
   
         YNORM=0.0
         DO 65 I=1,N
            YNORM=YNORM+ABS(WORK(I))
      65 CONTINUE
   
         CALL SOLVE(NDIM,N,A,WORK,IPVT)
   
         ZNORM=0.0
         DO 70 I=1,N
            ZNORM=ZNORM+ABS(WORK(I))
      70 CONTINUE
   

         COND=ANORM*ZNORM/YNORM
         IF(COND.LT.1.0)COND=1.0
         RETURN
   
      80 COND=1.0
         IF(A(1,1).NE.0.0)RETURN
   
      90 CONTINUE
         COND=1.0E+32
         RETURN
   END SUBROUTINE DECOMP
      
   SUBROUTINE SOLVE(NDIM,N,A,B,IPVT)
         INTEGER NDIM,N,IPVT(N)
         REAL(R_) A(NDIM,N),B(N)
         INTEGER KB,KM1,NM1,KP1,I,K,M
         REAL(R_) T
   
         IF(N.EQ.1) GO TO 50
         NM1=N-1
         DO 20 K=1,NM1
            KP1=K+1
            M=IPVT(K)
            T=B(M)
            B(M)=B(K)
            B(K)=T
            DO 10 I=KP1,N
               B(I)=B(I)+A(I,K)*T
      10   CONTINUE
      20 CONTINUE
   
         DO 40 KB=1,NM1
            KM1=N-KB
            K=KM1+1
            B(K)=B(K)/A(K,K)
            T=-B(K)
            DO 30 I=1,KM1
               B(I)=B(I)+A(I,K)*T
      30   CONTINUE
      40 CONTINUE
      50 B(1)=B(1)/A(1,1)
         RETURN
   END SUBROUTINE SOLVE
      
end program number_conditionaly
