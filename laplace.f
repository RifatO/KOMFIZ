       PROGRAM LAPLACE

!----------------------------------------------------------------------
!        PROGRAM ZA RJESAVANJE LAPLACEOVE JEDNACINE
!            U OBLASTI KVADRATA (PRVI ZADATAK)
!            Autor: R. Omerovic, septembar 2016.
!----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      DIMENSION V(1000,1000), Vn(1000,1000), ERROR(1000,1000)
      DIMENSION dX(1000), dY(1000)
      COMMON TOL, a, V0
      TOL=1.0E-8
      IT=0   ! BROJ ITERACIJA
      
      OPEN(UNIT=2,FILE='laplace1.dat',STATUS='UNKNOWN')
      OPEN(UNIT=3,FILE='laplace2.dat',STATUS='UNKNOWN')
      
      PRINT*,'---------------------------------------------------------'
      PRINT*,'PROGRAM RJESAVA LAPLACEOVU JEDNACINU'
      PRINT*,'U OBLASTI KVADRATA STRANICA a'
      PRINT*,'---------------------------------------------------------'
      PRINT*,'Unijeti duzinu stranice a, i broj tacaka u intervalu(0,a)'
      READ*,  a, N
      PRINT*, 'UNIJETI VRIJEDNOST Vo'
      READ*, V0
      PRINT*,'IZABRATI METODU RJESAVANJA LAPLACEOVE JEDNACINE'
      PRINT*, 'JAKOBIJEVA METODA - UNIJETI 1, GAUSS-SEIDEL 2, SOR 3'
      READ*, IZBOR
      
      DO I=1,N+2
         dX(I)=a*(I-1)/(N+1)
         dY(I)=a*(I-1)/(N+1)
      END DO
      
    !INICIJALIZACIJA POCETNIH I GRANICNIH VRIJEDNOSTI POTENCIJALA
      
      CALL initialV(N,V0,I,J,Vn)   
        ! Vn - varijabla iz posljednje iteracija
    !----- kraj inicijalizacije
    
    
    !---- KALKULACIJA NOVIH VRIJEDNOSTI ODABRANOM METODOM
    
   20 DO I=1,N+2
        DO J=1,N+2
        V(I,J)=Vn(I,J)
        END DO
      END DO
      
        
      SELECT CASE(IZBOR)
        CASE (1)
          CALL JACOBI(N,V,Vn)
        CASE(2)
          CALL GAUSSSEIDEL(N,V,Vn)
        CASE (3)
          CALL SOR(N,V,Vn,OMEGA)
        END SELECT
      ERR=0.
      DO I=2,N+1
         DO J=2,N+1
         ERR=ERR+((V(I,J)-Vn(I,J))**2)
         END DO
      END DO
         
      EERROR=sqrt(ERR)
      IT=IT+1

! Zaustavljanje postupka na zeljenom broju iteracija
!      IF (IT.LT.100) GO TO 20   
!      continue

       IF (EERROR .GT. TOL) GO TO 20  ! Usporedba greske i tolerancije
       CONTINUE
        
        
      ! UPISIVANJE PODATAK U MATRICNOJ FORMI U LAPLACE1.DAT'
       DO I=1,N+2
          WRITE (2,121), (Vn(I,J), J=1,N+2) 
       END DO
       WRITE(3,'(A4,A6,A12)') 'X','Y','V(X,Y)'
       WRITE(3,'(A24)') '------------------------------'
      
       !UPISIVANJE PODATAKA U FORMI TABELE U LAPLACE2.DAT
      DO I=1,N+2
          DO J=1,N+2
          WRITE(3,122) dX(I),dY(J),Vn(I,J)  
          END DO
       END DO
       
       PRINT*, '-----------------------------------------------------'
       WRITE(*,333) 'BROJ ITERACIJA:', IT
       WRITE(*,334) 'TEZINSKI FAKTOR (ZA SOR):', OMEGA
       WRITE(*,335) 'Postupak je konvergirao pri gresci:', EERROR
       PRINT*, '-----------------------------------------------------'
       PRINT*, 'PODACI SU ZAPISANI U laplace1.dat i laplace2.dat'
       PRINT*, '-----------------------------------------------------'       
  121  FORMAT(110F10.4)  
  122  FORMAT(2f6.2, f11.6)
  333  FORMAT(1X,A15,24X,I10)
  334  FORMAT(1X,A25,15X,f9.4)
  335  FORMAT(1X, A35,4x,e10.4)
      
      
      END PROGRAM LAPLACE
      
      
!----------------------------------------------------------------------      
   ! SUBROUTINA KOJA RACUNA POCENTE I GRANICNE VRIJEDNOSTI
   ! POTENCIJALA - PREMA USLOVIMA ZADATKA
   
      SUBROUTINE initialV(N,V0,K,M,U)
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      DIMENSION U(1000,1000)
         
      DO K=1,N+2
         DO M=1,N+2
         U(K,M)=0.       ! POPUNJVAMO MATRICU SA NULAMA
         END DO
      END DO
      
      DO K=1,N+1
        U(K,N+2)=V0*(K-1)/(N+1)   !POTENCIJAL NA GRANICI y=a
        U(N+2,K)=V0*(K-1)/(N+1)   ! POTENCIJAL NA GRANICI x=a
      END DO
      
      END 
!----------------------------------------------------------------------      
      
         
! JACOBIJEVA METODA ---------------------------------------------------
      SUBROUTINE JACOBI(N,U,Un)
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      DIMENSION U(1000,1000), Un(1000,1000)
      
      DO I=2,N+1
         DO J=2,N+1
         Un(I,J)=0.25*(U(I+1,J)+U(I-1,J)+U(I,J+1)+U(I,J-1))
         END DO
      END DO
      
      END 
! ---------------------------------------------------------------------
  
  
!  GAUSS-SEIDEL METODA ------------------------------------------------
      SUBROUTINE GAUSSSEIDEL(N,U,Un)
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      DIMENSION U(1000,1000), Un(1000,1000)
      
      DO I=2,N+1
         DO J=2,N+1
         Un(I,J)=0.25*(U(I+1,J)+Un(I-1,J)+U(I,J+1)+Un(I,J-1))
         END DO
      END DO
      
      END 
! ---------------------------------------------------------------------
  
  
  
!  SOR METODA ---------------------------------------------------------
      SUBROUTINE SOR(N,U,Un,OMEGA)
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      DIMENSION U(1000,1000), Un(1000,1000)
     
      PI=4.0*datan(1.0d0)
      OMEGA=2/(1+(PI/N))
!      OMEGA=1.96   ! Rucno unosenje vrijednosti tezisnkog faktora

      DO I=2,N+1
         DO J=2,N+1
      Un(I,J)=(1.-OMEGA)*U(I,J)+
     10.25*OMEGA*(U(I+1,J)+Un(I-1,J)+U(I,J+1)+Un(I,J-1))
         END DO
      END DO

      END 
! ---------------------------------------------------------------------
  
    
      
      
      
     
