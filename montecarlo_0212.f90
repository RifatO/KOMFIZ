PROGRAM Quantum_Monte_Carlo

! ------------------------PROGRAM ZA IZRAČUNAVANJE ENERGIJE OSNOVNOG STANJA ATOMA -------------------------------------!
! --------------------------VODIKA I HELIJUMA POMOĆU KVANTNOG MONTE-CARLO METODA---------------------------------------!

IMPLICIT NONE
INTEGER :: MC_koraci, nvar_param, termal, &
            atom_izbor
DOUBLE PRECISION ::  korak
DOUBLE PRECISION, DIMENSION(100) :: E_sum, E_sum2

    CALL input_data( atom_izbor, nvar_param, &       ! učitavanje ULAZNIH podataka
                    MC_koraci, termal, korak)

    E_sum=0.0d0 ;  E_sum2=0.0d0                      ! Početne vrijednosti suma energije i kvadrata energije

    CALL QMC( atom_izbor, nvar_param, &      ! Monte - Carlo uzorkovanje
                     termal, MC_koraci, korak, &
                     E_sum, E_sum2)

    CALL output(nvar_param, MC_koraci, E_sum, E_sum2) ! Ispisivanje rezultata
END PROGRAM Quantum_Monte_Carlo


!-----------------------------------------SUBRUTINE I FUNCKIJE---------------------------------------------------------!

SUBROUTINE QMC(atom_izbor,nvar_param,termal, MC_koraci, korak,   &
                                                E_sum, E_sum2)

!----- Subrutina QMC vrši uzorkovanje na zadatom broju tačaka u
!------kojima se izračunava vrijednost energije.
!------ Pri odabiru tačaka koristi se Metropolis algoritam.


IMPLICIT NONE
DOUBLE PRECISION, DIMENSION(nvar_param)  ::  E_sum, E_sum2
INTEGER  :: dimension, atom_izbor, nvar_param, &
                       MC_koraci, termal
DOUBLE PRECISION  ::   korak
INTEGER ::  ciklus, var_param, korak_ok, dim, i, j

DOUBLE PRECISION :: valna_0, valna_1, alpha, E, E2, &
                      delta_e, valna_funk, E_loc
DOUBLE PRECISION,  DIMENSION(3,2) :: r_0, r_1

    alpha = 1.69       ! Početna vrijednost varijacionog parametra


    DO var_param=1,nvar_param                ! Petlja po varijacionim parametrima

       alpha = alpha + 0.0d0                   ! Varijacionu parametar dobiva novu vrijednost
       E = 0.0d0; E2 = 0.0d0; korak_ok =0; delta_e=0.0d0;  ! Pocetne vrijednosti energija

       DO i = 1,  atom_izbor                   !  Pocetni probni polozaj cestice
          DO j=1, 3                            !  Uzima se u 3 dimenzije
             r_0(j,i) = korak*(rand()-0.5d0)/alpha
          ENDDO
       ENDDO
       valna_1 = valna_funk(r_0, alpha,  atom_izbor)

       DO ciklus = 1, MC_koraci+termal         ! Petlja po broju MC ciklusa

          DO i = 1,  atom_izbor                ! Novi polzaj elektrona
             DO j=1, 3
                r_1(j,i)=r_0(j,i) + korak*(rand()-0.5d0)/alpha
             ENDDO
          ENDDO
          valna_0 = valna_funk(r_1, alpha,  atom_izbor)

          IF (rand() <= valna_0*valna_0/valna_1/valna_1 ) THEN  !  METROPOLIS TEST
             r_0=r_1
             valna_1 = valna_0
             korak_ok = korak_ok+1
          ENDIF

          IF ( ciklus > termal ) THEN        ! Kada se zavrsi termalizacija, pocinje raucnanje lokalne energije
             delta_e = E_loc(r_0, alpha, valna_1, &
                                     atom_izbor)

             E = E+delta_e
             E2 = E2+delta_e*delta_e
          ENDIF
       ENDDO
       WRITE(*,*)'Varijacioni parametar i Broj prihvacenih koraka', &
                                                         alpha, korak_ok

       E_sum(var_param) = E/MC_koraci        ! Sumiranje energije
       E_sum2(var_param) = E2/MC_koraci      ! Sumiranje energije
    ENDDO

END SUBROUTINE  QMC



DOUBLE PRECISION FUNCTION valna_funk(r, alpha, atom_izbor)

!----- FUNKCIJA KOJA RAČUNA VRIJEDNOST VALNE FUNKCIJE-------
!------ARGUMENTI: polozaj elektrona, varijacioni parametar
!------           i broj elektrona -------------------------!

IMPLICIT NONE
INTEGER  :: dimension, atom_izbor
DOUBLE PRECISION, DIMENSION(3,2)  :: r
DOUBLE PRECISION  :: alpha, alpha3
INTEGER ::  i, j
DOUBLE PRECISION :: argument, r_electron
DOUBLE PRECISION :: PI

    PI=4.*DATAN(1.0d0)
    argument = 0.0d0
    DO i = 1, atom_izbor
       r_electron = 0.0d0
       DO j = 1, 3                ! Položaj pojedinačnog elektrona u 3D
          r_electron = r_electron+ r(j,i)*r(j,i)
       ENDDO
       argument = argument + SQRT(r_electron)
    ENDDO
    alpha3=alpha*alpha*alpha
    argument = argument*alpha
    IF (atom_izbor.EQ.1) THEN
       valna_funk=EXP(-argument)   !Valna funkcija za vodonik
    ELSE
     valna_funk = alpha3*EXP(-argument)/PI   !Valna funkcija za helijum
    END IF
END FUNCTION valna_funk



DOUBLE PRECISION FUNCTION E_loc(r, alpha, valna_1, &
                                       atom_izbor)
!----- Funkcija računa vrijednost lokalne energije------
!------ Argumenti: polozaj,varijacioni parametar,
!-------valna funkcija i broj elektrona----------------!
IMPLICIT NONE
INTEGER  :: dimension, atom_izbor
DOUBLE PRECISION, DIMENSION(3,2)  :: r
DOUBLE PRECISION  :: alpha, valna_1
INTEGER :: i, j , k
DOUBLE PRECISION :: e_local, valna_minus, valna_plus, E_k, E_p, &
                    r_12, r_electron, valna_funk
DOUBLE PRECISION, DIMENSION(3,2) :: r_plus, r_minus
DOUBLE PRECISION, PARAMETER  :: h = 0.001d0    !Korak za trazenje izvoda
DOUBLE PRECISION, PARAMETER  :: h2 = 1000000d0 !Inverzno kvadratno h

    r_plus = r; r_minus = r

    E_k = 0.0d0
    DO i = 1, atom_izbor              ! Izračunavanje kinetičke energije
       DO j = 1, 3
          r_plus(j,i) = r(j,i)+h
          r_minus(j,i) = r(j,i)-h
          valna_minus = valna_funk(r_minus, alpha, atom_izbor)
          valna_plus  = valna_funk(r_plus, alpha, atom_izbor)
          E_k=E_k-((valna_plus-valna_1)+(valna_minus-valna_1)) ! Izračunavanje Laplaciana metodom središnje razlike
          r_plus(j,i) = r(j,i)                                 ! prema jednačini
          r_minus(j,i) = r(j,i)
       ENDDO
    ENDDO

    E_k = 0.5d0*h2*E_k/valna_1 ! uključiti masu elektrona i (h/2pi)^2, pa podijeliti sa valnom funkcijom
    ! izračunavanje potencijala E
    E_p = 0.0d0
    ! doprinos potencijala interakcije elektron - proton
    DO i = 1, atom_izbor
       r_electron = 0.0d0
       DO j = 1, 3
          r_electron = r_electron+r(j,i)*r(j,i)
       ENDDO
       E_p = E_p-atom_izbor/SQRT(r_electron)
    ENDDO
    ! doprinos potencijala interakcije elektron - elektron
    DO i = 1, atom_izbor-1
       DO J= i+1, atom_izbor
	  r_12 = 0
          DO k = 1, 3
             r_12 = r_12+(r(k,i)-r(k,j))*(r(k,i)-r(k,j))
          ENDDO
	  E_p =  E_p+1/SQRT(r_12)
       ENDDO
    ENDDO
    E_loc = E_p+E_k

END FUNCTION E_loc



SUBROUTINE input_data( atom_izbor, nvar_param, &
                      MC_koraci, termal, korak)

!--------- Subroutina koja učitava ulazne podatke
!--------- sa konzole
!------------------------------------------------!

IMPLICIT NONE
INTEGER  :: dimension, atom_izbor, nvar_param, &
                          MC_koraci, termal
DOUBLE PRECISION  ::   korak
    WRITE(*,*)'Izbor atoma: za VODIK unijeti 1, za HELIJ 2:'
    READ(*,*) atom_izbor
    if (atom_izbor.NE.1 .AND. atom_izbor.NE. 2) then
     WRITE(*,*) "ERROR: Pogresan unos broja"
     STOP
    end if

    WRITE(*,*)'Broj Monte-Carlo ciklusa:'
    READ(*,*)MC_koraci
    WRITE(*,*)'Duzina koraka:'
    READ(*,*)korak
    WRITE(*,*)'Broj koraka do ekvilibrijuma:'
    READ(*,*) termal
    WRITE(*,*)'Broj varijacionih parametara:'
    READ(*,*) nvar_param
END SUBROUTINE input_data



SUBROUTINE output(nvar_param, MC_koraci, E_sum, E_sum2)

!---------- Subroutina koja zapisuje rezultate QMC simulacije
!----------- u file "QMC_rezultati.dat":
!----------- Varijacioni parametar, Energiju, Varijancu,
!----------- Standardnu devijaciju

DOUBLE PRECISION, DIMENSION(nvar_param)  ::  E_sum,E_sum2

INTEGER  :: nvar_param, MC_koraci
INTEGER :: i
DOUBLE PRECISION :: alpha, var, dev

    OPEN(UNIT=6,FILE='QMC_rezultati.dat')  !OUTPUT file
    alpha = 0.0
    WRITE(6,'(4A10)') 'VAR.PARAM',' ENERGIJA','VAR.','ST.DEV.'
    WRITE(6,*)'---------------------------------------'
    DO i=1, nvar_param
       alpha = alpha+0.1d0
       var = E_sum2(i)-E_sum(i)*E_sum(i)
       dev=SQRT(var/MC_koraci)
       WRITE(6,'(4F10.5)')alpha,E_sum(i),var, dev
    ENDDO
    WRITE(6,*)'---------------------------------------'
    CLOSE(6)

END SUBROUTINE  output
