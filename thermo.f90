        MODULE   THERMO

        IMPLICIT NONE

       real(kind=8), save :: runiv = 8.314d3
       integer, parameter :: nspeci = 5   
       REAL(KIND=8), SAVE :: MW(1:NSPECI) = (/27.0d0, 32.0d0, 43.0d0, 102.0d0, 28.0d0/) 

       REAL(KIND=8) :: THCOEF(1:7,1:2,1:NSPECI)
       REAL(KIND=8) :: TLOW(1:NSPECI), TUPP(1:NSPECI), TMID(1:NSPECI)
       CHARACTER (LEN = 10) :: SPNAME(1:NSPECI) 
       REAL(KIND=8) :: RGK(1:NSPECI) 


CONTAINS
         subroutine thermo_data() 

       INTEGER :: J, IOERROR, NS, NUNIT
       LOGICAL, DIMENSION(:), ALLOCATABLE :: EXISTDATA
       LOGICAL :: FLAG_EXIST
       REAL (KIND = 8) :: TLOWTEMP, TUPPTEMP, TMIDTEMP, TLOWDEF, TMIDDEF
       REAL (KIND = 8) :: TUPPDEF
       CHARACTER :: SPECIESNAME*18, DUMMY1*27


          NUNIT = 10

          SPNAME(1) = 'Al'
          SPNAME(2) = 'O2'
          SPNAME(3) = 'AlO'
          SPNAME(4) = 'AL2O3(L)'
          SPNAME(5) = 'N2'

          do ns = 1, nspeci
           RGK(ns) = runiv/MW(ns)
          enddo



          ALLOCATE (EXISTDATA(1:NSPECI))

      EXISTDATA = .FALSE.

      PRINT *,"Reading Chemkin thermodynamic database..."
      !OPEN(UNIT=NUNIT,FILE="sandiego20120907_therm.txt",FORM='FORMATTED')
      OPEN(UNIT=NUNIT,FILE="therm.dat",FORM='FORMATTED')
      READ(NUNIT,*)
      READ(NUNIT,'(3F10.0)') TLOWDEF, TMIDDEF, TUPPDEF
      DO
         READ(nunit,'(A18,A27,E10.0,E10.0,E8.0)',IOSTAT=IOERROR)    &
      SPECIESNAME, DUMMY1, TLOWTEMP, TUPPTEMP, TMIDTEMP
         IF (SPECIESNAME(1:3)=="END" .OR. IOERROR/=0) THEN
            EXIT
         END IF
         FLAG_EXIST = .FALSE.
         DO NS = 1, NSPECI
            IF (SPNAME(NS) == SPECIESNAME(1:10)) THEN
               IF (EXISTdATA(NS) .EQV. .TRUE.) THEN
                  WRITE(*,*)'SPECIES=',SPNAME(NS)
                  print*, 'SPECIES INFO FOUND TWICE'
                  CLOSE(NUNIT)
                  STOP
               END IF
               FLAG_EXIST = .TRUE.
               EXISTdATA(NS) = .TRUE.
               EXIT
            END IF
         END DO
         IF (FLAG_EXIST.EQV..TRUE.) THEN
            IF (TLOWTEMP == 0) THEN
               TLOW(NS) = TLOWDEF
            ELSE
               TLOW(NS) = TLOWTEMP
            END IF
            IF (TMIDTEMP == 0) THEN
               TMID(NS) = TMIDDEF
            ELSE
               TMID(NS) = TMIDTEMP
            END IF
            IF (TUPPTEMP == 0) THEN
               TUPP(NS) = TUPPDEF
            ELSE
               TUPP(NS) = TUPPTEMP
            END IF
            READ(NUNIT,'(5E15.0)') THCOEF(1:5,2,NS)
            READ(NUNIT,'(5E15.0)') THCOEF(6:7,2,NS), THCOEF(1:3,1,NS)
            READ(NUNIT,'(4E15.0)') THCOEF(4:7,1,NS)
         ELSE
            READ(NUNIT,*)
            READ(NUNIT,*)
            READ(NUNIT,*)
         END IF
      END DO
      CLOSE(NUNIT)
      DO NS = 1, NSPECI
         IF (EXISTdATA(NS).EQV..FALSE.) THEN
            WRITE(*,*)'SPECIES=',SPNAME(NS)
            print*, 'SPECIES INFO NOT FOUND IN THERMODYNAMIC DATA FILE'
         END IF
      END DO
      PRINT *," Done reading Chemkin database"


          DEALLOCATE (EXISTDATA)

        return
        end subroutine 


       END MODULE

!--------------------------------------------------------------------------------


