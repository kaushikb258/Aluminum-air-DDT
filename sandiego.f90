
             subroutine combustion(T,conc,omegadot)

          USE THERMO

      IMPLICIT NONE

          real(kind=8) :: T, conc(1:nspeci)
          real(kind=8), dimension(1:nspeci) :: omegadot
          real(kind=8) :: fa, SUMOMEGADOT
          
           !Al/O2 ratio
           fa = conc(1)/(conc(2)+1.0d-12)

           
              omegadot = 0.0d0


              if(T.ge.3750.0d0) then
               if(fa.ge.27.0d0/16.0d0) then
                ! fuel rich
                omegadot(1) = -conc(2)/16.0d0*27.0d0
                omegadot(2) = -conc(2) 
                omegadot(3) = conc(2)/16.0d0*43.0d0
               else
                ! fuel lean 
                omegadot(1) = -conc(1)
                omegadot(2) = -conc(1)/27.0d0*16.0d0 
                omegadot(3) = conc(1)/27.0d0*43.0d0
               endif
              else
               if(fa.ge.27.0d0/24.0d0) then
                ! fuel rich
                omegadot(1) = -conc(2)/24.0d0*27.0d0
                omegadot(2) = -conc(2) 
                omegadot(4) = conc(2)/24.0d0*51.0d0
               else
                ! fuel lean 
                omegadot(1) = -conc(1)
                omegadot(2) = -conc(1)/27.0d0*24.0d0 
                omegadot(4) = conc(1)/27.0d0*51.0d0
               endif
              endif  
     
           
                sumomegadot = sum(omegadot(:))
                if(abs(sumomegadot).gt.1.0d-6) then
                  print*, 'sum omegadot ', sumomegadot, omegadot
                  stop
                endif 


            return 
            end subroutine
          
!--------------------------------------------------------------------------------

           SUBROUTINE COMPUTE_CPK (TEMP, CPK)
      USE THERMO

      IMPLICIT NONE

      INTEGER :: NS
      REAL (KIND = 8), INTENT(IN) :: TEMP
      REAL (KIND = 8), INTENT(OUT) :: CPK(1:NSPECI)
      REAL (KIND = 8) :: CPP

       DO NS = 1, NSPECI   
         
         CPP = 0.0D0  

         IF(TEMP.LT.TLOW(NS)) THEN
            CPP = TLOW(NS) *        THCOEF(5,1,NS)
            CPP = TLOW(NS) *       (THCOEF(4,1,NS) + CPP)
            CPP = TLOW(NS) *       (THCOEF(3,1,NS) + CPP)
            CPP = TLOW(NS) *       (THCOEF(2,1,NS) + CPP)
            CPP =                  THCOEF(1,1,NS) + CPP
         ELSEIF ( TEMP .GT. TUPP(NS) ) THEN
            CPP = TUPP(NS) *        THCOEF(5,2,NS)
            CPP = TUPP(NS) *       (THCOEF(4,2,NS) + CPP)
            CPP = TUPP(NS) *       (THCOEF(3,2,NS) + CPP)
            CPP = TUPP(NS) *       (THCOEF(2,2,NS) + CPP)
            CPP =                  THCOEF(1,2,NS) + CPP
         ELSEIF ( TEMP > TMID(NS) ) THEN
            CPP = TEMP *  THCOEF(5,2,NS)
            CPP = TEMP * (THCOEF(4,2,NS) + CPP)
            CPP = TEMP * (THCOEF(3,2,NS) + CPP)
            CPP = TEMP * (THCOEF(2,2,NS) + CPP)
            CPP =         THCOEF(1,2,NS) + CPP
         ELSE
            CPP = TEMP *  THCOEF(5,1,NS)
            CPP = TEMP * (THCOEF(4,1,NS) + CPP)
            CPP = TEMP * (THCOEF(3,1,NS) + CPP)
            CPP = TEMP * (THCOEF(2,1,NS) + CPP)
            CPP =         THCOEF(1,1,NS) + CPP
         ENDIF
          CPK(NS) = CPP * RGK(NS)

       ENDDO 

      RETURN
      END SUBROUTINE
!--------------------------------------------------------------------

      SUBROUTINE COMPUTE_HFK (TEMP, HFK)

      USE THERMO

      IMPLICIT NONE

      INTEGER :: NS
      REAL (KIND = 8), INTENT(IN) :: TEMP
      REAL (KIND = 8), INTENT(OUT) :: HFK(1:NSPECI)
      REAL (KIND = 8) :: HFP, CPP

      REAL(KIND=8), SAVE :: R5I = 1.0D0/5.0D0  
      REAL(KIND=8), SAVE :: R4I = 1.0D0/4.0D0  
      REAL(KIND=8), SAVE :: R3I = 1.0D0/3.0D0  
      REAL(KIND=8), SAVE :: R2I = 1.0D0/2.0D0  

       DO NS = 1, NSPECI  

         HFP = 0.0D0
         CPP = 0.0D0

         IF(TEMP.LT.TLOW(NS)) THEN
            HFP = TLOW(NS) *  R5I * THCOEF(5,1,NS)
            HFP = TLOW(NS) * (R4I * THCOEF(4,1,NS) + HFP)
            HFP = TLOW(NS) * (R3I * THCOEF(3,1,NS) + HFP)
            HFP = TLOW(NS) * (R2I * THCOEF(2,1,NS) + HFP)
            HFP = TLOW(NS) * (      THCOEF(1,1,NS) + HFP)
            HFP =                   THCOEF(6,1,NS) + HFP

            CPP = TLOW(NS) *        THCOEF(5,1,NS)
            CPP = TLOW(NS) *       (THCOEF(4,1,NS) + CPP)
            CPP = TLOW(NS) *       (THCOEF(3,1,NS) + CPP)
            CPP = TLOW(NS) *       (THCOEF(2,1,NS) + CPP)
            CPP =                   THCOEF(1,1,NS) + CPP

            HFP = HFP + CPP * (TEMP - TLOW(NS))
         ELSE IF(TEMP .GT. TUPP(NS)) THEN
            CPP = TUPP(NS) *        THCOEF(5,2,NS)
            CPP = TUPP(NS) *       (THCOEF(4,2,NS) + CPP)
            CPP = TUPP(NS) *       (THCOEF(3,2,NS) + CPP)
            CPP = TUPP(NS) *       (THCOEF(2,2,NS) + CPP)
            CPP =                   THCOEF(1,2,NS) + CPP

            HFP = TUPP(NS) *  R5I * THCOEF(5,2,NS)
            HFP = TUPP(NS) * (R4I * THCOEF(4,2,NS) + HFP)
            HFP = TUPP(NS) * (R3I * THCOEF(3,2,NS) + HFP)
            HFP = TUPP(NS) * (R2I * THCOEF(2,2,NS) + HFP)
            HFP = TUPP(NS) * (      THCOEF(1,2,NS) + HFP)
            HFP =                   THCOEF(6,2,NS) + HFP

            HFP = HFP + CPP * (TEMP - TUPP(NS))
         ELSEIF ( TEMP > TMID(NS) ) THEN
            HFP = TEMP *  R5I * THCOEF(5,2,NS)
            HFP = TEMP * (R4I * THCOEF(4,2,NS) + HFP)
            HFP = TEMP * (R3I * THCOEF(3,2,NS) + HFP)
            HFP = TEMP * (R2I * THCOEF(2,2,NS) + HFP)
            HFP = TEMP * (      THCOEF(1,2,NS) + HFP)
            HFP =               THCOEF(6,2,NS) + HFP
         ELSE
            HFP = TEMP *  R5I * THCOEF(5,1,NS)
            HFP = TEMP * (R4I * THCOEF(4,1,NS) + HFP)
            HFP = TEMP * (R3I * THCOEF(3,1,NS) + HFP)
            HFP = TEMP * (R2I * THCOEF(2,1,NS) + HFP)
            HFP = TEMP * (      THCOEF(1,1,NS) + HFP)
            HFP =               THCOEF(6,1,NS) + HFP
         ENDIF
         HFK(NS) = HFP * RGK(NS)

      ENDDO


      RETURN
      END


!--------------------------------------------------------------------

      SUBROUTINE COMPUTE_SK (TEMP, SK)

      USE THERMO

      IMPLICIT NONE

      INTEGER :: NS
      REAL (KIND = 8), INTENT(IN) :: TEMP
      REAL (KIND = 8), INTENT(OUT) :: SK(1:NSPECI)
      REAL (KIND = 8) :: SP

      REAL(KIND=8), SAVE :: R5I = 1.0D0/5.0D0  
      REAL(KIND=8), SAVE :: R4I = 1.0D0/4.0D0  
      REAL(KIND=8), SAVE :: R3I = 1.0D0/3.0D0  
      REAL(KIND=8), SAVE :: R2I = 1.0D0/2.0D0  

       DO NS = 1, NSPECI  

         SP = 0.0D0

         IF(TEMP.LT.TLOW(NS)) THEN
            SP = R4I*THCOEF(5,1,NS)*TLOW(NS) + R3I*THCOEF(4,1,NS)
            SP = R2I*THCOEF(3,1,NS) + TLOW(NS)*SP 
            SP = THCOEF(2,1,NS) + TLOW(NS)*SP        
            SP = THCOEF(1,1,NS)*LOG(TLOW(NS)) + THCOEF(7,1,NS) + TLOW(NS)*SP
         ELSE IF(TEMP .GT. TUPP(NS)) THEN
            SP = R4I*THCOEF(5,2,NS)*TUPP(NS) + R3I*THCOEF(4,2,NS)
            SP = R2I*THCOEF(3,2,NS) + TUPP(NS)*SP 
            SP = THCOEF(2,2,NS) + TUPP(NS)*SP        
            SP = THCOEF(1,2,NS)*LOG(TUPP(NS)) + THCOEF(7,2,NS) + TUPP(NS)*SP
         ELSE IF ( TEMP > TMID(NS) ) THEN
            SP = R4I*THCOEF(5,2,NS)*TEMP + R3I*THCOEF(4,2,NS)
            SP = R2I*THCOEF(3,2,NS) + TEMP*SP 
            SP = THCOEF(2,2,NS) + TEMP*SP        
            SP = THCOEF(1,2,NS)*LOG(TEMP) + THCOEF(7,2,NS) + TEMP*SP
         ELSE
            SP = R4I*THCOEF(5,1,NS)*TEMP + R3I*THCOEF(4,1,NS)
            SP = R2I*THCOEF(3,1,NS) + TEMP*SP 
            SP = THCOEF(2,1,NS) + TEMP*SP        
            SP = THCOEF(1,1,NS)*LOG(TEMP) + THCOEF(7,1,NS) + TEMP*SP
         ENDIF
         SK(NS) = SP * RGK(NS)

      ENDDO


      RETURN
      END


!--------------------------------------------------------------------

        SUBROUTINE FROM_IE_TO_T(E_I, YK, TEMP)

      USE THERMO  

      IMPLICIT NONE

      INTEGER:: NS
      REAL (KIND = 8) :: E_I, TEMP, YK(1:NSPECI)
      REAL (KIND = 8) :: ERRMAX, ERROR, DTEMP, HFP, CPP,RGK_YK,    &
                                  CPMIX, HFMIX, FZERO, TEMP0
      REAL(KIND=8) :: RGL, HFL, CPL
      REAL(KIND=8), SAVE :: R5I = 1.0D0/5.0D0  
      REAL(KIND=8), SAVE :: R4I = 1.0D0/4.0D0  
      REAL(KIND=8), SAVE :: R3I = 1.0D0/3.0D0  
      REAL(KIND=8), SAVE :: R2I = 1.0D0/2.0D0  

      INTEGER ITER
      PARAMETER ( ERRMAX = 1.0D-5 )

      RGL = 0.0D0
      DO NS = 1, NSPECI
         RGL = RGL + YK(NS) * RGK(NS)
      ENDDO

         ITER = 0
         ERROR = 2.D0 * ERRMAX 
         DTEMP = 0.0d0
         TEMP0 = TEMP
         DO WHILE ( ERROR .GT. ERRMAX )
            TEMP = TEMP + DTEMP
            ITER = ITER + 1
            IF ( ITER .GT. 10 ) THEN
               PRINT*, E_I, TEMP, YK(:)
               PRINT*, 'TEMP NOT CONVERGING (ie2t) '
               stop
            END IF
            HFMIX = 0.0d0
            CPMIX = 0.0d0
            DO NS = 1, NSPECI
               RGK_YK = RGK(NS) * YK(NS)
               IF(TEMP.LT.TLOW(NS)) THEN
                  CPP = TLOW(NS) *        THCOEF(5,1,NS)
                  CPP = TLOW(NS) * (      THCOEF(4,1,NS) + CPP)
                  CPP = TLOW(NS) * (      THCOEF(3,1,NS) + CPP)
                  CPP = TLOW(NS) * (      THCOEF(2,1,NS) + CPP)
                  CPP =                  THCOEF(1,1,NS) + CPP

                  HFP = TLOW(NS) *  R5I * THCOEF(5,1,NS)
                  HFP = TLOW(NS) * (R4I * THCOEF(4,1,NS) + HFP)
                  HFP = TLOW(NS) * (R3I * THCOEF(3,1,NS) + HFP)
                  HFP = TLOW(NS) * (R2I * THCOEF(2,1,NS) + HFP)
                  HFP = TLOW(NS) * (      THCOEF(1,1,NS) + HFP)
                  HFP =                  THCOEF(6,1,NS) + HFP

                  HFP =  HFP + CPP * (TEMP - TLOW(NS))
               ELSEIF ( TEMP .GT. TUPP(NS) ) THEN
                  CPP = TUPP(NS) *        THCOEF(5,2,NS)
                  CPP = TUPP(NS) * (      THCOEF(4,2,NS) + CPP)
                  CPP = TUPP(NS) * (      THCOEF(3,2,NS) + CPP)
                  CPP = TUPP(NS) * (      THCOEF(2,2,NS) + CPP)
                  CPP =                  THCOEF(1,2,NS) + CPP

                  HFP = TUPP(NS) *  R5I * THCOEF(5,2,NS)
                  HFP = TUPP(NS) * (R4I * THCOEF(4,2,NS) + HFP)
                  HFP = TUPP(NS) * (R3I * THCOEF(3,2,NS) + HFP)
                  HFP = TUPP(NS) * (R2I * THCOEF(2,2,NS) + HFP)
                  HFP = TUPP(NS) * (      THCOEF(1,2,NS) + HFP)
                  HFP =                  THCOEF(6,2,NS) + HFP

                  HFP =  HFP + CPP * (TEMP - TUPP(NS))
               ELSEIF ( TEMP .GT. TMID(NS) ) THEN
                  CPP = TEMP *        THCOEF(5,2,NS)
                  CPP = TEMP * (      THCOEF(4,2,NS) + CPP)
                  CPP = TEMP * (      THCOEF(3,2,NS) + CPP)
                  CPP = TEMP * (      THCOEF(2,2,NS) + CPP)
                  CPP =                  THCOEF(1,2,NS) + CPP

                  HFP = TEMP *  R5I * THCOEF(5,2,NS)
                  HFP = TEMP * (R4I * THCOEF(4,2,NS) + HFP)
                  HFP = TEMP * (R3I * THCOEF(3,2,NS) + HFP)
                  HFP = TEMP * (R2I * THCOEF(2,2,NS) + HFP)
                  HFP = TEMP * (      THCOEF(1,2,NS) + HFP)
                  HFP =                  THCOEF(6,2,NS) + HFP
               ELSE
                  CPP = TEMP *        THCOEF(5,1,NS)
                  CPP = TEMP * (      THCOEF(4,1,NS) + CPP)
                  CPP = TEMP * (      THCOEF(3,1,NS) + CPP)
                  CPP = TEMP * (      THCOEF(2,1,NS) + CPP)
                  CPP =                  THCOEF(1,1,NS) + CPP

                  HFP = TEMP *  R5I * THCOEF(5,1,NS)
                  HFP = TEMP * (R4I * THCOEF(4,1,NS) + HFP)
                  HFP = TEMP * (R3I * THCOEF(3,1,NS) + HFP)
                  HFP = TEMP * (R2I * THCOEF(2,1,NS) + HFP)
                  HFP = TEMP * (      THCOEF(1,1,NS) + HFP)
                  HFP =                  THCOEF(6,1,NS) + HFP
               ENDIF
               HFMIX = HFMIX + HFP * RGK_YK
               CPMIX = CPMIX + CPP * RGK_YK
            ENDDO
            FZERO = HFMIX - RGL * TEMP - E_I
            DTEMP = -FZERO / (CPMIX - RGL)
            ERROR = ABS(FZERO) / RGL
         ENDDO
         HFL = HFMIX
         CPL = CPMIX


        RETURN
        END
     
!--------------------------------------------------------------------

          
        SUBROUTINE FROM_IH_TO_T(H, YK, TEMP)

      USE THERMO  

      IMPLICIT NONE

      INTEGER:: NS
      REAL (KIND = 8) :: H, TEMP, YK(1:NSPECI)
      REAL (KIND = 8) :: ERRMAX, ERROR, DTEMP, HFP, CPP,RGK_YK,    &
                                  CPMIX, HFMIX, FZERO, TEMP0
      REAL(KIND=8) :: RGL, HFL, CPL
      REAL(KIND=8), SAVE :: R5I = 1.0D0/5.0D0  
      REAL(KIND=8), SAVE :: R4I = 1.0D0/4.0D0  
      REAL(KIND=8), SAVE :: R3I = 1.0D0/3.0D0  
      REAL(KIND=8), SAVE :: R2I = 1.0D0/2.0D0  

      INTEGER ITER
      PARAMETER ( ERRMAX = 1.0D-5 )


     

      RGL = 0.0D0
      DO NS = 1, NSPECI
         RGL = RGL + YK(NS) * RGK(NS)
      ENDDO

         ITER = 0
         ERROR = 2.D0 * ERRMAX 
         DTEMP = 0.0d0
         TEMP0 = TEMP
         DO WHILE ( ERROR .GT. ERRMAX )
            TEMP = TEMP + DTEMP
            ITER = ITER + 1
            IF ( ITER .GT. 10 ) THEN
               PRINT*, H, TEMP, YK(:)
               PRINT*, 'TEMP NOT CONVERGING '
               stop
            END IF
            HFMIX = 0.0d0
            CPMIX = 0.0d0
            DO NS = 1, NSPECI
               RGK_YK = RGK(NS) * YK(NS)
               IF(TEMP.LT.TLOW(NS)) THEN
                  CPP = TLOW(NS) *        THCOEF(5,1,NS)
                  CPP = TLOW(NS) * (      THCOEF(4,1,NS) + CPP)
                  CPP = TLOW(NS) * (      THCOEF(3,1,NS) + CPP)
                  CPP = TLOW(NS) * (      THCOEF(2,1,NS) + CPP)
                  CPP =                  THCOEF(1,1,NS) + CPP

                  HFP = TLOW(NS) *  R5I * THCOEF(5,1,NS)
                  HFP = TLOW(NS) * (R4I * THCOEF(4,1,NS) + HFP)
                  HFP = TLOW(NS) * (R3I * THCOEF(3,1,NS) + HFP)
                  HFP = TLOW(NS) * (R2I * THCOEF(2,1,NS) + HFP)
                  HFP = TLOW(NS) * (      THCOEF(1,1,NS) + HFP)
                  HFP =                  THCOEF(6,1,NS) + HFP

                  HFP =  HFP + CPP * (TEMP - TLOW(NS))
               ELSEIF ( TEMP .GT. TUPP(NS) ) THEN
                  CPP = TUPP(NS) *        THCOEF(5,2,NS)
                  CPP = TUPP(NS) * (      THCOEF(4,2,NS) + CPP)
                  CPP = TUPP(NS) * (      THCOEF(3,2,NS) + CPP)
                  CPP = TUPP(NS) * (      THCOEF(2,2,NS) + CPP)
                  CPP =                  THCOEF(1,2,NS) + CPP

                  HFP = TUPP(NS) *  R5I * THCOEF(5,2,NS)
                  HFP = TUPP(NS) * (R4I * THCOEF(4,2,NS) + HFP)
                  HFP = TUPP(NS) * (R3I * THCOEF(3,2,NS) + HFP)
                  HFP = TUPP(NS) * (R2I * THCOEF(2,2,NS) + HFP)
                  HFP = TUPP(NS) * (      THCOEF(1,2,NS) + HFP)
                  HFP =                  THCOEF(6,2,NS) + HFP

                  HFP =  HFP + CPP * (TEMP - TUPP(NS))
               ELSEIF ( TEMP .GT. TMID(NS) ) THEN
                  CPP = TEMP *        THCOEF(5,2,NS)
                  CPP = TEMP * (      THCOEF(4,2,NS) + CPP)
                  CPP = TEMP * (      THCOEF(3,2,NS) + CPP)
                  CPP = TEMP * (      THCOEF(2,2,NS) + CPP)
                  CPP =                  THCOEF(1,2,NS) + CPP

                  HFP = TEMP *  R5I * THCOEF(5,2,NS)
                  HFP = TEMP * (R4I * THCOEF(4,2,NS) + HFP)
                  HFP = TEMP * (R3I * THCOEF(3,2,NS) + HFP)
                  HFP = TEMP * (R2I * THCOEF(2,2,NS) + HFP)
                  HFP = TEMP * (      THCOEF(1,2,NS) + HFP)
                  HFP =                  THCOEF(6,2,NS) + HFP
               ELSE
                  CPP = TEMP *        THCOEF(5,1,NS)
                  CPP = TEMP * (      THCOEF(4,1,NS) + CPP)
                  CPP = TEMP * (      THCOEF(3,1,NS) + CPP)
                  CPP = TEMP * (      THCOEF(2,1,NS) + CPP)
                  CPP =                  THCOEF(1,1,NS) + CPP

                  HFP = TEMP *  R5I * THCOEF(5,1,NS)
                  HFP = TEMP * (R4I * THCOEF(4,1,NS) + HFP)
                  HFP = TEMP * (R3I * THCOEF(3,1,NS) + HFP)
                  HFP = TEMP * (R2I * THCOEF(2,1,NS) + HFP)
                  HFP = TEMP * (      THCOEF(1,1,NS) + HFP)
                  HFP =                  THCOEF(6,1,NS) + HFP
               ENDIF
               HFMIX = HFMIX + HFP * RGK_YK
               CPMIX = CPMIX + CPP * RGK_YK
            ENDDO
            FZERO = HFMIX  - H
            DTEMP = -FZERO / (CPMIX)
            ERROR = ABS(FZERO) / RGL
         ENDDO
         HFL = HFMIX
         CPL = CPMIX


        RETURN
        END
     
!--------------------------------------------------------------------


              subroutine from_es_to_Ts(es,Ts)

          implicit none

             real(kind=8) :: es, e, Ts, Tm, Tb
             real(kind=8) :: Lm, Lb, e300, em, eb
               

              Tm = 933.0d0
              Tb = 2792.0d0  
              Lm = 10.7d3
              Lb = 290.0d3 

              e300 = 7.2d3 ! J/mol, es at 300 K    
              em = 24.924d3 ! J/mol, es at melt temperature 
              eb = em + Lm + 32.0d0*(Tb-Tm) ! J/mol, es at boil temperature 

              
                ! convert from J/kg to J/mol  
                  e = es*27.0d-3 


                  
                if(e.le.e300) then
                  Ts = e/24.0d0
                  return
                endif  


                if(e.gt.e300.and.e.lt.em) then
                  Ts = (-20.208d0 + sqrt(20.208d0**2.0d0 - 4.0d0*6.32d-3*(568.8d0-e)))/(2.0d0*6.32d-3)                    
                  return
                endif 
                  

                if(e.ge.em.and.e.le.em+Lm) then
                  Ts = Tm
                  return
                endif
                   

                if(e.gt.em+Lm.and.e.le.eb) then
                  Ts = Tm + (e-em-Lm)/32.0d0
                  return
                endif


                if(e.gt.eb.and.e.le.eb+Lb) then
                  Ts = Tb
                  return
                endif


                if(e.gt.eb+Lb) then
                  Ts = Tb + (e-eb-Lb)/32.0d0
                  return
                endif



                return
              end subroutine 

!-------------------------------------------------------



