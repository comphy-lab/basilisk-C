	SUBROUTINE CONFORMSYSTEM ( THETA,ntheta,
     &                             ReZmin,ImZmin,ReZmax,ImZmax,nslit,   
     &                             ReCn,ImCn,ncoef,Q,ReV0,ImV0,Phi0,
     &                             varout,ReZ,ImZ,
     &                             PHI,PSI,Vx,Vy,J_PHI,J_PSI,J_Vx,J_Vy,
! work arrays:
     &                             J_PHI1,J_PSI1,J_Vx1,J_Vy1,
     &                             Xdof,ndof)
	!DEC$ ATTRIBUTES DLLEXPORT,ALIAS:'_conformsystem_' :: ConformSystem 

	IMPLICIT NONE
	
	INTEGER nslit,ncoef,ntheta,ndof
        INTEGER varout(3)
        DOUBLE PRECISION THETA(ntheta)
        DOUBLE PRECISION ReZmin(nslit),ImZmin(nslit),
     &                   ReZmax(nslit),ImZmax(nslit)
	DOUBLE PRECISION ReCn(ncoef,nslit),ImCn(ncoef,nslit),Q(nslit)
	DOUBLE PRECISION ReV0,ImV0,Phi0

C       Output arrays (double, allocated in C wrapper):

	DOUBLE PRECISION ReZ(ntheta*nslit),ImZ(ntheta*nslit)
	DOUBLE PRECISION PHI(ntheta*nslit),PSI(ntheta*nslit),
     &                   Vx(ntheta*nslit),Vy(ntheta*nslit)
	DOUBLE PRECISION J_PHI(ntheta*nslit,ndof),
     &                   J_PSI(ntheta*nslit,ndof),
     &	                 J_Vx(ntheta*nslit,ndof),
     &                   J_Vy(ntheta*nslit,ndof)
	DOUBLE PRECISION J_PHI1(ndof),J_PSI1(ndof),
     &	                 J_Vx1(ndof),J_Vy1(ndof),Xdof(ndof)

C       Local variables
	COMPLEX*16 I / (0.0D0, 1.0D0) /
        INTEGER LOCAL_VAROUT(3)
        INTEGER n,e,dof
        INTEGER k,trg_elem,trg_point,row_offset,one
        DOUBLE PRECISION dummy1,dummy2,dummy3,dummy4,zero
	COMPLEX*16 dz,zmin,zmax,zmid,Omega,v0,Ztrg
        
        zero = 0.D0
        one = 1
         
        IF(ndof.NE.(nslit*(2*ncoef+1)))RETURN

        dof = 0
        DO e = 1,nslit
          DO n=1,ncoef 
          dof = dof+1
          Xdof(dof) = ReCn(n,e)
          dof = dof+1
          Xdof(dof) = ImCn(n,e)
          ENDDO
          dof=dof+1
          Xdof(dof) = Q(e) 
        ENDDO

C       Always compute Jacobian for local potential

        LOCAL_VAROUT(2)=1
        IF((varout(1).GT.0).OR.(varout(3).GT.0))THEN
          LOCAL_VAROUT(1)=1
          LOCAL_VAROUT(3)=1
        ELSE
          LOCAL_VAROUT(1)=0
          LOCAL_VAROUT(3)=0 
        ENDIF

        v0 = ReV0 + I*ImV0

C       Loop on all control points

        DO trg_elem = 1,nslit

          zmin = ReZmin(trg_elem) + I*ImZmin(trg_elem)
          zmax = ReZmax(trg_elem) + I*ImZmax(trg_elem)
          zmid = 0.5*(zmin+zmax)
          dz = zmax-zmin 
          
          DO k=1,ntheta

            trg_point = (trg_elem-1)*ntheta+k
!       zeta = cos(theta)
            Ztrg = 0.5*cos(THETA(k))*dz + zmid         
            ReZ(trg_point) = DBLE(Ztrg)
            ImZ(trg_point) = AIMAG(Ztrg)

            CALL SLIT( DBLE(Ztrg),AIMAG(Ztrg),one,
     &                 ReZmin,ImZmin,ReZmax,ImZmax,
     &                 nslit,ReCn,ImCn,ncoef,Q,zero,zero,zero,
     &                 LOCAL_VAROUT,dummy1,dummy2,dummy3,dummy4,
     &                 J_PHI1,J_PSI1,J_Vx1,J_Vy1)

!       Specific expression for potential & velocity created by the element (trg_elem)
!       on which the current target (control) point lies            
!       => update corresponding entries in the Jacobian matrix for Omega,
!          and in the Jacobian matrix for V if needed

            row_offset = (trg_elem-1)*(2*ncoef+1)+1

	    CALL SLITSIDE (THETA(k),one,
     &                     ReZmin(trg_elem),ImZmin(trg_elem),
     &                     ReZmax(trg_elem),ImZmax(trg_elem),
     &                     ReCn(1,trg_elem),ImCn(1,trg_elem),ncoef,
     &                     Q(trg_elem),LOCAL_VAROUT,
     &                     dummy1,dummy2,dummy3,dummy4,
     &                     J_PHI1(row_offset),J_PSI1(row_offset),
     &                     J_Vx1(row_offset),J_Vy1(row_offset) )
            
            Omega = Phi0 - CONJG(v0)*Ztrg
            PHI(trg_point) = DBLE(Omega)
            PSI(trg_point) = AIMAG(Omega)

            IF(varout(1).GT.0)THEN
              Vx(trg_point) = DBLE(v0)
              Vy(trg_point) = AIMAG(v0)
            ENDIF

            DO dof = 1,ndof

! Update complex potential:               
               PHI(trg_point) = PHI(trg_point)+J_PHI1(dof)*Xdof(dof)              
               PSI(trg_point) = PSI(trg_point)+J_PSI1(dof)*Xdof(dof)              

               IF(varout(1).GT.0)THEN
! Update velocity:     
                 Vx(trg_point) = Vx(trg_point)+J_Vx1(dof)*Xdof(dof)              
                 Vy(trg_point) = Vy(trg_point)+J_Vy1(dof)*Xdof(dof)              
               ENDIF

               IF(varout(2).GT.0)THEN
! Store row of Jacobian matrix for complex potential:
                 J_PHI(trg_point,dof) = J_PHI1(dof)              
                 J_PSI(trg_point,dof) = J_PSI1(dof)              
               ENDIF

               IF(varout(3).GT.0)THEN
! Store row of Jacobian matrix for velocity:
                 J_Vx(trg_point,dof) = J_Vx1(dof)              
                 J_Vy(trg_point,dof) = J_Vy1(dof)              
               ENDIF
               
            ENDDO

        ENDDO
        ENDDO

        END

!*******************************************************************************************************************

	SUBROUTINE SLITSIDE (THETA,ntheta,
     &                       ReZmin,ImZmin,ReZmax,ImZmax,
     &                       ReCn,ImCn,ncoef,Q,
     &                       varout,
     &                       PHI,PSI,Vx,Vy,J_PHI,J_PSI,J_Vx,J_Vy)
	!DEC$ ATTRIBUTES DLLEXPORT,ALIAS:'_slitside_' :: SlitSide 

	IMPLICIT NONE
	
	INTEGER ntheta,ncoef
        DOUBLE PRECISION THETA(ntheta)
        DOUBLE PRECISION ReZmin,ImZmin,ReZmax,ImZmax,Q
	DOUBLE PRECISION ReCn(ncoef),ImCn(ncoef)
        INTEGER varout(3)

C       Output arrays

	DOUBLE PRECISION PHI(ntheta),PSI(ntheta),Vx(ntheta),Vy(ntheta)
	DOUBLE PRECISION J_PHI(ntheta,2*ncoef+1),
     &                   J_PSI(ntheta,2*ncoef+1)
	DOUBLE PRECISION J_Vx(ntheta,2*ncoef+1),
     &                   J_Vy(ntheta,2*ncoef+1)

C       Local variables
	COMPLEX*16 I / (0.0D0, 1.0D0) /
        INTEGER k,m,n
        DOUBLE PRECISION L,pi
        COMPLEX*16 rot,dz,jv,jom
        
        pi = 2.D0*AIMAG(LOG(I))

!       rot = exp(i*alpha)
        dz = (ReZmax-ReZmin)+I*(ImZmax-ImZmin)
        L = 0.5*ABS(dz)
        rot = dz/(2.*L)

        DO m=1,ntheta

          PHI(m) = 0.
          PSI(m) = 0.
          Vx(m) = 0.
          Vy(m) = 0.

          DO n=1,ncoef

            jom = cos(DBLE(n)*THETA(m))-I*sin(DBLE(n)*THETA(m))

            J_PHI(m,2*n-1) = DBLE(jom)
            J_PSI(m,2*n-1) = AIMAG(jom)
            J_PHI(m,2*n) = DBLE(I*jom)
            J_PSI(m,2*n) = AIMAG(I*jom)

            PHI(m) = PHI(m)
     &             + J_PHI(m,2*n-1)*ReCn(n) 
     &             + J_PHI(m,2*n)  *ImCn(n) 

            PSI(m) = PSI(m) 
     &             + J_PSI(m,2*n-1)*ReCn(n) 
     &             + J_PSI(m,2*n)  *ImCn(n) 

            jv = DBLE(n)*(I*cos(DBLE(n)*THETA(m))
     &         -sin(DBLE(n)*THETA(m)))*rot/L/sin(THETA(m))

            J_Vx(m,2*n-1) = DBLE(jv)
            J_Vy(m,2*n-1) = AIMAG(jv)
            J_Vx(m,2*n) = DBLE(-I*jv)
            J_Vy(m,2*n) = AIMAG(-I*jv)

            Vx(m) = Vx(m) 
     &             + J_Vx(m,2*n-1)*ReCn(n) 
     &             + J_Vx(m,2*n)  *ImCn(n) 

            Vy(m) = Vy(m) 
     &             + J_Vy(m,2*n-1)*ReCn(n) 
     &             + J_Vy(m,2*n)  *ImCn(n) 

          ENDDO
!         w.r.t. Q:          
          J_PHI(m,2*ncoef+1) = -1.0
          J_PSI(m,2*ncoef+1) = THETA(m)/2./pi
          jv = -I*rot/(2.*pi*L*sin(THETA(m)))
          J_Vx(m,2*ncoef+1) = DBLE(jv)
          J_Vy(m,2*ncoef+1) = AIMAG(jv)

	  PHI(m) = PHI(m) + J_PHI(m,2*ncoef+1)*Q
	  PSI(m) = PSI(m) + J_PSI(m,2*ncoef+1)*Q
	  Vx(m) = Vx(m) + J_Vx(m,2*ncoef+1)*Q
	  Vy(m) = Vy(m) + J_Vy(m,2*ncoef+1)*Q

        ENDDO

 
        END   

!*******************************************************************************************************************

	SUBROUTINE SLIT (ReZ,ImZ,nz,
     &                   ReZmin,ImZmin,ReZmax,ImZmax,nslit,
     &                   ReCn,ImCn,ncoef,Q,ReV0,ImV0,Phi0,
     &                   varout,
     &                   PHI,PSI,Vx,Vy,J_PHI,J_PSI,J_Vx,J_Vy)
	!DEC$ ATTRIBUTES DLLEXPORT,ALIAS:'_slit_' :: Slit 

	IMPLICIT NONE
	
	INTEGER nz,nslit,ncoef
        INTEGER varout(3)
	DOUBLE PRECISION ReZ(nz),ImZ(nz)
        DOUBLE PRECISION ReZmin(nslit),ImZmin(nslit),
     &                   ReZmax(nslit),ImZmax(nslit)
	DOUBLE PRECISION ReCn(ncoef,nslit),ImCn(ncoef,nslit),Q(nslit)
	DOUBLE PRECISION ReV0,ImV0,Phi0

C       Output arrays (double, allocated in C wrapper):

	DOUBLE PRECISION PHI(nz),PSI(nz),Vx(nz),Vy(nz)
	DOUBLE PRECISION J_PHI(nz,nslit*(2*ncoef+1)),
     &                   J_PSI(nz,nslit*(2*ncoef+1))
	DOUBLE PRECISION J_Vx(nz,nslit*(2*ncoef+1)),
     &                   J_Vy(nz,nslit*(2*ncoef+1))
! Les coeff des jacobiens sont à calculer vis-à-vis de tous les dof : Re{Cn},Im{Cn},Q

C	Complex variables
		
	COMPLEX*16 I / (0.0D0, 1.0D0) /
	COMPLEX*16 z,dz,ZZ,zeta,zmin,zmax,Omega,v,jv,tq,v0,
     &             zmid,ZZn,fac,LOGZZ
	COMPLEX*16 rac1,rac2

	DOUBLE PRECISION alpha,pi,L,ReZZ,ImZZ
	INTEGER e,k,n,coljaco

        pi = 2.D0*AIMAG(LOG(I))
        v0 = ReV0 + I*ImV0 

        DO k=1,nz
C       Loop on query points
            
            z = ReZ(k) + I*ImZ(k)
            Omega = Phi0 - CONJG(v0)*z
            v = v0

            DO e=1,nslit	
C	    Loop on slit elements
	        
            	zmin = ReZmin(e) + I*ImZmin(e)
            	zmax = ReZmax(e) + I*ImZmax(e)
            	zmid = 0.5*(zmin+zmax)
            	dz = zmax-zmin
            	alpha = AIMAG(LOG(dz))
            	L = 0.5*ABS(dz)

            	zeta = 2.0D0*(z-zmid)/dz
		rac1 = SQRT(zeta+1.0D0)
            	rac2 = SQRT(zeta-1.0D0)
                
                IF((AIMAG(rac1)*AIMAG(rac2)).LT.(0.D0))THEN 
                	rac2 = CONJG(rac2)
	        ENDIF
                
                ZZ = zeta+rac1*rac2
		LOGZZ = LOG(ZZ)

C               Denominator of (A5) :
                fac = rac1*rac2*L*EXP(I*alpha)
                fac = CONJG(1.0D0/fac)
  
    	    	DO n=1,ncoef
C           	Loop on the coefficients of the Laurent serie
               		
                        ZZn = EXP(-DCMPLX(n)*LOGZZ) 
			Omega = Omega + (ReCn(n,e)+I*ImCn(n,e))*ZZn

                        jv = n*fac*CONJG(ZZn)
 
			IF (varout(1).GT.0) THEN
                              v = v + jv * (ReCn(n,e)-I*ImCn(n,e))
                        ENDIF
                        
                        IF (varout(2).GT.0) THEN
C                             Derivatives w.r.t. Re{Cn}  
    			      coljaco = (e-1)*(2*ncoef+1)+2*n-1
                              J_PHI(k,coljaco) = DBLE(ZZn)
                              J_PSI(k,coljaco) = AIMAG(ZZn)                    
C                             Derivatives w.r.t. Im{Cn}  
                              coljaco = (e-1)*(2*ncoef+1)+2*n
                              J_PHI(k,coljaco) = DBLE(I*ZZn)
                              J_PSI(k,coljaco) = AIMAG(I*ZZn)                    
                        ENDIF

                        IF (varout(3).GT.0) THEN
C                             Derivatives w.r.t. Re{Cn}  
    			      coljaco = (e-1)*(2*ncoef+1)+2*n-1
                              J_Vx(k,coljaco) = DBLE(jv)
                              J_Vy(k,coljaco) = AIMAG(jv)
C                             Derivatives w.r.t. Im{Cn}       
    			      coljaco = (e-1)*(2*ncoef+1)+2*n   
                              J_Vx(k,coljaco) = DBLE(-I*jv)
                              J_Vy(k,coljaco) = AIMAG(-I*jv)     
                        ENDIF
	
	    	ENDDO

C               Contribution of term Q
                
                tq = LOGZZ/2.0D0/pi - 1.0D0
                Omega = Omega + Q(e) * tq

C               Potential & Stream function
                PHI(k) = DBLE(Omega)		
                PSI(k) = AIMAG(Omega)		
                
		IF (varout(1).GT.0) THEN
                      v = v - fac*Q(e)/2.0D0/pi
                      Vx(k) = DBLE(v)		
                      Vy(k) = AIMAG(v)	
                ENDIF
        
		coljaco = (e-1)*(2*ncoef+1)+2*ncoef+1

                IF (varout(2).GT.0) THEN
                      J_PHI(k,coljaco) = DBLE(tq)
                      J_PSI(k,coljaco) = AIMAG(tq)                    
                ENDIF

                IF (varout(3).GT.0) THEN
                      J_Vx(k,coljaco) = -1.0D0*DBLE(fac)/2.0D0/pi
                      J_Vy(k,coljaco) = -1.0D0*AIMAG(fac)/2.0D0/pi      
                ENDIF

            ENDDO
                        
        ENDDO


	END