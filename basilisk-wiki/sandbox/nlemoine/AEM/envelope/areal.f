	SUBROUTINE AREALSINK ( Xvert,Yvert,nvert,   
     &                         Xquery,Yquery,nquery,gamma0,
! output arrays:
     &                         PHI,Vx,Vy,A)
	!DEC$ ATTRIBUTES DLLEXPORT,ALIAS:'_arealsink_' :: ArealSink 

	IMPLICIT NONE
	
	INTEGER nvert,nquery
        DOUBLE PRECISION Xvert(nvert),Yvert(nvert)
        DOUBLE PRECISION Xquery(nquery),Yquery(nquery)
        DOUBLE PRECISION gamma0

C       Output arrays (double, allocated in C wrapper):

	DOUBLE PRECISION PHI(nquery),Vx(nquery),Vy(nquery)
        DOUBLE PRECISION A

C       Local variables

	COMPLEX*16 I / (0.0D0, 1.0D0) /
        INTEGER q,j
        COMPLEX*16 z,W,nuj,Zj,Cj,smallfj,Fj
	DOUBLE PRECISION Lj2,Yj,pi
 
        pi = 2.D0*AIMAG(LOG(I))
        A = 0.D0

C       Loop on all query points

        DO q = 1,nquery

          z = Xquery(q)+I*Yquery(q)
          PHI(q) = 0.D0
          W = 0.D0 + I*0.D0
          
          DO j=1,nvert

            IF(j.lt.nvert)THEN
              nuj = Xvert(j+1)-Xvert(j)+I*(Yvert(j+1)-Yvert(j))
            ELSE
              nuj = Xvert(1)-Xvert(j)+I*(Yvert(1)-Yvert(j))
            ENDIF

            Lj2 = ABS(nuj)
            Lj2 = Lj2*Lj2

! Compute local coordinate of current query point z,
! using Eq. (8.562) p.336 defining (Zj)+1 hence Zj 
            Zj = 2.0*(z-Xvert(j)-I*Yvert(j))/nuj - 1.0
            Yj = AIMAG(Zj);
! Constant Cj Eq. (8.571) p.337
            Cj = -2.0*LOG(0.5*nuj)
! Function f(Zj) (small f) Eq. (8.564) p.336
            smallfj = LOG((Zj-1.0)/(Zj+1.0))
! Function F(Zj) (capital F) Eq. (8.585) p.341
            Fj = Zj*smallfj-LOG(Zj-1.0)-LOG(Zj+1.0)+2.0+Cj
! Update sum in Eq. (8.576) p.338 defining W, with shortcuts:
!    (a) conj(Zj)-Zj = -2*I*Im{Zj} = -2*I*Yj
!    (b) F(Zj)+F(conj(Zj)) = F(Zj) + conj(F(Zj)) = 2*Re{F(Zj)}
            W = W + CONJG(nuj)*(-2.0*I*Yj*smallfj-2.0*DBLE(Fj))
! Update sum in Eq. (8.581) p.338 defining Phi, real-valued potential
            PHI(q) = PHI(q) + 2.0*Yj*Lj2*DBLE(Fj)

            IF(q.EQ.1)THEN
              A = A + 0.25*Yj*Lj2
            ENDIF

          ENDDO

! Apply multiplicating factor for Phi, Eq. (8.581)
          PHI(q) = -(gamma0/(32.*pi))*PHI(q)
! Apply multiplicating factor for W, Eq. (8.576)
          W = -(gamma0/(16.*pi*I))*W
! Get real-valued velocity components of V = conj(W)
          Vx(q) = DBLE(W)
          Vy(q) = -AIMAG(W)

        ENDDO

        END


        SUBROUTINE INPOLYGON( Xvert,Yvert,nvert,   
     &                        Xquery,Yquery,nquery,
! Output array
     &                        RES)
	!DEC$ ATTRIBUTES DLLEXPORT,ALIAS:'_inpolygon_' :: InPolygon 

	IMPLICIT NONE
	
	INTEGER nvert,nquery
        DOUBLE PRECISION Xvert(nvert),Yvert(nvert)
        DOUBLE PRECISION Xquery(nquery),Yquery(nquery)
!        INTEGER RES(nquery)       
        DOUBLE PRECISION RES(nquery)       

C       Local variables

	COMPLEX*16 I / (0.0D0, 1.0D0) /
        INTEGER q,j
        COMPLEX*16 z,zv1,zvnext,dz,dznext
        DOUBLE PRECISION pi,sumArg

        pi = 2.D0*AIMAG(LOG(I))
        zv1 = Xvert(1)+I*Yvert(1)

        DO q=1,nquery

          z = Xquery(q)+I*Yquery(q)
          dz = z-zv1
          sumArg = 0.D0

          DO j=1,nvert
              
            IF(j.LT.nvert)THEN
              zvnext = Xvert(j+1)+I*Yvert(j+1)
            ELSE
              zvnext = Xvert(1)+I*Yvert(1)
            ENDIF

            dznext = z-zvnext
!           Use of Eq. (8.552) p.334 with shortcut
!           ln((z-z_{j+1})/(z-z_j)) - conj(ln((z-z_{j+1})/(z-z_j))) = 2*Im{ln(dznext/dz)}
            sumArg = sumArg + AIMAG(LOG(dznext/dz))  
            dz = dznext
    
          ENDDO   
   
!          RES(q) = IDINT(sumArg/2.D0/pi)        
          RES(q) = sumArg/2.D0/pi      

        ENDDO
  
        END
