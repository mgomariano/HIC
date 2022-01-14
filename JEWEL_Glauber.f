      program JEWEL_Glauber
      end program JEWEL_Glauber


      SUBROUTINE MEDINIT(FILE,id,etam,mass)
      IMPLICIT NONE
C--medium parameters
      COMMON/MEDPARAM/CENTRMIN,CENTRMAX,BREAL,CENTR,RAU,NF
      INTEGER NF
      DOUBLE PRECISION CENTRMIN,CENTRMAX,BREAL,CENTR,RAU
      COMMON/MEDPARAMINT/TAUI,TI,TC,D3,ZETA3,D,
     &N0,SIGMANN,A,WOODSSAXON
      DOUBLE PRECISION TAUI,TI,TC,ALPHA,BETA,GAMMA,D3,ZETA3,D,N0,
     &SIGMANN
      INTEGER A
      LOGICAL WOODSSAXON
C--max rapidity
	common/rapmax2/etamax2
	double precision etamax2
C--longitudinal boost of momentum distribution
	common/boostmed/boost
	logical boost
C--factor to vary Debye mass
	COMMON/MDFAC/MDFACTOR,MDSCALEFAC
	DOUBLE PRECISION MDFACTOR,MDSCALEFAC
C--nuclear thickness function
      COMMON /THICKFNC/ RMAX,TA(100,2)
      DOUBLE PRECISION RMAX,TA
C--geometrical cross section
      COMMON /CROSSSEC/ IMPMAX,CROSS(200,3)
      DOUBLE PRECISION IMPMAX,CROSS
C--identifier of log file
	common/logfile/logfid
	integer logfid

      DATA RAU/10./
      DATA D3/0.9d0/
      DATA ZETA3/1.2d0/
C--local variables
      INTEGER I,LUN,POS,IOS,id,mass
	double precision etam
      CHARACTER*100 BUFFER,LABEL,tempbuf
	CHARACTER*80 FILE
	character firstchar
	logical fileexist

	etamax2 = etam
	logfid = id

      IOS=0
      LUN=77

C--default settings
      TAUI=0.6d0
      TI=0.36d0
      TC=0.17d0
      WOODSSAXON=.TRUE.
      CENTRMIN=0.d0
      CENTRMAX=100.d0
      NF=3
      A=mass
      N0=0.17d0
      D=0.513d0
      SIGMANN=6.2
	MDFACTOR=0.45d0
	MDSCALEFAC=0.9d0
	boost = .true.

C--read settings from file
	write(logfid,*)
	inquire(file=FILE,exist=fileexist)
	if(fileexist)then
        write(logfid,*)'Reading medium parameters from ',FILE
        OPEN(unit=LUN,file=FILE,status='old',err=10)
	  do 20 i=1,1000
          READ(LUN, '(A)', iostat=ios) BUFFER
	    if (ios.ne.0) goto 30
	    firstchar = buffer(1:1)
	    if (firstchar.eq.'#') goto 20
          POS=SCAN(BUFFER,' ')
          LABEL=BUFFER(1:POS)
          BUFFER=BUFFER(POS+1:)
          IF (LABEL=="TAUI")THEN
            READ(BUFFER,*,IOSTAT=IOS) TAUI
          ELSE IF (LABEL=="TI") THEN
            READ(BUFFER,*,IOSTAT=IOS) TI
          ELSE IF (LABEL=="TC") THEN
            READ(BUFFER,*,IOSTAT=IOS) TC
          ELSE IF (LABEL=="WOODSSAXON") THEN
            READ(BUFFER,*,IOSTAT=IOS) WOODSSAXON
          ELSE IF (LABEL=="CENTRMIN") THEN
            READ(BUFFER,*,IOSTAT=IOS) CENTRMIN
          ELSE IF (LABEL=="CENTRMAX") THEN
            READ(BUFFER,*,IOSTAT=IOS) CENTRMAX
          ELSE IF (LABEL=="NF") THEN
            READ(BUFFER,*,IOSTAT=IOS) NF
          ELSE IF (LABEL=="N0") THEN
            READ(BUFFER,*,IOSTAT=IOS) N0
          ELSE IF (LABEL=="D") THEN
            READ(BUFFER,*,IOSTAT=IOS) D
          ELSE IF (LABEL=="SIGMANN") THEN
            READ(BUFFER,*,IOSTAT=IOS) SIGMANN
          ELSE IF (LABEL=="MDFACTOR") THEN
            READ(BUFFER,*,IOSTAT=IOS) MDFACTOR
          ELSE IF (LABEL=="MDSCALEFAC") THEN
            READ(BUFFER,*,IOSTAT=IOS) MDSCALEFAC
	    else
	      write(logfid,*)'unknown label ',label
	    endif
 20	  continue

 30	  close(LUN,status='keep')
	  write(logfid,*)'...done'
	  goto 40

 10     write(logfid,*)'Could not open medium parameter file, '//
     &	'will run with default settings.'

	else
	  write(logfid,*)'No medium parameter file found, '//
     &	'will run with default settings.'
	endif

 40   write(logfid,*)'using parameters:'
      write(logfid,*)'TAUI       =',TAUI
      write(logfid,*)'TI         =',TI
      write(logfid,*)'TC         =',TC
      write(logfid,*)'WOODSSAXON =',WOODSSAXON
      write(logfid,*)'CENTRMIN   =',CENTRMIN
      write(logfid,*)'CENTRMAX   =',CENTRMAX
      write(logfid,*)'NF         =',NF
      write(logfid,*)'A          =',A
      write(logfid,*)'N0         =',N0
      write(logfid,*)'D          =',D
      write(logfid,*)'SIGMANN    =',SIGMANN
      write(logfid,*)'MDFACTOR   =',MDFACTOR
      write(logfid,*)'MDSCALEFAC =',MDSCALEFAC
	write(logfid,*)
	write(logfid,*)
	write(logfid,*)

C--calculate T_A(x,y)
      CALL CALCTA
C--calculate geometrical cross section
      CALL CALCXSECTION

      END



      SUBROUTINE MEDNEXTEVT
      IMPLICIT NONE
C--medium parameters
      COMMON/MEDPARAM/CENTRMIN,CENTRMAX,BREAL,CENTR,RAU,NF
      INTEGER NF
      DOUBLE PRECISION CENTRMIN,CENTRMAX,BREAL,CENTR,RAU
      COMMON/MEDPARAMINT/TAUI,TI,TC,D3,ZETA3,D,
     &N0,SIGMANN,A,WOODSSAXON
      DOUBLE PRECISION TAUI,TI,TC,ALPHA,BETA,GAMMA,D3,ZETA3,D,N0,
     &SIGMANN
      INTEGER A
      LOGICAL WOODSSAXON
C--geometrical cross section
      COMMON /CROSSSEC/ IMPMAX,CROSS(200,3)
      DOUBLE PRECISION IMPMAX,CROSS
C--local variables
      integer i,j
      DOUBLE PRECISION PYR,R,b1,b2,gettemp

C--pick an impact parameter
      r=(pyr(0)*(centrmax-centrmin)+centrmin)/100.
      i=0
      do 130 j=1,200
       if ((r-cross(j,3)/cross(200,3)).ge.0.) then
        i=i+1
       else 
        goto 132
       endif
 130  continue
 132  continue
      b1 = (i-1)*0.1d0
      b2 = i*0.1d0
      breal = (b2*(cross(i,3)/cross(200,3)-r)
     &      +b1*(r-cross(i+1,3)/cross(200,3)))/
     &	(cross(i,3)/cross(200,3)-cross(i+1,3)/cross(200,3))
      centr = r;
      END

      SUBROUTINE SETCENTRALITY(MIN, MAX)
      IMPLICIT NONE
      COMMON/MEDPARAM/CENTRMIN,CENTRMAX,BREAL,CENTR,RAU,NF
      INTEGER NF
      DOUBLE PRECISION CENTRMIN,CENTRMAX,BREAL,CENTR,RAU
      DOUBLE PRECISION MIN,MAX
      CENTRMAX = MAX
      CENTRMIN = MIN
      END

      double precision function getb()
      implicit none
      COMMON/MEDPARAM/CENTRMIN,CENTRMAX,BREAL,CENTR,RAU,NF
      INTEGER NF
      DOUBLE PRECISION CENTRMIN,CENTRMAX,BREAL,CENTR,RAU
      getb=breal
      end


      double precision function getcentrality()
      implicit none
      COMMON/MEDPARAM/CENTRMIN,CENTRMAX,BREAL,CENTR,RAU,NF
      INTEGER NF
      DOUBLE PRECISION CENTRMIN,CENTRMAX,BREAL,CENTR,RAU
      getcentrality=centr
      end

      SUBROUTINE PICKVTX(X,Y)
      IMPLICIT NONE
      DOUBLE PRECISION X,Y
C--medium parameters
      COMMON/MEDPARAM/CENTRMIN,CENTRMAX,BREAL,CENTR,RAU,NF
      INTEGER NF
      DOUBLE PRECISION CENTRMIN,CENTRMAX,BREAL,CENTR,RAU
C--local variables
      DOUBLE PRECISION X1,X2,Y1,Y2,Z,XVAL,YVAL,ZVAL,NTHICK,PYR

      X1=BREAL/2.-RAU
      X2=RAU-BREAL/2.
      Y1=-SQRT(4*RAU**2-BREAL**2)/2.
      Y2=SQRT(4*RAU**2-BREAL**2)/2.
 131  XVAL=PYR(0)*(X2-X1)+X1
      YVAL=PYR(0)*(Y2-Y1)+Y1
      IF((NTHICK(XVAL-BREAL/2.,YVAL).EQ.0.d0).OR.
     &     NTHICK(XVAL+BREAL/2.,YVAL).EQ.0.d0) GOTO 131
      ZVAL=PYR(0)*NTHICK(-BREAL/2.,0d0)*NTHICK(BREAL/2.,0d0)
      Z=NTHICK(XVAL-BREAL/2.,YVAL)*NTHICK(XVAL+BREAL/2.,YVAL)
      IF(ZVAL.GT.Z) GOTO 131
      X=XVAL
      Y=YVAL
      END

	SUBROUTINE SETB(BVAL)
	IMPLICIT NONE
C--medium parameters
      COMMON/MEDPARAM/CENTRMIN,CENTRMAX,BREAL,CENTR,RAU,NF
      INTEGER NF
      DOUBLE PRECISION CENTRMIN,CENTRMAX,BREAL,CENTR,RAU
	DOUBLE PRECISION BVAL
	BREAL=BVAL
	END



      SUBROUTINE GETSCATTERER(X,Y,Z,T,TYPE,PX,PY,PZ,E,MS)
      IMPLICIT NONE
C--medium parameters
      COMMON/MEDPARAM/CENTRMIN,CENTRMAX,BREAL,CENTR,RAU,NF
      INTEGER NF
      DOUBLE PRECISION CENTRMIN,CENTRMAX,BREAL,CENTR,RAU
C--internal medium parameters
      COMMON/MEDPARAMINT/TAUI,TI,TC,D3,ZETA3,D,
     &N0,SIGMANN,A,WOODSSAXON
      DOUBLE PRECISION TAUI,TI,TC,ALPHA,BETA,GAMMA,D3,ZETA3,D,N0,
     &SIGMANN
      INTEGER A
      LOGICAL WOODSSAXON
C--longitudinal boost of momentum distribution
	common/boostmed/boost
	logical boost
C--function calls
      DOUBLE PRECISION GETTEMP,GETMD,GETMOM,GETMS
C--identifier of log file
	common/logfile/logfid
	integer logfid
C--local variables
      DOUBLE PRECISION X,Y,Z,T,MS,PX,PY,PZ,E,MD,TEMP
      INTEGER TYPE
      DOUBLE PRECISION R,PYR,pmax,wt,tau,theta,phi,pi,p,ys,pz2,e2
      DATA PI/3.141592653589793d0/

      R=PYR(0)
      IF(R.LT.(2.*12.*NF*D3/3.)/(2.*12.*NF*D3/3.+3.*16.*ZETA3/2.))THEN
         TYPE=2
      ELSE
         TYPE=21
      ENDIF
      MS=GETMS(X,Y,Z,T)
      MD=GETMD(X,Y,Z,T)
      TEMP=GETTEMP(X,Y,Z,T)
	tau=sqrt(t**2-z**2)
	if (boost) then
  	  ys = 0.5*log((t+z)/(t-z))
	else
	  ys = 0.d0
	endif
	pmax = 10.*temp

      IF(TEMP.LT.1.D-2)THEN
       write(logfid,*)'asking for a scattering centre without medium:'
       write(logfid,*)'at (x,y,z,t)=',X,Y,Z,T
       write(logfid,*)'making one up to continue but '//
     &	'something is wrong!'
       TYPE=21
       PX=0.d0
       PY=0.d0
       PZ=0.d0
       MS=GETMS(0.d0,0.d0,0.d0,0.d0)
       MD=GETMD(0.d0,0.d0,0.d0,0.d0)
       E=SQRT(PX**2+PY**2+PZ**2+MS**2)
       RETURN
      ENDIF

 10	p = pyr(0)**0.3333333*pmax
	E2 = sqrt(p**2+ms**2)
	if (type.eq.2) then
	  wt = (exp(ms/temp)-1.)/(exp(E2/temp)-1.)
	else
	  wt = (exp(ms/temp)+1.)/(exp(E2/temp)+1.)
	endif
	if (wt.gt.1.) write(logfid,*)'Error in getscatterer: weight = ',wt
	if (wt.lt.0.) write(logfid,*)'Error in getscatterer: weight = ',wt
	if (pyr(0).gt.wt) goto 10
	phi = pyr(0)*2.*pi
	theta = -acos(2.*pyr(0)-1.)+pi
	px  = p*sin(theta)*cos(phi)
	py  = p*sin(theta)*sin(phi)
	pz2 = p*cos(theta)
	E   = cosh(ys)*E2 + sinh(ys)*pz2
	pz  = sinh(ys)*E2 + cosh(ys)*pz2
      END


      SUBROUTINE AVSCATCEN(X,Y,Z,T,PX,PY,PZ,E,m)
      IMPLICIT NONE
C--longitudinal boost of momentum distribution
	common/boostmed/boost
	logical boost
C--max rapidity
	common/rapmax2/etamax2
	double precision etamax2
C--local variables
	double precision x,y,z,t,px,py,pz,e,getms,m,ys
	if (boost) then
  	  ys = 0.5*log((t+z)/(t-z))
	  if ((z.eq.0.d0).and.(t.eq.0.d0)) ys =0.d0
	  if (ys.gt.etamax2) ys=etamax2
	  if (ys.lt.-etamax2) ys=-etamax2
	else
	  ys = 0.d0
	endif
	m  = getms(x,y,z,t)
	e  = m*cosh(ys)
	px = 0.d0
	py = 0.d0
	pz = m*sinh(ys)
	end


      SUBROUTINE maxscatcen(PX,PY,PZ,E,m)
      IMPLICIT NONE
C--longitudinal boost of momentum distribution
	common/boostmed/boost
	logical boost
C--max rapidity
	common/rapmax2/etamax2
	double precision etamax2
C--local variables
	double precision px,py,pz,e,getmsmax,m,ys
	if (boost) then
  	  ys = etamax2
	else
	  ys = 0.d0
	endif
	m  = getmsmax()
	e  = m*cosh(ys)
	px = 0.d0
	py = 0.d0
	pz = m*sinh(ys)
	end
	


      DOUBLE PRECISION FUNCTION GETMD(X1,Y1,Z1,T1)
      IMPLICIT NONE
C--factor to vary Debye mass
	COMMON/MDFAC/MDFACTOR,MDSCALEFAC
	DOUBLE PRECISION MDFACTOR,MDSCALEFAC
      DOUBLE PRECISION X1,Y1,Z1,T1,GETTEMP
      GETMD=MDSCALEFAC*3.*GETTEMP(X1,Y1,Z1,T1)
      GETMD=MAX(GETMD,MDFACTOR)
      END



      DOUBLE PRECISION FUNCTION GETMS(X2,Y2,Z2,T2)
      IMPLICIT NONE
      DOUBLE PRECISION X2,Y2,Z2,T2,GETMD
      GETMS=GETMD(X2,Y2,Z2,T2)/SQRT(2.)
      END



      DOUBLE PRECISION FUNCTION GETNEFF(X3,Y3,Z3,T3)
      IMPLICIT NONE
      COMMON/MEDPARAM/CENTRMIN,CENTRMAX,BREAL,CENTR,RAU,NF
      INTEGER NF
      DOUBLE PRECISION CENTRMIN,CENTRMAX,BREAL,CENTR,RAU
      COMMON/MEDPARAMINT/TAUI,TI,TC,D3,ZETA3,D,
     &N0,SIGMANN,A,WOODSSAXON
      DOUBLE PRECISION TAUI,TI,TC,ALPHA,BETA,GAMMA,D3,ZETA3,D,N0,
     &SIGMANN
      INTEGER A
      LOGICAL WOODSSAXON
C--   local variables
      DOUBLE PRECISION X3,Y3,Z3,T3,PI,GETTEMP,tau,cosheta
      DATA PI/3.141592653589793d0/
	tau = sqrt(t3**2-z3**2)
	cosheta = t3/tau
      GETNEFF=(2.*6.*NF*D3*2./3. + 16.*ZETA3*3./2.)
     &     *GETTEMP(X3,Y3,Z3,T3)**3/PI**2
	getneff = getneff/cosheta
      END
      
      

      DOUBLE PRECISION FUNCTION GETTEMP(X4,Y4,Z4,T4)
      IMPLICIT NONE
C--medium parameters
      COMMON/MEDPARAM/CENTRMIN,CENTRMAX,BREAL,CENTR,RAU,NF
      INTEGER NF
      DOUBLE PRECISION CENTRMIN,CENTRMAX,BREAL,CENTR,RAU
      COMMON/MEDPARAMINT/TAUI,TI,TC,D3,ZETA3,D,
     &N0,SIGMANN,A,WOODSSAXON
      DOUBLE PRECISION TAUI,TI,TC,ALPHA,BETA,GAMMA,D3,ZETA3,D,N0,
     &SIGMANN
      INTEGER A
      LOGICAL WOODSSAXON
C--max rapidity
	common/rapmax2/etamax2
	double precision etamax2
C--local variables
      DOUBLE PRECISION X4,Y4,Z4,T4,TAU,NPART,EPS0,EPSIN,TEMPIN,PI,
     &NTHICK,ys
      DATA PI/3.141592653589793d0/

      GETTEMP=0.D0

      IF(ABS(Z4).GT.T4)RETURN

      TAU=SQRT(T4**2-Z4**2)
C--check for overlap region
      IF((NTHICK(X4-BREAL/2.,Y4).EQ.0.d0).OR.
     &NTHICK(X4+BREAL/2.,Y4).EQ.0.d0) RETURN

	ys = 0.5*log((t4+z4)/(t4-z4))
	if (abs(ys).gt.etamax2) return
C--determine initial temperature at transverse position
      IF(WOODSSAXON)THEN
         EPS0=(16.*8.+7.*2.*6.*NF)*PI**2*TI**4/240.
         EPSIN=EPS0*NPART(X4-BREAL/2.,Y4,X4+BREAL/2.,Y4)
     &        *PI*RAU**2/(2.*A)
         TEMPIN=(EPSIN*240./(PI**2*(16.*8.+7.*2.*6.*NF)))**0.25
      ELSE
         TEMPIN=TI
      ENDIF
C--calculate temperature if before initial time
      IF(TAU.LE.TAUI)THEN
	 GETTEMP=TEMPIN*TAU/TAUI
      ELSE
C--evolve temperature
       GETTEMP=TEMPIN*(TAUI/TAU)**0.3333
      ENDIF
      IF(GETTEMP.LT.TC) GETTEMP=0.d0
      END



      DOUBLE PRECISION FUNCTION GETTEMPMAX()
      IMPLICIT NONE
C--medium parameters
      COMMON/MEDPARAM/CENTRMIN,CENTRMAX,BREAL,CENTR,RAU,NF
      INTEGER NF
      DOUBLE PRECISION CENTRMIN,CENTRMAX,BREAL,CENTR,RAU
      COMMON/MEDPARAMINT/TAUI,TI,TC,D3,ZETA3,D,
     &N0,SIGMANN,A,WOODSSAXON
      DOUBLE PRECISION TAUI,TI,TC,ALPHA,BETA,GAMMA,D3,ZETA3,D,N0,
     &SIGMANN
      INTEGER A
      LOGICAL WOODSSAXON
C--function call
      DOUBLE PRECISION GETTEMP
      GETTEMPMAX=GETTEMP(0.D0,0.D0,0.D0,TAUI)
      END



      DOUBLE PRECISION FUNCTION GETMDMAX()
      IMPLICIT NONE
C--factor to vary Debye mass
	COMMON/MDFAC/MDFACTOR,MDSCALEFAC
	DOUBLE PRECISION MDFACTOR,MDSCALEFAC
      DOUBLE PRECISION GETTEMPMAX
      GETMDMAX=MDSCALEFAC*3.*GETTEMPMAX()
      GETMDMAX=MAX(GETMDMAX,MDFACTOR)
      END



      DOUBLE PRECISION FUNCTION GETMDMIN()
      IMPLICIT NONE
C--medium parameters
      COMMON/MEDPARAM/CENTRMIN,CENTRMAX,BREAL,CENTR,RAU,NF
      INTEGER NF
      DOUBLE PRECISION CENTRMIN,CENTRMAX,BREAL,CENTR,RAU
      COMMON/MEDPARAMINT/TAUI,TI,TC,D3,ZETA3,D,
     &N0,SIGMANN,A,WOODSSAXON
      DOUBLE PRECISION TAUI,TI,TC,ALPHA,BETA,GAMMA,D3,ZETA3,D,N0,
     &SIGMANN
      INTEGER A
      LOGICAL WOODSSAXON
C--factor to vary Debye mass
	COMMON/MDFAC/MDFACTOR,MDSCALEFAC
	DOUBLE PRECISION MDFACTOR,MDSCALEFAC
      DOUBLE PRECISION GETTEMPMAX
	GETMDMIN=MDSCALEFAC*3.*TC
      GETMDMIN=MAX(GETMDMIN,MDFACTOR)
      END



      DOUBLE PRECISION FUNCTION GETMSMAX()
      IMPLICIT NONE
      DOUBLE PRECISION GETMDMAX,SQRT
      GETMSMAX=GETMDMAX()/SQRT(2.D0)
      END



	DOUBLE PRECISION FUNCTION GETNATMDMIN()
	IMPLICIT NONE
C--medium parameters
      COMMON/MEDPARAM/CENTRMIN,CENTRMAX,BREAL,CENTR,RAU,NF
      INTEGER NF
      DOUBLE PRECISION CENTRMIN,CENTRMAX,BREAL,CENTR,RAU
      COMMON/MEDPARAMINT/TAUI,TI,TC,D3,ZETA3,D,
     &N0,SIGMANN,A,WOODSSAXON
      DOUBLE PRECISION TAUI,TI,TC,ALPHA,BETA,GAMMA,D3,ZETA3,D,N0,
     &SIGMANN
      INTEGER A
      LOGICAL WOODSSAXON
C--max rapidity
	common/rapmax2/etamax2
	double precision etamax2
C--factor to vary Debye mass
	COMMON/MDFAC/MDFACTOR,MDSCALEFAC
	DOUBLE PRECISION MDFACTOR,MDSCALEFAC,PI
      DATA PI/3.141592653589793d0/
C--local variables
	DOUBLE PRECISION T,GETMDMIN
	T=GETMDMIN()/(MDSCALEFAC*3.)
      GETNATMDMIN=(2.*6.*NF*D3*2./3. + 16.*ZETA3*3./2.)
     &     *T**3/PI**2
	END



	DOUBLE PRECISION FUNCTION GETLTIMEMAX()
	IMPLICIT NONE
C--medium parameters
      COMMON/MEDPARAM/CENTRMIN,CENTRMAX,BREAL,CENTR,RAU,NF
      INTEGER NF
      DOUBLE PRECISION CENTRMIN,CENTRMAX,BREAL,CENTR,RAU
      COMMON/MEDPARAMINT/TAUI,TI,TC,D3,ZETA3,D,
     &N0,SIGMANN,A,WOODSSAXON
      DOUBLE PRECISION TAUI,TI,TC,ALPHA,BETA,GAMMA,D3,ZETA3,D,N0,
     &SIGMANN
      INTEGER A
      LOGICAL WOODSSAXON
C--max rapidity
	common/rapmax2/etamax2
	double precision etamax2
C--function call
      DOUBLE PRECISION GETTEMPMAX
	GETLTIMEMAX=TAUI*(GETTEMPMAX()/TC)**3*cosh(etamax2)
	END



      DOUBLE PRECISION FUNCTION GETNEFFMAX()
      IMPLICIT NONE
      COMMON/MEDPARAM/CENTRMIN,CENTRMAX,BREAL,CENTR,RAU,NF
      INTEGER NF
      DOUBLE PRECISION CENTRMIN,CENTRMAX,BREAL,CENTR,RAU
      COMMON/MEDPARAMINT/TAUI,TI,TC,D3,ZETA3,D,
     &N0,SIGMANN,A,WOODSSAXON
      DOUBLE PRECISION TAUI,TI,TC,ALPHA,BETA,GAMMA,D3,ZETA3,D,N0,
     &SIGMANN
      INTEGER A
      LOGICAL WOODSSAXON
C--max rapidity
	common/rapmax2/etamax2
	double precision etamax2
C--   local variables
      DOUBLE PRECISION PI,GETTEMPMAX
      DATA PI/3.141592653589793d0/
      GETNEFFMAX=(2.*6.*NF*D3*2./3. + 16.*ZETA3*3./2.)
     &     *GETTEMPMAX()**3/PI**2
      END
      
      

      DOUBLE PRECISION FUNCTION NPART(XX1,YY1,XX2,YY2)
      IMPLICIT NONE
      COMMON/MEDPARAMINT/TAUI,TI,TC,D3,ZETA3,D,
     &N0,SIGMANN,A,WOODSSAXON
      DOUBLE PRECISION TAUI,TI,TC,ALPHA,BETA,GAMMA,D3,ZETA3,D,N0,
     &SIGMANN
      INTEGER A
      LOGICAL WOODSSAXON
C--local variables
      DOUBLE PRECISION XX1,YY1,XX2,YY2,NTHICK
      
      NPART = NTHICK(XX1,YY1)*(1.-EXP(-SIGMANN*NTHICK(XX2,YY2))) +
     &        NTHICK(XX2,YY2)*(1.-EXP(-SIGMANN*NTHICK(XX1,YY1)))
      END
      DOUBLE PRECISION FUNCTION NTHICK(X1,Y1)
      IMPLICIT NONE
C--medium parameters
      COMMON/MEDPARAM/CENTRMIN,CENTRMAX,BREAL,CENTR,RAU,NF
      INTEGER NF
      DOUBLE PRECISION CENTRMIN,CENTRMAX,BREAL,CENTR,RAU
      COMMON/MEDPARAMINT/TAUI,TI,TC,D3,ZETA3,D,
     &N0,SIGMANN,A,WOODSSAXON
      DOUBLE PRECISION TAUI,TI,TC,ALPHA,BETA,GAMMA,D3,ZETA3,D,N0,
     &SIGMANN
      INTEGER A
      LOGICAL WOODSSAXON
C--identifier of log file
	common/logfile/logfid
	integer logfid
C--nuclear thickness function
      COMMON /THICKFNC/ RMAX,TA(100,2)
      DOUBLE PRECISION RMAX,TA
      INTEGER LINE,LMIN,LMAX,I
      DOUBLE PRECISION X1,Y1,XA(4),YA(4),Y,DY,R,C,B,DELTA
  
      R=SQRT(X1**2+Y1**2)
      IF(R.GT.TA(100,1))THEN
	 NTHICK=0.
      ELSE
	 LINE=INT(R*99.d0/TA(100,1)+1)
	 LMIN=MAX(LINE,1)
	 LMIN=MIN(LMIN,99)
	 IF((R.LT.TA(LMIN,1)).OR.(R.GT.TA(LMIN+1,1)))
     &	write(logfid,*)LINE,LMIN,R,TA(LMIN,1),TA(LMIN+1,1)
	 XA(1)=TA(LMIN,1)
	 XA(2)=TA(LMIN+1,1)
	 YA(1)=TA(LMIN,2)
	 YA(2)=TA(LMIN+1,2)
	 C=(YA(2)-YA(1))/(XA(2)-XA(1))
	 B=YA(1)-C*XA(1)
	 NTHICK=C*R+B
      ENDIF
      END



      SUBROUTINE CALCTA()
      IMPLICIT NONE
C--medium parameters
      COMMON/MEDPARAM/CENTRMIN,CENTRMAX,BREAL,CENTR,RAU,NF
      INTEGER NF
      DOUBLE PRECISION CENTRMIN,CENTRMAX,BREAL,CENTR,RAU
      COMMON/MEDPARAMINT/TAUI,TI,TC,D3,ZETA3,D,
     &N0,SIGMANN,A,WOODSSAXON
      DOUBLE PRECISION TAUI,TI,TC,ALPHA,BETA,GAMMA,D3,ZETA3,D,N0,
     &SIGMANN
      INTEGER A
      LOGICAL WOODSSAXON
C--   nuclear thickness function
      COMMON /THICKFNC/ RMAX,TA(100,2)
      DOUBLE PRECISION RMAX,TA
C--variables for integration
      COMMON/INTEG/B,R
      DOUBLE PRECISION B,R
C--local variables
      INTEGER NSTEPS,I
      DOUBLE PRECISION EPS,HFIRST,Y

      NSTEPS=100
      EPS=1.E-4
      HFIRST=0.1D0

      R=1.12*A**(0.33333)-0.86*A**(-0.33333)
c      R=2.608;
      RMAX=2.*R

      DO 10 I=1,NSTEPS
C--set transverse position
       B=(I-1)*2.D0*R/NSTEPS
       Y=0.D0
C--integrate along longitudinal line
       CALL ODEINT(Y,-2*R,2*R,EPS,HFIRST,0.d0,101)
       TA(I,1)=B
       TA(I,2)=Y
 10   CONTINUE
      END



      SUBROUTINE CALCXSECTION()
      IMPLICIT NONE
C--medium parameters
      COMMON/MEDPARAM/CENTRMIN,CENTRMAX,BREAL,CENTR,RAU,NF
      INTEGER NF
      DOUBLE PRECISION CENTRMIN,CENTRMAX,BREAL,CENTR,RAU
      COMMON/MEDPARAMINT/TAUI,TI,TC,D3,ZETA3,D,
     &N0,SIGMANN,A,WOODSSAXON
      DOUBLE PRECISION TAUI,TI,TC,ALPHA,BETA,GAMMA,D3,ZETA3,D,N0,
     &SIGMANN
      INTEGER A
      LOGICAL WOODSSAXON
C--   geometrical cross section
      COMMON /CROSSSEC/ IMPMAX,CROSS(200,3)
      DOUBLE PRECISION IMPMAX,CROSS
C--local variables
      INTEGER IX,IY,IB
      DOUBLE PRECISION B,P,PROD,X,Y,NTHICK,NPART,pprev

      pprev=0.
      DO 30 IB=1,200
       B=0.1d0*IB
       PROD=1.d0
       DO 10 IX=1,100
        DO 20 IY=1,100
         X=-20.d0+IX*0.4d0
         Y=-20.d0+IY*0.4d0
         PROD=PROD*
     &EXP(-NTHICK(X+B/2.D0,Y)*SIGMANN)**(0.16d0*NTHICK(X-B/2.D0,Y))
 20     CONTINUE
 10    CONTINUE
       P=(1.D0-PROD)*8.8D0/14.D0*B
       CROSS(IB,1)=B
       CROSS(IB,2)=P
       if (ib.eq.1) then
        cross(ib,3)=0.
       else
        cross(ib,3)=cross(ib-1,3)+(p+pprev)/2.*0.1
       endif
       pprev=p
 30   CONTINUE
      IMPMAX=19.95
      END



      DOUBLE PRECISION FUNCTION MEDDERIV(XVAL,W)
      IMPLICIT NONE
      DOUBLE PRECISION XVAL
      INTEGER W
C--medium parameters
      COMMON/MEDPARAMINT/TAUI,TI,TC,D3,ZETA3,D,
     &N0,SIGMANN,A,WOODSSAXON
      DOUBLE PRECISION TAUI,TI,TC,ALPHA,BETA,GAMMA,D3,ZETA3,D,N0,
     &SIGMANN
      INTEGER A
      LOGICAL WOODSSAXON
C--variables for integration
      COMMON/INTEG/B,R
      DOUBLE PRECISION B,R

      IF (W.EQ.1) THEN
C--XVAL corresponds to z-coordinate
       MEDDERIV=N0/(1+EXP((SQRT(B**2+XVAL**2)-R)/D))
      ELSE 
       MEDDERIV=0.D0
      ENDIF
      END

***********************************************************************
***	  function odeint
***********************************************************************
      subroutine odeint(ystart,a,b,eps,h1,hmin,w1)
	implicit none
C--identifier of file for hepmc output and logfile
	common/hepmcid/hpmcfid,logfid
	integer hpmcfid,logfid
C--local variables
	integer nmax,nstep,w1
	double precision ystart,a,b,eps,h1,hmin,x,h,y,dydx,
     &deriv,yscale,hdid,hnew
	data nmax/100000/

	x = a
	y = ystart
	h = sign(h1,b-a)
	do 20 nstep=1,nmax
	  dydx = deriv(x,w1)
	  yscale = abs(y) + abs(h*dydx) + 1.e-25
	  if (((x + h - b)*h).gt.0.) h = b-x
	  call rkstepper(x,y,dydx,h,hdid,hnew,yscale,eps,w1)
	  if ((x - b)*h.ge.0) then
	    ystart = y
	    return
	  endif
	  h = hnew
	  if (abs(h).lt.abs(hmin)) then
	    write(logfid,*)'Error in odeint: stepsize too small',w1
     &	,ystart,a,b,h1
	    return
	  endif	  
 20	continue
	write(logfid,*)'Error in odeint: too many steps',w1
     &	,ystart,a,b,h1
	end



***********************************************************************
***	  function rkstepper
***********************************************************************
	subroutine rkstepper(x,y,dydx,htest,hdid,hnew,yscale,eps,w1)
	implicit none
C--identifier of file for hepmc output and logfile
	common/hepmcid/hpmcfid,logfid
	integer hpmcfid,logfid
C--local variables
	integer w1
	double precision x,y,dydx,htest,hdid,hnew,yscale,eps,
     &yhalf,y1,y2,rk4step,dydxhalf,xnew,delta,err,h,safety, powerdown,
     &powerup,maxup,maxdown,deriv,fac
	logical reject
	data powerdown/0.25/
	data powerup/0.2/
	data safety/0.9/
	data maxdown/10./
	data maxup/5./

	reject = .false.
	h = htest
 10	xnew = x + h
	if (x.eq.xnew) then
	  write(logfid,*)'Error in rkstepper: step size not significant'
	  return
	endif
	yhalf = rk4step(x,y,dydx,h/2.,w1)
	dydxhalf = deriv(x+h/2.,w1)
	y2 = rk4step(x+h/2.,yhalf,dydxhalf,h/2.,w1)
	y1 = rk4step(x,y,dydx,h,w1)
	delta = y2-y1
	err = abs(delta)/(yscale*eps)
	if (err.gt.1.) then
	  reject = .true.
	  fac = max(1./maxdown,safety/err**powerdown)
	  h = h*fac
	  goto 10 
	else
	  if (reject) then
	    hnew = h
	  else
	    fac = min(maxup,safety/err**powerup)
	    hnew = fac*h
	  endif
	  x = xnew
	  y = y2 + delta/15.
	  hdid = h
	endif
	end



***********************************************************************
***	  function rk4step
***********************************************************************
	double precision function rk4step(x,y,dydx,h,w1)
	implicit none
	integer w1
	double precision x,y,dydx,h,k1,k2,k4,yout,deriv
	k1 = h*dydx
	k2 = h*deriv(x+h/2.,w1)
	k4 = h*deriv(x+h,w1)
	yout = y+k1/6.+2.*k2/3.+k4/6.
	rk4step = yout
	end

***********************************************************************
***	  function deriv
***********************************************************************
      DOUBLE PRECISION FUNCTION DERIV(XVAL,W4)
      IMPLICIT NONE
C--Parameter common block
	COMMON/PARAM/Q0,LPS,LQCD,LTIME,SCALEFACM,ANGORD,SCATRECOIL,
     &ALLHAD,compress,NF
      INTEGER NF
	DOUBLE PRECISION Q0,LQCD,LTIME,LPS,SCALEFACM
      LOGICAL ANGORD,SCATRECOIL,ALLHAD,compress
C--variables for splitting function integration
	COMMON/INTSPLITF/QQUAD,FM
	DOUBLE PRECISION QQUAD,FM
C--variables for Sudakov integration
	COMMON/SUDAINT/QA,ZA2,EB,T,INSTATE,TYP
	DOUBLE PRECISION QA,ZA2,EB,T
	CHARACTER*2 TYP
	LOGICAL INSTATE
C--variables for pdf integration
	COMMON/PDFINTV/XMAX,Z
	DOUBLE PRECISION XMAX,Z
C--variables for cross section integration 
	COMMON/XSECV/QLOW,MDX
	DOUBLE PRECISION QLOW,MDX
C--local variables
	INTEGER W4
      DOUBLE PRECISION XVAL,GETSPLITI,PI,ALPHAS,GETINSPLITI,
     &GETINSUDAFAST,SCATPRIMFUNC,PQQ,PQG,PGG,PGQ,
     &MEDDERIV
	DATA PI/3.141592653589793d0/

      DERIV=MEDDERIV(XVAL,W4-100)
      END

C...PYR
C...Generates random numbers uniformly distributed between
C...0 and 1, excluding the endpoints.
 
      FUNCTION PYR(IDUMMY)
 
C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
C...Commonblocks.
      COMMON/PYDATR/MRPY(6),RRPY(100)
      SAVE /PYDATR/
C...Equivalence between commonblock and local variables.
      EQUIVALENCE (MRPY1,MRPY(1)),(MRPY2,MRPY(2)),(MRPY3,MRPY(3)),
     &(MRPY4,MRPY(4)),(MRPY5,MRPY(5)),(MRPY6,MRPY(6)),
     &(RRPY98,RRPY(98)),(RRPY99,RRPY(99)),(RRPY00,RRPY(100))
 
C...Initialize generation from given seed.
      IF(MRPY2.EQ.0) THEN
        IJ=MOD(MRPY1/30082,31329)
        KL=MOD(MRPY1,30082)
        I=MOD(IJ/177,177)+2
        J=MOD(IJ,177)+2
        K=MOD(KL/169,178)+1
        L=MOD(KL,169)
        DO 110 II=1,97
          S=0D0
          T=0.5D0
          DO 100 JJ=1,48
            M=MOD(MOD(I*J,179)*K,179)
            I=J
            J=K
            K=M
            L=MOD(53*L+1,169)
            IF(MOD(L*M,64).GE.32) S=S+T
            T=0.5D0*T
  100     CONTINUE
          RRPY(II)=S
  110   CONTINUE
        TWOM24=1D0
        DO 120 I24=1,24
          TWOM24=0.5D0*TWOM24
  120   CONTINUE
        RRPY98=362436D0*TWOM24
        RRPY99=7654321D0*TWOM24
        RRPY00=16777213D0*TWOM24
        MRPY2=1
        MRPY3=0
        MRPY4=97
        MRPY5=33
      ENDIF
 
C...Generate next random number.
  130 RUNI=RRPY(MRPY4)-RRPY(MRPY5)
      IF(RUNI.LT.0D0) RUNI=RUNI+1D0
      RRPY(MRPY4)=RUNI
      MRPY4=MRPY4-1
      IF(MRPY4.EQ.0) MRPY4=97
      MRPY5=MRPY5-1
      IF(MRPY5.EQ.0) MRPY5=97
      RRPY98=RRPY98-RRPY99
      IF(RRPY98.LT.0D0) RRPY98=RRPY98+RRPY00
      RUNI=RUNI-RRPY98
      IF(RUNI.LT.0D0) RUNI=RUNI+1D0
      IF(RUNI.LE.0D0.OR.RUNI.GE.1D0) GOTO 130
 
C...Update counters. Random number to output.
      MRPY3=MRPY3+1
      IF(MRPY3.EQ.1000000000) THEN
        MRPY2=MRPY2+1
        MRPY3=0
      ENDIF
      PYR=RUNI
 
      RETURN
      END
 
