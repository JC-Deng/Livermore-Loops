C
C***********************************************
      SUBROUTINE KERNEL( TK)
C***********************************************************************
C                                                                      *
C            KERNEL     executes 24 samples of Fortran computation     *
c               TK(1) - total cpu time to execute only the 24 kernels. *
c               TK(2) - total Flops executed by the 24 Kernels         *
C***********************************************************************
C                                                                      *
C     L. L. N. L.   F O R T R A N   K E R N E L S:   M F L O P S       *
C                                                                      *
C   These kernels measure  Fortran  numerical  computation rates for a *
C   spectrum of  CPU-limited  computational  structures.  Mathematical *
C   through-put is measured  in  units  of  millions of floating-point *
C   operations executed per Second, called Mega-Flops/Sec.             *
C                                                                      *
C   This program  measures  a realistic  CPU performance range for the *
C   Fortran programming system  on  a  given day.  The CPU performance *
C   rates depend  strongly  on  the maturity of the Fortran compiler's *
C   ability to translate Fortran code into efficient machine code.     *
C   [ The CPU hardware  capability  apart  from  compiler maturity (or *
C   availability), could be measured (or simulated) by programming the *
C   kernels in assembly  or machine code directly.  These measurements *
C   can also  serve  as a framework for tracking the maturation of the *
C   Fortran compiler during system development.]                       *
C                                                                      *
C     Fonzi's Law: There is not now and there never will be a language *
C                  in which it is the least bit difficult to write     *
C                  bad programs.                                       *
C                                                    F.H.MCMAHON  1972 *
C***********************************************************************
C
C     l1 :=  param-dimension governs the size of most 1-d arrays
C     l2 :=  param-dimension governs the size of most 2-d arrays
C
C     Loop :=  multiple pass control to execute kernel long enough to time.
C     n  :=  DO loop control for each kernel.  Controls are set in subr. SIZES
C
C     ******************************************************************
C
cANSI IMPLICIT  DOUBLE PRECISION (A-H,O-Z)
cIBM  IMPLICIT  REAL*8           (A-H,O-Z)
C
C/      PARAMETER( l1= 1001, l2=  101, l1d= 2*1001 )
C/      PARAMETER( l13=  64, l13h= l13/2, l213= l13+l13h, l813= 8*l13 )
C/      PARAMETER( l14=2048, l16=  75, l416= 4*l16 , l21= 25 )
C/      PARAMETER( kn= 47, kn2= 95, np= 3, ls= 3*47, krs= 24)
C
C
C/      PARAMETER( nk= 47, nl= 3, nr= 8 )
      INTEGER TEST
C
      COMMON /ALPHA/ mk,ik,im,ml,il,Mruns,Nruns,jr,iovec,NPFS(8,3,47)
      COMMON /BETA / tic, TIMES(8,3,47), SEE(5,3,8,3),
     1              TERRS(8,3,47), CSUMS(8,3,47),
     2              FOPN(8,3,47), DOS(8,3,47)
C
      COMMON /SPACES/ ion,j5,k2,k3,MULTI,laps,Loop,m,kr,LP,n13h,ibuf,nx,
     1 L,npass,nfail,n,n1,n2,n13,n213,n813,n14,n16,n416,n21,nt1,nt2,
     2 last,idebug,mpy,mpy100, intbuf(16)
C
      COMMON /SPACER/ A11,A12,A13,A21,A22,A23,A31,A32,A33,
     2                AR,BR,C0,CR,DI,DK,
     3  DM22,DM23,DM24,DM25,DM26,DM27,DM28,DN,E3,E6,EXPMAX,FLX,
     4  Q,QA,R,RI,S,SCALE,SIG,STB5,T,XNC,XNEI,XNM
C
CPFM  COMMON /KAPPA/ iflag1, ikern, statis(100,20), istats(100,20)
C
      COMMON /SPACE0/ TIME(47), CSUM(47), WW(47), WT(47), ticks,
     1                FR(9), TERR1(47), SUMW(7), START,
     2              SKALE(47), BIAS(47), WS(95), TOTAL(47), FLOPN(47),
     3                IQ(7), NPF, NPFS1(47)
C
      COMMON /SPACEI/ WTP(3), MUL(3), ISPAN(47,3), IPASS(47,3)
C
C/      INTEGER    E,F,ZONE
C/      COMMON /ISPACE/ E(l213), F(l213),
C/     1  IX(l1), IR(l1), ZONE(l416)
C/C
C/      COMMON /SPACE1/ U(l1), V(l1), W(l1),
C/     1  X(l1), Y(l1), Z(l1), G(l1),
C/     2  DU1(l2), DU2(l2), DU3(l2), GRD(l1), DEX(l1),
C/     3  XI(l1), EX(l1), EX1(l1), DEX1(l1),
C/     4  VX(l14), XX(l14), RX(l14), RH(l14),
C/     5  VSP(l2), VSTP(l2), VXNE(l2), VXND(l2),
C/     6  VE3(l2), VLR(l2), VLIN(l2), B5(l2),
C/     7  PLAN(l416), D(l416), SA(l2), SB(l2)
C/C
C/      COMMON /SPACE2/ P(4,l813), PX(l21,l2), CX(l21,l2),
C/     1  VY(l2,l21), VH(l2,7), VF(l2,7), VG(l2,7), VS(l2,7),
C/     2  ZA(l2,7)  , ZP(l2,7), ZQ(l2,7), ZR(l2,7), ZM(l2,7),
C/     3  ZB(l2,7)  , ZU(l2,7), ZV(l2,7), ZZ(l2,7),
C/     4  B(l13,l13), C(l13,l13), H(l13,l13),
C/     5  U1(5,l2,2),  U2(5,l2,2),  U3(5,l2,2)
C
C     ******************************************************************
C
C
C/      PARAMETER( l1=   1001, l2=   101, l1d= 2*1001 )
C/      PARAMETER( l13= 64, l13h= 64/2, l213= 64+32, l813= 8*64 )
C/      PARAMETER( l14= 2048, l16= 75, l416= 4*75 , l21= 25)
C
C
care
C
      INTEGER    E,F,ZONE
      COMMON /ISPACE/ E(96), F(96),
     1  IX(1001), IR(1001), ZONE(300)
C
      COMMON /SPACE1/ U(1001), V(1001), W(1001),
     1  X(1001), Y(1001), Z(1001), G(1001),
     2  DU1(101), DU2(101), DU3(101), GRD(1001), DEX(1001),
     3  XI(1001), EX(1001), EX1(1001), DEX1(1001),
     4  VX(1001), XX(1001), RX(1001), RH(2048),
     5  VSP(101), VSTP(101), VXNE(101), VXND(101),
     6  VE3(101), VLR(101), VLIN(101), B5(101),
     7  PLAN(300), D(300), SA(101), SB(101)
C
      COMMON /SPACE2/ P(4,512), PX(25,101), CX(25,101),
     1  VY(101,25), VH(101,7), VF(101,7), VG(101,7), VS(101,7),
     2  ZA(101,7)  , ZP(101,7), ZQ(101,7), ZR(101,7), ZM(101,7),
     3  ZB(101,7)  , ZU(101,7), ZV(101,7), ZZ(101,7),
     4  B(64,64), C(64,64), H(64,64),
     5  U1(5,101,2),  U2(5,101,2),  U3(5,101,2)
C
C     ******************************************************************
C
      DIMENSION     ZX(1023), XZ(1500), TK(6)
      EQUIVALENCE ( ZX(1), Z(1)), ( XZ(1), X(1))
C
C
C//      DIMENSION       E(96), F(96), U(1001), V(1001), W(1001),
C//     1  X(1001), Y(1001), Z(1001), G(1001),
C//     2  DU1(101), DU2(101), DU3(101), GRD(1001), DEX(1001),
C//     3  IX(1001), XI(1001), EX(1001), EX1(1001), DEX1(1001),
C//     4  VX(1001), XX(1001), IR(1001), RX(1001), RH(2048),
C//     5  VSP(101), VSTP(101), VXNE(101), VXND(101),
C//     6  VE3(101), VLR(101), VLIN(101), B5(101),
C//     7  PLAN(300), ZONE(300), D(300), SA(101), SB(101)
C//C
C//      DIMENSION       P(4,512), PX(25,101), CX(25,101),
C//     1  VY(101,25), VH(101,7), VF(101,7), VG(101,7), VS(101,7),
C//     2  ZA(101,7)  , ZP(101,7), ZQ(101,7), ZR(101,7), ZM(101,7),
C//     3  ZB(101,7)  , ZU(101,7), ZV(101,7), ZZ(101,7),
C//     4  B(64,64), C(64,64), H(64,64),
C//     5  U1(5,101,2),  U2(5,101,2),  U3(5,101,2)
C//C
C//C     ******************************************************************
C//C
C//      COMMON /POINT/ ME,MF,MU,MV,MW,MX,MY,MZ,MG,MDU1,MDU2,MDU3,MGRD,
C//     1  MDEX,MIX,MXI,MEX,MEX1,MDEX1,MVX,MXX,MIR,MRX,MRH,MVSP,MVSTP,
C//     2  MVXNE,MVXND,MVE3,MVLR,MVLIN,MB5,MPLAN,MZONE,MD,MSA,MSB,
C//     3  MP,MPX,MCX,MVY,MVH,MVF,MVG,MVS,MZA,MZP,MZQ,MZR,MZM,MZB,MZU,
C//     4  MZV,MZZ,MB,MC,MH,MU1,MU2,MU3
C//C
C//      POINTER  (ME,E), (MF,F), (MU,U), (MV,V), (MW,W),
C//     1         (MX,X), (MY,Y), (MZ,Z), (MG,G),
C//     2         (MDU1,DU1),(MDU2,DU2),(MDU3,DU3),(MGRD,GRD),(MDEX,DEX),
C//     3         (MIX,IX), (MXI,XI), (MEX,EX), (MEX1,EX1), (MDEX1,DEX1),
C//     4         (MVX,VX), (MXX,XX), (MIR,IR), (MRX,RX), (MRH,RH),
C//     5         (MVSP,VSP), (MVSTP,VSTP), (MVXNE,VXNE), (MVXND,VXND),
C//     6         (MVE3,VE3), (MVLR,VLR), (MVLIN,VLIN), (MB5,B5),
C//     7         (MPLAN,PLAN), (MZONE,ZONE), (MD,D), (MSA,SA), (MSB,SB)
C//C
C//      POINTER  (MP,P), (MPX,PX), (MCX,CX),
C//     1         (MVY,VY), (MVH,VH), (MVF,VF), (MVG,VG), (MVS,VS),
C//     2         (MZA,ZA), (MZP,ZP), (MZQ,ZQ), (MZR,ZR), (MZM,ZM),
C//     3         (MZB,ZB), (MZU,ZU), (MZV,ZV), (MZZ,ZZ),
C//     4         (MB,B), (MC,C), (MH,H),
C//     5         (MU1,U1), (MU2,U2), (MU3,U3)
C..      COMMON DUMMY(2000)
C..      LOC(X)  =.LOC.X
C..      IQ8QDSP = 64*LOC(DUMMY)
C
C     ******************************************************************
C
C     STANDARD PRODUCT COMPILER DIRECTIVES MAY BE USED FOR OPTIMIZATION
C
CDIR$ VECTOR
CLLL. OPTIMIZE LEVEL i
CLLL. OPTION INTEGER (7)
CLLL. OPTION ASSERT (NO HAZARD)
CLLL. OPTION NODYNEQV
C
C     ******************************************************************
C       BINARY MACHINES MAY USE THE  AND(P,Q)  FUNCTION IF AVAILABLE
C       IN PLACE OF THE FOLLOWING CONGRUENCE FUNCTION (SEE KERNEL 13, 14)
C                                 IFF:   j= 2**N
c      INTEGER AND
c     IAND(j,k) = AND(j,k)
CLLL. IAND(j,k) = j.INT.k
      MOD2N(i,j)= MOD(i,j)
c      MOD2N(i,j)= IAND(i,j-1)
C                             i  is Congruent to  MOD2N(i,j)   mod(j)
C     ******************************************************************
C
C
C
C
C
      CALL TRACE ('KERNEL  ')
C
      CALL SPACE
C
CPFM       call  OUTPFM( 0, ion)
      mpy   = 1
      mpy100= 1
      L     = 1
      Loop  = 1
      LP    = Loop
      it0   = TEST(0)
CPFM  iflag1= 13579
C
C*******************************************************************************
C***  KERNEL 1      HYDRO FRAGMENT
C*******************************************************************************
C
cdir$ ivdep
 1001    DO 1 k = 1,n
    1       X(k)= Q + Y(k) * (R * ZX(k+10) + T * ZX(k+11))
C
C...................
      IF( TEST(1) .GT. 0) GO TO 1001
C                   we must execute    DO k= 1,n  repeatedly for accurate timing
C
C*******************************************************************************
C***  KERNEL 2      ICCG EXCERPT (INCOMPLETE CHOLESKY - CONJUGATE GRADIENT)
C*******************************************************************************
C
C
 1002     II= n
       IPNTP= 0
  222   IPNT= IPNTP
       IPNTP= IPNTP+II
          II= II/2
           i= IPNTP+1
cdir$ ivdep
c:ibm_dir:ignore recrdeps (x)
C
      DO 2 k= IPNT+2,IPNTP,2
           i= i+1
    2   X(i)= X(k) - V(k) * X(k-1) - V(k+1) * X(k+1)
          IF( II.GT.1) GO TO 222
C
C...................
      IF( TEST(2) .GT. 0) GO TO 1002
C
C*******************************************************************************
C***  KERNEL 3      INNER PRODUCT
C*******************************************************************************
C
C
 1003      Q= 0.000d0
      DO 3 k= 1,n
    3      Q= Q + Z(k) * X(k)
C
C...................
      IF( TEST(3) .GT. 0) GO TO 1003
C
C*******************************************************************************
C***  KERNEL 4      BANDED LINEAR EQUATIONS
C*******************************************************************************
C
              m= (1001-7)/2
             fw= 1.000d-25
C
 1004 DO 404  k= 7,1001,m
             lw= k-6
           temp= XZ(k-1)
cdir$ ivdep
      DO   4  j= 5,n,5
         temp  = temp   - XZ(lw) * Y(j)
    4        lw= lw+1
        XZ(k-1)= Y(5) * temp
 404  CONTINUE
C
C...................
      IF( TEST(4) .GT. 0) GO TO 1004
C
C*******************************************************************************
C***  KERNEL 5      TRI-DIAGONAL ELIMINATION, BELOW DIAGONAL (NO VECTORS)
C*******************************************************************************
C
C
cdir$ novector
 1005 DO 5 i = 2,n
    5    X(i)= Z(i) * (Y(i) - X(i-1))
cdir$ vector
C
C...................
      IF( TEST(5) .GT. 0) GO TO 1005
C
C*******************************************************************************
C***  KERNEL 6      GENERAL LINEAR RECURRENCE EQUATIONS
C*******************************************************************************
C
C
 1006 DO  6  i= 2,n
          W(i)= 0.0100d0
      DO  6  k= 1,i-1
          W(i)= W(i)  + B(i,k) * W(i-k)
    6 CONTINUE
C
C...................
      IF( TEST(6) .GT. 0) GO TO 1006
C
C*******************************************************************************
C***  KERNEL 7      EQUATION OF STATE FRAGMENT
C*******************************************************************************
C
C
cdir$ ivdep
 1007 DO 7 k= 1,n
        X(k)=     U(k  ) + R*( Z(k  ) + R*Y(k  )) +
     .        T*( U(k+3) + R*( U(k+2) + R*U(k+1)) +
     .        T*( U(k+6) + Q*( U(k+5) + Q*U(k+4))))
    7 CONTINUE
C
C...................
      IF( TEST(7) .GT. 0) GO TO 1007
C
C
C*******************************************************************************
C***  KERNEL 8      A.D.I. INTEGRATION
C*******************************************************************************
C
C
 1008          nl1 = 1
               nl2 = 2
                fw= 2.000d0
      DO  8     kx = 2,3
cdir$ ivdep
      DO  8     ky = 2,n
            DU1(ky)=U1(kx,ky+1,nl1)  -  U1(kx,ky-1,nl1)
            DU2(ky)=U2(kx,ky+1,nl1)  -  U2(kx,ky-1,nl1)
            DU3(ky)=U3(kx,ky+1,nl1)  -  U3(kx,ky-1,nl1)
      U1(kx,ky,nl2)=U1(kx,ky,nl1) +A11*DU1(ky) +A12*DU2(ky) +A13*DU3(ky)
     .       + SIG*(U1(kx+1,ky,nl1) -fw*U1(kx,ky,nl1) +U1(kx-1,ky,nl1))
      U2(kx,ky,nl2)=U2(kx,ky,nl1) +A21*DU1(ky) +A22*DU2(ky) +A23*DU3(ky)
     .       + SIG*(U2(kx+1,ky,nl1) -fw*U2(kx,ky,nl1) +U2(kx-1,ky,nl1))
      U3(kx,ky,nl2)=U3(kx,ky,nl1) +A31*DU1(ky) +A32*DU2(ky) +A33*DU3(ky)
     .       + SIG*(U3(kx+1,ky,nl1) -fw*U3(kx,ky,nl1) +U3(kx-1,ky,nl1))
    8 CONTINUE
C
C...................
      IF( TEST(8) .GT. 0) GO TO 1008
C
C*******************************************************************************
C***  KERNEL 9      INTEGRATE PREDICTORS
C*******************************************************************************
C
C
 1009 DO 9  k = 1,n
      PX( 1,k)= DM28*PX(13,k) + DM27*PX(12,k) + DM26*PX(11,k) +
     .          DM25*PX(10,k) + DM24*PX( 9,k) + DM23*PX( 8,k) +
     .          DM22*PX( 7,k) +  C0*(PX( 5,k) +      PX( 6,k))+ PX( 3,k)
    9 CONTINUE
C
C...................
      IF( TEST(9) .GT. 0) GO TO 1009
C
C*******************************************************************************
C***  KERNEL 10     DIFFERENCE PREDICTORS
C*******************************************************************************
C
C
 1010 DO 10  k= 1,n
      AR      =      CX(5,k)
      BR      = AR - PX(5,k)
      PX(5,k) = AR
      CR      = BR - PX(6,k)
      PX(6,k) = BR
      AR      = CR - PX(7,k)
      PX(7,k) = CR
      BR      = AR - PX(8,k)
      PX(8,k) = AR
      CR      = BR - PX(9,k)
      PX(9,k) = BR
      AR      = CR - PX(10,k)
      PX(10,k)= CR
      BR      = AR - PX(11,k)
      PX(11,k)= AR
      CR      = BR - PX(12,k)
      PX(12,k)= BR
      PX(14,k)= CR - PX(13,k)
      PX(13,k)= CR
   10 CONTINUE
C
C...................
      IF( TEST(10) .GT. 0) GO TO 1010
C
C*******************************************************************************
C***  KERNEL 11     FIRST SUM.   PARTIAL SUMS.              (NO VECTORS)
C*******************************************************************************
C
C
 1011     X(1)= Y(1)
cdir$ novector
      DO 11 k = 2,n
   11     X(k)= X(k-1) + Y(k)
cdir$ vector
C
C...................
      IF( TEST(11) .GT. 0) GO TO 1011
C
C*******************************************************************************
C***  KERNEL 12     FIRST DIFF.
C*******************************************************************************
C
C
cdir$ ivdep
 1012 DO 12 k = 1,n
   12     X(k)= Y(k+1) - Y(k)
C
C...................
      IF( TEST(12) .GT. 0) GO TO 1012
C
C*******************************************************************************
C***  KERNEL 13      2-D PIC   Particle In Cell
C*******************************************************************************
C
                fw= 1.000d0
C
 1013 DO  13     k= 1,n
                i1= P(1,k)
                j1= P(2,k)
                i1=       1 + MOD2N(i1,64)
                j1=       1 + MOD2N(j1,64)
            P(3,k)= P(3,k)  + B(i1,j1)
            P(4,k)= P(4,k)  + C(i1,j1)
            P(1,k)= P(1,k)  + P(3,k)
            P(2,k)= P(2,k)  + P(4,k)
                i2= P(1,k)
                j2= P(2,k)
                i2=            MOD2N(i2,64)
                j2=            MOD2N(j2,64)
            P(1,k)= P(1,k)  + Y(i2+32)
            P(2,k)= P(2,k)  + Z(j2+32)
                i2= i2      + E(i2+32)
                j2= j2      + F(j2+32)
          H(i2,j2)= H(i2,j2) + fw
   13 CONTINUE
C
C...................
      IF( TEST(13) .GT. 0) GO TO 1013
C
C*******************************************************************************
C***  KERNEL 14      1-D PIC   Particle In Cell
C*******************************************************************************
C
C
               fw= 1.000d0
C
 1014 DO   141  k= 1,n
            VX(k)= 0.0d0
            XX(k)= 0.0d0
            IX(k)= INT(  GRD(k))
            XI(k)= REAL( IX(k))
           EX1(k)= EX   ( IX(k))
          DEX1(k)= DEX  ( IX(k))
 141  CONTINUE
C
      DO   142  k= 1,n
            VX(k)= VX(k) + EX1(k) + (XX(k) - XI(k))*DEX1(k)
            XX(k)= XX(k) + VX(k)  + FLX
            IR(k)= XX(k)
            RX(k)= XX(k) - IR(k)
            IR(k)= MOD2N(  IR(k),2048) + 1
            XX(k)= RX(k) + IR(k)
 142  CONTINUE
C
      DO  14    k= 1,n
      RH(IR(k)  )= RH(IR(k)  ) + fw - RX(k)
      RH(IR(k)+1)= RH(IR(k)+1) + RX(k)
  14  CONTINUE
C
C...................
      IF( TEST(14) .GT. 0) GO TO 1014
C
C
C
C
C
C
C
C
C
C
C
C
C
C
C
C
C
C
C
C*******************************************************************************
C***  KERNEL 15     CASUAL FORTRAN.  DEVELOPMENT VERSION.
C*******************************************************************************
C
C
C       CASUAL ORDERING OF SCALAR OPERATIONS IS TYPICAL PRACTICE.
C       THIS EXAMPLE DEMONSTRATES THE NON-TRIVIAL TRANSFORMATION
C       REQUIRED TO MAP INTO AN EFFICIENT MACHINE IMPLEMENTATION.
C
C
 1015          NG= 7
               NZ= n
               AR= 0.05300d0
               BR= 0.07300d0
        DO 45  j = 2,NG
        DO 45  k = 2,NZ
               IF( j-NG) 31,30,30
   30     VY(k,j)= 0.0d0
                   GO TO 45
   31          IF( VH(k,j+1) -VH(k,j)) 33,33,32
   32           T= AR
                   GO TO 34
   33           T= BR
   34          IF( VF(k,j) -VF(k-1,j)) 35,36,36
   35           R= MAX( VH(k-1,j), VH(k-1,j+1))
                S= VF(k-1,j)
                   GO TO 37
   36           R= MAX( VH(k,j),   VH(k,j+1))
                S= VF(k,j)
   37     VY(k,j)= SQRT( VG(k,j)**2 +R*R)*T/S
               IF( k-NZ) 40,39,39
   39     VS(k,j)= 0.0d0
                   GO TO 45
   40          IF( VF(k,j) -VF(k,j-1)) 41,42,42
   41           R= MAX( VG(k,j-1), VG(k+1,j-1))
                S= VF(k,j-1)
                T= BR
                   GO TO 43
   42           R= MAX( VG(k,j),   VG(k+1,j))
                S= VF(k,j)
                T= AR
   43     VS(k,j)= SQRT( VH(k,j)**2 +R*R)*T/S
   45    CONTINUE
C
C...................
      IF( TEST(15) .GT. 0) GO TO 1015
C
C
C
C
C
C
C
C
C
C
C
C
C
C
C*******************************************************************************
C***  KERNEL 16     MONTE CARLO SEARCH LOOP
C*******************************************************************************
C
            II= n/3
            LB= II+II
            k2= 0
            k3= 0
C
C
 1016        m= 1
            i1= m
  410       j2= (n+n)*(m-1)+1
      DO 470 k= 1,n
            k2= k2+1
            j4= j2+k+k
            j5= ZONE(j4)
            IF( j5-n      ) 420,475,450
  415       IF( j5-n+II   ) 430,425,425
  420       IF( j5-n+LB   ) 435,415,415
  425       IF( PLAN(j5)-R) 445,480,440
  430       IF( PLAN(j5)-S) 445,480,440
  435       IF( PLAN(j5)-T) 445,480,440
  440       IF( ZONE(j4-1)) 455,485,470
  445       IF( ZONE(j4-1)) 470,485,455
  450       k3= k3+1
            IF( D(j5)-(D(j5-1)*(T-D(j5-2))**2+(S-D(j5-3))**2
     .                        +(R-D(j5-4))**2)) 445,480,440
  455        m= m+1
            IF( m-ZONE(1) ) 465,465,460
  460        m= 1
  465       IF( i1-m) 410,480,410
  470 CONTINUE
  475 CONTINUE
  480 CONTINUE
  485 CONTINUE
C
C...................
      IF( TEST(16) .GT. 0) GO TO 1016
C
C*******************************************************************************
C***  KERNEL 17     IMPLICIT, CONDITIONAL COMPUTATION       (NO VECTORS)
C*******************************************************************************
C
C         RECURSIVE-DOUBLING VECTOR TECHNIQUES CAN NOT BE USED
C         BECAUSE CONDITIONAL OPERATIONS APPLY TO EACH ELEMENT.
C
                 dw= 5.0000d0/3.0000d0
                 fw= 1.0000d0/3.0000d0
                 tw= 1.0300d0/3.0700d0
cdir$ novector
C
 1017             k= n
                  j= 1
                ink= -1
              SCALE= dw
                XNM= fw
                 E6= tw
                     GO TO 61
C                                            STEP MODEL
  60             E6= XNM*VSP(k)+VSTP(k)
            VXNE(k)= E6
                XNM= E6
             VE3(k)= E6
                  k= k+ink
                 IF( k.EQ.j) GO TO  62
  61             E3= XNM*VLR(k) +VLIN(k)
               XNEI= VXNE(k)
            VXND(k)= E6
                XNC= SCALE*E3
C                                            SELECT MODEL
                 IF( XNM .GT.XNC) GO TO  60
                 IF( XNEI.GT.XNC) GO TO  60
C                                            LINEAR MODEL
             VE3(k)= E3
                 E6= E3+E3-XNM
            VXNE(k)= E3+E3-XNEI
                XNM= E6
                  k= k+ink
                 IF( k.NE.j) GO TO 61
   62 CONTINUE
cdir$ vector
C
C...................
      IF( TEST(17) .GT. 0) GO TO 1017
C
C*******************************************************************************
C***  KERNEL 18     2-D EXPLICIT HYDRODYNAMICS FRAGMENT
C*******************************************************************************
C
C
 1018           T= 0.003700d0
                S= 0.004100d0
               KN= 6
               JN= n
         DO 70  k= 2,KN
         DO 70  j= 2,JN
          ZA(j,k)= (ZP(j-1,k+1)+ZQ(j-1,k+1)-ZP(j-1,k)-ZQ(j-1,k))
     .            *(ZR(j,k)+ZR(j-1,k))/(ZM(j-1,k)+ZM(j-1,k+1))
          ZB(j,k)= (ZP(j-1,k)+ZQ(j-1,k)-ZP(j,k)-ZQ(j,k))
     .            *(ZR(j,k)+ZR(j,k-1))/(ZM(j,k)+ZM(j-1,k))
   70    CONTINUE
C
         DO 72  k= 2,KN
         DO 72  j= 2,JN
          ZU(j,k)= ZU(j,k)+S*(ZA(j,k)*(ZZ(j,k)-ZZ(j+1,k))
     .                    -ZA(j-1,k) *(ZZ(j,k)-ZZ(j-1,k))
     .                    -ZB(j,k)   *(ZZ(j,k)-ZZ(j,k-1))
     .                    +ZB(j,k+1) *(ZZ(j,k)-ZZ(j,k+1)))
          ZV(j,k)= ZV(j,k)+S*(ZA(j,k)*(ZR(j,k)-ZR(j+1,k))
     .                    -ZA(j-1,k) *(ZR(j,k)-ZR(j-1,k))
     .                    -ZB(j,k)   *(ZR(j,k)-ZR(j,k-1))
     .                    +ZB(j,k+1) *(ZR(j,k)-ZR(j,k+1)))
   72    CONTINUE
C
         DO 75  k= 2,KN
         DO 75  j= 2,JN
          ZR(j,k)= ZR(j,k)+T*ZU(j,k)
          ZZ(j,k)= ZZ(j,k)+T*ZV(j,k)
   75    CONTINUE
C
C...................
      IF( TEST(18) .GT. 0) GO TO 1018
C
C*******************************************************************************
C***  KERNEL 19      GENERAL LINEAR RECURRENCE EQUATIONS    (NO VECTORS)
C*******************************************************************************
C
 1019            KB5I= 0
C
C     IF( JR.LE.1 )  THEN
cdir$ novector
             DO 191 k= 1,n
           B5(k+KB5I)= SA(k) +STB5*SB(k)
                 STB5= B5(k+KB5I) -STB5
  191        CONTINUE
C     ELSE
C
             DO 193 i= 1,n
                    k= n-i+1
           B5(k+KB5I)= SA(k) +STB5*SB(k)
                 STB5= B5(k+KB5I) -STB5
  193        CONTINUE
C     ENDIF
cdir$ vector
C
C...................
      IF( TEST(19) .GT. 0) GO TO 1019
C
C*******************************************************************************
C***  KERNEL 20     DISCRETE ORDINATES TRANSPORT: RECURRENCE (NO VECTORS)
C*******************************************************************************
C
           dw= 0.200d0
cdir$ novector
C
 1020 DO 20 k= 1,n
           DI= Y(k)-G(k)/( XX(k)+DK)
           DN= dw
           IF( DI.NE.0.0) DN= MAX( S,MIN( Z(k)/DI, T))
         X(k)= ((W(k)+V(k)*DN)* XX(k)+U(k))/(VX(k)+V(k)*DN)
      XX(k+1)= (X(k)- XX(k))*DN+ XX(k)
   20 CONTINUE
cdir$ vector
C
C...................
      IF( TEST(20) .GT. 0) GO TO 1020
C
C*******************************************************************************
C***  KERNEL 21     MATRIX*MATRIX PRODUCT
C*******************************************************************************
C
C
 1021 DO 21 k= 1,25
      DO 21 i= 1,25
      DO 21 j= 1,n
      PX(i,j)= PX(i,j) +VY(i,k) * CX(k,j)
   21 CONTINUE
C
C...................
      IF( TEST(21) .GT. 0) GO TO 1021
C
C
C
C
C
C
C
C*******************************************************************************
C***  KERNEL 22     PLANCKIAN DISTRIBUTION
C*******************************************************************************
C
C
C      EXPMAX= 234.500d0
       EXPMAX= 20.0000d0
           fw= 1.00000d0
         U(n)= 0.99000d0*EXPMAX*V(n)
C
 1022 DO 22 k= 1,n
care       IF( U(k) .LT. EXPMAX*V(k))  THEN
                                            Y(k)= U(k)/V(k)
care                                   ELSE
care                                        Y(k)= EXPMAX
care    ENDIF
         W(k)= X(k)/( EXP( Y(k)) -fw)
   22 CONTINUE
C...................
      IF( TEST(22) .GT. 0) GO TO 1022
C
C*******************************************************************************
C***  KERNEL 23     2-D IMPLICIT HYDRODYNAMICS FRAGMENT
C*******************************************************************************
C
            fw= 0.17500d0
C
 1023 DO 23  j= 2,6
      DO 23  k= 2,n
            QA= ZA(k,j+1)*ZR(k,j) +ZA(k,j-1)*ZB(k,j) +
     .          ZA(k+1,j)*ZU(k,j) +ZA(k-1,j)*ZV(k,j) +ZZ(k,j)
   23  ZA(k,j)= ZA(k,j) +fw*(QA -ZA(k,j))
C
C...................
      IF( TEST(23) .GT. 0) GO TO 1023
C
C*******************************************************************************
C***  KERNEL 24     FIND LOCATION OF FIRST MINIMUM IN ARRAY
C*******************************************************************************
C
C      X( n/2)= -1.000d+50
       X( n/2)= -1.000d+10
C
 1024        m= 1
      DO 24  k= 2,n
            IF( X(k).LT.X(m))  m= k
   24 CONTINUE
C
C            m= imin1( n,x,1)        35 nanosec./element STACKLIBE/CRAY
C...................
      IF( TEST(24) .NE. 0) GO TO 1024
C
C*******************************************************************************
C
CPFM  iflag1= 0
           sum= 0.00d0
           som= 0.00d0
      DO 999 k= 1,mk
           sum= sum + TIME (k)
      TIMES(jr,il,k)= TIME (k)
      TERRS(jr,il,k)= TERR1(k)
      NPFS (jr,il,k)= NPFS1(k)
      CSUMS(jr,il,k)= CSUM (k)
      DOS  (jr,il,k)= TOTAL(k)
      FOPN (jr,il,k)= FLOPN(k)
           som= som + FLOPN(k) * TOTAL(k)
  999 continue
C
         TK(1)= TK(1) + sum
         TK(2)= TK(2) + som
C                                  Dumpout Checksums
c     WRITE ( 7,706) jr, il
c 706 FORMAT(1X,2I3)
c     WRITE ( 7,707) ( CSUM(k), k= 1,mk)
c 707 FORMAT(5X,'&',1PE21.15,',',1PE21.15,',',1PE21.15,',')
C
      CALL TRACK ('KERNEL  ')
      RETURN
      END

      SUBROUTINE KERN12(lo, n, x, sum)
cANSI IMPLICIT  DOUBLE PRECISION (A-H,O-Z)
      dimension x(n+1)
      integer k, m
c     This function executes kernel 12. It is placed in this file
c     so that it will be optimized the same way as the kernels used 
c     in the benchmark timings. In the original version of the loops,
c     this loop was inline in the timing verification routine and 
c     was run at reduced optimization when there were difficulties with
c     the I/O routines. The result was a serious underestimate of
c     the number of passes needed to get reasonable elapsed times.
c
c     the sum accumulated below is intended to prevent
c     over-clever optimization.
      sum= 0.0
      do 100 m=1,lo
c       the following loop has the minus ones to mimic the
c       macro used elsewhere in the program to get Fortran-like
c       indexing
        do 50 k=1,n
          x(k)= x(k+1) - x(k)
 50     continue
        sum = sum+x(n)
 100  continue
      return
      end
