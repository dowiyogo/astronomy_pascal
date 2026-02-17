(*-----------------------------------------------------------------------*)
(* NUTEQU: transformation of mean to true coordinates                    *)
(*         (including terms >0.1" according to IAU 1980)                 *)
(*         T = (JD-2451545.0)/36525.0                                    *)
(*-----------------------------------------------------------------------*)
PROCEDURE NUTEQU(T:REAL;VAR X,Y,Z:REAL);
  CONST ARC=206264.8062;          (* arcseconds per radian = 3600*180/pi *)
        P2 =6.283185307;          (* 2*pi                                *)
  VAR   LS,D,F,N,EPS : REAL;
        DPSI,DEPS,C,S: REAL;
        DX,DY,DZ     : REAL;
  FUNCTION FRAC(X:REAL):REAL;
    (* with several compilers it may be necessary to replace TRUNC *)
    (* by LONG_TRUNC or INT if T<-24!                              *)
    BEGIN  FRAC:=X-TRUNC(X) END;
  BEGIN
    LS  := P2*FRAC(0.993133+  99.997306*T); (* mean anomaly Sun          *)
    D   := P2*FRAC(0.827362+1236.853087*T); (* diff. longitude Moon-Sun  *)
    F   := P2*FRAC(0.259089+1342.227826*T); (* mean argument of latitude *)
    N   := P2*FRAC(0.347346-   5.372447*T); (* longit. ascending node    *)
    EPS := 0.4090928-2.2696E-4*T;           (* obliquity of the ecliptic *)
    DPSI := ( -17.200*SIN(N)   - 1.319*SIN(2*(F-D+N)) - 0.227*SIN(2*(F+N))
              + 0.206*SIN(2*N) + 0.143*SIN(LS) ) / ARC;
    DEPS := ( + 9.203*COS(N)   + 0.574*COS(2*(F-D+N)) + 0.098*COS(2*(F+N))
              - 0.090*COS(2*N)                 ) / ARC;
    C := DPSI*COS(EPS);  S := DPSI*SIN(EPS);
    DX := -(C*Y+S*Z); DY := (C*X-DEPS*Z); DZ := (S*X+DEPS*Y);
    X  :=  X + DX;        Y  := Y + DY;       Z  := Z + DZ;
  END;

