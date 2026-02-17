(*-----------------------------------------------------------------------*)
(* ELEMENT: calculates orbital elements from two positions               *)
(*          for elliptic and hyperbolic orbits                           *)
(*                                                                       *)
(*   JDA,JDB: time of passage of points A and B (Julian Date)            *)
(*   RA, RB : position vectors of points A and B                         *)
(*   TP     : perihelion time (in Julian centuries since J2000)          *)
(*   Q      : perihelion distance                                        *)
(*   ECC    : eccentricity                                               *)
(*   INC    : inclination (in deg)                                       *)
(*   LAN    : longitude of the ascending node (in deg)                   *)
(*   AOP    : argument of perihelion (in deg)                            *)
(*-----------------------------------------------------------------------*)
PROCEDURE ELEMENT ( JDA,JDB: REAL; RA,RB: VECTOR;
                    VAR TP,Q,ECC,INC,LAN,AOP: REAL);

  CONST KGAUSS = 0.01720209895;
        RAD    = 0.01745329252; (* 180/pi *)

  VAR   TAU,ETA,P,AX,N,NY,E,M,U         : REAL;
        SA,SB,S0,FAC,DUMMY,SHH          : REAL;
        COS_DNY,SIN_DNY,ECOS_NY,ESIN_NY : REAL;
        EA,R0,E0,R                      : VECTOR;
        I                               : INDEX;

  BEGIN

    (* calculate vector R0 (fraction of RB perpendicular to RA)   *)
    (* and the magnitudes of RA, RB and R0                        *)
    SA := NORM(RA);     FOR I:=X TO Z DO  EA[I]:=RA[I]/SA;
    SB := NORM(RB);
    FAC := DOT(RB,EA);  FOR I:=X TO Z DO  R0[I]:=RB[I]-FAC*EA[I];
    S0 := NORM(R0);     FOR I:=X TO Z DO  E0[I]:=R0[I]/S0;

    (* inclination and ascending node *)
    CROSS (EA,E0,R);
    POLAR ( -R[Y],R[X],R[Z], DUMMY,INC,LAN );  INC := 90.0-INC;
    U  :=  ATN2 ( (+E0[X]*R[Y]-E0[Y]*R[X]) , (-EA[X]*R[Y]+EA[Y]*R[X]) );
    IF INC=0.0 THEN U:=ATN2(RA[Y],RA[X]);

    (* semilatus rectum p *)
    TAU := KGAUSS * ABS(JDB-JDA);   ETA := FIND_ETA(RA,RB,TAU);
    P := SA*S0*ETA / TAU;   P := P*P;

    (* eccentricity, true anomaly and longitude of perihelion *)
    COS_DNY := FAC/SB;    SIN_DNY := S0/SB;
    ECOS_NY := P/SA-1.0;  ESIN_NY := (ECOS_NY*COS_DNY-(P/SB-1.0))/SIN_DNY;
    POLAR ( ECOS_NY,ESIN_NY,0.0, ECC,DUMMY,NY );
    AOP := U-NY;  WHILE (AOP<0.0) DO AOP:=AOP+360.0;

    (* perihelion distance, semimajor axis and mean daily motion *)
    Q  := P/(1.0+ECC);     AX := Q/(1.0-ECC);
    N  := KGAUSS / SQRT(ABS(AX*AX*AX));

    (* mean anomaly and time of perihelion passage *)
    IF (ECC<1.0)
      THEN
        BEGIN
          E := ATN2 ( SQRT((1.0-ECC)*(1.0+ECC))*ESIN_NY, ECOS_NY+ECC*ECC );
          E := RAD*E;  M := E-ECC*SIN(E); ;
        END
      ELSE
        BEGIN
          SHH := SQRT((ECC-1.0)*(ECC+1.0))*ESIN_NY / (ECC+ECC*ECOS_NY) ;
          M   := ECC*SHH - LN(SHH+SQRT(1.0+SHH*SHH))
        END;
    TP := ( (JDA-M/N) - 2451545.0 ) / 36525.0;

  END;

