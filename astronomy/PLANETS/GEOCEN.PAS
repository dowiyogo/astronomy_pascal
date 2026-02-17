(*-----------------------------------------------------------------------*)
(*                                                                       *)
(* GEOCEN: geocentric coordinates (geometric and light-time corrected)   *)
(*                                                                       *)
(*   T:        time in Julian centuries since J2000                      *)
(*             T=(JD-2451545.0)/36525.0                                  *)
(*   LP,BP,RP: ecliptic heliocentric coordinates of the planet           *)
(*   LS,BS,RS: ecliptic geocentric   coordinates of the sun              *)
(*                                                                       *)
(*   IPLAN:    planet (0=Sun,1=Mercury,2=Venus,3=Earth,...,9=Pluto       *)
(*   IMODE:    desired type of coordinates (see description of X,Y,Z)    *)
(*             (0=geometric,1=astrometric,2=apparent)                    *)
(*   XP,YP,ZP: ecliptic heliocentric coordinates of the planet           *)
(*   XS,YS,ZS: ecliptic geocentric coordinates of the Sun                *)
(*   X, Y, Z : ecliptic geocentric cordinates of the planet (geometric   *)
(*             if IMODE=0, astrometric if IMODE=1, apparent if IMODE=2)  *)
(*   DELTA0:   geocentric distance (geometric)                           *)
(*                                                                       *)
(*   (all angles in degrees, distances in AU)                            *)
(*                                                                       *)
(*-----------------------------------------------------------------------*)

PROCEDURE GEOCEN(T, LP,BP,RP, LS,BS,RS: REAL; IPLAN,IMODE: INTEGER;
                 VAR XP,YP,ZP, XS,YS,ZS, X,Y,Z,DELTA0: REAL);

  CONST P2=6.283185307;

  VAR DL,DB,DR, DLS,DBS,DRS, FAC: REAL;
      VX,VY,VZ, VXS,VYS,VZS, M  : REAL;

  FUNCTION FRAC(X:REAL):REAL;
    BEGIN  X:=X-TRUNC(X); IF (X<0) THEN X:=X+1; FRAC:=X  END;

  PROCEDURE POSVEL(L,B,R,DL,DB,DR: REAL; VAR X,Y,Z,VX,VY,VZ:REAL);
    VAR  CL,SL,CB,SB: REAL;
    BEGIN
      CL:=CS(L); SL:=SN(L); CB:=CS(B); SB:=SN(B);
      X := R*CL*CB;  VX := DR*CL*CB-DL*R*SL*CB-DB*R*CL*SB;
      Y := R*SL*CB;  VY := DR*SL*CB+DL*R*CL*CB-DB*R*SL*SB;
      Z := R*SB;     VZ := DR*SB              +DB*R*CB;
    END;

  BEGIN

    DL:=0.0; DB:=0.0; DR:=0.0; DLS:=0.0; DBS:=0.0; DRS:=0.0;

    IF (IMODE>0) THEN

      BEGIN

        M := P2*FRAC(0.9931266+ 99.9973604*T);              (* Sun     *)
        DLS := 172.00+5.75*SIN(M);  DRS := 2.87*COS(M); DBS := 0.0;

        (* dl,db in 1e-4 rad/d, dr in 1e-4 AU/d *)
        CASE IPLAN OF
          0: BEGIN  DL:=0.0; DB:=0.0; DR:=0.0;  END;        (* Sun     *)
          1: BEGIN                                          (* Mercury *)
               M := P2*FRAC(0.4855407+415.2014314*T);
               DL := 714.00+292.66*COS(M)+71.96*COS(2*M)+18.16*COS(3*M)
                     +4.61*COS(4*M)+3.81*SIN(2*M)+2.43*SIN(3*M)
                     +1.08*SIN(4*M);
               DR := 55.94*SIN(M)+11.36*SIN(2*M)+2.60*SIN(3*M);
               DB := 73.40*COS(M)+29.82*COS(2*M)+10.22*COS(3*M)
                     +3.28*COS(4*M)-40.44*SIN(M)-16.55*SIN(2*M)
                     -5.56*SIN(3*M)-1.72*SIN(4*M);
             END;
          2: BEGIN                                          (* Venus   *)
               M := P2*FRAC(0.1400197+162.5494552*T);
               DL := 280.00+3.79*COS(M);   DR := 1.37*SIN(M);
               DB := 9.54*COS(M)-13.57*SIN(M);
             END;
          3: BEGIN  DL:=DLS; DR:=DRS; DB:=-DBS;  END;       (* Earth   *)
          4: BEGIN                                          (* Mars    *)
               M := P2*FRAC(0.0538553+53.1662736*T);
               DL := 91.50+17.07*COS(M)+2.03*COS(2*M);
               DR := 12.98*SIN(M)+1.21*COS(2*M);
               DB :=  0.83*COS(M)+2.80*SIN(M);
             END;
          5: BEGIN                                          (* Jupiter *)
               M := P2*FRAC(0.0565314+8.4302963*T);
               DL := 14.50+1.41*COS(M); DR:=3.66*SIN(M); DB:=0.33*SIN(M);
             END;
          6: BEGIN                                          (* Saturn  *)
               M := P2*FRAC(0.8829867+3.3947688*T);
               DL := 5.84+0.65*COS(M); DR:=3.09*SIN(M); DB:=0.24*COS(M);
             END;
          7: BEGIN                                          (* Uranus  *)
               M := P2*FRAC(0.3967117+1.1902849*T);
               DL := 2.05+0.19*COS(M); DR:=1.86*SIN(M); DB:=-0.03*SIN(M);
             END;
          8: BEGIN                                          (* Neptune *)
               M := P2*FRAC(0.7214906+0.6068526*T);
               DL := 1.04+0.02*COS(M); DR:=0.27*SIN(M); DB:=0.03*SIN(M);
             END;
          9: BEGIN                                          (* Pluto   *)
               M := P2*FRAC(0.0385795+0.4026667*T);
               DL := 0.69+0.34*COS(M)+0.12*COS(2*M)+0.05*COS(3*M);
               DR := 6.66*SIN(M)+1.64*SIN(2*M);
               DB := -0.08*COS(M)-0.17*SIN(M)-0.09*SIN(2*M);
             END;
        END;

      END;

    POSVEL (LS,BS,RS,DLS,DBS,DRS,XS,YS,ZS,VXS,VYS,VZS);
    POSVEL (LP,BP,RP,DL ,DB ,DR, XP,YP,ZP,VX ,VY ,VZ );
    X:=XP+XS;  Y:=YP+YS;  Z:=ZP+ZS;   DELTA0 := SQRT(X*X+Y*Y+Z*Z);
    IF IPLAN=3 THEN BEGIN X:=0.0; Y:=0.0; Z:=0.0; DELTA0:=0.0 END;

    FAC := 0.00578 * DELTA0 * 1E-4;
    CASE IMODE OF
      1: BEGIN X:=X-FAC*VX;  Y:=Y-FAC*VY;  Z:=Z-FAC*VZ; END;
      2: BEGIN X:=X-FAC*(VX+VXS); Y:=Y-FAC*(VY+VYS); Z:=Z-FAC*(VZ+VZS);END;
     END;

  END;

