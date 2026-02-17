(*-----------------------------------------------------------------------*)
(*                                                                       *)
(* ILLUM: Computes the illumination parameters of a planet               *)
(*                                                                       *)
(*   X,Y,Z      Heliocentric coordinates of the planet                   *)
(*   XE,YE,ZE   Heliocentric coordinates of the Earth                    *)
(*   R          Heliocentric distance of the planet                      *)
(*   D          Geocentric distance of the planet                        *)
(*   ELONG      Elongation (deg)                                         *)
(*   PHI        Phase angle (deg)                                        *)
(*   K          Phase                                                    *)
(*                                                                       *)
(* Note: All coordinates must refer to the same coordinate system.       *)
(*                                                                       *)
(*-----------------------------------------------------------------------*)

PROCEDURE ILLUM ( X,Y,Z, XE,YE,ZE: REAL; VAR R,D,ELONG,PHI,K: REAL );

  VAR   XP,YP,ZP, RE, C_PHI :  REAL;

  BEGIN

    (* Compute the planet's geocentric position *)

    XP:=X-XE; YP:=Y-YE; ZP:=Z-ZE;

    (* Compute the distances in the Sun-Earth-planet triangle *)

    R  := SQRT (  X*X  +  Y*Y  +  Z*Z  );   (* Sun-planet distance   *)
    RE := SQRT ( XE*XE + YE*YE + ZE*ZE );   (* Sun-Earth  distance   *)
    D  := SQRT ( XP*XP + YP*YP + ZP*ZP );   (* Earth-planet distance *)

    (* Compute elongation, phase angle and phase *)

    ELONG := ACS ( ( D*D + RE*RE - R*R ) / ( 2.0*D*RE ) );

    C_PHI := ( D*D + R*R - RE*RE ) / ( 2.0*D*R ) ;
    PHI   := ACS ( C_PHI );
    K     := 0.5*(1.0+C_PHI);

  END; 

