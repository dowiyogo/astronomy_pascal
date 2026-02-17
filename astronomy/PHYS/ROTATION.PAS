(*-----------------------------------------------------------------------*)
(*                                                                       *)
(* ROTATION: computes the rotation parameters of a planet                *)
(*                                                                       *)
(*   X,Y,Z    Geocentric equatorial coordinates of the planet (AU)       *)
(*   A        Right ascension of the axis of rotation (J2000)            *)
(*   D        Declination of the axis of rotation (J2000)                *)
(*   W        Orientation of the prime meridian with respect to the      *)
(*            intersection of the Earth's equator of J2000 and the       *)
(*            planetary equator of date                                  *)
(*   SENSE    Sense of rotation (DIRECT or RETROGRADE)                   *)
(*   FLATT    Geometric flattening of the planet                         *)
(*   AX,AY,AZ Rotation axis unit vector (J2000)                          *)
(*   LONG     Planetographic longitude of the Earth (deg)                *)
(*   LAT      Planetographic latitude of the Earth (deg)                 *)
(*   DEC      Planetocentric latitude of the Earth (deg)                 *)
(*                                                                       *)
(* The heliocentric coordinates of the planet may be substituted for     *)
(* X,Y,Z to obtain the planetographic coordinates of the Sun.            *)
(*                                                                       *)
(*-----------------------------------------------------------------------*)

PROCEDURE ROTATION ( X,Y,Z, A,D,W: REAL; SENSE: ROTATION_TYPE; FLATT: REAL; 
                     VAR AX,AY,AZ, LONG,LAT,DEC: REAL );

  VAR   SX,SY,SZ,T: REAL;
        E         : REAL33;

  BEGIN

    (* Compute unit vectors E(*,1) (intersection of prime meridian and *)
    (* planetary equator), E(*,3) (parallel to the rotation axis) and  *)
    (* E(*,2) perpendicular to E(*,1) and E(*,3))                      *)

    GAUSVEC ( 90.0+A,90.0-D,W, E );

    (* Copy rotation axis unit vector *)

    AX := E[1,3];  AY := E[2,3];  AZ := E[3,3];

    (* Compute planetocentric latitude and longitude *)

    SX := - ( E[1,1]*X + E[2,1]*Y + E[3,1]*Z );
    SY := - ( E[1,2]*X + E[2,2]*Y + E[3,2]*Z );
    SZ := - ( E[1,3]*X + E[2,3]*Y + E[3,3]*Z );

    T  := SZ / SQRT(SX*SX+SY*SY);

    DEC  := ATN ( T );
    LONG := ATN2 ( SY, SX );

    (* Compute planetographic latitude and longitude *)

    IF (SENSE=DIRECT) THEN LONG:=-LONG;
    IF LONG<0.0 THEN LONG:=LONG+360.0;
    LAT := ATN ( T / ((1.0-FLATT)*(1.0-FLATT)) );

  END; 

