(*---------------------------------------------------------------------------*)
(* Unit PHYLIB: physical ephemerides                                         *)
(*---------------------------------------------------------------------------*)

UNIT PHYLIB;

INTERFACE

  USES MATLIB,SPHLIB,PLALIB;

  TYPE  ROTATION_TYPE = ( DIRECT, RETROGRADE );
        SYSTEM_TYPE   = ( SYS_I, SYS_II, SYS_III );

  FUNCTION  POSANG   ( X,Y,Z, DX,DY,DZ : REAL ): REAL;
  PROCEDURE SHAPE    ( PLANET: PLANET_TYPE; VAR R_EQU, FL: REAL );
  PROCEDURE ORIENT   ( PLANET: PLANET_TYPE; SYSTEM: SYSTEM_TYPE; T: REAL;
                       VAR A,D,W: REAL; VAR SENSE: ROTATION_TYPE );
  PROCEDURE ROTATION ( X,Y,Z, A,D,W: REAL;  SENSE: ROTATION_TYPE; FLATT: REAL; 
                       VAR AX,AY,AZ, LONG,LAT,DEC: REAL );
  FUNCTION  BRIGHT   ( PLANET: PLANET_TYPE; 
                       R,DELTA,PHI,DEC,DLONG: REAL ): REAL;
  PROCEDURE ILLUM    ( X,Y,Z, XE,YE,ZE : REAL; VAR R,D,ELONG,PHI,K : REAL  );
  
IMPLEMENTATION

  (*$I PHYS\POSANG.PAS   *)
  (*$I PHYS\SHAPE.PAS    *)
  (*$I PHYS\ORIENT.PAS   *)
  (*$I PHYS\ROTATION.PAS *)
  (*$I PHYS\ILLUM.PAS    *)
  (*$I PHYS\BRIGHT.PAS   *)

END.
