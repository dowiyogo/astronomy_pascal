(*-----------------------------------------------------------------------*)
(* Unit PLALIB : Planetary position subroutines                          *)
(*-----------------------------------------------------------------------*)

UNIT PLALIB;

INTERFACE

  TYPE  PLANET_TYPE = ( MERCURY, VENUS,  EARTH,  MARS,
                        JUPITER, SATURN, URANUS, NEPTUNE, PLUTO );

  PROCEDURE POSITION ( PLANET: PLANET_TYPE; T: REAL;  VAR X,Y,Z: REAL );


IMPLEMENTATION

  USES MATLIB,SPHLIB,KEPLIB;

  (*$I \PLANETS\POSITION.PAS *)

END.
