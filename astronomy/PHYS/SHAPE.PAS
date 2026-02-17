(*-----------------------------------------------------------------------*)
(*                                                                       *)
(* SHAPE:  returns the equatorial radii and flattening of planets        *)
(*                                                                       *)
(*   PLANET   Name of the planet                                         *)
(*   R_EQU    Equatorial radius (km)                                     *)
(*   FL       Geometric flattening                                       *)
(*                                                                       *)
(* Source: M.E. Davies et al., Report of the IAU Working Group on        *)
(*         Cartographic Coordinates and Rotational Elements, 1982        *)
(*                                                                       *)
(*-----------------------------------------------------------------------*)
PROCEDURE SHAPE ( PLANET: PLANET_TYPE; VAR R_EQU, FL: REAL );
  BEGIN
    CASE PLANET OF
      MERCURY:  BEGIN  R_EQU :=  2439.0;  FL := 0.0         END;
      VENUS  :  BEGIN  R_EQU :=  6051.0;  FL := 0.0         END;
      EARTH  :  BEGIN  R_EQU :=  6378.14; FL := 0.00335281  END;
      MARS   :  BEGIN  R_EQU :=  3393.4;  FL := 0.0051865   END;
      JUPITER:  BEGIN  R_EQU := 71398.0;  FL := 0.0648088   END;
      SATURN :  BEGIN  R_EQU := 60000.0;  FL := 0.1076209   END;
      URANUS :  BEGIN  R_EQU := 25400.0;  FL := 0.030       END;
      NEPTUNE:  BEGIN  R_EQU := 24300.0;  FL := 0.0259      END;
      PLUTO  :  BEGIN  R_EQU :=  1500.0;  FL := 0.0         END;
    END;
  END; 

