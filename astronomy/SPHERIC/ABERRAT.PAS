(*-----------------------------------------------------------------------*)
(* ABERRAT: velocity vector of the Earth in equatorial coordinates       *)
(*          (in units of the velocity of light)                          *)
(*-----------------------------------------------------------------------*)
PROCEDURE ABERRAT(T: REAL; VAR VX,VY,VZ: REAL);
  CONST P2=6.283185307;
  VAR L,CL: REAL;
  FUNCTION FRAC(X:REAL):REAL;
    BEGIN  X:=X-TRUNC(X); IF (X<0) THEN X:=X+1; FRAC:=X  END;
  BEGIN
    L := P2*FRAC(0.27908+100.00214*T);  CL:=COS(L); 
    VX := -0.994E-4*SIN(L); VY := +0.912E-4*CL; VZ := +0.395E-4*CL;
  END;

