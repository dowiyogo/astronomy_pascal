(*-----------------------------------------------------------------------*)
(* SITE:  calculates geocentric from geographic coordinates              *)
(*        RCPHI:  r * cos(phi') (geocentric; in earth radii)             *)
(*        RSPHI:  r * sin(phi') (geocentric; in earth radii)             *)
(*        PHI:    geographic latitude (deg)                              *)
(*-----------------------------------------------------------------------*)
PROCEDURE SITE ( PHI: REAL; VAR RCPHI,RSPHI: REAL );
  CONST E2 = 0.006694;         (* e**2=f(2-f) for flattening f=1/298.257 *)
  VAR   N,SNPHI: REAL;
  BEGIN
    SNPHI := SN(PHI);    N := 1.0/SQRT(1.0-E2*SNPHI*SNPHI);
    RCPHI := N*CS(PHI);  RSPHI := (1.0-E2)*N*SNPHI;
  END;

