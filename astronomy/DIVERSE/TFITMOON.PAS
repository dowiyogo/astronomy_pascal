(*-----------------------------------------------------------------------*)
(* T_FIT_MOON: approximates the equatorial lunar coordinates             *)
(*             of the moon by Chebyshev expansions for a                 *)
(*             given period of time of at most one month                 *)
(*                                                                       *)
(*  TA     : first date (in julian centuries since J2000)                *)
(*  TB     : last date ( TB < TA + 1 month )                             *)
(*  N      : highest order                                               *)
(*  RA_POLY: coefficients for right ascencion                            *)
(*  DE_POLY: coefficients for declination                                *)
(*  R_POLY : coefficients for geocentric distance                        *)
(*-----------------------------------------------------------------------*)
PROCEDURE T_FIT_MOON ( TA,TB: REAL; N: INTEGER;
                       VAR RA_POLY,DE_POLY,R_POLY: TPOLYNOM);
  CONST PI = 3.1415926535898;
        NDIM = 27;
  VAR   I,J,K          : INTEGER;
        FAC,BPA,BMA,PHI: REAL;
        T,H,RA,DE,R    : ARRAY[0..NDIM] OF REAL;
  BEGIN
    IF (NDIM<2*MAX_TP_DEG+1) THEN WRITELN(' NDIM too small in T_FIT_MOON');
    IF (N>MAX_TP_DEG) THEN WRITELN(' N too large in T_FIT_MOON');
    RA_POLY.M := N;    DE_POLY.M := N;    R_POLY.M := N;
    RA_POLY.A := TA;   DE_POLY.A := TA;   R_POLY.A := TA;
    RA_POLY.B := TB;   DE_POLY.B := TB;   R_POLY.B := TB;
    BMA := (TB-TA)/2.0;   BPA := (TB+TA)/2.0;
    FAC := 2.0/(N+1);
    PHI:=PI/(2*N+2);                              (* h(k)=cos(pi*k/N/2)  *)
    H[0]:=1.0; H[1]:=COS(PHI);
    FOR I:=2 TO (2*N+1) DO H[I]:=2*H[1]*H[I-1]-H[I-2];
    FOR K:=1 TO N+1 DO T[K] := H[2*K-1]*BMA+BPA;  (* subdivision points  *)
    FOR K:=1 TO N+1 DO MOONEQU(T[K],RA[K],DE[K],R[K]);
    FOR K := 2 TO N+1 DO                          (* make RA continuous  *)
      IF (RA[K-1]<RA[K]) THEN RA[K]:=RA[K]-360.0; (* in [-360,+360] deg  *)
    FOR J := 0 TO N DO                            (* calculate Chebyshev *)
      BEGIN                                       (* coefficients C[J]   *)
        PHI:=PI*J/(2*N+2); H[1]:=COS(PHI);
        FOR I:=2 TO (2*N+1) DO H[I] := 2*H[1]*H[I-1]-H[I-2];
        RA_POLY.C[J]:=0.0; DE_POLY.C[J]:=0.0; R_POLY.C[J]:=0.0;
        FOR K:=1 TO N+1 DO
          BEGIN
            RA_POLY.C[J] := RA_POLY.C[J] + H[2*K-1]*RA[K];
            DE_POLY.C[J] := DE_POLY.C[J] + H[2*K-1]*DE[K];
            R_POLY.C[J]  := R_POLY.C[J]  + H[2*K-1]*R[K];
          END;
        RA_POLY.C[J]:=RA_POLY.C[J]*FAC; DE_POLY.C[J]:=DE_POLY.C[J]*FAC;
        R_POLY.C[J] :=R_POLY.C[J]*FAC;
      END;
  END;

