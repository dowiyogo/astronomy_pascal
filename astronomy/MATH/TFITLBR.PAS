(*-----------------------------------------------------------------------*)
(* T_FIT_LBR: expands lunar or planetary coordinates into series of      *)
(*            Chebyshev polynomials for longitude, latitude and radius   *)
(*            that are valid for a specified period of time              *)
(*                                                                       *)
(*  POSITION: routine for calculating the coordinates L,B,R              *)
(*  TA      : first date of desired period of time                       *)
(*  TB      : last date                                                  *)
(*  N       : highest order of Chebyshev polynomials (N<=MAX_TP_DEG)     *)
(*  L_POLY  : coefficients for longitude                                 *)
(*  B_POLY  : coefficients for latitude                                  *)
(*  R_POLY  : coefficients for radius                                    *)
(*                                                                       *)
(* note:                                                                 *)
(*  . the interval [TA,TB] must be shorter than one revolution!          *)
(*  . the routine will only work for heliocentric planetary or           *)
(*    geocentric lunar but not for geocentric planetary coordinates!     *)
(*-----------------------------------------------------------------------*)
PROCEDURE T_FIT_LBR ( PROCEDURE POSITION (T:REAL; VAR LL,BB,RR: REAL);
                      TA,TB: REAL; N: INTEGER;
                      VAR L_POLY,B_POLY,R_POLY: TPOLYNOM);
  CONST PI = 3.1415926535898;
        NDIM = 27;
  VAR   I,J,K          : INTEGER;
        FAC,BPA,BMA,PHI: REAL;
        T,H,L,B,R      : ARRAY[0..NDIM] OF REAL;
  BEGIN
    IF (NDIM<2*MAX_TP_DEG+1) THEN WRITELN(' NDIM too small in T_FIT_LBR');
    IF (N>MAX_TP_DEG) THEN WRITELN(' N too large in T_FIT_LBR');
    L_POLY.M := N;    B_POLY.M := N;    R_POLY.M := N;
    L_POLY.A := TA;   B_POLY.A := TA;   R_POLY.A := TA;
    L_POLY.B := TB;   B_POLY.B := TB;   R_POLY.B := TB;
    BMA := (TB-TA)/2.0;   BPA := (TB+TA)/2.0;
    FAC := 2.0/(N+1);
    PHI:=PI/(2*N+2);                              (* h(k)=cos(pi*k/N/2)  *)
    H[0]:=1.0; H[1]:=COS(PHI);
    FOR I:=2 TO (2*N+1) DO H[I]:=2*H[1]*H[I-1]-H[I-2];
    FOR K:=1 TO N+1 DO T[K] := H[2*K-1]*BMA+BPA;  (* subdivison points   *)
    FOR K:=1 TO N+1 DO POSITION(T[K],L[K],B[K],R[K]);
    FOR K := 2 TO N+1 DO                          (* make L continuous   *)
      IF (L[K-1]<L[K]) THEN L[K]:=L[K]-360.0;     (* in [-360,+360] !    *)
    FOR J := 0 TO N DO                            (* calculate Chebyshev *)
      BEGIN                                       (* coefficients C(j)   *)
        PHI:=PI*J/(2*N+2); H[1]:=COS(PHI);
        FOR I:=2 TO (2*N+1) DO H[I] := 2*H[1]*H[I-1]-H[I-2];
        L_POLY.C[J]:=0.0; B_POLY.C[J]:=0.0; R_POLY.C[J]:=0.0;
        FOR K:=1 TO N+1 DO
          BEGIN
            L_POLY.C[J] := L_POLY.C[J] + H[2*K-1]*L[K];
            B_POLY.C[J] := B_POLY.C[J] + H[2*K-1]*B[K];
            R_POLY.C[J] := R_POLY.C[J] + H[2*K-1]*R[K];
          END;
        L_POLY.C[J]:=L_POLY.C[J]*FAC; B_POLY.C[J]:=B_POLY.C[J]*FAC;
        R_POLY.C[J]:=R_POLY.C[J]*FAC;
      END;
  END;

