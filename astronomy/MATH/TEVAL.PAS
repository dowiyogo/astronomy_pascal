(*-----------------------------------------------------------------------*)
(* T_EVAL: evaluates the approximation of a function by Chebyshev        *)
(*         polynomials of maximum order F.M over the interval [F.A,F.B]  *)
(*  F : record containing the Chebyshev coefficients                     *)  
(*  X : argument                                                         *)
(*-----------------------------------------------------------------------*)
 FUNCTION T_EVAL(F: TPOLYNOM; X: REAL): REAL;
  VAR F1,F2,OLD_F1,XX,XX2 : REAL;
      I                : INTEGER;
  BEGIN
    IF ( (X<F.A) OR (F.B<X) ) THEN
      BEGIN WRITELN(' T_EVAL : x not within [a,b]'); END;
    F1 := 0.0;  F2 := 0.0;
    XX := (2.0*X-F.A-F.B)/(F.B-F.A);  XX2 := 2.0*XX;
    FOR I := F.M DOWNTO 1 DO
      BEGIN OLD_F1 := F1; F1 := XX2*F1-F2+F.C[I];  F2 := OLD_F1;  END;
    T_EVAL := XX*F1-F2+0.5*F.C[0]
  END;

