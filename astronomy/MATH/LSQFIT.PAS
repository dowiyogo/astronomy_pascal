(*-----------------------------------------------------------------------*)
(* LSQFIT:                                                               *)
(*   solution of an overdetermined system of linear equations            *)
(*   A[i,1]*s[1]+...A[i,m]*s[m] - A[i,m+1] = 0   (i=1,..,n)              *)
(*   according to the method of least squares using Givens rotations     *)
(*   A: matrix of coefficients                                           *)
(*   N: number of equations  (rows of A)                                 *)
(*   M: number of unknowns   (M+1=columns of A, M=elements of S)         *)
(*   S: solution vector                                                  *)
(*-----------------------------------------------------------------------*)
PROCEDURE LSQFIT ( A: LSQMAT; N, M: INTEGER; VAR S: LSQVEC );

  CONST EPS = 1.0E-10;  (* machine accuracy *)

  VAR I,J,K: INTEGER;
      P,Q,H: REAL;

  BEGIN

    FOR J:=1 TO M DO  (* loop over columns 1...M of A *)

      (* eliminate matrix elements A[i,j] with i>j from column j *)

      FOR I:=J+1 TO N DO
        IF A[I,J]<>0.0 THEN
          BEGIN
            (* calculate p, q and new A[j,j]; set A[i,j]=0 *)
            IF ( ABS(A[J,J])<EPS*ABS(A[I,J]) )
              THEN
                BEGIN
                  P:=0.0; Q:=1.0; A[J,J]:=-A[I,J]; A[I,J]:=0.0;
                END
              ELSE
                BEGIN
                  H:=SQRT(A[J,J]*A[J,J]+A[I,J]*A[I,J]);
                  IF A[J,J]<0.0 THEN H:=-H;
                  P:=A[J,J]/H; Q:=-A[I,J]/H; A[J,J]:=H; A[I,J]:=0.0;
                END;
            (*  calculate rest of the line *)
            FOR K:=J+1 TO M+1 DO
              BEGIN
                H      := P*A[J,K] - Q*A[I,K];
                A[I,K] := Q*A[J,K] + P*A[I,K];
                A[J,K] := H;
              END;
          END;

    (* backsubstitution *)

    FOR I:=M DOWNTO 1 DO
      BEGIN
        H:=A[I,M+1];
        FOR K:=I+1 TO M DO H:=H-A[I,K]*S[K];
        S[I] := H/A[I,I];
      END;

  END;  (* LSQFIT *)

