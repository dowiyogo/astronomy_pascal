(*-----------------------------------------------------------------------*)
(* QUAD: finds a parabola through 3 points                               *)
(*       (-1,Y_MINUS), (0,Y_0) und (1,Y_PLUS),                           *)
(*       that do not lie on a straight line.                             *)
(*                                                                       *)
(*      Y_MINUS,Y_0,Y_PLUS: three y-values                               *)
(*      XE,YE   : x and y of the extreme value of the parabola           *)
(*      ZERO1   : first root within [-1,+1] (for NZ=1,2)                 *)
(*      ZERO2   : second root within [-1,+1] (only for NZ=2)             *)
(*      NZ      : number of roots within the interval [-1,+1]            *)
(*-----------------------------------------------------------------------*)
PROCEDURE QUAD(Y_MINUS,Y_0,Y_PLUS: REAL; 
               VAR XE,YE,ZERO1,ZERO2: REAL; VAR NZ: INTEGER);
  VAR A,B,C,DIS,DX: REAL;
  BEGIN
    NZ := 0;
    A  := 0.5*(Y_MINUS+Y_PLUS)-Y_0; B := 0.5*(Y_PLUS-Y_MINUS); C := Y_0;
    XE := -B/(2.0*A); YE := (A*XE + B) * XE + C;
    DIS := B*B - 4.0*A*C; (* discriminant of y = axx+bx+c *)
    IF (DIS >= 0) THEN    (* parabola intersects x-axis   *)
      BEGIN
        DX := 0.5*SQRT(DIS)/ABS(A); ZERO1 := XE-DX; ZERO2 := XE+DX;
        IF (ABS(ZERO1) <= +1.0) THEN NZ := NZ + 1;  
        IF (ABS(ZERO2) <= +1.0) THEN NZ := NZ + 1;
        IF (ZERO1<-1.0) THEN ZERO1:=ZERO2;
      END;
    END;
