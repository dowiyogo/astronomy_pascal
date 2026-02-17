(*-----------------------------------------------------------------------*)
(* EPSMACH computes the machine accuracy u (1.0+u>1.0,1.0+u/2=1.0)       *)
(*-----------------------------------------------------------------------*)

FUNCTION EPSMACH: REAL;
  VAR ONE,TWO,U: REAL;
  BEGIN
    ONE:=1.0; TWO:=2.0; U:=1.0;
    REPEAT
      U:=U/TWO;
    UNTIL ( (ONE+U)=ONE );
    EPSMACH := TWO*U;
  END;

