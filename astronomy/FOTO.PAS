(*-----------------------------------------------------------------------*)
(*                                 FOTO                                  *)
(*              astrometric analysis of photographic plates              *)
(*                           version 93/07/01                            *)
(*-----------------------------------------------------------------------*)
PROGRAM FOTO (INPUT,OUTPUT,FOTINP);

  USES  MATLIB, SPHLIB;

  CONST MAXDIM      = 30;      (* maximum number of objects on the photo *)
        NAME_LENGTH = 12;
        ARC         = 206264.8;                 (* arcseconds per radian *)

  TYPE  NAME_TYPE  =  ARRAY[1..NAME_LENGTH] OF CHAR;
        REAL_ARRAY =  ARRAY[1..MAXDIM] OF REAL;
        NAME_ARRAY =  ARRAY[1..MAXDIM] OF NAME_TYPE;

  VAR I,J,K, NREF,NOBJ, DEG,MIN  : INTEGER;
      RA0,DEC0, A,B,C,D,E,F, SEC : REAL;
      RA_OBS,DEC_OBS, D_RA,D_DEC : REAL;
      DET, FOC_LEN, SCALE        : REAL;
      RA,DEC, X,Y, XX,YY, DELTA  : REAL_ARRAY;
      S                          : LSQVEC;
      AA                         : LSQMAT;
      NAME                       : NAME_ARRAY;
      FOTINP                     : TEXT;


(*-----------------------------------------------------------------------*)
(* GETINP: read input data from file FOTINP                              *)
(*-----------------------------------------------------------------------*)
PROCEDURE GETINP ( VAR RA0,DEC0: REAL;  VAR NOBJ: INTEGER;
                   VAR NAME: NAME_ARRAY;  VAR RA,DEC,X,Y: REAL_ARRAY );

  VAR I,K,H,M: INTEGER;
      S      : REAL;
      C      : CHAR;

  BEGIN

    WRITELN;
    WRITELN ('     FOTO: astrometric analysis of photographic plates  ');
    WRITELN ('                      version 93/07/01                  ');
    WRITELN ('        (c) 1993 Thomas Pfleger, Oliver Montenbruck     ');
    WRITELN;

    (* open file for reading *
    (* RESET(FOTINP); *)                              (* Standard Pascal *)
    ASSIGN(FOTINP,'FOTINP.DAT'); RESET(FOTINP);       (* TURBO Pascal    *)

    WRITELN (' Input data file: FOTINP'); WRITELN;

    (*read coordinates of the plate center *)
    FOR K:=1 TO NAME_LENGTH DO READ(FOTINP,C);
    READ  (FOTINP,H,M,S); DDD(H,M,S,RA0); RA0:=15.0*RA0;
    READLN(FOTINP,H,M,S); DDD(H,M,S,DEC0);

    (* read name, plate coordinates (and equatorial coordinates) *)
    I := 0;
    REPEAT
      I := I+1;
      FOR K:=1 TO NAME_LENGTH DO READ(FOTINP,NAME[I][K]);   (* name *)
      IF NAME[I][1]='*'
        THEN (* reference star *)
          BEGIN
            READ  (FOTINP,X[I],Y[I]);
            READ  (FOTINP,H,M,S); DDD(H,M,S,RA[I]); RA[I]:=15.0*RA[I];
            READLN(FOTINP,H,M,S); DDD(H,M,S,DEC[I]);
          END
        ELSE (* unknown object *)
          BEGIN
            READLN (FOTINP,X[I],Y[I]);  RA[I]:=0.0; DEC[I]:=0.0;
          END;
    UNTIL EOF(FOTINP);

    NOBJ := I;

  END;

(*-----------------------------------------------------------------------*)

BEGIN  (* FOTO *)

  (* read input from file *)

  GETINP ( RA0,DEC0, NOBJ, NAME, RA,DEC,X,Y );

  (* calculate standard coordinates of reference stars; *)
  (* fill elements of matrix AA (column AA[*,5] serves  *)
  (* as intermediate storage)                           *)

  J:=0;
  FOR I:=1 TO NOBJ DO
    IF NAME[I][1]='*' THEN
      BEGIN
        J := J+1;
        EQUSTD ( RA0,DEC0, RA[I],DEC[I], XX[I],YY[I] );
        AA[J,1]:=  X[I]; AA[J,2]:=  Y[I]; AA[J,3]:=1.0;
        AA[J,4]:=+XX[I]; AA[J,5]:=+YY[I];
      END;
  NREF := J;    (* number of reference stars *)

  (* calculate plate constants a,b,c *)

  LSQFIT ( AA, NREF, 3, S );  A:=S[1]; B:=S[2]; C:=S[3];

  (* calculate plate constants d,e,f *)

  FOR I:=1 TO NREF DO  AA[I,4]:=AA[I,5]; (* copy column A[*,5]->A[*,4] *)
  LSQFIT ( AA, NREF, 3, S );  D:=S[1]; E:=S[2]; F:=S[3];

  (* calculate equatorial coordinates (and errors for reference stars) *)

  FOR I:=1 TO NOBJ DO
    BEGIN
      XX[I] := A*X[I]+B*Y[I]+C;
      YY[I] := D*X[I]+E*Y[I]+F;
      STDEQU ( RA0,DEC0, XX[I],YY[I], RA_OBS,DEC_OBS );
      IF NAME[I][1]='*' THEN (* find error in arcseconds *)
        BEGIN
          D_RA  := (RA_OBS-RA[I])*CS(DEC[I]);
          D_DEC := (DEC_OBS-DEC[I]);
          DELTA[I] := 3600.0 * SQRT ( D_RA*D_RA + D_DEC*D_DEC );
        END;
      RA[I] := RA_OBS;  DEC[I] := DEC_OBS;
    END;

  (* focal length *)

  DET := A*E-D*B;
  FOC_LEN := 1.0/SQRT(ABS(DET));
  SCALE := ARC / FOC_LEN;

  (* output *)

  WRITELN (' Plate constants:' );
  WRITELN;
  WRITELN ('  a =',a:12:8,'  b =',b:12:8,'  c =',c:12:8);
  WRITELN ('  d =',d:12:8,'  e =',e:12:8,'  f =',f:12:8);
  WRITELN;
  WRITELN (' Effective focal length and image scale:');
  WRITELN;
  WRITELN ('  F =',FOC_LEN:9:2,' mm');
  WRITELN ('  m =',  SCALE:7:2,' "/mm');
  WRITELN;
  WRITELN (' Coordinates:');
  WRITELN;
  WRITELN (' Name':11, 'x':9,'y':7,'X':8,'Y':8,
           'RA':12, 'Dec':13, 'Error':9 );
  WRITELN ('mm':20, 'mm':7, ' ':23,
           'h  m  s  ', 'o  ''  " ':12, ' " ':6 );

  FOR I:=1 TO NOBJ DO
    BEGIN
      WRITE('  '); FOR K:=1 TO NAME_LENGTH DO WRITE(NAME[I][K]);
      WRITE(X[I]:7:1,Y[I]:7:1,XX[I]:9:4,YY[I]:8:4);
      DMS(RA[I]/15.0,DEG,MIN,SEC); WRITE(DEG:5,MIN:3,SEC:6:2);
      DMS(DEC[I],DEG,MIN,SEC); WRITE(DEG:4,MIN:3,SEC:5:1);
      IF NAME[I][1]='*' THEN WRITE(DELTA[I]:6:1);
      WRITELN;
    END;
  WRITELN;

END. (* FOTO *)

(*-----------------------------------------------------------------------*)
