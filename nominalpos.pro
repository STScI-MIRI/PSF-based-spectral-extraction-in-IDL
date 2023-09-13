FUNCTION nominalpos,fovid,SLT=slt,XSLT=xslt

;  6 Dec 11 created
;
; given a FOVID, return the nominal position of the trace
; position given as the fractional column where the trace crosses row 0
; this information is gathered here in one routine to be called by many
;   since it is hardcoded but could change
;
; INPUT
;   fovid = the FOVID (field-of-view ID number) from the FITS header
;   SLT   = keyword to force the return of the position for the slit center
;   XSLT  = keyword to force position of other slit center
; set tfovid to slit center (if slt set), other slit center (xslt set),
;    or fovid (neither set)

if (keyword_set(slt) ne 0) then begin
  case fovid of
    26 : tfovid=28
    27 : tfovid=28
    32 : tfovid=34
    33 : tfovid=34
    38 : tfovid=40
    39 : tfovid=40
    44 : tfovid=46
    45 : tfovid=46
    else : tfovid=fovid
  endcase
endif else if (keyword_set(xslt) ne 0) then begin
  case fovid of
    26 : tfovid=34
    27 : tfovid=34
    32 : tfovid=28
    33 : tfovid=28
    38 : tfovid=46
    39 : tfovid=46
    44 : tfovid=40
    45 : tfovid=40
    else : tfovid=fovid
  endcase
endif else tfovid=fovid

; return the nominal position - see notes 28 Nov 11

case tfovid of  ; nominal position in pixels on row 0 of the array
  26 : nominal =  9.67 ; SL1 nod 1
  27 : nominal = 20.29 ; SL1 nod 2
  28 : nominal = 14.98 ; SL1 center 
  32 : nominal = 52.81 ; SL2 nod 1
  33 : nominal = 62.92 ; SL2 nod 2
  34 : nominal = 57.87 ; SL2 center 
  38 : nominal = 47.44 ; LL1 nod 1
  39 : nominal = 36.64 ; LL1 nod 2
  40 : nominal = 42.04 ; LL1 center 
  44 : nominal = 85.58 ; LL2 nod 1
  45 : nominal = 74.47 ; LL2 nod 2
  46 : nominal = 80.03 ; LL2 center 
  else : nominal = -9.99
endcase

RETURN,nominal
END
