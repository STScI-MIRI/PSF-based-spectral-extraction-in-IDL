PRO mk_pos,extpos,header,XAP=xap,VERBOSE=verbose

; 18 Feb 12 tested to see that xap=1 works, made verbose a keyword
;  6 Dec 11 moving nominal pos'n to separate routine, allowing xap ne 0 now
; 28 Nov 11 created
;
; given an extracted position (extpos) and an IRS FITS header
;   find the extracted RA and DEC
;
; INPUT
;   extpos  - the extracted position along the slit
;             defined to be the position of the spectral trace in row 0 (y_int)
;   header  - the FITS header - keywords will be read from here and added
;   xap     - keyword - if 1, the source is in the non-targeted aperture
;   verbose - keyword - if 1, a lot of information printed to screen


if (keyword_set(xap) eq 0) then xap=0
if (keyword_set(verbose) eq 0) then verbose=0 

; extract needed keywords from FITS header

module   =  sxpar(header,'CHNLNUM') ; 0=SL, 1=SH, 2=LL, 3=LH
fovid    =  sxpar(header,'FOVID')   ; FOV ID number
ra_fov   =  sxpar(header,'RA_FOV')  ; RA at FOV
dec_fov  =  sxpar(header,'DEC_FOV') ; Dec at FOV
pa_fov   =  sxpar(header,'PA_FOV')  ; position angle of slit at FOV
ra_slt   =  sxpar(header,'RA_SLT')  ; RA at slit center
dec_slt  =  sxpar(header,'DEC_SLT') ; Dec at slit center
pa_slt   =  sxpar(header,'PA_SLT')  ; position angle at slit center
ra_xslt  =  sxpar(header,'RA_XSLT') ; RA at center of other slit 
dec_xslt =  sxpar(header,'DEC_XSLT'); Dec at center of other slit
pa_xslt  =  sxpar(header,'PA_XSLT') ; position angle at other slit center

if (verbose eq 1) then begin 
  print,'RA  (FOV, SLT, XSLT)',ra_fov, ra_slt, ra_xslt
  print,'DEC (FOV, SLT, XSLT)',dec_fov,dec_slt,dec_xslt
  print,'PA  (FOV, SLT, XSLT)',pa_fov, pa_slt, pa_xslt
endif

; convert ra, dec, pa to generic variables
; use SLT values, unless xap set, in which case, use XSLT values

if (xap eq 0) then begin
  ra  = ra_fov   & dec = dec_fov   & pa  = pa_fov
endif else begin
  ra  = ra_xslt  & dec = dec_xslt  & pa  = pa_xslt
endelse
; convert extpos to an offset in arcsec
; pixel scale based on module, nominal nod pos'n based on fovid
;   see notes 28 Nov 11 (and 21 Aug 09) for nominal nod pos'ns

case module of ; scale is in arcsec per pixel
  0 : scale = 1.80 
  2 : scale = 5.10
endcase

if (xap eq 0) then nominal=nominalpos(fovid) else $
                   nominal=nominalpos(fovid,/XSLT)

offset = double ( scale * (extpos-nominal) ) / 3600.0 ; offset now in degrees

; find the extracted RA and DEC by finding deltas
;   and then shifting from RA and DEC (which is usually from FOV)

rascale = cos(!PI * dec / 180.0)

ra_delta  = offset * cos(!PI * (90.0-pa) / 180.0) / rascale
dec_delta = offset * sin(!PI * (90.0-pa) / 180.0)

ra_ext  = ra  + ra_delta
dec_ext = dec + dec_delta

; write EXTPOS and the new extracted coordinates to the FITS header

sxaddpar,header,'EXTPOS',extpos,$                    ; extracted pos'n in pix
  ' [pix] Position of extracted spectrum in row 0',before='DATE   '
sxaddpar,header,'RA_EXT',ra_ext,$                    ; extracted RA in degrees
  ' [deg] RA of extracted position of spectrum',before='DATE   '
sxaddpar,header,'DEC_EXT',dec_ext,$                  ; extracted Dec in degrees
  ' [deg] DEC of extracted position of spectrum',before='DATE   '

if (verbose eq 1) then begin
  print,'NOMINAL ',nominal,' EXTRACTED ',extpos,' OFFSET ',offset
  print,'RA  (FOV/XSLT, DELTA, EXT): ',ra,ra_delta, ra_ext
  print,'DEC (FOV/XSLT, DELTA, EXT): ',dec,dec_delta, dec_ext
endif

END
