FUNCTION exspecrow,datarow,errow,bmrow,psfrow,lamrow,order,DEGREE=degree,SKY=sky,WTMETHOD=wtmethod,ROW=row,DIAG=diag,_extra=e

; 16 Feb 12 enabled sky removal by using n_elements to check degree keyword
; 25 Aug 09 now dividing weighted wavelength sum by total of profile
; 24 Aug 09 now returning sky polynomial values
; 20 Aug 09 added row keyword, default degree now -1, stopped normalizing PSF
; 22 Jul 09 modifications to run properly for multiple sources
; 13 Jul 09 disabling diagnostics
; 29 Jun 09 modified to receive and use a row of errors and bmask values
; 22 Jun 09 heavily modified to use Vianney's new fitting kernel
; 17 Jun 08 implementing some improvements from SMART version of this routine
;  2 May 07 don't load data unless sum PSF > 0.1, added _extra keyword
; 26 Apr 07 renamed algorithm keyword to exmethod, added diag keyword
; 21 Apr 07 created in its original form
;
; exspecrow measures the flux from ONE spectrum in one row of an image
;
; INPUT
;   datarow  - 1D - relevant row of image, out-of-order information is masked
;   errow    - 1D - relevant row of errors, out-of-order information is masked
;   bmrow    - 1D - relevant row of bmasks, out-of-order information is masked
;   psfrow   - 2D - matching PSF laid down to match position of source in data
;            - one PSF is supplied for each source to be extracted
;   lamrow   - 1D - corresponding row of wavelength map
;   resrow   - 1D - residuals after subtraction of fitted PSFs and background
;   order    - the order being extracted
;   degree   - the degree of the polynomial to fit to the background
;   bmval    - keyword to set bmval, defaults to 4096
;   wtmethod - keyword to specify the weighting method for extraction
;            - method 0 = all weights are 1 
;            - method 1 = simple weights taken from errow and bmrow - DEFAULT
;   row      - keyword to pass row number from calling routine for diagnostics
;   diag     - keyword to specify diagnostic mode - LIMITED SUPPORT
;
; OUTPUT    - an array with five elements for each source to be extracted
;             col 0 - wavelength (um), col 1 - flux (DN), col 2 - error (DN)
;             col 3 - order, col 4 - mask information
;             for now, mask information is set to zero 
;
; PLANNED IMPROVEMENTS
;   pass bmval, allow wtmethod=1
;
; algorithm
; calls opt_mregress, based on sm_mregress
; ensure that PSF is normalized so that sum=1
; generate arrays to pass to opt_mregress:  psfpass, datapass, weight
; weights will ultimately be adjustable, but not in this test version
; assemble results of opt_mgregress call, return results to calling procedure
;
; caveats
; no error analysis is performed, errors are set to zero
; no masking is performed - invalid rows are set to NaN for later filtering

; check keywords, reset wtmethod and bmval if necessary

if (n_elements(degree) eq 0) then degree=-1           ; 0 means flat background
if (keyword_set(diag) eq 0) then diag=0 
if (keyword_set(wtmethod) eq 0) then wtmethod=0       ; set if not set
if (wtmethod lt 0 or wtmethod gt 1) then wtmethod=0   ; reset if illegally set
if (keyword_set(bmval) eq 0) then bmval=4096
if (n_elements(row) eq 0) then row=-1

; check input data to ensure they are nonzero

if (total(datarow) eq 0) then goflag=0

; check input PSF row for dimensionality, set goflag=0 if illegal
; create psfpass and load PSF row into it

sz=size(psfrow)
case sz[0] of
  1 : begin
    incol=n_elements(psfrow)
    npsf=1
    psfpass=dblarr(incol,npsf+degree+1)
    psfpass[*,0]=psfrow
    goflag=1
  end
  2 : begin
    incol=n_elements(psfrow[*,0])
    npsf=n_elements(psfrow[0,*])
    psfpass=dblarr(incol,npsf+degree+1) 
    psfpass[*,0:npsf-1]=psfrow
    goflag=1
  end
  else : goflag=0
endcase

; define columns for output spectral vector (hardcoded)

outcol=5 & lcol=0 & fcol=1 & ecol=2 & ocol=3 & mcol=4

specdata=dblarr(outcol,npsf) ; output array
specdata[ocol,*]=order       ; load order information immediately

if (goflag eq 1) then begin

; clean psfpass array - check that PSF rows are nonzero
; npsf rows already in psfpass, one for each PSF - normalize each row to 1
; load additional rows with col^degree for background fit, but only if

;  for n=0,npsf-1 do begin                        ; go through PSF(s)
;    if (total(psfpass[*,n]) eq 0) then goflag=0  ; goflag meaningless here
;  endfor
  if (degree gt -1) then $                        ; prepare psfpass for sky
    for m=0,degree do psfpass[*,m+npsf] = dindgen(incol)^m  

; set weights based on wtmethod

  case wtmethod of
    0 : weight = dblarr(incol)+1                  ; uniformly set to unity
    1 : weight = double(bmrow lt bmval) / errow 
  endcase

; create datapass vector - modified datarow to pass to opt_mregress
; need to zero first and last columns of data, and all NaN data

  datapass=double(datarow)
  idx=where(finite(datapass)-1)           ; will flag NaNs
  if (idx[0] gt -1) then idx=[0,idx,incol-1] else idx=[0,incol-1]
  datapass[idx] = 0                       ; set relevant data to zero
  weight[idx]   = 0                       ; set corresponding weights to zero

; call opt_mregress, returns the flux of the PSFs and background coefficients

  coeff=opt_mregress(psfpass,datapass,weights=weight,status=stat,sigma=err)

; load specdata array, one PSF at a time (order column already loaded)

  for n=0,npsf-1 do begin

; find the wavelength, PSF-weighted average of the wavelengths in the row
; have to divide by total of weights because there is no guarantee that sum = 1

    specdata[lcol,n] = total(psfpass[*,n]*lamrow) / total(psfpass[*,n])

; load data into flux and error files

    specdata[fcol,n] = coeff[n]
    specdata[ecol,n] = err[n]

  endfor

; if background set, then load sky polynomial

  if (degree gt -1) then begin
    sky=fltarr(degree+1)
    for m=0,degree do sky[m] = coeff[m+npsf]
  endif

; error analysis to be added later

  plotflag=1

endif else begin ; invalid data or PSF, set wavelength and flux to NaN

  plotflag=0
  specdata[lcol,*] = !VALUES.F_NAN
  specdata[fcol,*] = !VALUES.F_NAN
  resrow=datarow   ; residuals = input since nothing found or subtracted

endelse

if (diag eq 4 and row eq 64) then begin
  model=fltarr(n_elements(datapass))
  for n=0,npsf-1 do model = model + coeff[n]*psfpass[*,n]
  plot,datapass-model,psym=10,th=2
endif

;if (diag eq 5 and plotflag eq 1) then begin ; this block NOT set for npsf > 1
;  plot,datarow,xsty=1,title='order '+string(order,format='(I2)'),_extra=e
;  oplot,flux*psf,th=2
;  wait,0.15
;endif

RETURN,specdata
END
