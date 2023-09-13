FUNCTION sum_profile,image,tracecoeff,medianrow=medianrow,rowrange=rowrange,diag=diag,plan=plan

; 16 Feb 12 improved comments
; 17 Aug 09 modifying default plan back to 2
; 10 Aug 09 new plan keyword - can change how returned profile is calculated
; 25 Jul 09 now returning the spatial profile with the medianrow keyword
; 20 Jul 09 added rowrange keyword to limit range of rows for chi-squared sum
; 13 Jul 09 created
;
; given an image and polynomial coefficients for a spectral trace
;   find the error (sum of absolute value or chi-squared error)
; called by pos_faint.pro
;
; INPUT
;   image - 2D residual image
;   tracecoeff - polynomial coefficients defining the spectral trace
;   rowrange   - keyword to limit range over which chi-squared error is summed
;   medianrow  - keyword to return rectified profile to calling routine
;   plan       - keyword to set how returned profile is calculated
;                1 = sum of abs. value, 2 = sum of profile squared (default)
;   diag       - diagnostic keyword - not supported
; OUTPUT
;   summed squared profile after rectifying and smoothing residual image

; check keywords and image parameters

if (keyword_set(plan) eq 0) then plan=2
if (keyword_set(diag) eq 0) then diag=0
ncol=n_elements(image[*,0])
nrow=n_elements(image[0,*])

; colrange set to far L and far R of array
;   other routines (such as pos_faint) require this

col0=0 & col1=ncol-1 

; check rowrange, can be an input, otherwise set to top and bottom of array

if (n_elements(rowrange) eq 2) then begin
  row0=rowrange[0] & row1=rowrange[1]
endif else begin
  row0=0 & row1=nrow-1
endelse

; if (diag eq 1) then print,'sum_profile rowrange ',rowrange

; rectify image
; shift each row by an amount determined from trace coefficients
; this shifts the profile to row 0

shifted=fltarr(ncol,nrow)
for j=row0,row1 do begin ; need to skip NaNs
  row = reform(image[*,j])
  x = findgen(ncol)+poly(float(j),tracecoeff)
  shifted[*,j] = interpol(row,findgen(ncol),x,/spline)
endfor

; median smooth each column and find the median value - save in finalrow vector

medianrow=fltarr(ncol)

for i=col0,col1 do begin
  col = reform(shifted[i,row0:row1])
  medianrow[i] = median(medsmooth(col,3))
endfor

; return the sum of the median-smoothed finalrow vector
; DEFAULT = plan=2 = total of medianrow squared
; if plan=1, then return total of absolute value of medianrow

if (plan eq 1) then checkrow=abs(medianrow) else checkrow=medianrow^2

RETURN,total(checkrow,/nan)
END
