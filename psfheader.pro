FUNCTION psfheader,header,subpix,centcol,minorder

; 22 Jun 09 modified to account for new FITS header construction
;  1 May 07 created
;
; psfheader reads a FITS header for a PSF file and loads needed values
; returns a 1D tracecoeff array if the PSF file is 2D, otherwise returns 2D
;
; INPUT
;   header - the FITS header
;
; OUTPUT
;   tracecoeffs - return value from the function
;                 a 2D or 3D array of polynomial trace coefficients
;                 separate set of coefficients for each nod in each order
;   subpix      - the number of subpixels per pixel
;   centcol     - the central column
;   minorder    - the order number of the plane 0 in the image

; first determine if the PSF image is 2D or 3D

naxis=sxpar(header,'NAXIS')
if (naxis eq 2) then nplanes=1 else nplanes=sxpar(header,'NAXIS3')

; read the polynomial degree for each plane

for j=0,nplanes-1 do begin
  label='TRACED'+string(j,format='(I2.2)') ; can have over 10 sets of coeff.
  deg = sxpar(header,label)
  if (j eq 0) then degree=deg else degree=[degree,deg]
endfor

; create the tracecoeffs array - needs to be big enough for max degree

maxdegree=max(degree)
tracecoeffs=dblarr(maxdegree+1,2,nplanes)

; read the polynomial coefficients for each plane into tracecoeffs
; trace coeff keywords have format TRACEOOP
; OO=order-minorder, P=power of coefficient

for j=0,nplanes-1 do for i=0,maxdegree do begin
  label='TRACEA'+string(j,format='(I1)')+string(i,format='(I1)')
  tracecoeffs[i,0,j] = sxpar(header,label)
  label='TRACEB'+string(j,format='(I1)')+string(i,format='(I1)')
  tracecoeffs[i,1,j] = sxpar(header,label)
endfor

; read the other needed keywords

subpix     = sxpar(header,'SUBPIX')
centcol    = sxpar(header,'CENTCOL')
minorder   = sxpar(header,'MINORDER',COUNT=mincount)

RETURN,reform(tracecoeffs)
END
