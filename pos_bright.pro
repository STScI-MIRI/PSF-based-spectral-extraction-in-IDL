FUNCTION pos_bright,inimage,psfimage,subpix,centcol,tracecoeff,omap,order,posvector=posvector,sigpos=sigpos,range=range

;  8 Jul 09 added posrange keyword to limit the searchable range of positions
;  4 Jul 09 optimizing iteration to hone in on source position
; 22 Jun 09 created
;
; find the position of a source in an image for a given order
;
; INPUT
;   inimage    - 2D image with the source in one order
;   psfimage   - 2D image with the PSF at subpixel resolution in each valid row
;   subpix     - number of subpixels pwer pixel in PSF image
;   centcol    - center column of PSF in psfimage
;   tracecoeff - polynomial trace coefficients defining the spectral trace
;   omap       - order map
;   order      - spectral order containing the spectral trace of the source
;   range      - 2-element vector defining the allowable range of positions
; OUTPUT
;   FUNCTION RETURNS the position of the shifted spectral trace in row 0
;   posvector  - optional, return full vector of positions
;   sigpos     - optional, return std dev of final estimate of position

; determine dimensions of input image

image=inimage
ncol=n_elements(image[*,0])
nrow=n_elements(image[0,*])

; set range if not already defined

if (n_elements(range) ne 2) then range=[0,nrow]

; mask image with order mask, then zero negative data

idx=where(omap ne order)
if (max(idx) gt -1) then image[idx]=0.0
idx=where(image lt 0.0)
if (max(idx) gt -1) then image[idx]=0.0

; setup

tracecoeff[0]=0.0 ; zero y_intercept of trace coefficients (just in case!)
posvector=fltarr(nrow)
chivector=fltarr(nrow)

; iterate through rows (dropping top and bottom rows)
; normalize each row to one
; find position of spectral trace in each row iteratively
; save position and chi_squared value of best fit

for j=1,nrow-2 do begin

  imrow=image[*,j]

  if (total(imrow) ne 0) then begin

;   normalize row

    imrow=imrow/total(imrow,/nan)

;   take starting position to be the maximum column

;    gfit=gaussfit(findgen(ncol),imrow,cc,nterms=3)
;    startpos=cc[1]+0.5

    maxval  = max(imrow,maxpos)
    startpos=maxpos-0.5 
    delta=2.5

;   begin recursive iteration, each will be nine steps over +/- delta
;   now resetting all position values outside range to nearest edge

    goflag=1 & istart=-4 & istop=4 & count=0
    while (goflag eq 1) do begin

      for i=istart,istop do begin

        position = startpos + i*delta/istop
        if (position gt max(range)) then position=max(range)
        if (position lt min(range)) then position=min(range)
      
        psf=build_psfrow( reform(psfimage[*,j]),subpix,centcol, $
                          reform(omap[*,j]),order,position )
        residual=imrow-psf
        chisq=total(residual^2,/nan)

        if (i eq istart and count eq 0) then begin
          bestchi=chisq & bestpos=position
        endif else if (chisq lt bestchi) then begin
          bestchi=chisq & bestpos=position
        endif

        count=count+1
      endfor
     
      startpos=bestpos 
      delta=0.11*delta
      if (delta lt 0.001) then goflag=0 
    endwhile

;   load results

    posvector[j]=bestpos
    chivector[j]=bestchi

  endif
endfor

; use trace coefficients to convert positions to estimates of source position

spectrace = poly(findgen(nrow),tracecoeff)
posvector = posvector - spectrace

; renormalize chi vector so that the mean=1

chivector = chivector/mean(chivector,/nan)

; first pass - find mean position for data where position is non-zero

idx     = where(posvector gt 0.0)
meanpos = mean(posvector[idx],/nan)
sigpos  = stddev(posvector[idx],/nan)

; second pass - find mean position for data where position within a std dev

idx     = where(posvector gt 0.0 and abs(posvector-meanpos) lt sigpos)
meanpos = mean(posvector[idx],/nan)
sigpos  = stddev(posvector[idx],/nan)


RETURN,meanpos
END
