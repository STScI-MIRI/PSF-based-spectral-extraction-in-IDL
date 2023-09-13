FUNCTION pos_faint,inimage,error,bmask,psfimage,subpix,centcol,tracecoeff,omap,order,lmap,position=position,colrange=colrange,nfind=nfind,rowrange=rowrange,diag=diag,residual=residual,fat=fat,slow=slow,fast=fast,_extra=e

; NOTE - can't pass degree as _extra=e in the current configuration
; 19 Aug 12 changed initial conditions for bestchi and improved diagnostics
; 17 Feb 12 modified to extract two sources
; 22 Aug 09 added fast and slow keywords, sped iteration back up
; 17 Aug 09 playing with the iteration method
; 10 Aug 09 setting _extra=e in definition and sum_profile call
; 26 Jul 09 can now fattening the PSF with the fat keyword
; 25 Jul 09 using a modified sum_profile to calculate starting profile
; 22 Jul 09 improving ability to use colrange to pick out a faint target
; 21 Jul 09 added rowrange, colrange, residual, diag keywords
; 16 Jul 09 modified to call sum_profile to help find position with min residual
; 10 Jul 09 created
;
; finds the position of a source in an image for a given order
; based on Vianney's algorithm, which recursively iterates through
;   possible positions, seeking the residual image with the minimum flux
; uses the residual image now returned by exspecim.pro
; we're trying to minimize the sum of the residual squared (or abs. value)
;
; INPUT
;   inimage    - 2-D image
;   error      - 2-D error plane
;   bmask      - 2-D bmask
;   psfimage   - 2-D PSF image
;   subpix     - number of subpix/pix in PSF image
;   centcol    - central column of PSF image
;   tracecoeff - polynomial trace coefficients
;   omap       - 2-D order map
;   order      - order under consideration
;   lmap       - 2-D wavelength map
;   position   - keyword for sources with known positions
;   nfind      - keyword - number of sources to find, DEFAULT=1
;   colrange   - keyword to set range of columns to consider (per source)
;   rowrange   - keyword to set range of rows to consider
;   residual   - keyword to return residual image
;   fat        - keyword - ratio to fatten PSF by convolving with a gaussian
;                need one element for each source, to find or fixed
;                values for sources to find come BEFORE fixed sources
;                any extra unspecified elements padded to 1
;   slow       - keyword - force slow build_psf algorithm all the time
;   fast       - keyword - force fast build_psf algorithm all the time
;   diag       - keyword, diagnostic flag
; OUTPUT - returns position(s)

; STEP ONE - SETUP, CHECKING INPUTS, ADJUSTING KEYWORDS

; determine dimensions of input image

image=inimage
ncol=n_elements(image[*,0])
nrow=n_elements(image[0,*])

; check diag keywords

if (keyword_set(diag) eq 0) then diag=0

; set fast1 and fast2, used for search and final residual image, respectively
; default:  call build_psf with fast=1 while searching
;   with fast=2 for final residual image
; adjust these if fast or slow keywords set, note that fast overrides slow

fast1=1 & fast2=0 
if (keyword_set(slow) ne 0) then fast1=0                           ; all slow
if (keyword_set(fast) ne 0) then begin & fast1=1 & fast2=1 & endif ; all fast

; check position, nfind, and colrange keywords to determine nfix and nfind
;   nfix = number of sources with known positions
;   nfind = number of sources to find, can be either 1 or 2 - DEFAULT=1
; if colrange gives more ranges than nfind specifies, it takes priority
; if colrange has fewer rows than nfind, it needs to be padded

if (keyword_set(position) ne 0) then nfix=n_elements(position) else nfix=0
if (keyword_set(nfind) eq 0) then nfind=1
if (n_elements(colrange) gt 1) then begin
  nranges=n_elements(colrange[0,*])
  if (nranges gt nfind ) then nfind=nranges ; bump nfind to reflect colrange
endif else nranges=0
if (nranges lt nfind) then begin ; need to build or pad colrange array
  if (keyword_set(colrange) eq 0) then colrange=[0,ncol]
  for n=0,nfind-nranges-1 do colrange=[[colrange],[0,ncol]]
endif

if (nfind gt 2) then begin
  print,'Warning.  Will only find two sources of the ',nfind,' requested.'
  if (nfix gt 0) then print,'POS ',position else print,'NFIX ',nfixcolrange
  for n=0,nfind-1 do print,colrange[n,*],format='(2f8.4)'
  nfind=2
endif

case nfind of
  0 : print,'Oops!  Nothing to find, but we still called pos_faint.'
  1 : fmt='(i5,i4,i4,4(f9.4),2(e12.4))'
  2 : fmt='(i5,i4,i4,5(f9.4),2(e12.4))'
endcase

; check fat keyword, if fewer elements than npos+nfind, fill it out with 1's

if (n_elements(fat) eq 0) then fat=1.0 
while (n_elements(fat) lt nfix+nfind) do fat=[fat,1.0]

; mask image with order mask, then zero negative data
; NOTE - MASKING IS REDUNDANT CODE - ALREADY DONE

idx=where(omap ne order)
if (idx[0] ne -1) then image[idx]=0.0
idx=where(image lt 0.0)
if (idx[0] ne -1) then image[idx]=0.0

; check rowrange, if not set, set to min and max of non-zero rows

if (n_elements(rowrange) ne 2) then begin
  testcol=fltarr(nrow)
  for i=0,nrow-1 do testcol[i]=total(image[*,i],/nan)
  idx=where(testcol ne 0)
  rowrange=[min(idx),max(idx)]
endif

; zero y_intercept of trace coefficients (this should already be done)

tracecoeff[0]=0.0 

; STEP TWO - PREPARE bestchi, COLRANGE, POSARRAY, PSFPASS, PSF for iteration

; call sum_profile to determine initial bestchi (with no source removed)
;   18 Aug 12 - bestchi now immediately replaced below
;   and generate a collapsed profile used to set up searching ranges

bestchi=sum_profile(image,tracecoeff,medianrow=row,rowrange=rowrange,$,
  diag=diag,_extra=e) 

; shift colrange values so that they do not exceed the range of nonzero data

idx=where(row ne 0)
for n=0,nfind-1 do colrange[*,n] = $
  [ max([colrange[0,n],min(idx)]), min([colrange[1,n],max(idx)]) ]

; initialize several parameters and arrays for the search below
; note that nsteps must be at least 2 elements, even for one source

delta = 0.80 ; step size, fixed initially to 0.80 pixels (from 22 Aug 09)
len = 150    ; 150 should be sufficient, if len > 150 in next loop, we crash!
nsteps   = intarr(max([nfind,2]))+1 ; will hold len for each source to find
posarray = fltarr(len,nfind)      ; will hold search positions for each source
bestpos  = fltarr(nfind)          ; array to hold best position for each source
goflag=1 & count=0                ; to start array

; set tolerance based on fall-off to either side of peak in median row
; tolerance limits the columns considered when searching for a source
; tolerance will be 0.05 unless the peak is a noise spike

for n=0,nfind-1 do begin

  maxp = max(row[colrange[0,n]:colrange[1,n]],maxpos,/nan) ; max within colrange
  maxpos = maxpos+colrange[0,n]   ; reset maxpos to correspond to full row array
  if (maxpos eq 0)      then maxpos = 1
  if (maxpos eq ncol-1) then maxpos = ncol-2

  if (row[maxpos-1] gt 0.05*maxp and row[maxpos+1] gt 0.05*maxp) then $
    tolerance=0.05 else tolerance=0.01

; use tolerance to further limit colrange, which sets the initial search range

  idx = where(row gt tolerance*maxp)
  colrange[*,n] = [max([colrange[0,n],min(idx)]), min([colrange[1,n],max(idx)])]

; use colrange to determine startpos, the min position
; use median row to adjust colrange and build initial array of positions
; posarray covers colrange, first and last elements are colrange[0,n] and [1,n]

  nsteps[n]                 = fix ( (colrange[1,n]-colrange[0,n])/delta ) + 1
  posarray[0:nsteps[n]-1,n] = colrange[0,n] + findgen(nsteps[n])*delta
  posarray[nsteps[n]-1,n]   = float(colrange[1,n])

  bestpos[n] = -99.0
endfor

inrow=row ; for plotting diagnostics - row was generated by sum_profile above
if (diag ge 1) then begin
  plot,row,psym=10,th=2
  print,'tolerance ',tolerance
  print,'iterating over columns ',posarray[0],posarray[len-1]
endif

; prepare psfpass and psf arrays - psfpass used to generate individual PSFs
; one plane per source, either to find or fixed

psfpass = fltarr(n_elements(psfimage[*,0]),n_elements(psfimage[0,*]),nfind+nfix)
psf     = fltarr(ncol,nrow,nfind+nfix)

for n=0,nfind+nfix-1 do begin

; load psfpass for each source, fattening if necessary

  if (fat[n] gt 1.0) then psfpass[*,*,n] = fatten(psfimage,fat[n]) $
    else                  psfpass[*,*,n] = psfimage

; load psf for fixed sources here, since these won't change
; default for call to build_psf here is fast1=1 for less accuracy and speed

  if (n gt nfind-1) then begin
    psfone     = build_psf(reform(psfpass[*,*,n]),subpix,centcol,omap,order,$
                 position[n-nfind],tracecoeff,fast=fast1,_extra=e)
    psf[*,*,n] = psfone
  endif

endfor

; STEP THREE - iterate to find the sources we are looking for
;   currently can search for up to two sources

while (goflag eq 1) do begin

  for i=0,nsteps[0]-1 do begin   ; use i to loop through first source
    for j=0,nsteps[1]-1 do begin ; use j to loop through second source


;   load the planes of the psf array for the sources we're trying to find
;   calling build_psf with fast option for less accurate and faster PSF

      for n=0,nfind-1 do begin
         case n of
           0 : checkpos = posarray[i,n]
           1 : checkpos = [ checkpos , posarray[j,n] ]
         endcase
         psfone = build_psf(reform(psfpass[*,*,n]),subpix,centcol,omap,order,$
                  checkpos[n],tracecoeff,fast=fast1,_extra=e) ; default: fast1=1
         psf[*,*,n] = psfone
      endfor

;   extract the spectra with the positions as set

      spectrum=exspecim(image,error,bmask,psf,lmap,order,$
        residual=residual,diag=diag,_extra=e)

;   call sum_profile to determine the error in the residual image
;   sum_profile rectifies the residual image, generates a smoothed profile
;   and returns the summed squared error

      chisq=sum_profile(residual,tracecoeff,$
        medianrow=row,rowrange=rowrange,diag=diag,_extra=e) 

;   save the positions corresponding to the least error in the bestpos array

      if (chisq lt bestchi or count eq 0) then begin
        bestchi = chisq & bestpos = checkpos
      endif

      if (diag ge 2) then print,$
        'pos_faint - count,delta,checkpos,bestpos,chisq,bestchi ',$
         count,delta,checkpos,bestpos,chisq,bestchi,$
         format='(a55,i3,f7.4,2(f7.2),2(e10.3))'
;      print,count,i,j,delta,checkpos,bestpos,chisq,bestchi,format=fmt
      if (diag eq 3) then begin
        plot,row,th=2,psym=10
        zq=get_kbrd(1)
      endif
      count=count+1

    endfor 
  endfor

; reset step size and range for iterating

  if (delta lt 0.002) then goflag=0
  delta = 0.20*delta     ; 22 Aug 09 - was 0.20

; load new positions and update nstpeps for next iteration
; check posarray for each source and reduce if out of bounds of colrange

  for n=0,nfind-1 do begin
    len=11               ; for 11 steps ; 22 Aug 09 - was 10 for 11 steps
    startpos            = bestpos[n] - len*delta/2
    posarray[0:len-1,n] = startpos + findgen(len)*delta
    nsteps[n]           = len

    idx=where(posarray[0:nsteps[n]-1,n] ge colrange[0,n] $
          and posarray[0:nsteps[n]-1,n] le colrange[1,n])

    if (idx[0] ne -1) then begin
      len                 = n_elements(posarray[idx,n])
      posarray[0:len-1,n] = posarray[idx,n] 
      nsteps[n]           = len
    endif else begin
      posarray[0,n]       = bestpos[n]
      nsteps[n]           = 1
    endelse
  endfor

endwhile

; STEP FOUR - finish up

; diagnostic block

if (diag ge 1) then plot,inrow,psym=10,th=2
if (diag ge 1) then oplot,row,psym=10,th=1.5,col=3
if (diag ge 1) then print,'pos_faint - count,bestpos,bestchi ',$
                    count,bestpos,bestchi

; repeat extraction at best positions in order to return residual
; have to call build_psf for each source we think we found and rebuild psf
; the default for build_psf is now fast2=0, slower and more accurate

for n=0,nfind-1 do begin
  psfone     = build_psf(reform(psfpass[*,*,n]),subpix,centcol,omap,order,$
               bestpos[n],tracecoeff,fast=fast2,_extra=e)
  psf[*,*,n] = psfone
endfor

spectrum   = exspecim(image,error,bmask,psf,lmap,order,$
             residual=residual,diag=0,_extra=e) ; diag turned off

; build returnpos, the array to pass back to the calling procedure
; return the positions of the found sources first, followed by the fixed

returnpos=bestpos
if (nfix gt 0) then returnpos=[bestpos,position]

RETURN,returnpos
END
