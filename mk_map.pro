PRO mk_map,module,ordermap,lambdamap,posnmap,filename

;  1 Mar 11 updated all wavsamp files to S19 (duplicating S18.18 files)
; 23 Jun 09 updated all wavsamp files to S18, updated pixel scales
; 22 Apr 07 updated to add filename parameter
; 17 Apr 07 updated all wavsamp files to S15
;  7 Apr 06 updated all wavsamp files to S13
; 12 Nov 05 updated LL to S13
; 18 Aug 05 ensured use of latest lores wavesamps
; 27 Aug 04 moved hardcoded path information to top of file for ease of editing
; 29 Jan 04 changing wavsamp file version numbers
;           this really shouldn't be hardcoded!
; 13 Jan 04 modified module so that the numbers run 0-3 and not 1-4
; 28 Nov 03 modified to use rd_ipac_tbl and read new wavsamp files
; 23 Jan 03 playing with width of order
; 16 Jan 03 new paths to wave-samp files
; 15 Jan 03 fixed bug in sign of position

; mk_map makes maps of order number, wavelength, and slit position for the IRS
; 
; INPUT 
;   module   - the module number, 0 SL, 1 SH, 2 LL, 3 LH
;   filename - optional - the name of the wavesamp file to read
;              if not supplied, then one of the hardcoded files below is used
; OUTPUT
;   ordermap  - 2D image - order number for each pixel
;   lambdamap - 2D image - wavelength (in um) for each pixel
;   posnmapp  - 2D image - wavelength (in um) for each pixel

; read and set general parameters

fn=['b0_wavsamp_19.tbl','b1_wavsamp_19.tbl',$
    'b2_wavsamp_19.tbl','b3_wavsamp_19.tbl']

if (keyword_set(filename) eq 0) then path='/home/sloan/irs/calfiles/'$
else path='./'
if (keyword_set(filename) eq 0) then filename = fn[module]

side=128
deg=2

; assign relevant values for reading data
; pixel scales based on IRS pocket guide at SSC website examined 23 Jun 09

case module of
  0 : begin ; SL
;        len=253
        scale=1.80 & minorder=1 & maxorder=3
;        dumstop=3
      end
  1 : begin ; SH
;        len=998
        scale=2.30 & minorder=11 & maxorder=20
;        dumstop=3
      end
  2 : begin ; LL
;        len=260
        scale=5.10 & minorder=1 & maxorder=3
        dumstop=3
      end
  3 : begin ; LH
;        len=1147
        scale=4.50 & minorder=11 & maxorder=20
;        dumstop=4
      end
endcase

; read in wavesamp file

a=rd_ipac_tbl(path+filename)

; zero arrays

ordermap = intarr(side,side)
lambdamap= dblarr(side,side)
posnmap  = dblarr(side,side)

c_yx = dblarr(deg+1,maxorder+1)
c_xy = dblarr(deg+1,maxorder+1)
c_sl = dblarr(deg+1,maxorder+1)
c_lx = dblarr(deg+1,maxorder+1)
c_ly = dblarr(deg+1,maxorder+1)

; for each wavesamp file
; find top and bottom of box and slope of line of constant wavelength

xtop =(a[5,*]+a[11,*])*0.5
ytop =(a[6,*]+a[12,*])*0.5
xbtm =(a[7,*]+a[ 9,*])*0.5
ybtm =(a[8,*]+a[10,*])*0.5
slope=((ytop-ybtm)/(xtop-xbtm))

; increment through orders and find coefficients of various fits
; for future use

for i=minorder,maxorder do begin
  ind_m=where(a[0,*] eq float(i))                   ; indices for order
  c_yx[*,i] = poly_fit(a[1,ind_m],a[2,ind_m],deg)   ; y as a f'n of x
  c_xy[*,i] = poly_fit(a[2,ind_m],a[1,ind_m],deg)   ; x as a f'n of y
  c_sl[*,i] = poly_fit(a[3,ind_m],slope[ind_m],deg) ; slope vs. lambda
  c_lx[*,i] = poly_fit(a[1,ind_m],a[3,ind_m],deg)   ; lambda vs. x
  c_ly[*,i] = poly_fit(a[2,ind_m],a[3,ind_m],deg)   ; lambda vs. y
endfor

; increment through orders and fill ordermap array

for i=minorder,maxorder do begin

  ind_m=where(a[0,*] eq float(i))                   ; indices for order
;
; fit x(y) for top and bottom of slit
; step through all rows (y)
; find x for bottom and top of slit, round down and up respectively
; fill ordermap with all x's in between with value of order

  cxyb=poly_fit(ybtm[ind_m],xbtm[ind_m],deg,fxyb)
  cxyt=poly_fit(ytop[ind_m],xtop[ind_m],deg,fxyt)

  jstart=fix(min([min(ybtm[ind_m]),min(ytop[ind_m])])-0.5)
  if (jstart lt 0) then jstart=0
  jstop=fix(max([max(ybtm[ind_m]),max(ytop[ind_m])])+0.5)
  if (jstop gt side-1) then jstop=side-1

  for j=jstart,jstop do begin
    y=double(j)
    x0= cxyb[0] + cxyb[1]*y + cxyb[2]*y*y
    xstart = fix(x0 + 0.0) ; -0.5
    if (xstart lt 0) then xstart=0
    x1 = cxyt[0] + cxyt[1]*y + cxyt[2]*y*y
    xstop  = fix(x1 + 0.0) ; +0.5
    if (xstop gt side-1) then xstop=side-1
    ordermap[xstart:xstop,j]=i
  endfor

endfor

; pass through all pixels, find order for that pixel, then iterate to lambda
for i=0,side-1 do begin
  x=double(i)+0.5

  for j=0,side-1 do begin
    y=double(j)+0.5
    m=ordermap[i,j]

    iter=0
    if (m gt 0) then begin ; pixel exposed to valid order, continue

;     the initial guess for lambda is the lambda at same y on spectrum

      l_0 = poly(y,c_ly[*,m])
      l_old = l_0 + 1.0

      while (abs(l_0-l_old) gt 0.001) do begin ; iterate until delta is small
        l_old = l_0
        slope_0 = poly(l_0,c_sl[*,m])         ; find slope for test lambda
        if (slope_0 ne 0.0) then begin
          islope_0 = 1.0/slope_0              ; find inverse slope
          x_int   = x - islope_0*y            ; fit islope through fixed x and y
        endif
        if (deg eq 2 and slope_0 ne 0.0) then begin

;         find intersections of line and quadratic if slope not 0

          aa = c_xy[2,m]
          bb = c_xy[1,m] - islope_0
          cc = c_xy[0,m] - x_int

          operand=bb*bb-4.0*aa*cc
          if (operand gt 0.0) then operand=sqrt(operand) else begin
            lambdamap[i,j]=0.0
            repeat_flag=0
          endelse
          y_1 = (-bb + operand)/(2.0*aa)
          y_2 = (-bb - operand)/(2.0*aa)

          l_1 = poly(y_1,c_ly[*,m])
          l_2 = poly(y_2,c_ly[*,m])

;         pick solution with the wavelength closest to l_0  

          if (abs(l_0-l_1) lt abs(l_0-l_2)) then begin
            l_0 = l_1
            y_0 = y_1
          endif else begin
            l_0 = l_2
            y_0 = y_2
          endelse

        endif else begin

          l_0 = poly(y,c_ly[*,m])
          y_0 = y

        endelse

      endwhile

      lambdamap[i,j] = l_0
      x_0 = poly(y_0,c_xy[*,m])
      if (x ge x_0) then psign=1.0 else psign=-1.0 ; above/below slit center
      posnmap[i,j]   = psign * scale * sqrt ((x-x_0)^2 + (y-y_0)^2)

    endif else begin ; pixel not illuminated by any order
      lambdamap[i,j] = 0.0
      posnmap[i,j]   = 0.0
    endelse

  endfor
endfor

end
