PRO mk_script,infile,outfile,option,SUBOPTION=suboption,CAMPAIGN=campaign,MANY=many,CLEAN=clean,ZERO=zero,POS=pos,DIFF=diff,CMASK=cmask,ROGUE=rogue,CLEANOUT=cleanout,CROSS=cross,SKYIN=skyin,CLUSTERSKY=clustersky,FULL=full,EXKEY=exkey,PRKEY=prkey,OFFLINE=offline,COADD=coadd

;  7 Dec 11 added coadd keyword for imspex
; 10 Feb 11 modified extract option to print script heeders needed for CUPID
;           and added offline keyword to avoid these
;  6 Oct 10 added clustersky suboption for diff mode - use for clusters
; 10 Mar 10 improved handling of skies in diff/for mode so that skies can be
;           used in some EXPIDs and not others within an AOR
;  3 Dec 09 repaired a bug with SL nod diff
; 25 Aug 09 repairing minor bug in imspex block (for when pos'n defined)
; 22 Jul 09 updated imfspex calls since the syntax has changed
;  8 Jul 09 added imfspex code
;  7 Apr 09 corrected a bug for 'diff' and 'for' modes (for diff modes for nods)
; 18 Mar 09 cleaning up modes and output for 'diff' and 'for'
;  3 Feb 09 finally got the 'extract' option running
; 27 Jan 09 created
;           based on mk_bat_file_v4.pro dated 24 Apr 08
;
; INPUT
;   infile     - input module list file - space delimited with columns:
;                target, pointing, aorkey, expid, dcenum, module, aperture/order
;                optional columns:  difference mode
;   outfile    - name of output script file
;   option     - string which determines the nature of the script
;   suboption  - string to determine how the script will be written
;   campaign   - IS THIS USED?
;   many       - to warn mk_script that the EXPID is 3 digits, not 2
;   clean      - keyword for option 'diff' or 'for' - diff cleaned images - 'ic'
;   zero       - keyword for option 'diff' or 'for' - set DCE 00 for neg images
;   diff       - keyword for option 'clean' - clean differenced images - 'id'
;   rogue      - keyword for option 'clean' - specifies rogue file
;   cmask      - keyword for option 'clean' - specifies cmask file
;   cleanout   - keyword for option 'clean' - set cleaned suffix, default='c'
;   pos        - keyword for option 'diff' or 'for' - MORE ...
;   exkey      - keyword for option 'extract' - set im_ex and bm_ex file names
;   prkey      - keyword for option 'extract' - set im_pr and bm_pr file names
;   offline    - keyword for option 'extract' - suppress CUPID script headers
;   skyin      - keyword for option 'imfsky' or 'diff'
;   clustersky - keyword for option 'diff'
;   full       - keyword for option 'exfirs' - full slit extraction
;   coadd      - keyword for option 'imspex' - changes imspex input files
;
; options and valid suboptions and additional keywords
;   diff     - generate a script of imdiff calls
;              generally, the difference mode is read from the module list file
;              if the diff mode column is not present, the default is nod diff
;              the diff mode can be overridden with the following suboptions:
;              'none', 'nod', 'ap', 'cross', 'sky' (corresponds to 0-4)
;              note that the software assumes that the DCEs are available in
;                the relevant EXPIDs for the various diff modes
;              other allowed keywords:
;              '/clean' to diff cleaned images (ic.fits instead of im.fits)
;              '/zero' to set DCE of negative image to zero
;              'skyin' sets the sky image manually
;              'clustersky' - text from "_" is replaced with "sky"
;              pos = no. of pos'ns per nod to count EXPIDs in nonstandard AORs 
;   for      - generate a script of imfor calls - see 'diff' above for details
;   clean    - generate a script of imclean calls
;              '/diff' to clean differenced images (id.fits instead of im.fits)
;              '/cmask' for code to specify a campaign mask
;              '/rogue' for code to use the input file in a call to imrogue
;              cleanout to set last character of output suffix
;                (default is 'c' in 'ic')
;   extract  - writes commands to run offline pipeline at Cornell
;                echo file root (and again after ridge call below)
;                copy image and bmask files into correct dummy file names
;                calls to profile, ridge, and extract
;                copy dummy output files up a directory and with file root
;              exkey - keyword to modify extract file suffix - single letter
;              prkey - keyword to modify profile file suffix - single letter
;              in both cases, letter could be "c" for "ic.fits", "bc.fits", etc.
;                using "m", "d", "x" typically
;   convert  - converts IPAC table files to spectral FITS files
;   offset   - calls find_offset to locate position of source in slit
;   imsky    - generate a script of imfsky calls
;              skyin keyword to reset last character of input suffix
;                (default is 'n' in 'in')
;   exfirs   - generate a script of exfirs calls - extraction software for IOC
;              ',/full' for full-slit extraction
;   findpos  - 0 - NOT YET WRITTEN!!
;   imspex   - generate a script of imspex calls
;              '/coadd' to extract from one image per EXPID, DCE = 'cc'
;
; module list file format
; col 0 - target name
; col 1 - pointing - either a single letter or campaign number and a letter
; col 2 - AOR key
; col 3 - 4-digit EXPID
; col 4 - 4-digit DCENUM
; col 5 - module number (SL 0, SH 1, LL 2, LH 3)
; col 6 - order/aperture number (1-2 for SL and LL, 0 for SH and LH)
; col 7 - optional, diff mode (none 0, nod 1, ap 2, cross 3, sky 4)
; col 8 - optional, position of source in row 0 (in pixels) for imspex calls
;   will change to string indicating position(s) for imspex calls
; additional columns are for additional sources - to be modified
; a position of -N means that imspex will be called with a N-order background
;   - to be modified (but backgrounds not yet tested within imspex)
; note that if the position is included, col 7 must be filled with a diff mode

if (keyword_set(pos) eq 0) then pos=1 ; set no. of positions per nod, default=1

; open input files

openr,fi,infile,/get_lun
openw,fo,outfile,/get_lun

; set module based on infile, the module input file

modname=strmid(infile,0,2)
if (strcmp(modname,'sl',2) eq 1) then module=0
if (strcmp(modname,'sh',2) eq 1) then module=1
if (strcmp(modname,'ll',2) eq 1) then module=2
if (strcmp(modname,'lh',2) eq 1) then module=3

if (module eq 1 or module eq 3) then hiresflag=1 else hiresflag=0
if (module eq 2) then gflag=1 else gflag=0 ; flip glo,ghi in ridge calls for LL

; set wavelength grid file and module string based on module

gridfile = 'cal/'+modname+'.lam.fits'
modstr   = string(module,format='(I1)')

; options 'diff' and 'for' use same code, so set 'for' to 'diff' with forflag=1 
; if option 'diff' (or 'for' on input), set default diff modes for modules

if (option eq 'for') then begin
  option='diff' & forflag=1
endif else forflag=0
if (option eq 'diff') then diffdefault=[2,4,1,4] ; SL ap, SH sky, LL nod, LH sky

; various initializations before reading input list

line=' ' & testline=' '
s=' '    & p='.'        & q=string(byte(39))
c=','    & f='.fits'
cq=','+q & fq='.fits'+q
old_expid='99'
locroot  ='extract/'
target_old='_____'
expid_old='9999'
foundflag=0 
count=0

; parse input list, build elements to complete file names, and print

while (not eof(fi)) do begin

; parse input line

  readf,fi,line
  spl_line=strsplit(line,' ',/extract)
  target  =spl_line[0]
  ptg     =spl_line[1] ; pointing
  aorkey  =spl_line[2]
  expid   =spl_line[3]
  if (keyword_set(many) eq 0) then expid_sh=strmid(spl_line[3],2,2) else $
    expid_sh=strmid(spl_line[3],1,3)
  expnum=fix(expid_sh)
  dce     =spl_line[4]
  dce_sh  =strmid(spl_line[4],2,2)
  module  =spl_line[5]
  order   =spl_line[6]

; input line may also contain spl_line[7], check in the diff option block
; input line may have additional position fields, check in imspex option block

; option blocks

  case option of

    'diff' : begin ; suboptions 0,1,2,3, 8 (4-7,9 = 0-3,8 with cleanflag=1)

;     check for presence of spl_line[7] to determine diff mode
;     if not present, set to default for module SL-ap SH-sky LL-nod LH-sky
;     if the suboption keyword contains a legal diff mode, override other modes

      if (n_elements(spl_line) gt 7) then diff=fix(spl_line[7]) $
        else diff=diffdefault[fix(module)]
      if (keyword_set(suboption) ne 0) then case suboption of
        'none'      : diff=0
        'nod'       : diff=1
        'ap'        : diff=2
        'aper'      : diff=2 ; alternative spelling
        'cross'     : diff=3
        'sky'       : diff=4
        else        :        ; do nothing
      endcase

;     setup for new target
;     set foundflag=0 and and record first expnum as expstart

      if (target ne target_old) then begin
        foundflag=0
        expstart=expnum
      endif

;     various checks
;     check forflag and set iroot and command accordingly
;     check clean keyword and modify iroot accordingly
;     check zero keyword and set dce_neg accordingly (either dce_sh or '00')
;     if diff=4 ('sky') check target, if it's a sky, then reset to nod diff

      if (forflag eq 0) then begin
        iroot='i' & command='imdiff'
      endif else begin
        iroot='b' & command='imfor'
      endelse
      if (keyword_set(clean) eq 0) then filetype = iroot+'m' $
        else filetype=iroot+'c'
      if (keyword_set(zero) eq 0)  then dce_neg=dce_sh else dce_neg='00'
      if (diff eq 4 and strpos(target,'sky') ne -1) then diff=1

;     set defaults

      postype=filetype & negtype=filetype
      posfile=locroot+target+p+ptg+order+p+expid_sh+dce_sh+p+postype+fq
      ordernum=fix(order)

;     prepare script depending on diff mode

      case diff of 
        0 : begin ; no diff
	  negfile='nosubtraction'+q
          difftype=iroot+'d'
        end
        1 : begin ; nod diff
          dex=pos*1
          if ((expnum/2)-(expnum/2.0) eq 0) then dex=dex else dex=-dex
          expnum=expnum+dex
          if (expnum lt 10) then expid_ne='0'+string(expnum,format='(I1)') $
            else begin
              if (expnum lt 100) then expid_ne=string(expnum,format='(I2)') $
              else expid_ne=string(expnum,format='(I3)')
            endelse
          negfile=locroot+target+p+ptg+order+p+expid_ne+dce_neg+p+negtype+fq
          difftype=iroot+'d'
        end 
        2 : begin ; aper diff
          dex=pos*2
          ordernum=fix(order)
          if (order eq 2) then begin 
             expnum=expnum+dex
             ordernum=ordernum-1
          endif else begin 
             expnum=expnum-dex
             ordernum=ordernum+1
          endelse
          if (expnum lt 10) then expid_ne='0'+string(expnum,format='(I1)') $
            else expid_ne=string(expnum,format='(I2)')
          order_ne=string(ordernum,format='(I1)')
          negfile=locroot+target+p+ptg+order_ne+p+expid_ne+dce_neg+p+negtype+fq 
          difftype=iroot+'d'
        end
        3 : begin ; cross diff
          dex1=pos*3 & dex2=pos
          if (order eq 2) then begin 
            if ((expnum/2)-(expnum/2.0) eq 0) then delex=dex1 else delex=dex2
            expnum=expnum+delex
            ordernum=ordernum-1
          endif else begin 
            if ((expnum/2)-(expnum/2.0) eq 0) then delex=-dex2 else delex=-dex1
            expnum=expnum+delex
            ordernum=ordernum+1
          endelse
          if (expnum lt 10) then expid_ne='0'+string(expnum,format='(I1)') $
            else expid_ne=string(expnum,format='(I2)')
          order_ne=string(ordernum,format='(I1)')
          negfile=locroot+target+p+ptg+order_ne+p+expid_ne+dce_neg+p+negtype+fq
          difftype=iroot+'d' ; was "iroot+'x'"
        end

	4 : begin ; sky diff

;         open infile again to search for "sky" target

	  if (foundflag eq 0) then begin 
            if (keyword_set(skyin) eq 0) then begin
              if (keyword_set(clustersky) ne 0) then begin
                split_tgt=strsplit(target,'_',/extract)
                targetsky=split_tgt[0]+'sky'   ; skyin not set, clustersky set
              endif else targetsky=target+'sky' ; skyin and clustersky not set
            endif else targetsky=skyin          ; skyin set

            openr,fsky,infile,/get_lun
            while (foundflag eq 0 and not eof(fsky)) do begin
              readf,fsky,testline
              spl_test=strsplit(testline,' ',/extract)
              targettest  =spl_test[0]
              if (targettest eq targetsky) then begin
                expidtest   =spl_test[3]
                delex=fix(expidtest)-expstart
                foundflag=1
              endif
            endwhile
	    free_lun,fsky
          endif 
	  if (foundflag eq 1) then expsky=expnum+delex else expsky=0
          if (expsky lt 10) then expid_ne='0'+string(expsky,format='(I1)') $
            else expid_ne=string(expsky,format='(I2)')
          order_ne=string(ordernum,format='(I1)')
          if (foundflag eq 1) then negfile=locroot+targetsky+p+ptg+$
	    order_ne+p+expid_ne+dce_neg+p+negtype+fq $
	    else negfile='nosubtraction'+q
          difftype=iroot+'d'
	end
      endcase

      diffile=locroot+target+p+ptg+order+p+expid_sh+dce_sh+p+difftype+fq
      printf,fo,command+cq+posfile+cq+negfile+cq+diffile
      old_expid=expid_sh

    end

    'clean' : begin 

;     set input and output file suffixes based on keywords

      if (keyword_set(diff) eq 0) then begin
        in='im' & bin='bm' 
      endif else begin
        in='id' & bin='bd'
      endelse
      if (keyword_set(cleanout) eq 0) then begin
        out='ic' & bout='bc'
      endif else begin
        out='i'+cleanout & bout='b'+cleanout
      endelse

;     set file names

      imfile = locroot+target+p+ptg+order+p+expid_sh+dce_sh+p+ in +fq
      outfile= locroot+target+p+ptg+order+p+expid_sh+dce_sh+p+ out+fq
      inmask = locroot+target+p+ptg+order+p+expid_sh+dce_sh+p+ bin+fq
      outmask= locroot+target+p+ptg+order+p+expid_sh+dce_sh+p+bout+fq

;     set strings of output line

      string_1='imfclean,'+q+imfile+cq+outfile+c+module+c
      string_2='bmask='+q+inmask+c+'outmask='+q+outmask+c
      string_3=''                  ; modified if cmask or rogue keywords set
      string_4='maskval=4096,/nan'
      string_5=''                  ; modified in some cases if rogue set

;     check cmask keyword and update string_3 if needed to set a campaign mask 

      if (keyword_set(cmask) ne 0) then begin
        if (strlen(cmask) ge 1 and strlen(cmask) le 2) then $
          cstr='_'+string(cmask,format='(i2.2)') else cstr=''
        cmaskfile='superrogue_'+string(module,format='(i1)')+cstr+fq
        string_3='cmask='+q+cmaskfile+c
      endif

;     if cmask not set and rogue is set, then update string_3 to call imrogue
;     if both are set, the block above is executed and rogue is ignored
;     there is a bit of logic required here to make the imrogue call correctly

      if (keyword_set(cmask) eq 0 and keyword_set(rogue) eq 1) then begin
        if (expnum/2 eq expnum/2.0) then begin
          expnod=expnum+1 
        endif else begin
          expnod=expnum-1
          string_5=',/flip'
        endelse
       if (expnod lt 10) then begin
         expnod_sh='0'+string(expnod,format='(I1)')
       endif else begin
         if (expnod lt 100) then expnod_sh=string(expnod,format='(I2)') $
         else expnod_sh=string(expnod,format='(I3)')
       endelse
       rogfile=locroot+target+p+ptg+order+p+expnod_sh+dce_sh+p+'im'+fq
       string_3='rogue='+q+rogfile+c
      endif
     
;     generate output

      printf,fo,string_1+string_2+string_3+string_4+string_5

    end

    'extract' : begin 

;     if count=0 and offline not set, print header for CUPID script

      if (count eq 0 and keyword_set(offline) eq 0) then begin
        printf,fo,'#!/bin/csh'
        printf,fo,'source ~/irs/cupidv2.csh'
      endif

;     set input files for image, bmask, fu for calls to extract and profile

      root=target+p+ptg+order+p+expid_sh+dce_sh+p
      if (keyword_set(exkey) eq 0) then exkey='c'
      if (keyword_set(prkey) eq 0) then prkey='c'
      imfile_pr = root + 'i' + prkey + '.fits'
      bmfile_pr = root + 'b' + prkey + '.fits'
      imfile_ex = root + 'i' + exkey + '.fits'
      bmfile_ex = root + 'b' + exkey + '.fits'
      fufile=root+'fu.fits'

;     determine nod position and thus values of glo and ghi

      if (hiresflag eq 0) then begin  ; lores only
        if ((expnum/2)-(expnum/2.0) eq 0) then begin
          if (module eq 0) then begin ; SL nod 1
            glo=20 & ghi=50
          endif else begin            ; LL nod 1
            glo=50 & ghi=80
          endelse
        endif else begin 
          if (module eq 0) then begin ; SL nod 2
            glo=50 & ghi=80
          endif else begin            ; LL nod 2
            glo=20 & ghi=50
          endelse
        endelse
      endif

;     write extraction script for one DCE

      printf,fo,'echo ',root
      printf,fo,'/bin/cp ',imfile_pr,s,'im_pr.fits'
      printf,fo,'/bin/cp ',bmfile_pr,s,'bm_pr.fits'
      printf,fo,'/bin/cp ',imfile_ex,s,'im_ex.fits'
      printf,fo,'/bin/cp ',bmfile_ex,s,'bm_ex.fits'
      printf,fo,'/bin/cp ',fufile,s,'fu.fits'
      printf,fo,'profile -n profile.nl'
      if (hiresflag eq 0 and gflag eq 1) then $
        printf,fo,'ridge -glo',glo,' -ghi',ghi,' -n ridge.nl' $
      else printf,fo,'ridge -n ridge.nl'
      printf,fo,'echo ',root
      printf,fo,'extract -n extract.nl'
      printf,fo,'/bin/cp extract.tbl ','../'+root+'sp.tbl'
      printf,fo,'/bin/cp profile.tbl ','../'+root+'pro.tbl'

;     increment old EXPID

      old_expid=expid_sh

    end

    'convert' : begin
      tsuffix='sp.tbl' & fsuffix='sp.fits'
      tfile=target+p+ptg+order+p+expid_sh+dce_sh+p+tsuffix+q
      ffile=target+p+ptg+order+p+expid_sh+dce_sh+p+fsuffix+q
      printf,fo,'ipac2fits'+cq+tfile+cq+ffile
    end

    'offset' : begin ; this block of code is outdated
      fsuffix='im.fits'
      sfile=target+p+ptg+order+p+expid_sh+dce_sh
      ffile=locroot+sfile+p+fsuffix
      printf,fo,'print',cq+sfile+q+c+'find_offset('+q+ffile+q+')'
    end

    'imfsky' : begin ; this block of code is outdated
      if (keyword_set(skyin) eq 0) then insuffix='in.fits' $
        else insuffix='i'+skyin+'.fits'
      outsuffix='is.fits'
      infile =locroot+target+p+ptg+order+p+expid_sh+dce_sh+p+insuffix+q
      outfile=locroot+target+p+ptg+order+p+expid_sh+dce_sh+p+outsuffix+q
      printf,fo,'imfsky,'+q+infile+cq+outfile+c+module+c+order
    end

    'exfirs' : begin  ; this block of code is outdated
      if (keyword_set(full) eq 0) then extra='' else extra=',/full'
      prefix1=locroot
      prefix2='./'
      if (hiresflag eq 0) then suffix1='id.fits' else suffix1='im.fits'
      suffix2='ex.fits'
      imfile=prefix1+target+p+ptg+order+p+expid_sh+dce_sh+p+suffix1+q
      spfile=prefix2+target+p+ptg+order+p+expid_sh+dce_sh+p+suffix2+q
      printf,fo,'exfirs,'+q+imfile+cq+spfile+c+module+c+'0,gridfile='+q+$
        gridfile+q+extra
    end

    'imspex' : begin

;     keywords to check for:  bm=bm suffix, fu=fu suffix

      imsuffix='ic.fits' ; this needs to be adjustable, but not yet
      spsuffix='sp.fits' ; fixed

;     create order string ("1" for order 1, "[2,3]" for order 2)

      case order of
        '1' :  ordstr = ',1,'
        '2' :  ordstr = ',[2,3],'
        else : ordstr = ',[11,12,13,14,15,16,17,18,19,20],'
      endcase

;     build posstr to specify source positions if necessary
;     for now, positions are intercepts of spectral trace in row 0 (in pixels)

      if (n_elements(spl_line) le 8) then posstr='' else begin
        posstr=',position='
        npos=n_elements(spl_line)-8
        if (npos eq 1) then posstr=posstr+spl_line[8] else begin
          posstr=posstr+'['
          for i=1,npos do begin
            pos=spl_line[7+i]
            posstr=posstr+pos
            if (i ne npos) then posstr=posstr+',' else posstr=posstr+']'
          endfor
        endelse
      endelse

;     if coadd keyword set, then set order='' and dce_sh = 'cc' 

      if (keyword_set(coadd) ne 0) then begin 
        order='' & dce_sh='cc'
      endif

;     load output strings for files - NOTE HARDCODING for psfile and wavfile

      imfile=target+p+ptg+order+p+expid_sh+dce_sh+p+imsuffix+q
      spfile='../'+target+p+ptg+order+p+expid_sh+dce_sh+p+spsuffix+q
      psffile='cal/psf_'+modstr+'_18v1.fits'+q      ; needs to be adjustable
      wavfile='cal/b'+modstr+'_wavsamp_18.tbl'+q    ; needs to be adjustable

;     do not print output line if in coadd mode and expid hasn't changed

      if (keyword_set(coadd) eq 0 or expid ne expid_old) then $
        printf,fo,'imfspex,'+q+imfile+cq+spfile+ordstr+q+psffile+cq+gridfile+$
        q+c+'wavsamp='+q+wavfile+posstr

    end

  endcase

  expid_old=expid
  target_old=target
  count=count+1

endwhile

free_lun,fi
free_lun,fo
END
