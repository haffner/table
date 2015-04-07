; pro simrad
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
;
; author
;  colin seftor
;
; purpose
;  pull nr's out of the master table and feed them to the omps
;  tc algorithm for simulations
;
; purpose
;
;    changes
;      23-december-2003 NGST OMPS Drop 2.1 TC EDR Software Package
;                       SPCR ALG00000270 James P. Done
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
; set OMPS instrument wavelengths
;
  pointer = [20,25,29,30,33,36,38,41,43,48,50,52,60,67,69,74,76,86,152,160,171,183]
; pix     = [9,179,349,519,689]
  pixel   = fltarr(699)
  tmp  = fltarr(193,699)
  openr,20,'tables/ntcbandcntrs4-12-01.txt'
  readf,20,tmp
  close,20
  pixel = tmp(0,*)
  twv = tmp(1:192,where(pixel eq 0))
  bcw = twv(pointer)
;
; set up parameters
  prf = ''
  atm = ''
  read,'What year? ',yr
  read,'What day of year? ',doy
  xlon = 0.
  read,'What profile? ',prf
  if strmid(prf,0,1) eq 'l' then xlat = 0.
  if strmid(prf,0,1) eq 'm' then xlat = 45.
  if strmid(prf,0,1) eq 'h' then xlat = 80.
  read,'What pressure (1.0, 0.7, 0.4, 0.2 atm)? ',prs
  atm = strmid(strcompress(string(fix(prs*10.)),/remove_all),0,2)
  if prs lt 1.0 then atm = '0'+atm
  read,'What solar zenith angle? ',sza
  read,'What satellite zenith angle? ',satza
  read,'What azimuth angle? ',phi
  read,'What reflectivity (0.0-1.0)? ',refl
;
; get part of OMPS master table
;
  sza_node    = [0.,30.,45.,60.,70.,77.,81.,84.,86.,88]
  satza_node  = [0.,15.,30.,45.,60.,70.]
  lsec_sza   = alog10(1./cos(!dtor*sza_node))
  lsec_satza = alog10(1./cos(!dtor*satza_node))
;
; get trop ozone from data base
;
  edoy    = [31,59,90,120,151,181,212,243,273,304,334,365]
  tropoz  = fltarr(4,360,180,12)
  tozclim = fltarr(4)
  openr,20,'tables/tropoz.database'
  readu,20,tropoz
  close,20
  ilon = fix(181. + xlon)
  if ilon gt 360 then ilon = 360
  ilat = fix(91.  + xlat)
  if ilat gt 180 then ilat = 180
  if yr mod 4 eq 0 then edoy(1:11) = edoy(1:11) + 1
  pointer = where(doy lt edoy)
  ptr = pointer(0) - 1
  if ptr ge 0 then begin
    mfrac = (doy - edoy(ptr)) /  (edoy(ptr+1) - edoy(ptr))
    tozclim(0) = (1.-mfrac) * tropoz(0,ilon,ilat,ptr) + $
                   mfrac    * tropoz(0,ilon,ilat,ptr+1)
    tozclim(1) = (1.-mfrac) * tropoz(1,ilon,ilat,ptr) + $
                   mfrac    * tropoz(1,ilon,ilat,ptr+1)
    tozclim(2) = (1.-mfrac) * tropoz(2,ilon,ilat,ptr) + $
                   mfrac    * tropoz(2,ilon,ilat,ptr+1)
    tozclim(3) = (1.-mfrac) * tropoz(3,ilon,ilat,ptr) + $
                   mfrac    * tropoz(3,ilon,ilat,ptr+1)
  endif else begin
    mfrac = (doy - edoy(11)) /  (edoy(0) - edoy(11))
    tozclim(0) = (1.-mfrac) * tropoz(0,ilon,ilat,11) + $
                   mfrac    * tropoz(0,ilon,ilat,0)
    tozclim(1) = (1.-mfrac) * tropoz(1,ilon,ilat,11) + $
                   mfrac    * tropoz(1,ilon,ilat,0)
    tozclim(2) = (1.-mfrac) * tropoz(2,ilon,ilat,11) + $
                   mfrac    * tropoz(2,ilon,ilat,0)
    tozclim(3) = (1.-mfrac) * tropoz(3,ilon,ilat,11) + $
                   mfrac    * tropoz(3,ilon,ilat,0)
  endelse
;
  nwav = 3047
  dum  = lonarr(3)
  wavl = dblarr(nwav)
  i0   = dblarr(nwav,10,6)
  z1   = dblarr(nwav,10,6)
  z2   = dblarr(nwav,10,6)
  i1   = dblarr(nwav,10,6)
  i2   = dblarr(nwav,10,6)
  t    = dblarr(nwav,10,6)
  sb   = dblarr(nwav)
  sbn  = dblarr(nwav,10,6)
  q1   = dblarr(nwav,10,6)
  q2   = dblarr(nwav,10,6)
  rad  = dblarr(nwav,10,6)
  openr,20,'MASTER/atm'+atm+'/'+prf+'_atm'+atm+'.dat',/f77_unformatted
  readu,20,dum
  readu,20,wavl
  readu,20,i0
  readu,20,z1
  readu,20,z2
  readu,20,t
  readu,20,sb
  close,20
;
; calculate q1 and q2 and pad out sb
  for i = 0,9 do begin
    for j = 0,5 do begin
      q1(*,i,j)  = (-3./8.) * cos(sza_node(i)*!dtor) * sin(sza_node(i)*!dtor) * sin(satza_node(j)*!dtor)
      q2(*,i,j)  = (3./32.) * (sin(sza_node(i)*!dtor))^2 * (sin(satza_node(j)*!dtor))^2 / cos(satza_node(j)*!dtor)
      sbn(*,i,j) = sb(*)
    endfor
  endfor
;
  wavl = wavl/10.
; get rid of pesky pi flux normalization
  i0 = i0/!pi
  z1 = z1/!pi
  z2 = z2/!pi
  t  = t/!pi
; convert z1 and z2 into i1 and i2
  i1 = z1*q1
  i2 = z2*q2
; set phi
  phi = 90.
;
; calculate normalized radiances
  rad = i0 + i1*cos(!dtor*phi) + i2*cos(2.*!dtor*phi) + refl*t/(1.-refl*sbn)
;
; get solar flux
  swvl = fltarr(4228)
  srad = fltarr(4228)
  tmp  = fltarr(2)
  openr,20,'tables/solar_spectrum_calib.dat'
  for i = 0,4227 do begin
    readf,20,tmp
    swvl(i) = tmp(0)
    srad(i) = tmp(1)
  endfor
  close,20
;
; get spectral functions
;
  tlam = fltarr(994)
  twt  = fltarr(994)
  dlam = fltarr(22,994)
  wt   = fltarr(22,994)
  openr,40,'tables/omps_spec_function.dat'
  header = ''
  readf,40,header
  readf,40,header
  readf,40,header
  readf,40,header
  readf,40,header
  readf,40,header
  for iw = 0,21 do begin
    readf,40,header
    readf,40,header
    readf,40,tlam
    dlam(iw,*) = tlam
    readf,40,header
    readf,40,twt
    wt(iw,*) = twt
  endfor
  close,40
; dlam = 1.12*dlam
; 
; interpolate solar flux to OMPS master table wavelengths
  sflux = interpol(srad,swvl,wavl)
;
; now perform convolution of OMPS spectral functions
  nr        = fltarr(22,10,6)
  nr_satza  = fltarr(6)
  nr_sza    = fltarr(10)
  nval      = fltarr(22)
  sfl_cnv   = fltarr(22)
  for i = 0,21 do begin
    d=wavl-bcw(i)
    wght=interpol(wt(i,*),dlam(i,*),d)
    wght(where(wght eq 0)) = 0.
;   convolve solar flux
    sfl_cnv(i) = total(sflux*wght) / total(wght)
;   print,sfl_cnv(i)
;   print,sfl_cnv(i)
;   calculate flux weighted integral over omps slit function
    for isza = 0,9 do begin
      for isatza = 0,5 do begin
        rr = rad(*,isza,isatza)
        nr_satza(isatza) = total(rr*sflux*wght)/total(sflux*wght)
      endfor
;     if i eq 0 then print,nr_satza
      a = lsec_satza
      b = nr_satza
      c = alog10(1./cos(!dtor*satza))
      polint,a,b,c,d
      nr_sza(isza) = d
    endfor
;   if i eq 0 then print,nr_sza
    a = lsec_sza
    b = nr_sza
    c = alog10(1./cos(!dtor*sza))
    polint,a,b,c,d
    nval(i) = d
  endfor
  print,nval
;
; 10% cloud
  nnv = [0.00935516,0.0168803,0.0244937,0.0271614,0.0299210,0.0387965,0.0406650,0.0475966,0.0506375,0.0589570, $
         0.0710417,0.0679831,0.0735384,0.0790805,0.0891372,0.0839624,0.0897140,0.0915818,0.0804197,0.0789980, $
         0.0771316,0.0751480]
; 10% ground
; nnv = [0.0140698,0.0268035,0.0399594,0.0445938,0.0496428,0.0653398,0.0688838,0.0815241,0.0872173,0.103145, $
;        0.125159,0.120278,0.132835,0.145398,0.164869,0.157111,0.168750,0.176453,0.180499,0.180394, $
;        0.180259,0.180116]
  nval(18:21) = nnv(18:21)
  print,nval
;
  nval = -100.*alog10(nval)
  print,nval
  dd = ''
; read,dd
;
  getprf,prf,prfoz,prftemp
  if prs lt 1. then tozclim(0:3) = prfoz(0:3)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
; build "sdr"
;
  openw,42,'omps.sdr',/f77_unformatted
; file id
  id = 'omps.sdr                         '
  id = strmid(id,0,20)
;
; fill solar flux, limb profile,  pres profile
; (we're only doing build 1 so far)
;
  sflux      = fltarr(6)
  presprf    = fltarr(11)
  presprf(*) = -99.
  pcld       = prs
  pgrd       = 1.0
  ozprf      = prfoz
  tempprf    = prftemp
  clfrac     = 1.
  sdev       = 0.
; print,'pgrd,pcld = ',pgrd,pcld
; viirs cloud fraction
  vcf        = 0.
; set snow to no snow, activate with flag later
  snow       = 0l
  tozclim    = ozprf(0:3)
; fill solar azimuth angle (because I don't want to calculate it and it's not needed right now)
  writeu,42,id,long(yr),long(doy),xlat,xlon,sza,.253,satza,nval,sfl_cnv,phi, $
            tozclim,sdev,clfrac,vcf,pcld,pgrd,tempprf,presprf,snow,ozprf,bcw,0.08
  close,42
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
  dum = ''
  residue = fltarr(5)
  sens    = fltarr(5)
  refsens = fltarr(6)
;
;    nim,id,ozone,ref,algflg,residue,sens,refsens,ozcld,clfrac,prfrac,so2ind
;
  omps_nadir,/xprf,/tprf,/ptrop,id,ozone,ref,algflg,residue,sens,refsens,ozcld,clfrac,prfrac,so2ind
; print,ozone,ref,algflg,ozcld,clfrac,prfrac,so2ind
; print,'Ozone = ',ozone
; print,'Reflec = ',ref
  end

