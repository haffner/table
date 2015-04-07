  nwav = 1521
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
  sza  = dblarr(10)
  satza= dblarr(6)
  openr,20,'/home/cseftor/dev/tomrad/OMTO3LUT/tomrad/v2.24/atm10/m425_atm10.dat',/f77_unformatted
  readu,20,dum
  readu,20,satza
  readu,20,sza
  readu,20,wavl
  readu,20,i0
  readu,20,z1
  readu,20,z2
  readu,20,t
  readu,20,sb
  close,20
  print,wavl(0),wavl(1122)
 
; calculate q1 and q2 and pad out sb
  sza_node    = [0.,30.,45.,60.,70.,77.,81.,84.,86.,88]
  satza_node  = [0.,15.,30.,45.,60.,70.]
  for i = 0,9 do begin
    for j = 0,5 do begin
      q1(*,i,j)  = (-3./8.) * cos(sza_node(i)*!dtor) * sin(sza_node(i)*!dtor) * sin(satza_node(j)*!dtor)
      q2(*,i,j)  = (3./32.) * (sin(sza_node(i)*!dtor))^2 * (sin(satza_node(j)*!dtor))^2 / cos(satza_node(j)*!dtor)
      sbn(*,i,j) = sb(*)
    endfor
  endfor
 
  wavl = wavl/10.
; get rid of pesky pi flux normalization
  i0 = i0/!pi
  z1 = z1/!pi
  z2 = z2/!pi
  t  = t/!pi
; convert z1 and z2 into i1 and i2
  i1 = z1*q1
  i2 = z2*q2

; grab the values for a solar zenith angle of 45 degrees
  i0 = i0(*,2,0)
  i1 = i1(*,2,0)
  i2 = i2(*,2,0)
  t  = t(*,2,0)

  plot,wavl,t
  
  save,filename='m425_bass.sav',wavl,i0,i1,i1,i2,t,sb

  end

