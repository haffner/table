  logi0 = fltarr(6,10,26,22,4,7)
  z1i0  = fltarr(6,10,26,22,4,7)
  z2i0  = fltarr(6,10,26,22,4,7)
  ti0   = fltarr(6,10,26,22,4,7)
  ssb   = fltarr(26,22,4,7)
  sb    = fltarr(6,10,26,22,4,7)

; calculate q1 and q2 and pad out sb
  sza_node    = [0.,30.,45.,60.,70.,77.,81.,84.,86.,88]
  satza_node  = [0.,15.,30.,45.,60.,70.]

  openr,10,'data/TABLE.OMPS',/f77_unformatted
  readu,10, logi0, z1i0, z2i0, ti0, ssb

  q1 = fltarr(6,10,26,22,4,7)
  q2 = fltarr(6,10,26,22,4,7)

  for i = 0,9 do begin
    for j = 0,5 do begin
      for k = 0,25 do begin
        for m = 0,21 do begin
          for n = 0,3 do begin
            for o = 0,6 do begin
              q1(j,i,k,m,n,o)  = (-3./8.) * cos(sza_node(i)*!dtor) * sin(sza_node(i)*!dtor) * sin(satza_node(j)*!dtor)
              q2(j,i,k,m,n,o)  = (3./32.) * (sin(sza_node(i)*!dtor))^2 * (sin(satza_node(j)*!dtor))^2 / cos(satza_node(j)*!dtor)
              sb(j,i,k,m,n,o)  = ssb(k,m,n,o)
            endfor
          endfor
        endfor
      endfor
    endfor
  endfor
  i0 = 10.^logi0
  i1 = z1i0 * i0 * q1
  i2 = z2i0 * i0 * q2
  t  = ti0 * i0

  rad = i0 + i1*cos(!pi/2.) + i2*cos(!pi/2.) + 0.08*t / (1.-0.08*sb)

  nval = -100.*alog10(rad)

  for i = 0,4 do begin

    print,transpose(t(0,0,10,*,0,i))

    print,''

    print,transpose(i0(0,0,10,*,0,i))

    print,''

    print,transpose(rad(0,0,10,*,0,i))
  
    print,''

    print,transpose(nval(0,0,10,*,0,i))

    dd = ''
    read,dd
    print,''
    print,''
    print,''

  endfor

  end
