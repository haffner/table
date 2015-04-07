  line = ''
  coe = fltarr(6,2031)
  openr,10,'/omi/home/cseftor/dev/tomrad/OMTO3LUT/tomrad/v2.24/test/coe.dat'
  readf,10,line
  readf,10,coe
  close,10
  openw,10,'wavl.dat'
  printf,10,transpose(coe(0,*))
  close,10
  end
