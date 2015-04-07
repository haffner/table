; openr,10,'TABLE.one',/f77_unformatted
  openr,10,'/home/cseftor/dev/TC_EDR_NASA/tbl/TABLE.OMPS',/f77_unformatted
  logi0 = fltarr(137280,7)
  z1 = fltarr(137280,7)
  z2 = fltarr(137280,7)
  t = fltarr(137280,7)
  sb = fltarr(2288,7)
  readu,10,logi0,z1,z2,t,sb
  close,10

  openr,10,'TABLE.two',/f77_unformatted
  nlogi0 = fltarr(137280,7)
  nz1 = fltarr(137280,7)
  nz2 = fltarr(137280,7)
  nt = fltarr(137280,7)
  nsb = fltarr(2288,7)
  readu,10,nlogi0,nz1,nz2,nt,nsb
  close,10

  plot,100.*(logi0(*,0)-nlogi0(*,0))/nlogi0(*,0)
  dd = ''
  read,dd

  plot,100.*(logi0(*,3)-nlogi0(*,3))/nlogi0(*,3)
  dd = ''
  read,dd

  plot,100.*(logi0(*,6)-nlogi0(*,6))/nlogi0(*,6)
  dd = ''
  read,dd

  plot,100.*(z1(*,0)-nz1(*,0))/nz1(*,0)
  dd = ''
  read,dd

  plot,100.*(z1(*,3)-nz1(*,3))/nz1(*,3)
  dd = ''
  read,dd

  plot,100.*(z1(*,6)-nz1(*,6))/nz1(*,6)
  dd = ''
  read,dd

  plot,100.*(z2(*,3)-nz2(*,3))/nz2(*,3)
  dd = ''
  read,dd

  plot,100.*(t(*,3)-nt(*,3))/nt(*,3)
  dd = ''
  read,dd

  plot,100.*(sb(*,3)-nsb(*,3))/nsb(*,3)
  dd = ''
  read,dd

  end
