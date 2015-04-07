  fid = H5F_OPEN('O3_CRS_BDM.h5')

  ind = H5D_OPEN(fid,'coefficients/wavelength')
  wvl = H5D_READ(ind)

  ind = H5D_OPEN(fid,'coefficients/c0')
  c0  = H5D_READ(ind)

  ind = H5D_OPEN(fid,'coefficients/c1')
  c1  = H5D_READ(ind)

  ind = H5D_OPEN(fid,'coefficients/c2')
  c2  = H5D_READ(ind)

  openw,10,'bdm_coef.dat'
  for i = 0,4999 do begin
  printf,10,wvl/10.,c0,c1,c2,format='(4f15.8)'
  endfor
  close,10

  end
