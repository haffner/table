  logi0n = fltarr(11,10,26,4,4)
  z1i0n  = fltarr(11,10,26,4,4)
  z2i0n  = fltarr(11,10,26,4,4)
  ti0n   = fltarr(11,10,26,4,4)
  sbn    = fltarr(26,4,4)

  openr,10,'TABLE.D',/f77_unformatted
  readu,10,logi0n,z1i0n,z2i0n,ti0n,sbn
  close,10

  nlogi0n = logi0n(*,*,*,0:2,*)
  nz1i0n  = z1i0n(*,*,*,0:2,*)
  nz2i0n  = z2i0n(*,*,*,0:2,*)
  nti0n   = ti0n(*,*,*,0:2,*)
  nsbn    = sbn(*,0:2,*)

; lgi0rn  = logi0n(*,*,*,3,*)
; z1i0rn  = z1i0n(*,*,*,3,*)
; z2i0rn  = z2i0n(*,*,*,3,*)
; ti0rn   = ti0n(*,*,*,3,*)
; sbrn    = sbn(*,3,*)

  lgi0rn  = logi0n(*,*,7,3,*)
  z1i0rn  = z1i0n(*,*,7,3,*)
  z2i0rn  = z2i0n(*,*,7,3,*)
  ti0rn   = ti0n(*,*,7,3,*)
  sbrn    = sbn(7,3,*)

  openw,10,'TABLE.DSCOVR',/f77_unformatted
  writeu,10,nlogi0n,nz1i0n,nz2i0n,nti0n,nsbn
  writeu,10,lgi0rn,z1i0rn,z2i0rn,ti0rn,sbrn
  close,10

  end
