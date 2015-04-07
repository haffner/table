      program tg
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c  date
c    october 2001
c
c  author
c    colin seftor 
c
c  purpose
c    make sensor table from the master table:
c     1) preforms convolution over spectral wavelength functions
c
c  variables
c
c    name          type                description
c    ----          ----                -----------
c
c    spectral weighting function
c    dlam        real(101)             distance away from band-centered 
c                                      wavelength (in nm) 
c    weight      real(101)             weighting function value at dlam
c
c    solar flux
c    wlflux      real(2047)            wavelength of solar flux value
c    flux        real(2047)            solar flux value at flux
c
c    input from radiative transfer code
c    i0_in       real(1521,10,11)      i0 at 1521 wavelengths, 10 sza, 
c                                      and 10 satza
c    z1_in       real(1521,10,11)      z1 at same node points
c    z2_in       real(1521,10,11)      z2 at same node points
c    t_in        real(1521,10,11)      t  at same node points
c    sb_in       real(1521)            sb at same node points
c
c    ouput to table
c    logi0       real(11,10,26,10,4,5)log of i0 (after convolution)
c                                     at 10 satza, 10 sza, 26 oz profiles,
c                                     22 wavelengths, and 4 pressures
c    z1i0        real(11,10,26,10,4)  z1 / i0 at same node points 
c    z2i0        real(11,10,26,10,4)  z2 / i0 at same node points
c    tti0        real(11,10,26,10,4)  z2 / i0 at same node points
c    sb          real(26,4,4)         sb at same node points
c
c
c    character strings used to open files
c    prf         char(26)
c    atm         char(4)   
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      real*8 i0_in(1521,10,11), z1_in(1521,10,11), z2_in(1521,10,11),
     1 sb_in(1521), t_in(1521,10,11)
      real wlflux(2047),flux(2047)
      character*4 prf(26)
      character*5 atm(4)
      real logi0(11,10,26,4,4), z1i0(11,10,26,4,4)
      real z2i0(11,10,26,4,4),tti0(11,10,26,4,4),sb(26,4,4)
      real i0
      integer lun
      real dlam, weight
c
      common/array/wl(1521),fitflux(1521),cntrwl(4),output(6,5),
     1 nmb(4),stiwl(1521,5),dlam(61),weight(61,4)
c

c     -- names of ozone profiles
      data prf/'l225','l275','l325','l375','l425','l475',
     1 'm125','m175',
     1 'm225','m275','m325','m375','m425','m475','m525','m575',
     2 'h125','h175','h225','h275','h325','h375','h425','h475',
     3 'h525','h575'/
c     -- define pi
      data pi/3.1415926/
      data cntrwl/3175.0,3250.0,3400.0,3880.0/
c
c     ------------------------
c
c     -- read in spectral functions (wavelengths and weights)
      open (40,file='dscovr_bps.dat',status='old')
      read (40,*) dlam
      read (40,*) weight
      close (40)
      print *,dlam
c     read (5,*)
c     -- read in solar flux data
      open (30,
     1 file='/home/cseftor/dev/tomrad/OMTO3LUT/SolarFlux/SOLSTICE',
     2 status='old')
      open (32,file='dcoe.dat',status='old')
c     -- read solar flux data
      i = 0
3     continue
        read (30,*,end=400) aa,bb
        i = i+1
        wlflux(i) = aa*10.
        flux(i) = bb
        goto 3
400   continue
      close (30)
c     -- read wl from 'ABS_COEF.DAT'
c     -- this is the same file (and, therefore, these are the same 
c     -- wavelengths) used to generate the master table
      j = 0
      read (32,*)
1     continue
      read(32,*,end=800)a,b,x,d,e
      j = j+1
      wl(j) = a 
      goto 1
800   continue
      print *,wl(1),wl(1521)
      close (32)
c     -- interpolate solar flux data to instrument wavelengths
      do i=1,1521
        do k=1,2047
          if (wlflux(k).ge.wl(i)) then
            dl = wlflux(k) - wlflux(k-1)
            fitflux(i) = flux(k-1)+(flux(k)-flux(k-1)) *
     1       (wl(i)-wlflux(k-1))/dl
            goto 801
          endif
        enddo
801     continue
      enddo
c     -- get bandcenters for the different parts of the CCD array
c     -- (the first set correspond to the edge, the last set 
c     -- correspond to the middle)
      do jj = 1,1
c       -- open instrument table and set up wavelengths for this
c       -- part of the CCD array
	open (31,file='data/TABLE.D',status='unknown',
     1   form='unformatted')
c       -- figure out where each instrument wavelengths is
c       -- in master table
	do k=1,4
	  if (cntrwl(k).ge.300.0) then
	    do j=1,1521
	      if (wl(j).ge.cntrwl(k)) then
		nmb(k)=j
                goto 111
	      endif
	    enddo
111         continue
	  endif
	  write(6,*)'table. nmb=',nmb(k)
	enddo
c       stop
c       -- read from the master tables
	do ipres=1,4
	  do iozprf=1,26
c           -- open file corresponding to correct ozone profile
            if (ipres.eq.1) then
              open (33,
     1      file='/home/cseftor/dev/tomrad/OMTO3LUT/tomrad/v2.24/atm10/'
     1        //prf(iozprf)//'.atm10.dat',status='old',
     2        form='unformatted')
            else if (ipres.eq.2) then
              open (33,
     1      file='/home/cseftor/dev/tomrad/OMTO3LUT/tomrad/v2.24/atm07/'
     1        //prf(iozprf)//'.atm07.dat',status='old',
     2        form='unformatted')
            else if (ipres.eq.3) then
              open (33,
     1      file='/home/cseftor/dev/tomrad/OMTO3LUT/tomrad/v2.24/atm04/'
     1        //prf(iozprf)//'.atm04.dat',status='old',
     2        form='unformatted')
            else if (ipres.eq.4) then
              open (33,
     1      file='/home/cseftor/dev/tomrad/OMTO3LUT/tomrad/v2.24/atm02/'
     1        //prf(iozprf)//'.atm02.dat',status='old',
     2        form='unformatted')
            endif
	    print *,prf(iozprf)
c           -- read table values
	    call rdtabl(33,i0_in,z1_in,z2_in,t_in,sb_in)
	    close(33)
	    do isza=1,10
	      do iscn=1,11
c               -- set up array values to pass to fsavg
		do iwav=1,1521
		  stiwl(iwav,1) = i0_in(iwav,isza,iscn)
		  stiwl(iwav,2) = z1_in(iwav,isza,iscn)
		  stiwl(iwav,3) = z2_in(iwav,isza,iscn)
		  stiwl(iwav,4) = t_in(iwav,isza,iscn)
		  stiwl(iwav,5) = sb_in(iwav)
		enddo
c               -- do the actual convolution
		call fsavg
c               -- fill storage arrays with actual variables to
c               -- store in instrument array
		do icntrwl=1,4
		  if (nmb(icntrwl).ne.0) then
		    logi0(iscn,isza,iozprf,icntrwl,ipres)=
     1               alog10(output(icntrwl,1)/pi)
		    i0 = output(icntrwl,1)
		    z1i0(iscn,isza,iozprf,icntrwl,ipres)=
     1               output(icntrwl,2) / i0
		    z2i0(iscn,isza,iozprf,icntrwl,ipres)=
     1               output(icntrwl,3) / i0 
		    tti0(iscn,isza,iozprf,icntrwl,ipres)=
     1               output(icntrwl,4) / i0
		    if (isza.eq.1.and.iscn.eq.1) 
     1               sb(iozprf,icntrwl,ipres)=output(icntrwl,5)
		  endif
		enddo   
	      enddo
	    enddo
	  enddo
	enddo
c       -- output arrays logi0, z1, z2, t, sb
	write (31) logi0,z1i0,z2i0,tti0,sb
	close (31)
      enddo
      stop
      end
      subroutine fsavg
c     perform convolution over the slit function
      common/array/wl(1521),fitflux(1521),cntrwl(4),output(6,5),
     1 nmb(4),stiwl(1521,5),dlam(61),weight(61,4)
      real dlam, weight
      real sum(5), sumw(5)
c     
c     -------------
c
c     -- loop over all 4 instrument wavelengths
      do mm=1,4
	if (cntrwl(mm).ne.0.) then
c         -- start somewhere far enough away from instrument wavelength
c         -- to include entire slit function, but close enough that we
c         -- don't waste too much time
	  k1 = nmb(mm)-25
	  k2 = nmb(mm)+25
	  do m = 1,5
	    sum(m) = 0.
	    sumw(m) = 0.
	  enddo
c         -- do the convolution
	  do j = k1,k2
	    if (fitflux(j).ne.0.0) then
	      d = (wl(j)-cntrwl(mm)) / 10.
	      if (d.ge.dlam(1) .and. d.le.dlam(61)) then
		call locate(dlam,41,d,i)
		frac = (d-dlam(i)) /
     1                 (dlam(i+1)-dlam(i))
		wt   = weight(i,mm)*(1.-frac) +
     1                 weight(i+1,mm)*frac
		do m = 1,4
		  sum(m)  = sum(m)+wt*fitflux(j)*stiwl(j,m)
		  sumw(m) = sumw(m)+wt*fitflux(j)
		enddo
c               -- sb should not be weighted with solar flux
		sum(5)  = sum(5)  + wt*stiwl(j,5)
		sumw(5) = sumw(5) + wt
	      endif
	    endif
	  enddo
	  do m = 1,5
	    if (sumw(m).ne.0.00000) output(mm,m)=sum(m)/sumw(m)
	  enddo
	endif
      enddo
      end
      subroutine rdtabl(lun,i0_in,z1_in,z2_in,t_in,sb_in)
c     -- subroutine to read master table
      real*8 i0_in(1521,10,11), z1_in(1521,10,11), z2_in(1521,10,11), 
     1 sb_in(1521), t_in(1521,10,11)
      real*8 wavl(1521), scan(11), sza(10)
      integer nwavl, nscan, nsza
      read (lun) nwavl, nscan, nsza
      read (lun) scan
      read (lun) sza
      read (lun) wavl
      read (lun) i0_in
      read (lun) z1_in
      read (lun) z2_in
      read (lun) t_in
      read (lun) sb_in
      return
      end
      subroutine locate(xx,n,x,j)
      integer j,n
      real x,xx(n)
      integer jl,jm,ju
      jl=0
      ju=n+1
10    if(ju-jl.gt.1)then
        jm=(ju+jl)/2
        if((xx(n).ge.xx(1)).eqv.(x.ge.xx(jm)))then
          jl=jm
        else
          ju=jm
        endif
      goto 10
      endif
      if(x.eq.xx(1))then
        j=1
      else if(x.eq.xx(n))then
        j=n-1
      else
        j=jl
      endif
      return
      end
