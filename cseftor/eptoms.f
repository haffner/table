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
c    i0_in       real(3040,10,6)       i0 at 3040 wavelengths, 10 sza, 
c                                      and 6 satza
c    z1_in       real(3040,10,6)       z1 at same node points
c    z2_in       real(3040,10,6)       z2 at same node points
c    t_in        real(3040,10,6)       t  at same node points
c    sb_in       real(3040)            sb at same node points
c
c    ouput to table
c    logi0       real(6,10,21,6,4,5)  log of i0 (after convolution)
c                                      at 6 satza, 10 sza, 21 oz profiles,
c                                      22 wavelengths, and 4 pressures
c    z1i0        real(6,10,21,6,4)    z1 / i0 at same node points 
c    z2i0        real(6,10,21,6,4)    z2 / i0 at same node points
c    tti0        real(6,10,21,6,4)    z2 / i0 at same node points
c    sb          real(21,6,4)         sb at same node points
c
c    position    char(4)               character string to denote
c                                      which set of bandcenters
c                                      were used for the convolution    
c
c    character strings used to open files
c    prf         char(21)
c    atm         char(4)   
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      real*8 i0_in(3040,10,6), z1_in(3040,10,6), z2_in(3040,10,6),
     1 sb_in(3040), t_in(3040,10,6)
      real wlflux(2047),flux(2047)
      character*4 prf(21)
      character*5 atm(4)
      real logi0(6,10,21,6,4), z1i0(6,10,21,6,4)
      real z2i0(6,10,21,6,4),tti0(6,10,21,6,4),sb(21,6,4)
      real i0
      character*1 position(7)
      integer lun
      real dlam, weight
c
      common/array/wl(3040),fitflux(3040),cntrwl(6),output(6,5),
     1 nmb(6),stiwl(3040,5),dlam(41),weight(41)
c
c     -- Master table path
      data position/'0','1','2','3','4','5','6'/

c     -- names of ozone profiles
      data prf/'L225','L275','L325',
     1 'M225','M275','M325','M375','M425','M475','M525','M575',
     2 'H125','H175','H225','H275','H325','H375','H425','H475',
     3 'H525','H575'/
c     -- define pi
      data pi/3.1415926/
c
c     ------------------------
c
c     -- read in spectral functions (wavelengths and weights)
      open (40,file='toms_slit.dat',status='old')
      read (40,*) dlam
      print *,dlam
      read (40,*) weight
      close (40)
c     read (5,*)
c     -- read in solar flux data
      open (30,
     1 file='/home/cseftor/dev/tomrad/OMTO3LUT/SolarFlux/SOLSTICE',
     2 status='old')
      open (32,file='coe_daum_05b.dat',status='old')
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
      print *,wl(1),wl(3040)
      close (32)
c     -- interpolate solar flux data to instrument wavelengths
      do i=1,3040
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
      open (20,file='eptoms_wav.dat',status='old')
      read (20,*) cntrwl
      read (20,*) 
      close (20)
      do jj = 1,1
c       -- open instrument table and set up wavelengths for this
c       -- part of the CCD array
	open (31,file='data/TABLE.EP',status='unknown',
     1   form='unformatted')
c       -- figure out where each instrument wavelengths is
c       -- in master table
	do k=1,6
	  if (cntrwl(k).ge.3000.) then
	    do j=1,3040
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
	  do iozprf=1,21
c           -- open file corresponding to correct ozone profile
            if (ipres.eq.1) then
              open (33,file='/home/dhaffner/tables/master_daum2/rad/'
     1        //prf(iozprf)//'.atm10.rad.dat',status='old',
     2        form='unformatted')
            else if (ipres.eq.2) then
              open (33,file='/home/dhaffner/tables/master_daum2/rad/'
     1        //prf(iozprf)//'.atm07.rad.dat',status='old',
     2        form='unformatted')
            else if (ipres.eq.3) then
              open (33,file='/home/dhaffner/tables/master_daum2/rad/'
     1        //prf(iozprf)//'.atm04.rad.dat',status='old',
     2        form='unformatted')
            else if (ipres.eq.4) then
              open (33,file='/home/dhaffner/tables/master_daum2/rad/'
     1        //prf(iozprf)//'.atm025.rad.dat',status='old',
     2        form='unformatted')
            endif
	    print *,prf(iozprf)
c           -- read table values
	    call rdtabl(33,i0_in,z1_in,z2_in,t_in,sb_in)
	    close(33)
	    do isza=1,10
	      do iscn=1,6
c               -- set up array values to pass to fsavg
		do iwav=1,3040
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
		do icntrwl=1,6
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
      common/array/wl(3040),fitflux(3040),cntrwl(6),output(6,5),
     1 nmb(6),stiwl(3040,5),dlam(41),weight(41)
      real dlam, weight
      real sum(5), sumw(5)
c     
c     -------------
c
c     -- loop over all 6 instrument wavelengths
      do mm=1,6
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
	      if (d.ge.dlam(1) .and. d.le.dlam(41)) then
		call locate(dlam,41,d,i)
		frac = (d-dlam(i)) /
     1                 (dlam(i+1)-dlam(i))
		wt   = weight(i)*(1.-frac) +
     1                 weight(i+1)*frac
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
      real*8 i0_in(3040,10,6), z1_in(3040,10,6), z2_in(3040,10,6), 
     1 sb_in(3040), t_in(3040,10,6)
      real*8 wavl(3040), scan(6), sza(10)
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
