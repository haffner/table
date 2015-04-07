      program tab
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
c    dlam        real(201)             distance away from band-centered 
c                                      wavelength (in nm) 
c    weight      real(201)             weighting function value at dlam
c
c    solar flux
c    wlflux      real(4228)            wavelength of solar flux value
c    flux        real(4228)            solar flux value at flux
c
c    input from radiative transfer code
c    i0_in       real(1521,10,6)       i0 at 1521 wavelengths, 10 sza, 
c                                      and 6 satza
c    z1_in       real(1521,10,6)       z1 at same node points
c    z2_in       real(1521,10,6)       z2 at same node points
c    t_in        real(1521,10,6)       t  at same node points
c    sb_in       real(1521)            sb at same node points
c
c    ouput to table
c    logi0       real(6,10,26,5,2)     log of i0 (after convolution)
c                                      at 6 satza, 10 sza, 26 oz profiles,
c                                      6 wavelengths, and 4 pressures
c    z1i0        real(6,10,26,5,2)     z1 / i0 at same node points 
c    z2i0        real(6,10,26,5,2)     z2 / i0 at same node points
c    tti0        real(6,10,26,5,2)     z2 / i0 at same node points
c    sb          real(26,5,2)          sb at same node points
c    logi0r      real(6,10,26,2)       log of i0 (after convolution)
c                                      at 6 satza, 10 sza, 26 oz profiles,
c                                      6 wavelengths, and 4 pressures
c    z1i0r       real(6,10,26,2)       z1 / i0 at same node points 
c    z2i0r       real(6,10,26,2)       z2 / i0 at same node points
c    tti0r       real(6,10,26,2)       z2 / i0 at same node points
c    sbr         real(26,2)            sb at same node points
c
c
c    character strings used to open files
c    prf         char(26)
c    direc       char(4)   
c    atm         char(4)   
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      real*8 i0_in(1521,10,6), z1_in(1521,10,6), z2_in(1521,10,6),
     1 sb_in(1521), t_in(1521,10,6)
      real wlflux(4228),flux(4228)
      character*4 prf(26)
      character*5 atm(2)
      real logi0(6,10,26,5,2), z1i0(6,10,26,5,2)
      real z2i0(6,10,26,5,2),tti0(6,10,26,5,2),sb(26,5,2)
      real logi0r(6,10,2), z1i0r(6,10,2)
      real z2i0r(6,10,2),tti0r(6,10,2),sbr(2)
      real i0
      character*37  direc(2)
      integer lun
      real dlam, weight
      real cntrwl
c
      common/array/wl(1521),fitflux(1521),cntrwl(6),output(6,5),
     1 nmb(6),stiwl(1521,5),dlam(201),weight(201)
c
c     -- Master table path
      data direc/'/home/cseftor/dev/tomrad/v2.24/atm10/',
     &           '/home/cseftor/dev/tomrad/v2.24/atm04/'/
      data atm/'atm10','atm04'/

c     -- names of ozone profiles
      data prf/'L225','L275','L325','L375','L425','L475','M125','M175',
     1 'M225','M275','M325','M375','M425','M475','M525','M575',
     2 'H125','H175','H225','H275','H325','H375','H425','H475',
     3 'H525','H575'/
      data cntrwl/3087.0,3126.1,3176.2,3224.2,3313.4,3601.5/

c     -- define pi
      data pi/3.1415926/
c
c     ------------------------
c
c     -- read in spectral functions (wavelengths and weights)
      open (40,file='new_spec.dat',
     1 status='old')
      read (40,*) dlam
      read (40,*)
      read (40,*) weight
      close (40)
c     -- read in solar flux data
      open (30,file='data/solar_spectrum.dat',status='old')
      open (32,file='wavl.dat',status='old')
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
c     -- read wl from 'wavl.dat'
c     -- this is the same file (and, therefore, these are the same 
c     -- wavelengths) used to generate the master table
      read (32,*) wl
      close (32)
c     -- interpolate solar flux data to instrument wavelengths
      do i=1,1521
        jj = 0
        do k=1,4228
          if (jj.ne.1 .and. wlflux(k).ge.wl(i)) then
            jj = 1
            dl = wlflux(k) - wlflux(k-1)
            fitflux(i) = flux(k-1)+(flux(k)-flux(k-1)) *
     1       (wl(i)-wlflux(k-1))/dl
          endif
        enddo
      enddo
c     -- figure out where each instrument wavelengths is
c     -- in master table
      do k=1,6
        kk=0
        if (cntrwl(k).ge.3000.) then
          do j=1,1521
            if (kk.ne.1 .and. wl(j).ge.cntrwl(k)) then
              kk=1
              nmb(k)=j
              goto 111
            endif
          enddo
111       continue
        endif
        write(6,*)'table. nmb=',nmb(k)
      enddo
c       -- read from the master tables
      do ipres=1,2
        do iozprf=1,26
c           -- open file corresponding to correct ozone profile
          open (33,file=direc(ipres)//prf(iozprf)//'.'
     1       //atm(ipres)//'.dat',status='old',form='unformatted')
          print *,direc(ipres)//prf(iozprf)//'.'
     1       //atm(ipres)//'.dat'
c           -- read table values
          call rdtabl(33,i0_in,z1_in,z2_in,t_in,sb_in)
          close(33)
          do isza=1,10
            do iscn=1,6
c             -- set up array values to pass to fsavg
              do iwav=1,1521
                stiwl(iwav,1) = i0_in(iwav,isza,iscn)
                stiwl(iwav,2) = z1_in(iwav,isza,iscn)
                stiwl(iwav,3) = z2_in(iwav,isza,iscn)
                stiwl(iwav,4) = t_in(iwav,isza,iscn)
                stiwl(iwav,5) = sb_in(iwav)
              enddo
c             -- do the actual convolution
              call fsavg
c             -- fill storage arrays with actual variables to
c             -- store in instrument array
              do icntrwl=1,5
                if (nmb(icntrwl).ne.0) then
                  logi0(iscn,isza,iozprf,icntrwl,ipres) =
     1                 alog10(output(icntrwl,1)/pi)
                  i0 = output(icntrwl,1)
                  z1i0(iscn,isza,iozprf,icntrwl,ipres) =
     1                 output(icntrwl,2) / i0
                  z2i0(iscn,isza,iozprf,icntrwl,ipres) =
     1                 output(icntrwl,3) / i0 
                  tti0(iscn,isza,iozprf,icntrwl,ipres) =
     1                 output(icntrwl,4) / i0
                  if (isza.eq.1.and.iscn.eq.1) 
     1                 sb(iozprf,icntrwl,ipres)=output(icntrwl,5)
                endif
              enddo
              if (iozprf.eq.1 .and. nmb(6).ne.0) then
                logi0r(iscn,isza,ipres) =
     1               alog10(output(6,1)/pi)
                i0 = output(6,1)
                z1i0r(iscn,isza,ipres) =
     1               output(6,2) / i0
                z2i0r(iscn,isza,ipres) =
     1               output(6,3) / i0 
                tti0r(iscn,isza,ipres) =
     1               output(6,4) / i0
                if (isza.eq.1.and.iscn.eq.1) 
     1               sbr(ipres)=output(6,5)
              endif
            enddo
          enddo
        enddo
      enddo
c       -- output arrays logi0, z1, z2, t, sb
      write (31) logi0,z1i0,z2i0,tti0,sb
      write (31) logi0r,z1i0r,z2i0r,tti0r,sbr
      close (31)
      stop
      end
      subroutine fsavg
c     perform convolution over the slit function
      common/array/wl(1521),fitflux(1521),cntrwl(6),output(6,5),
     1 nmb(6),stiwl(1521,5),dlam(201),weight(201)
      real dlam, weight
      real sum(5), sumw(5)
c     
c     -------------
c
c     -- loop over all 6 instrument wavelengths
      do l=1,6
      if (cntrwl(l).ne.0.) then
c         -- start somewhere far enough away from instrument wavelength
c         -- to include entire slit function, but close enough that we
c         -- don't waste too much time
        k1 = nmb(l)-100
        k2 = nmb(l)+100
        do m = 1,5
          sum(m) = 0.
          sumw(m) = 0.
        enddo
c         -- do the convolution
        do j = k1,k2
          if (fitflux(j).ne.0.0) then
            d = (wl(j)-cntrwl(l)) / 10.
            if (d.ge.dlam(1) .and. d.le.dlam(201)) then
            call locate(dlam,201,d,i)
            frac = (d-dlam(i)) /
     1                 (dlam(i+1)-dlam(i))
            wt   = weight(i)*(1.-frac) +
     1                 weight(i+1)*frac
            do m = 1,4
              sum(m) = sum(m)+wt*fitflux(j)*stiwl(j,m)
              sumw(m) = sumw(m)+wt*fitflux(j)
            enddo
c           -- sb should not be weighted with solar flux
            sum(5)  = sum(5)  + wt*stiwl(j,5)
            sumw(5) = sumw(5) + wt
            endif
          endif
        enddo
        do m = 1,5
          if (sumw(m).ne.0.00000) output(l,m)=sum(m)/sumw(m)
        enddo
      endif
      enddo
      end
      subroutine rdtabl(lun,i0_in,z1_in,z2_in,t_in,sb_in)
c     -- subroutine to read master table
      real*8 i0_in(1521,10,6), z1_in(1521,10,6), z2_in(1521,10,6), 
     1 sb_in(1521), t_in(1521,10,6)
      real*8 wavl(1521), satza(6), sza(10)
      integer nwavl, nscan, nsza
      read (lun) nwavl, nscan, nsza
      read (lun) satza
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
