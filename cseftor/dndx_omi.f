      program table
c
c     purpose: read i0,z1,z2,t and s from the master tables, take slit 
c     and solar flux average for given lamda0 & lamda0+dwl, & write outputs.  
c
c input
      real*8 i0_in(1521,10,6), z1_in(1521,10,6), z2_in(1521,10,6),
     1 sb_in(1521), t_in(1521,10,6)
      real wlflux(4228),flux(4228)
      character*4 prf(26)
      character*5 atm(4)
c
c output
      real*4 logi0(6,10,26,22,4,12), z1i0(6,10,26,22,4,12)
      real*4 z2i0(6,10,26,22,4,12),tti0(6,10,26,22,4,12),sb(26,22,4,12)
c
      real i0
      character*30 direc
      character*60 sunfile,outfile,specfile
      character satellite*5
      integer isatellite,lun
      real*4 dlam, weight
c
      common/array/wl(1521),fitflux(1521),cntrwl(22),output(22,5),
     1 nmb(22),stiwl(1521,5),dlam(17),weight(17)
      data cntrwl
     1 /308.418,310.518,312.199,312.619,313.879,315.139,
     2  315.979,317.239,318.079,320.180,321.020,321.860,
     3  325.220,328.160,329.000,331.101,331.941,336.141,
     4  363.859,367.219,371.838,376.877/
c    4  363.859,363.859,363.859,363.859/
c     -- Master table path
      data direc/'/home/cseftor/dev/tomrad/dndx/'/
      data pi/3.1415926/
      data prf/'l225','l275','l325','l375','l425','l475','m125','m175',
     1 'm225','m275','m325','m375','m425','m475','m525','m575',
     2 'h125','h175','h225','h275','h325','h375','h425','h475',
     3 'h525','h575'/
      data atm/'atm10','atm07','atm04','atm01'/
c
      open (31,file='DNDX.OMPS_OMI',status='unknown',form='unformatted')
      do i = 1,22
        cntrwl(i) = 10.*cntrwl(i)
      enddo
c
c     -- read in spectral functions (wavelengths and weights)
      open (40,file='data/omi_slit.new',status='old')
      read (40,*) dlam
      read (40,*)
      read (40,*) weight
      close (40)
c
      open (30,file='data/solar_spectrum.dat',status='old')
      open (32,file='wavl.dat',status='old')
c     open (32,file='data/DNDX_COEF.DAT',status='old')
c     -- parameters needed for calculation of itotal & partial derivatives
c     -- (phi=90 degrees)
      t = -50.
      phi = pi/2.
      refl = 0.3
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
      read (32,*) wl
      close (32)
      do 10 i=1,1521
         jj = 0
         do 15 k=1,4228
            if (jj.ne.1 .and. wlflux(k).ge.wl(i)) then
               jj = 1
               dl = wlflux(k) - wlflux(k-1)
               fitflux(i) = flux(k-1)+(flux(k)-flux(k-1)) *
     1          (wl(i)-wlflux(k-1))/dl
            endif
15       continue
10    continue
      do 25 k=1,22
         kk=0
         if (cntrwl(k).ge.3000.) then
            do 20 j=1,1521
               if (kk.ne.1 .and. wl(j).ge.cntrwl(k)) then
                  kk=1
                  nmb(k)=j
               endif
20          continue
         endif
         write(6,*)'table, nmb=',nmb(k)
25    continue
c     -- read from the master tables
      do 85 ipres=1,4
         lun=33
         do 75 iozprf=1,26
            open (lun,file=direc//atm(ipres)//'/dndx.'//prf(iozprf),
     1       status='old',form='unformatted')
            print *,direc//atm(ipres)//'/dndx.'//prf(iozprf)
c           -- read in unperturbed layer and throw away
            do ilyr = 1,12
            iilyr = ilyr
            if (ilyr.ne.1) iilyr = 14-ilyr
c           -- read table values for pertubed layer ilyr
            call rdtabl(lun,i0_in,z1_in,z2_in,t_in,sb_in)
            do 65 isza=1,10
               do 55 iscn=1,6
c                 -- keep track of record count
                  do 40 iwav=1,1521
                     stiwl(iwav,1) = i0_in(iwav,isza,iscn)
                     stiwl(iwav,2) = z1_in(iwav,isza,iscn)
                     stiwl(iwav,3) = z2_in(iwav,isza,iscn)
                     stiwl(iwav,4) = t_in(iwav,isza,iscn)
                     stiwl(iwav,5) = sb_in(iwav)
40                continue
                  call fsavg
                  do 45 icntrwl=1,22
                     if (nmb(icntrwl).ne.0) then
                        logi0(iscn,isza,iozprf,icntrwl,ipres,iilyr)=
     1                   alog10(output(icntrwl,1)/pi)
                        i0 = output(icntrwl,1)
                        z1i0(iscn,isza,iozprf,icntrwl,ipres,iilyr)=
     1                   output(icntrwl,2) / i0
                        z2i0(iscn,isza,iozprf,icntrwl,ipres,iilyr)=
     1                   output(icntrwl,3) / i0 
                        tti0(iscn,isza,iozprf,icntrwl,ipres,iilyr)=
     1                   output(icntrwl,4) / i0
                        if (isza.eq.1.and.iscn.eq.1) 
     1                  sb(iozprf,icntrwl,ipres,iilyr)=output(icntrwl,5)
                     endif
45                continue
55             continue
65          continue
            enddo
            close(lun)
75       continue
85    continue
c        -- output arrays i0, z1, z2, t, sb; & n-values for conparison with
c        -- current table
      write (31) logi0,z1i0,z2i0,tti0,sb
      close (31)
      stop
1000  format(1x,f8.2,3x,f10.2)
1100  format(1x,5e12.6)
c1100 format(1p,6e12.6)
2000  format('Enter the satellite index (integer):',/,
     &'OMPS1=1',/,'OMPS2=2',/,'OMPS3=3',/,'OMPS4=4')
2010  format('Wrong number.  Try again')
3100  format(1x,6e13.6)
      end
      subroutine fsavg
c
c     calculate the flux & slit average of i0,z1,z2,t,s; & c0,c1,c2 at
c     lamda0 
c
      common/array/wl(1521),fitflux(1521),cntrwl(22),output(22,5),
     1 nmb(22),stiwl(1521,5),dlam(17),weight(17)
      real sum(5), sumw(5)
c     -- initialize array csf(12,3) and sumc(12,3)
c     -- index for lamda:1 to 12; index for parameter:1 to 3
      do l=1,22
        if (cntrwl(l).ne.0.) then
          k1 = nmb(l)-100
          k2 = nmb(l)+100
          do m = 1,5
             sum(m) = 0.
             sumw(m) = 0.
          enddo
          do j = k1,k2
            if (fitflux(j).ne.0.0) then
              d = (wl(j)-cntrwl(l)) / 10.
              if (d.ge.dlam(1) .and. d.le.dlam(17)) then
                call locate(dlam,17,d,i)
                frac = (d-dlam(i)) / 
     1                 (dlam(i+1)-dlam(i))
                wt = weight(i)*(1.-frac) +
     1               weight(i+1)*frac
                do m = 1,4
                  sum(m) =
     1             sum(m)+wt*fitflux(j)*stiwl(j,m)
                  sumw(m) = sumw(m)+wt*fitflux(j)
                enddo
c                -- sb should not be weighted with solar flux
                sum(5)  = sum(5) + wt*stiwl(j,5)
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
      real*8 wavl(1521), scan(6), sza(10)
      integer nwavl, nscan, nsza
      read (lun) nwavl, nscan, nsza
c     read (lun) scan
c     read (lun) sza
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
