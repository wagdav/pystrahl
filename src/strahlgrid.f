      program strahlgrid


c     define rho_volume and flux surface averages
c     and write results into STRAHL geometry file
c
c
c     INPUT VARIABLES:
c     shot	     shot  number
c     time	     time in shot
c     r0	     R-coordinate of magnetic axis
c     z0	     Z-coordinate of magnetic axis
c     uloop	     loop voltage
c     rho(nr)	     array of normalized flux labels
c     rmaj_in(nr)    array of R-coordinates for z=z0 at HFS (inboard)
c       	     for the given array of rho values
c     rmaj_out(nr)   array of R-coordinates for z=z0 at LFS (outboard)
c       	     for the given array of rho values
c     r(nr,NANGLE)   R-coordinate of flux surfaces for NANGLE poloidal
c       	     angles (equidistant grid)
c     z(nr,NANGLE)   Z-coordinate of flux surfaces
c     bp(nsep,NANGLE) poloidal field strength of flux surfaces for
c                     NANGLE poloidal angles (equidistant grid, 0-360)
c     bt(nsep,NANGLE) toroidal field strength of flux surfaces for
c                     NANGLE poloidal angles (equidistant grid, 0-360)
c
c     OUTPUT VARIABLES:
c     fvol(nsep)     array of volumes encolosed by flux surface
c     rho_vol(nr)         volume rho
c     b(nsep)             < total magnetic field>
c     q(nsep)             safety factor
c     r_bt(nsep)          < R B_T>
c     r_bp_b(nsep)        < R**2 B_p**2/B**2>
c     r2i(nsep)           < 1/R**2 >
c     fc(nsep)            fraction of circulating particles
c     cos_b_sq(nsep,NF)   cosine fourier coefficient of B^2
c     sin_b_sq(nsep,NF)   sine fourier coefficient of B^2
c     cos_b_lnB(nsep,NF)  cosine fourier coefficient of BlnB
c     sin_b_lnB(nsep,NF)  sine fourier coefficient of BlnB
c

      implicit none
c
c     number of rho values
c
      integer NRMAX
      parameter (NRMAX = 60)
c
c     number of rho values up to separatrix
c
      integer nsepMAX
      parameter (NSEPMAX = 60)
c
c     number of theta values on flux_tube (must be odd)
c
      integer NANGLE
      parameter (NANGLE = 103)
c
c     number of lambda values (must be odd)
c
      integer NLAMBDA
      parameter (NLAMBDA = 73)
c
c     number of fourier coeficients
c
      integer NF
      parameter (NF = 3)
c
c     number of thetata values for dense grid
c
      integer NANDENS
      parameter (NANDENS = NANGLE*(2*NF+1))
c
c     constants
c
      real pi,radian
      parameter(pi=3.14159)
      parameter(radian=180./pi)

c     ************************************
c     global variables
c     ************************************
      integer
     >     shot, nr, nsep

      real
     >     time,                !time in shot in s
     >     r0,z0,               !r0,z0 of magnetic axis
     >     uloop,               !loop voltage
     >     rho(NRMAX),             !rho_norm grid
     >     rho_vol(NRMAX),         !rho_volume grid
     >     fvol(NRMAX),            !volume inside flux surfaces
     >     rmaj_out(NRMAX),        !major radius low field side
     >     rmaj_in(NRMAX),         !major radius high field side
     >     r(NRMAX,NANGLE),        !R co-ordinate of flux surfaces
     >     z(NRMAX,NANGLE),        !Z co-ordinate of flux surfaces
     >     bp(NSEPMAX, NANGLE),    !poloidal field on flux surfaces
     >     bt(NSEPMAX, NANGLE)     !toroidal field on flux surfaces

c     magnetic fields and averages

      real
     >     b(NSEPMAX),             !< total magnetic field>
     >     r_bt(NSEPMAX),          !< R B_T>
     >     r_bp_b(NSEPMAX),        !< R**2 B_p**2/B**2>
     >     r2i(NSEPMAX),           !< 1/R**2 >
     >     fc(NSEPMAX),            !fraction of circulating particles
     >     cos_b_sq(NSEPMAX,NF),   !cosine fourier coefficient of B^2
     >     sin_b_sq(NSEPMAX,NF),   !sine fourier coefficient of B^2
     >     cos_b_lnB(NSEPMAX,NF),  !cosine fourier coefficient of BlnB
     >     sin_b_lnB(NSEPMAX,NF)   !sine fourier coefficient of BlnB

c     q_profile

      real
     >     q(NSEPMAX)

c     ************************************
c     local variables
c     ************************************
      integer
     >     i,j,k,npoly

      real
     >     b_max

c     functions
      real
     >     rz_vol, interp, simpson_int

c     rho ....
      real
     >     grad(NRMAX),            ! < gradient(rho_vol) >
     >     grad_sq(NRMAX)          ! < gradient(rho_vol)**2 >

c     dummy variables
      real
     >     dum,
     >     dum1(NANGLE),dum2(NANGLE),
     >     y(NANGLE), yy(NANDENS)

c     poloidal angle
      real
     >     dtheta,theta1,theta2,
     >     theta(NANGLE),
     >     dth_dens,
     >     th_dens(NANDENS)

c     poloidal length
      real
     >     dlp_dth(NSEPMAX, NANGLE)

c     integration variable lambda for trapped particle fraction
      real
     >     lambda(NLAMBDA), dlambda, ff(NLAMBDA)

c     magnetic fields and averages

      real
     >     b_sq(NSEPMAX),          !< B**2>
     >     b_sq_inv(NSEPMAX),      !< 1./B**2>
     >     gen_theta(NSEPMAX,NANDENS), !generalized theta
     >     bp_invers(NSEPMAX)     !< 1/B_poloidal >

c     save file
      character*12 file

c     home directory
      character*40 home

      integer nstrlen, index


      open(unit=10,file='strahlgrid.in',status='old')
      read(10,*)  shot
      read(10,*)  time
      read(10,*)  r0, z0, uloop
      read(10,*)  nr, nsep
      read(10,*) (rho(i),i=1,nr)
      read(10,*) (rmaj_in(i),i=1,nr)
      read(10,*) (rmaj_out(i),i=1,nr)
      read(10,*) ((r(i,j),i=1,nr),j=1,NANGLE)
      read(10,*) ((z(i,j),i=1,nr),j=1,NANGLE)
      read(10,*) ((bp(i,j),i=1,nsep),j=1,NANGLE)
      read(10,*) ((bt(i,j),i=1,nsep),j=1,NANGLE)
      close(10)
c*************************************************************
c     rho_vol
c*************************************************************
c
c
c     define rho_vol from the volume enclosed by the flux
c     surface rho_vol = sqrt(Vol/2/Pi**2/Ro)
c     unit = m
c
      do j=2, nr
         do i=1,NANGLE
            dum1(i) = r(j,i)
            dum2(i) = z(j,i)
         enddo
         fvol(j) = rz_vol(dum1,dum2,NANGLE)
         rho_vol(j) = sqrt(fvol(j)/2./pi**2/r0)
      enddo
c
c*************************************************************
c     end of rho_vol grid
c*************************************************************
c
c     large radius in midplane (Low Field Side)
c
      do j=1,nr
         rmaj_out(j) = rmaj_out(j)/r0
      enddo
c
c     large radius in midplane (High Field Side)
c
      do j=1,nr
         rmaj_in(j) = rmaj_in(j)/r0
      enddo
c*************************************************************
c     The following quantities are just computed inside the
c     spearatrix
c*************************************************************
c
c     define grid for poloidal angle theta
c
      theta1 = 0.0
      theta2 = 360.0
      dtheta= (theta2-theta1)/(NANGLE-1)
      do i=1,NANGLE
         theta(i) = theta1 + dtheta*(i-1)
      enddo
c
c     poloidal angle to radian
c
      dtheta = dtheta/radian
      do i=1,NANGLE
         theta(i)=theta(i)/radian
      enddo
c
c     change of poloidal flux surface length with theta
c
      do j=2,nsep-1
         do i=1,NANGLE
            if (i.eq.1.or.i.eq.NANGLE) then
               dlp_dth(j,i) = sqrt( (r(j,2)-r(j,NANGLE-1))**2 +
     >              (z(j,2)-z(j,NANGLE-1))**2 ) /2. / dtheta
            else
               dlp_dth(j,i) = sqrt( (r(j,i+1)-r(j,i-1))**2+
     >              (z(j,i+1)-z(j,i-1))**2 )/ 2. / dtheta
            endif
         enddo
      enddo
c
c     q = 1/2/pi * < Bt/ R B_poloidal >
c
      do j=2,nsep-1
         do i=1,NANGLE
            y(i) = bt(j,i)*dlp_dth(j,i)/bp(j,i)/r(j,i)
         enddo
         q(j) = simpson_int(y, NANGLE, dtheta)/2./Pi
      enddo
      q(1) = interp(rho(2),q(2),4,1,0.0)
      q(nsep) = interp(rho,q,nsep-1,1,rho(nsep))
c
c     < 1/ B_poloidal >  =
c     Normalization Integral for all averages
c
      do j=2,nsep-1
         do i=1,NANGLE
            y(i) = dlp_dth(j,i)/bp(j,i)
         enddo
         bp_invers(j) = simpson_int(y, NANGLE, dtheta)
      enddo
      bp_invers(1) = interp(rho(2),bp_invers(2),4,1,0.0)
      bp_invers(nsep) = interp(rho,bp_invers,nsep-1,1,rho(nsep))
c
c    < grad(rho_vol) > = <R B_p> dPsi/drho_vol
c                      = <R B_p> * Int(dl_p/B_p) / (2 pi R0 rho_vol)
      grad(1) = 1.
      do j=2,nsep-1
         do i=1,NANGLE
            y(i) = r(j,i)*dlp_dth(j,i)
         enddo
         grad(j) = simpson_int(y, NANGLE, dtheta) /
     >        ( 2.*pi*rho_vol(j)*r0 )
      enddo
      grad(nsep) = interp(rho,grad,nsep-1,1,rho(nsep))
      grad(2) = grad(1) + (grad(3)-grad(1))*rho_vol(2)/rho_vol(3)
c
c    < grad(rho_vol) >
c       = <(R B_p)**2> (dPsi/drho_vol)**2
c       = <(R B_p)**2> * (Int(dl_p/B_p) / (2 pi R0 rho_vol) )**2
c
      grad_sq(1) = 1.0
      do j=2,nsep-1
         do i=1,NANGLE
            y(i) = r(j,i)**2*bp(j,i)*dlp_dth(j,i)
         enddo
         grad_sq(j) = simpson_int(y, NANGLE, dtheta) * bp_invers(j) /
     >        (2.*pi *rho_vol(j)*r0)**2
      enddo
      grad_sq(1) = 1.
      grad_sq(nsep) = interp(rho,grad_sq,nsep-1,1,rho(nsep))
      grad_sq(2) = grad_sq(1) +
     >     (grad_sq(3)-grad_sq(1))*rho_vol(2)/rho_vol(3)

c*************************************************************
c     Now Surface Averages for NEOART-Code
c*************************************************************
c
c     < B_total >
c
      b(1) = abs(bt(1,1))
      do j=2,nsep-1
         do i=1,NANGLE
            y(i) = sqrt(bp(j,i)**2+bt(j,i)**2)*dlp_dth(j,i)/bp(j,i)
         enddo
         b(j) = simpson_int(y, NANGLE, dtheta) / bp_invers(j)
      enddo
      b(nsep) = interp(rho,b,nsep-1,1,rho(nsep))
c
c     < B_total**2 >
c
      b_sq(1) = bt(1,1)**2
      do j=2,nsep-1
         do i=1,NANGLE
            y(i) = (bp(j,i)**2+bt(j,i)**2)*dlp_dth(j,i)/bp(j,i)
         enddo
         b_sq(j) = simpson_int(y, NANGLE, dtheta) / bp_invers(j)
      enddo
      b_sq(nsep) = interp(rho,b_sq,nsep-1,1,rho(nsep))
c
c     < 1./B_total**2 >
c
      b_sq_inv(1) = 1./bt(1,1)**2
      do j=2,nsep-1
         do i=1,NANGLE
            y(i) = 1./(bp(j,i)**2+bt(j,i)**2)*dlp_dth(j,i)/bp(j,i)
         enddo
         b_sq_inv(j) = simpson_int(y, NANGLE, dtheta) / bp_invers(j)
      enddo
      b_sq_inv(nsep) = interp(rho,b_sq_inv,nsep-1,1,rho(nsep))
c
c     < R B_T >
c
      r_bt(1) = r0*abs(bt(1,1))
      do j=2,nsep-1
         do i=1,NANGLE
            y(i) = r(j,i)*abs(bt(j,i))*dlp_dth(j,i)/bp(j,i)
         enddo
         r_bt(j) = simpson_int(y, NANGLE, dtheta) / bp_invers(j)
      enddo
      r_bt(nsep) = interp(rho,r_bt,nsep-1,1,rho(nsep))
c
c     < R**2 B_p**2/ B**2 >
c
      r_bp_b(1) = 0.
      do j=2,nsep-1
         do i=1,NANGLE
            y(i) = r(j,i)**2 * bp(j,i)**2 / (bt(j,i)**2+bp(j,i)**2)
     >           *     dlp_dth(j,i)/bp(j,i)
         enddo
         r_bp_b(j) = simpson_int(y, NANGLE, dtheta) / bp_invers(j)
      enddo
      r_bp_b(nsep) = interp(rho,r_bp_b,nsep-1,1,rho(nsep))
c
c     < 1/R**2 >
c
      r2i(1) = 1./r0**2
      do j=2,nsep-1
         do i=1,NANGLE
            y(i) = 1./ r(j,i)**2
     >           *     dlp_dth(j,i)/bp(j,i)
         enddo
         r2i(j) = simpson_int(y, NANGLE, dtheta) / bp_invers(j)
      enddo
      r2i(nsep) = interp(rho,r2i,nsep-1,1,rho(nsep))
c
c     Int_0^1/Bmax (lambda dlambda / < sqrt(1 - lambda B) >)
c
      fc(1) = 1.
      ff(1) = 0.
      do j=2,nsep-1
         b_max = 0.0
         do i=1,NANGLE
            dum1(i) = sqrt(bp(j,i)**2+bt(j,i)**2)
            b_max = max(b_max, dum1(i))  ! find maximum field
         enddo
         dlambda = 1./b_max / (NLAMBDA-1.)
         do k=2,NLAMBDA
            lambda(k) =dlambda * (k-1)
c
c     < sqrt(1 - lambda B) >
c
            do i=1,NANGLE
               y(i) = sqrt(max(0.,1.-dum1(i)*lambda(k)))
               y(i) = y(i)*dlp_dth(j,i)/bp(j,i)
            enddo
            ff(k) = lambda(k)*bp_invers(j) /
     >              simpson_int(y, NANGLE, dtheta)
         enddo
         fc(j) = 3./4.*b_sq(j)*simpson_int(ff, NLAMBDA, dlambda)
      enddo
      fc(nsep) = interp(rho,fc,nsep-1,1,rho(nsep))

C    *****************************************************************
C    Now Fourier Coefficients for NEOART up to m=NF
c    *****************************************************************
c
c     fine theta grid
c
      dth_dens= dtheta*(NANGLE-1.)/(NANDENS-1.)
      do i=1,NANDENS
         th_dens(i) = dth_dens *(i-1)
      enddo
c
c    definde generalized poloidal angle on grid
c
      npoly=2
      dum1(1) = 0.
      do i=3, NANGLE, 2
         dum1((i+1)/2) = theta(i)
      enddo
      do j=2,nsep-1
         dum2(1) = 0.
         do i=3, NANGLE, 2    !integration for all odd # of points
            do k=1,i
               y(k) =  sqrt(bp(j,k)**2+bt(j,k)**2)*dlp_dth(j,k)/bp(j,k)
            enddo
            dum2((i+1)/2)  = 2.*pi*simpson_int(y, i, dtheta)
     >           /b(j)/bp_invers(j)
         enddo
         gen_theta(j,1) = 0.   !interpolation
         gen_theta(j,NANDENS) = 2.*pi
         do i=2,NANDENS-1
            gen_theta(j,i) = interp(dum1,dum2,(NANGLE+1)/2,
     >           npoly,th_dens(i))
         enddo
         do i=1,NANDENS
            yy(i) = gen_theta(j,i)
         enddo
      enddo
c
c    < cos (m theta_gen) * B_total^2 >
c    < sin (m theta_gen) * B_total^2 >
c
      do k=1, NF
         cos_b_sq(1,k) = 0.0
         sin_b_sq(1,k) = 0.0
      enddo
      do j=2,nsep-1
         do i=1,NANGLE
            dum1(i) = (bp(j,i)**2+bt(j,i)**2)*dlp_dth(j,i)/bp(j,i)
         enddo
         do k=1,NF
            do i=1,NANDENS       !cosine terms
               yy(i) = interp(theta,dum1,NANGLE,npoly,th_dens(i))
     >              * cos( k * gen_theta(j,i))
            enddo
            cos_b_sq(j,k) =
     >           simpson_int(yy, NANDENS, dth_dens) / bp_invers(j)
            do i=1,NANDENS      !sine terms
               yy(i) = interp(theta,dum1,NANGLE,npoly,th_dens(i))
     >              * sin( k * gen_theta(j,i))
            enddo
            sin_b_sq(j,k) =
     >           simpson_int(yy, NANDENS, dth_dens) / bp_invers(j)
         enddo
      enddo
      do k=1, NF
         cos_b_sq(nsep,k) =
     >        interp(rho,cos_b_sq(1,k),nsep-1,1,rho(nsep))
         sin_b_sq(nsep,k) =
     >        interp(rho,sin_b_sq(1,k),nsep-1,1,rho(nsep))
      enddo
c
c    < cos (m theta_gen) * BlnB >
c    < sin (m theta_gen) * BlnB >
c
      do k=1, NF
         cos_b_lnB(1,k) = 0.0
         sin_b_lnB(1,k) = 0.0
      enddo
      do j=2,nsep-1
         do i=1,NANGLE
            dum = sqrt(bp(j,i)**2+bt(j,i)**2)
            dum1(i) = dum*alog(dum)*dlp_dth(j,i)/bp(j,i)
         enddo
         do k=1,NF
            do i=1,NANDENS       !cosine terms
               yy(i) = interp(theta,dum1,NANGLE,npoly,th_dens(i))
     >              * cos( k * gen_theta(j,i))
            enddo
            cos_b_lnB(j,k) =
     >           simpson_int(yy, NANDENS, dth_dens) / bp_invers(j)
            do i=1,NANDENS      !sine terms
               yy(i) = interp(theta,dum1,NANGLE,npoly,th_dens(i))
     >              * sin( k * gen_theta(j,i))
            enddo
            sin_b_lnB(j,k) =
     >           simpson_int(yy, NANDENS, dth_dens) / bp_invers(j)
         enddo
      enddo
      do k=1, NF
         cos_b_lnb(nsep,k) =
     >        interp(rho,cos_b_lnb(1,k),nsep-1,1,rho(nsep))
         sin_b_lnb(nsep,k) =
     >        interp(rho,sin_b_lnb(1,k),nsep-1,1,rho(nsep))
      enddo

c*************************************************************
c     Make Plots
c*************************************************************



c*************************************************************
c     Save everything in a file
c*************************************************************
c
      write(*,*) ' '
          write(*,'(a)') 'cv  rho volume(LCFS)[cm]   R_axis[cm]  '//
     >                ' U_loop[V]    time[s] '
      write(*,'(9x,f7.1,10x,f7.1,6x,f7.4,6x,f7.4)')
     >            rho_vol(nsep)*100., r0*100., uloop, time

      write(*,*) ' '
      write(*,*) ' '
      write(*,'(a)') 'cv  number of grid points  '
     >     //'points up to separtrix'
     >     //'  fourier coefficients'
      write(*,'(i15,i22,i22)') nr,nsep,NF

      write(*,*) ' '
      write(*,*) ' '
      write(*,'(a)') 'cv  sqrt( (Psi-Psi_ax) / (Psi_sep - Psi_ax) )   '
      write(*,1000) (rho(i),i=1,nr)

      write(*,*) ' '
      write(*,*) ' '
      write(*,'(a)') 'cv     rho volume / rho_volume(LCFS)   '
      write(*,1000) (rho_vol(i)/rho_vol(nsep),i=1,nr)

      write(*,*) ' '
      write(*,*) ' '
      write(*,'(a)') 'cv   large radius low field side / R_axis'
      write(*,1000) (rmaj_out(i),i=1,nr)

      write(*,*) ' '
      write(*,*) ' '
      write(*,'(a)') 'cv   large radius high field side / R_axis'
      write(*,1000) (rmaj_in(i),i=1,nr)

      write(*,*) ' '
      write(*,*) ' '
      write(*,'(a)') 'cv  safety factor'
      write(*,1000) (q(i),i=1,nsep)

      write(*,*) ' '
      write(*,*) ' '
      write(*,'(a)') 'cv  fraction of circulating particles'
      write(*,1000) (fc(i),i=1,nsep)

      write(*,*) ' '
      write(*,*) ' '
      write(*,'(a)') 'cv  Integral( dl_p / B_p) [m/T]'
      write(*,1000) (bp_invers(i),i=1,nsep)

      write(*,*) ' '
      write(*,*) ' '
      write(*,'(a)') 'cv  < B_total > [T]'
      write(*,1000) (b(i),i=1,nsep)

      write(*,*) ' '
      write(*,*) ' '
      write(*,'(a)') 'cv  < B_total**2 > [T**2]'
      write(*,1000) (b_sq(i),i=1,nsep)

      write(*,*) ' '
      write(*,*) ' '
      write(*,'(a)') 'cv  < 1./B_total**2 > [1/T**2]'
      write(*,1000) (b_sq_inv(i),i=1,nsep)

      write(*,*) ' '
      write(*,*) ' '
      write(*,'(a)') 'cv  < R B_T > [m*T]'
      write(*,1000) (r_bt(i),i=1,nsep)

      write(*,*) ' '
      write(*,*) ' '
      write(*,'(a)') 'cv  < R**2 B_p**2/B**2 > [m**2]'
      write(*,1000) (r_bp_b(i),i=1,nsep)

      write(*,*) ' '
      write(*,*) ' '
      write(*,'(a)') 'cv  < 1/R**2 > [1/m**2]'
      write(*,1000) (r2i(i),i=1,nsep)

      write(*,*) ' '
      write(*,*) ' '
      write(*,'(a)') 'cv  <cos (m theta) B_total**2> [T**2]'
      do k=1,NF
         write(*,2000) (cos_b_sq(i,k),i=1,nsep)
      enddo

      write(*,*) ' '
      write(*,*) ' '
      write(*,'(a)') 'cv  <sin (m theta) B_total**2> [T**2]'
      do k=1,NF
         write(*,2000) (sin_b_sq(i,k),i=1,nsep)
      enddo

      write(*,*) ' '
      write(*,*) ' '
      write(*,'(a)') 'cv  <cos (m theta) B_total ln(B_total)> [Tln(T)]'
      do k=1,NF
         write(*,2000) (cos_b_lnb(i,k),i=1,nsep)
      enddo

      write(*,*) ' '
      write(*,*) ' '
      write(*,'(a)') 'cv  <sin (m theta) B_total ln(B_total)> [Tln(T)]'
      do k=1,NF
         write(*,2000) (sin_b_lnb(i,k),i=1,nsep)
      enddo

      write(*,*) ' '
      write(*,*) ' '
      write(*,'(a)') 'cv  < grad(rho_volume) > '
      write(*,1000) (grad(i),i=1,nsep)

      write(*,*) ' '
      write(*,*) ' '
      write(*,'(a)') 'cv  < grad(rho_volume**2) > '
      write(*,1000) (grad_sq(i),i=1,nsep)


      write(*,*) ' '
      write(*,*) ' '
      write(*,'(a)') 'cv  Bp at LFS equator [T] '
      write(*,1000) (bp(i,1),i=1,nsep)


 1000 format(7(f10.5))
 2000 format(7(1PE10.2))

c     save data

c     open(unit=10,file='strahlgrid.out')
c     write(10,*)  nr, nsep, nf
c     write(10,*) (fvol(i),i=1,nr)
c     write(10,*) (rho_vol(i),i=1,nr)
c     write(10,*) (b(i),i=1,nsep)
c     write(10,*) (q(i),i=1,nsep)
c     write(10,*) (r_bt(i),i=1,nsep)
c     write(10,*) (r_bp_b(i),i=1,nsep)
c     write(10,*) (r2i(i),i=1,nsep)
c     write(10,*) (fc(i),i=1,nsep)
c     write(10,*) ((cos_b_sq(i,j),i=1,nsep),j=1,NF)
c     write(10,*) ((sin_b_sq(i,j),i=1,nsep),j=1,NF)
c     write(10,*) ((cos_b_lnb(i,j),i=1,nsep),j=1,NF)
c     write(10,*) ((sin_b_lnb(i,j),i=1,nsep),j=1,NF)
c     close(10)

      stop

      end
c
c     ***********************************************************
c
      real function simpson_int(y,n,h)

      implicit none

      integer i,n
      real y(n),h

      if (mod(n,2).ne.1) then
       write(*,*) 'error during Simpson Integration! '
       write(*,*) ' Number of Points must be odd !!'
       simpson_int=0.0
       return
      endif

      simpson_int = y(1)+y(n)
      do i=2,n-1,2
       simpson_int = simpson_int+4.*y(i)
      enddo
      do i=3,n-2,2
       simpson_int = simpson_int+2.*y(i)
      enddo
      simpson_int = simpson_int*h/3.

      return
      end
c
c     ***********************************************************
c
      real function interp(x,y,npts,npoly,xin)
c
c     interpolation and extrapolation routine from BEVINGTON
c     for montonically increasing x-grid
c
c     x - independent variable of data array
c     y - dependent varaible of data array
c     npts - number of data points
c     npoly - degree of fitting polynomial (up to 9)
c     xin - input value of x
c

      implicit none

      integer
     >npts,nterms,npoly,
     >i,k,j,i1,i2,ix,
     >imax,ixmax

      real
     >x(npts),y(npts),
     >xin,denom,
     >delta(10),a(10),deltax,prod,sum

      nterms = npoly+1
c
c     search for appropriate value of x(npts)
c

11    do 19 i=1,npts
      if (xin - x(i)) 13,17,19
13    i1 = i - nterms/2
      if (i1) 15,15,21
15    i1=1
      goto 21
17    interp=y(i)
18    goto 61
19    continue

      i1 = npts - nterms + 1
21    i2 = i1 + nterms -1
      if (npts - i2) 23,31,31
23    i2 = npts
      i1 = i2 - nterms + 1
25    if (i1) 26,26,31
26    i1 = 1
27    nterms = i2 - i1 + 1

c
c     evaluate deviations delta
c

31    denom = x(i1+1) - x(i1)
      deltax =  (xin - x(i1)) / denom
      do 35 i=1,nterms
      ix = i1 + i - 1
35    delta(i) = (x(ix) - x(i1)) / denom

c
c     accumulate coefficients
c

      a(1) =  y(i1)
      do 50 k=2,nterms
      prod=1.
      sum=0.
      imax = k - 1
      ixmax = i1 + imax
      do 49 i=1,imax
      j = k - i
      prod = prod * (delta(k) - delta(j))
49    sum = sum - a(j)/prod
50    a(k) = sum + y(ixmax)/prod

c
c     accumulate sum of expansion
c

51    sum= a(1)
      do 57 j=2,nterms
      prod=1.
      imax = j-1
      do 56 i=1,imax
56    prod = prod*(deltax - delta(i))
57    sum = sum + a(j)*prod
60    interp = sum
61    return
      end

      real function rz_vol(r,z,npts)
c
c     define volume of the torus that is described by the
c     given r,z-contour. The first and the last point are
c     connected by a horizontal line.
c
c     r - array of npts large radii
c     z - array of npts z-coordinates
c     npts - number of data points
c

      implicit none

      real Pi
      parameter (Pi = 3.14159)

      integer
     >i, npts

      real
     >r(npts), z(npts),
     >dz, r_mean

      rz_vol = 0.
      do i=2, npts
       dz = z(i) - z(i-1)
       r_mean = .5*(r(i) + r(i-1))
       rz_vol = rz_vol + r_mean**2*dz
      enddo
      rz_vol = abs(rz_vol) * Pi

      return
      end
