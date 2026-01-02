program main
  use M_kracken
  use plplot
  use utils
  use cosmology
  use structure
  use survey2, setenv2 => setenv
  implicit none
  integer, parameter :: N = 100
  real(kind=plflt) :: z1, z2, dOmega, flux_lim
  real(kind=plflt) :: w, Omega_Lambda_0, Omega_M_0, sigma_8, tau_0, lambda_0, alpha, beta, sig_tau
  real(kind=plflt) :: sig_lambda, phi, psi, rho, chi, gamma
  integer :: Nz, Nl
  integer :: plparseopts_rc

  call kracken ('cmd', '-w -1.0 -om 0.30 -s8 0.9 -g 0.55 -h .false.')
  if (lget('cmd_h')) then
     print *, ''
     print *, 'Usage:'
     print *, ''
     print *, '  cosmotest [[OPTION]] [[PLPLOT OPTIONS]]*'
     print *, ''
     print *, 'Options (Default):'
     print *, '  -h                Print help information'
     print *, '  -w   (-1.0)       Dark Energy equation of state parameter w'
     print *, '  -om  (0.30)       Matter density parameter Omega_M'
     print *, '  -s8  (0.80)       Power spectrum normalization sigma_8'
     print *, '  -g   (0.55)       Linder growth index gamma'
     print *, ''
     print *, '*NOTE:  See http://plplot.sourceforge.net/docbook-manual/plplot-html-5.9.5/advanced.html#arguments'
     print *, ''
  else
     theta_G = f_lcdm()
     z1 = 0.01
     z2 = 2.5
     flux_lim = 1.25e-13
     dOmega = 4.0*pi
     theta_G%w = rget('cmd_w')
     theta_G%Omega_M_0 = rget('cmd_om')
     theta_G%sigma_8 = rget('cmd_s8')
     theta_G%gamma = rget('cmd_g')
     plparseopts_rc = plparseopts(PL_PARSE_FULL)
     call plinit()
     call plfontld(1)
     call plot1()
     call plot2()
     call plot3()
     call plot4()
     call plot5()
     call plot6()
     call plot7()
     call plot8()
     call plot9()
     call plend()
  end if
contains
  subroutine plot1 ()
    integer :: i
    real(kind=plflt) :: z(N), dv_dOmega_dz(N)

    do i = 1, N
       z(i) = z1 + (i-1)*(z2 - z1)/real(N - 1)
       dv_domega_dz(i) = f_dV_dOmega_dz(real(z(i)))
    end do
    call pladv(0)
    call plvasp(1.0_plflt)
    call plwind(minval(z), maxval(z), minval(dv_domega_dz), 8d0)
    call plcol0(1)
    call pllsty(1)
    call plschr(0.0_plflt, 1.2_plflt)
    call plbox('abcnts', 0._plflt, 0, 'abcnts', 0._plflt, 0)
    call plcol0(2)
    call pllab( '#frz', '#frdV#dco#u/dz    (H#d0#u/c)#u3#d', &
         &'#frComoving Volume' )
    call plcol0(3)
    call plline(z, 4*3.14159*dv_dOmega_dz)
  end subroutine plot1

  subroutine plot2 ()
    real(kind=plflt) :: z(N), G(N)
    integer :: i

    do i = 1, N
       z(i) = z1 + (i-1)*(z2 - z1)/real(N - 1)
       G(i) = f_G(real(z(i)))/f_G(0.0)
    end do
    call plschr(0.0_plflt, 1.5_plflt)
    call pladv(0)
    call plvasp(1.0_plflt)
    call plwind(z1, z2, 0d0, 1d0)
    call plcol0(1)
    call pllsty(1)
    call plschr(0.0_plflt, 1.2_plflt)
    call plbox('abcnts', 0._plflt, 0, 'abcnts', 0._plflt, 0)
    call plcol0(2)
    call pllab('#frz', '#frG(z)', '#frGrowth Factor' )
    call plcol0(3)
    call plline(z, G)
  end subroutine plot2

  subroutine plot3 ()
    real(kind=plflt), parameter :: mu0 = 13.0, mu1 = 15.0
    real(kind=plflt) :: mu(N), ln_inv_sigma(N)
    integer :: i

    do i = 1, N
       mu(i) = mu0 + (i-1)*(mu1 - mu0)/real(N - 1)
       ln_inv_sigma(i) = f_ln_inv_sigma(10.0**real(mu(i)), 0.0)
    end do
    call plschr(0.0_plflt, 1.5_plflt)
    call pladv(0)
    call plvasp(1.0_plflt)
    call plwind(minval(mu), maxval(mu), minval(ln_inv_sigma), maxval(ln_inv_sigma))
    call plcol0(1)
    call pllsty(1)
    call plschr(0.0_plflt, 1.2_plflt)
    call plbox('abclnts', 0._plflt, 0, 'abcnts', 0._plflt, 0)
    call plcol0(2)
    call pllab( '#frM', '#frln#gs#fr#u-1', '#frMass-Function Factor #gs#fr(M)')
    call plcol0(3)
    call plline(mu, ln_inv_sigma) 
  end subroutine plot3

  subroutine plot4 ()
    real(kind=plflt), parameter :: mu0 = 13.0, mu1 = 15.0
    real(kind=plflt) :: mu(N), alpha_eff(N)
    integer :: i

    do i = 1, N
       mu(i) = mu0 + (i-1)*(mu1 - mu0)/real(N - 1)
       alpha_eff(i) = f_alpha_eff(10.0**real(mu(i)))
    end do
    call plschr(0.0_plflt, 1.5_plflt)
    call pladv(0)
    call plvasp(1.0_plflt)
    call plwind(minval(mu), maxval(mu), minval(alpha_eff), maxval(alpha_eff))
    call plcol0(1)
    call pllsty(1)
    call plschr(0.0_plflt, 1.2_plflt)
    call plbox('abclnts', 0._plflt, 0, 'abcnts', 0._plflt, 0)
    call plcol0(2)
    call pllab( '#frM', '#ga#fr#deff#fr#u(M)', '#ga#fr#deff')
    call plcol0(3)
    call plline(mu, alpha_eff) 
  end subroutine plot4

  subroutine plot5 ()
    real(kind=plflt), parameter :: mu0 = 13.0, mu1 = 15.0
    real(kind=plflt) :: mu(N), mass_fraction(N)
    integer :: i

    do i = 1, N
       mu(i) = mu0 + (i-1)*(mu1 - mu0)/real(N - 1)
       mass_fraction(i) = f_mass_fraction(f_ln_inv_sigma(10.0**real(mu(i)), 0.0))
    end do
    call plschr(0.0_plflt, 1.5_plflt)
    call pladv(0)
    call plvasp(1.0_plflt)
    call plwind(minval(mu), maxval(mu), minval(mass_fraction), maxval(mass_fraction))
    call plcol0(1)
    call pllsty(1)
    call plschr(0.0_plflt, 1.2_plflt)
    call plbox('abclnts', 0._plflt, 0, 'abcnts', 0._plflt, 0)
    call plcol0(2)
    call pllab( '#frM', '#frMass Fraction (M)', '#frMass Fraction')
    call plcol0(3)
    call plline(mu, mass_fraction) 
  end subroutine plot5

  subroutine plot6 ()
    real(kind=plflt) :: mu0 = log10(8e12), mu1 = log10(5e15)
    real(kind=plflt) :: mu(N), dN_dV_dmu1(N), dN_dV_dmu2(N), dN_dV_dmu3(N)
    integer :: i

    do i = 1, N
       mu(i) = mu0 + i*(mu1 - mu0)/real(N - 1)
       dN_dV_dmu1(i) = f_dN_dV_dmu(10.0**real(mu(i)), 0.0)/log(10.)
       dN_dV_dmu2(i) = f_dN_dV_dmu(10.0**real(mu(i)), 0.5)/log(10.)
       dN_dV_dmu3(i) = f_dN_dV_dmu(10.0**real(mu(i)), 1.0)/log(10.)
    end do
    call plschr(0.0_plflt, 1.5_plflt)
    call pladv(0)
    call plvasp(1.0_plflt)
    call plwind(mu0, mu1, -10.0d0, -3.0d0)
    call plcol0(1)
    call pllsty(1)
    call plschr(0.0_plflt, 1.2_plflt)
    call plbox('abclnts', 0._plflt, 0, 'abclnts', 0._plflt, 0)
    call plcol0(2)
    call pllab( '#frM', '#frdn/dlnM (M)', '#frMass Function')
    call plcol0(1)
    call plline(mu, log10(dN_dV_dmu1)) 
    call plptex(15.0_plflt, -3.5_plflt, 0.0_plflt, 0.0_plflt, 1.0_plflt, '#frz = 0.0')
    call plline([15.1_plflt, 15.5_plflt], [-3.5_plflt, -3.5_plflt])
    call plcol0(2)
    call pllsty(2)
    call plline(mu, log10(dN_dV_dmu2)) 
    call plptex(15.0_plflt, -4.0_plflt, 0.0_plflt, 0.0_plflt, 1.0_plflt, '#frz = 0.5')
    call plline([15.1_plflt, 15.5_plflt], [-4.0_plflt, -4.0_plflt])
    call plcol0(3)
    call pllsty(3)
    call plline(mu, log10(dN_dV_dmu3)) 
    call plptex(15.0_plflt, -4.5_plflt, 0.0_plflt, 0.0_plflt, 1.0_plflt, '#frz = 1.0')
    call plline([15.1_plflt, 15.5_plflt], [-4.5_plflt, -4.5_plflt])
  end subroutine plot6

  subroutine plot7 ()
    real(kind=plflt) :: z(N), dN_dV(N)
    integer :: i

    call setenv2(Mlim = 3e14)
    do i = 1, N
       z(i) = z1 + (i-1)*(z2 - z1)/real(N - 1)
       dN_dV(i) = f_dN_dV(real(z(i)))
    end do
    call plschr(0.0_plflt, 1.5_plflt)
    call pladv(0)
    call plvasp(1.0_plflt)
    call plwind(z1, z2, minval(log10(dN_dV)), maxval(log10(dN_dV)))
    call plcol0(1)
    call pllsty(1)
    call plschr(0.0_plflt, 1.2_plflt)
    call plbox('abcnts', 0.0_plflt, 0, 'abclnts', 0.0_plflt, 0)
    call plcol0(2)
    call pllab('', '', '')
    call plcol0(3)
    call plline(z, log10(dN_dV))
  end subroutine plot7

  subroutine plot8 ()
    real(kind=plflt) :: z(N), dN_dz(N)
    integer :: i

    call setenv2(Mlim = 3e14)
    do i = 1, N
       z(i) = z1 + (i-1)*(z2 - z1)/real(N - 1)
       dN_dz(i) = f_dN_dz(real(z(i)), real(dOmega))
    end do
    call plschr(0.0_plflt, 1.5_plflt)
    call pladv(0)
    call plvasp(1.0_plflt)
    call plwind(z1, z2, 0.0d0, 1d5)
    call plcol0(1)
    call pllsty(1)
    call plschr(0.0_plflt, 1.2_plflt)
    call plbox('abcnts', 0._plflt, 0, 'abcnts', 0._plflt, 0)
    call plcol0(2)
    call pllab('#frz', '#frdN/dV#dCO#u (z)', '#frCluster Counts')
    call plptex(1.0_plflt, 0.9e5_plflt, 0.0_plflt, 0.0_plflt, 0.5_plflt, '#frM > 3 #(0727) 10#u14#d M#dSun')
    call plcol0(3)
    call plline(z, dN_dz)
  end subroutine plot8

  subroutine plot9 ()
    real(kind=plflt) :: z(N), lum(N)
    integer :: i

    do i = 1, N
       z(i) = z1 + i*(z2 - z1)/real(N - 1)
       lum(i) = f_luminosity(real(flux_lim), real(z(i)))
    end do
    call plschr(0.0_plflt, 1.5_plflt)
    call pladv(0)
    call plvasp(1.0_plflt)
    call plwind(z1, z2, minval(lum), maxval(lum))
    call plcol0(1)
    call pllsty(1)
    call plschr(0.0_plflt, 1.2_plflt)
    call plbox('abcnts', 0.0_plflt,  0, 'abclnts', 0.0_plflt, 0)
    call plcol0(2)
    call pllab('#frz', '#frL (z)', '#frFlux Limit')
    call plschr(0.0_plflt, 0.9_plflt)
    ! call plptex(1.0_plflt, minval(lum) + 0.1*(maxval(lum)-minval(lum))&
    !      &, 0.0_plflt, 0.0_plflt, 0.5_plflt, '#frflux = 1 #(0727) 10#u-13#d ergs/s/cm#u2')
    call plcol0(3)
    call plline(z, lum)
  end subroutine plot9
end program main

