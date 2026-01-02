program main
  use M_kracken
  use plplot
  use utils
  use cosmology
  use structure
  use survey3
  implicit none
  integer :: Nz, Nl
  real(kind=plflt) :: z1, z2, l1, l2, dOmega, flux_lim
  real(kind=plflt) :: w, Omega_Lambda_0, Omega_M_0, sigma_8, tau_0, lambda_0, alpha, beta, sig_tau
  real(kind=plflt) :: sig_lambda, phi, psi, rho, chi, gamma
  character(len=255) :: cmd
  integer :: plparseopts_rc

  plparseopts_rc = plparseopts(PL_PARSE_FULL)

  cmd = &
       &' -h .false. -z1 0.01 -z2 2.5 -nz 5 -l1 44.0 -l2 46.0 -nl 5 '//&
       &' -do 12.0 -fl 1.25e-13 -w -1.0 -om 0.30 -s8 0.8 -g 0.55 '
  call kracken ('cmd', cmd)
  if (lget('cmd_h')) then
     print *, ''
     print *, 'Usage:'
     print *, ''
     print *, '  cosmosurvey [[OPTION]] [[PLPLOT OPTIONS]]*'
     print *, ''
     print *, 'Options (Default):'
     print *, '  -h                Print help information'
     print *, '  -z1  (0.01)       Redshift start'
     print *, '  -z2  (2.50)       Redshift end'
     print *, '  -nz  (5)          Number of redshift bins'
     print *, '  -l1  (44.0)       log10(L1) luminosity start, [L1] = ergs/s'
     print *, '  -l2  (46.0)       log10(L2) luminosity end, [L2] = ergs/s'
     print *, '  -nl  (5)          Number of luminosity bins'
     print *, '  -do  (12.0)       Survey solid angle dOmega, [dOmega] = steradians'
     print *, '  -fl  (1.25e-13)   Survey flux limit, [fl] = ergs/s/cm^2'
     print *, '  -w   (-1.0)       Dark Energy equation of state parameter w'
     print *, '  -om  (0.30)       Matter density parameter Omega_M'
     print *, '  -s8  (0.80)       Power spectrum normalization sigma_8'
     print *, '  -g   (0.55)       Linder growth index gamma'
     print *, ''
     print *, '*NOTE:  See http://plplot.sourceforge.net/docbook-manual/plplot-html-5.9.5/advanced.html#arguments'
     print *, ''
  else
     z1 = rget('cmd_z1')
     z2 = rget('cmd_z2')
     Nz = iget('cmd_nz')
     l1 = rget('cmd_l1')
     l2 = rget('cmd_l2')
     Nl = iget('cmd_nl')
     dOmega = rget('cmd_do')
     flux_lim = rget('cmd_fl')

     theta_G = f_lcdm()
     theta_G%w = rget('cmd_w')
     theta_G%Omega_M_0 = rget('cmd_om')
     theta_G%sigma_8 = rget('cmd_s8')
     theta_G%gamma = rget('cmd_g')

     ! call plsdev('psc')
     ! call plsetopt('-o', '-')
     ! call plsetopt('-ori', '1')
     ! call plsetopt('-geometry', '400x300')
     ! call plsetopt('-width', '1')
     call plinit()
     call plfontld(1)
     call plot_survey()
     call plend()
  end if
contains
  subroutine plot_survey ()
    real(kind=plflt) :: dz, dl
    real, allocatable :: dN(:, :), z(:), l(:)
    real, dimension(:), pointer :: N
    integer :: i, j
    real(kind=plflt), allocatable :: xx(:), yy(:), zz(:, :)

    allocate(dN(Nz, Nl))
    allocate(z(Nz))
    allocate(l(Nl))

    dz = (z2 - z1)/real(Nz)
    dl = (l2 - l1)/real(Nl)
    z = [(z1 + (real(i-1)+0.5)*(z2 - z1)/real(Nz), i = 1, Nz)]
    l = [(l1 + (real(i-1)+0.5)*(l2 - l1)/real(Nl), i = 1, Nl)]
    N => f_N(theta_G, real(z1), real(z2), real(l1), real(l2), Nz, Nl, real(dOmega), real(flux_lim))
    dN = reshape(source=N, shape=[Nz, Nl], order=[2, 1])
    where (dN.lt.1e-1) dN = 1e-1
    dN = log10(dN)
    call hist_3d(z, l, dN)

    allocate(xx(size(z, 1)))
    allocate(yy(size(l, 1)))
    allocate(zz(size(dN, 1), size(dN, 2)))
    xx = z
    yy = l
    zz = dN
    call plschr(0.0_plflt, 1.0_plflt)
    call pladv(0)
    call plvpor(0.0_plflt, 1.0_plflt, 0.0_plflt, 1.0_plflt )
    call plwind(-2.0_plflt, 2.0_plflt, -2.0_plflt, 4.0_plflt)
    call plw3d(2.0_plflt, 2.0_plflt, 4.0_plflt, minval(xx), maxval(xx), &
         &minval(yy), maxval(yy), minval(zz), maxval(zz), 45.0_plflt, -60.0_plflt)
    call plcol0(3)
    call plot3d(xx, yy, zz, 3, .TRUE.)
    call plcol0(1)
    call plbox3('bnstu','', 0.0_plflt, 0, 'blnstu', '', 0.0_plflt, 0, 'bmlnstu', '', 0.0_plflt, 0)
    call plcol0(2)
    call plschr(0.0_plflt, 1.0_plflt)
    call plbox3('u','#frz', 0.0_plflt, 0, 'u', '#frL#dObs#u (ergs/s/cm#u2#d)', 0.0_plflt, 0, 'u', '#frdN', 0.0_plflt, 0)
    call plcol0(0)
    call plcol0(2)
    call plschr(0.0_plflt, 1.0_plflt)
    call plmtex('t', -0.9_plflt, 0.5_plflt, 0.5_plflt, '#frCluster Number Counts')
  end subroutine plot_survey
end program main

