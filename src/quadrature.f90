! *****************************************************************************
! * WRAPPER FUNCTIONS FOR THE INTEGRATION ROUTINES
! *****************************************************************************
module quadrature
  use quadpack
  use iso_c_binding
  implicit none
  double precision, parameter :: eps_G = 1d-3 !desired integration accuracy
  integer, parameter :: maxpts = 1000000 !maximum allowed invocations of integrand

  ! Cuba 4.x C interface
  interface
     subroutine Cuhre(ndim, ncomp, integrand, userdata, nvec, &
          epsrel, epsabs, flags, mineval, maxeval, key, &
          statefile, spin, nregions, neval, fail, &
          integral, error, prob) bind(C, name="Cuhre")
       import :: c_int, c_double, c_ptr, c_funptr, c_char
       integer(c_int), value :: ndim, ncomp, nvec, flags, mineval, maxeval, key
       real(c_double), value :: epsrel, epsabs
       type(c_funptr), value :: integrand
       type(c_ptr), value :: userdata, spin
       character(c_char), intent(in) :: statefile(*)
       integer(c_int), intent(out) :: nregions, neval, fail
       real(c_double), intent(out) :: integral(*), error(*), prob(*)
     end subroutine Cuhre

     subroutine Divonne(ndim, ncomp, integrand, userdata, nvec, &
          epsrel, epsabs, flags, seed, mineval, maxeval, &
          key1, key2, key3, maxpass, border, maxchisq, mindeviation, &
          ngiven, ldxgiven, xgiven, nextra, peakfinder, &
          statefile, spin, nregions, neval, fail, &
          integral, error, prob) bind(C, name="Divonne")
       import :: c_int, c_double, c_ptr, c_funptr, c_char
       integer(c_int), value :: ndim, ncomp, nvec, flags, seed
       integer(c_int), value :: mineval, maxeval
       integer(c_int), value :: key1, key2, key3, maxpass
       real(c_double), value :: epsrel, epsabs
       real(c_double), value :: border, maxchisq, mindeviation
       integer(c_int), value :: ngiven, ldxgiven, nextra
       type(c_funptr), value :: integrand, peakfinder
       type(c_ptr), value :: userdata, spin
       real(c_double), intent(in) :: xgiven(*)
       character(c_char), intent(in) :: statefile(*)
       integer(c_int), intent(out) :: nregions, neval, fail
       real(c_double), intent(out) :: integral(*), error(*), prob(*)
     end subroutine Divonne
  end interface

  ! Abstract interface for user integrands (old style)
  abstract interface
     subroutine user_integrand(ndim, x, ncomp, f)
       integer, intent(in) :: ndim, ncomp
       double precision :: x(ndim), f(ncomp)
     end subroutine user_integrand
  end interface

  ! Module variable to store integrand procedure pointer
  procedure(user_integrand), pointer, save :: current_integrand => null()

contains
! *****************************************************************************
! * MULTI-DIMENSIONAL INTEGRATORS
! *****************************************************************************

  ! C-compatible wrapper for the integrand
  function cuba_integrand_wrapper(ndim, x, ncomp, f, userdata) result(ret) bind(C)
    integer(c_int), intent(in) :: ndim, ncomp
    real(c_double), intent(in) :: x(ndim)
    real(c_double), intent(out) :: f(ncomp)
    type(c_ptr), value :: userdata
    integer(c_int) :: ret

    call current_integrand(int(ndim), x, int(ncomp), f)
    ret = 0  ! 0 = success
  end function cuba_integrand_wrapper

  ! Calls CUBA's fastest non-deterministic (Monte Carlo) integrator
  subroutine cuba_divonne(Integrand, ndim, val, err, n)
    procedure(user_integrand) :: Integrand
    integer, intent(in) :: ndim
    real, intent(out) :: val, err
    integer, intent(out) :: n
    integer(c_int) :: nregions, neval, fail
    real(c_double) :: integral(1), error(1), prob(1)
    real(c_double) :: xgiven(1)
    character(c_char) :: statefile(1)

    current_integrand => Integrand
    statefile(1) = c_null_char

    call Divonne(int(ndim, c_int), 1_c_int, &
         c_funloc(cuba_integrand_wrapper), c_null_ptr, 1_c_int, &
         real(eps_G, c_double), 1d-12, 0_c_int, 0_c_int, &
         0_c_int, int(maxpts, c_int), &
         47_c_int, 1_c_int, 1_c_int, 5_c_int, &
         0d0, 10d0, 0.25d0, &
         0_c_int, int(ndim, c_int), xgiven, 0_c_int, c_null_funptr, &
         statefile, c_null_ptr, nregions, neval, fail, &
         integral, error, prob)

    val = real(integral(1))
    err = real(error(1))
    n = int(neval)
  end subroutine cuba_divonne

  ! Calls CUBA's deterministic integrator
  subroutine cuba_cuhre(Integrand, ndim, val, tol, err, n)
    procedure(user_integrand) :: Integrand
    integer, intent(in) :: ndim
    real, intent(out) :: val, err
    real, intent(in) :: tol
    integer, intent(out) :: n
    integer(c_int) :: nregions, neval, fail
    real(c_double) :: integral(1), error(1), prob(1)
    character(c_char) :: statefile(1)

    current_integrand => Integrand
    statefile(1) = c_null_char

    call Cuhre(int(ndim, c_int), 1_c_int, &
         c_funloc(cuba_integrand_wrapper), c_null_ptr, 1_c_int, &
         real(tol, c_double), 1d-12, 0_c_int, &
         0_c_int, int(maxpts, c_int), 0_c_int, &
         statefile, c_null_ptr, nregions, neval, fail, &
         integral, error, prob)

    val = real(integral(1))
    err = real(error(1))
    n = int(neval)
  end subroutine cuba_cuhre

! *****************************************************************************
! * ONE-DIMENSIONAL INTEGRATORS
! *****************************************************************************
  ! Manual quadrature, which doesn't use Netlib.  It Just adds up boxes.
  function f_quad(f_dy_dx, x1, x2, N, argv) result(val)
    real, intent(in) :: x1, x2
    integer, intent(in) :: N
    real, intent(in), dimension(:), optional :: argv
    interface
       real function f_dy_dx(x, argv)
         real, intent(in) :: x
         real, intent(in), dimension(:), optional :: argv
       end function f_dy_dx
    end interface
    real abserr,epsabs,epsrel,f,val,work
    integer ier,neval
    real :: dx
    integer :: i
    real, dimension(N) :: dy_dx
    epsabs = 0.0
    epsrel = eps_G
    dx = (x2 - x1)/(N - 1)
    dy_dx(1) = 0.5*f_dy_dx(x1, argv)
    dy_dx(N) = 0.5*f_dy_dx(x2, argv)
    do i = 2, N - 1
       dy_dx(i) = f_dy_dx(x1 + (i - 1)*dx, argv)
    end do
    val = sum(dy_dx)*dx
  end function f_quad

 ! Uses Netlib's adaptive quadrature algorithm 'qags'
 function f_qags(f_dy_dx, x1, x2, N, argv) result(val)
   real, intent(in) :: x1, x2
   integer, intent(in) :: N
   real, intent(in), dimension(:), optional :: argv
   interface
      real function f_dy_dx(x, argv)
        real, intent(in) :: x
        real, intent(in), dimension(:), optional :: argv
      end function f_dy_dx
   end interface
   real abserr,epsabs,epsrel,f,val,work
   integer ier,iwork,last,lenw,limit,neval
   dimension iwork(100),work(400)
   epsabs = 0.0
   epsrel = eps_G
   limit = N
   lenw = limit*4
   call qags(f_dy_dx,x1,x2,epsabs,epsrel,val,abserr,neval,ier,argv)
 end function f_qags

 ! Uses Netlib's adaptive quadrature algorithm 'qag'
 function f_qag(f_dy_dx, x1, x2, N, argv) result(val)
   real, intent(in) :: x1, x2
   integer, intent(in) :: N
   real, intent(in), dimension(:), optional :: argv
   interface
      real function f_dy_dx(x, argv)
        real, intent(in) :: x
        real, intent(in), dimension(:), optional :: argv
      end function f_dy_dx
   end interface
   real abserr,epsabs,epsrel,f,val,work
   integer ier,iwork,last,lenw,limit,neval
   dimension iwork(100),work(400)
   epsabs = 0.0
   epsrel = eps_G
   limit = N
   lenw = limit*4
   call qag(f_dy_dx,x1,x2,epsabs,epsrel,1,val,abserr,neval,ier,argv)
 end function f_qag

 ! Uses Netlib's quadrature algorithm for open (infinite) intervals, 'qagi'
 function f_qagi(f_dy_dx, x1, inf, argv) result(val)
   real, intent(in) :: x1
   integer, intent(in) :: inf
   real, intent(in), dimension(:), optional :: argv
   interface
      real function f_dy_dx(x, argv)
        real, intent(in) :: x
        real, intent(in), dimension(:), optional :: argv
      end function f_dy_dx
   end interface
   real abserr,epsabs,epsrel,f,val,work
   integer ier,iwork,last,lenw,limit,neval
   dimension iwork(100),work(400)
   epsabs = 0.0e0
   epsrel = eps_G
   limit = 100
   lenw = limit*4
   call qagi(f_dy_dx,x1,inf,epsabs,epsrel,val,abserr,neval,ier,argv)
 end function f_qagi

 ! Uses's Netlib's non-adaptive integrator 'qng'
 function f_qng(f_dy_dx, x1, x2, N, argv) result(val)
   real, intent(in) :: x1, x2
   integer, intent(in) :: N
   real, intent(in), dimension(:), optional :: argv
   interface
      real function f_dy_dx(x, argv)
        real, intent(in) :: x
        real, intent(in), dimension(:), optional :: argv
      end function f_dy_dx
   end interface
   real abserr,epsabs,epsrel,f,val,work
   integer ier,neval
   epsabs = 0.0
   epsrel = eps_G
   call qng(f_dy_dx,x1,x2,epsabs,epsrel,val,abserr,neval,ier,argv)
 end function f_qng
end module quadrature
