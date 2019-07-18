module time_mod
  use healpix_types
  use params
  use spline_1D_mod
  use ode_solver
  implicit none

  integer(i4b)                           :: n_t                ! Number of x-values
  real(dp),    allocatable, dimension(:) :: x_t                ! Grid of relevant x-values
  real(dp),    allocatable, dimension(:) :: a_t                ! Grid of relevant a-values
  real(dp),    allocatable, dimension(:) :: z_t                ! Grid of z values

  integer(i4b)                           :: n_eta              ! Number of eta grid poins
  real(dp),    allocatable, dimension(:) :: x_eta              ! Grid points for eta
  real(dp),    allocatable, dimension(:) :: eta, eta2          ! Eta and eta'' at each grid point
  real(dp),    allocatable, dimension(:) :: z_eta              ! redshift grid for eta grid

  ! Density parameters
  real(dp),    allocatable, dimension(:) :: Omega_Radiation, Omega_Mass, Omega_Baryons, Omega_L

contains

  subroutine initialize_time_mod
    implicit none

    integer(i4b) :: i, n, n1, n2
    real(dp)     :: z_start_rec, z_end_rec, z_0, x_start_rec, x_end_rec, x_0, dx, x_eta1, x_eta2, a_init

    ! Define two epochs, 1) during and 2) after recombination.
    n1          = 200                       ! Number of grid points during recombination
    n2          = 300                       ! Number of grid points after recombination
    n_t         = n1 + n2                   ! Total number of grid points
    z_start_rec = 1630.4d0                  ! Redshift of start of recombination
    z_end_rec   = 614.2d0                   ! Redshift of end of recombination
    z_0         = 0.d0                      ! Redshift today
    x_start_rec = -log(1.d0 + z_start_rec)  ! x of start of recombination
    x_end_rec   = -log(1.d0 + z_end_rec)    ! x of end of recombination
    x_0         = 0.d0                      ! x today

    n_eta       = 1000                      ! Number of eta grid points (for spline)
    a_init      = 1.d-10                    ! Start value of a for eta evaluation
    x_eta1      = log(a_init)               ! Start value of x for eta evaluation
    x_eta2      = 0.d0                      ! End value of x for eta evaluation


    ! ==========================================================
    ! Making x and a grids
    allocate(x_t(n_t))
    allocate(a_t(n_t))
    allocate(z_t(n_t))

    z_t(1) = z_start_rec
    x_t(1) = x_start_rec
    a_t(1) = exp(x_t(1))
    do i=1, n_t-1

       if (i <= n1) then
          z_t(i+1) = z_t(i) + (z_end_rec - z_start_rec)/ (n1)
          x_t(i+1) = -log(1.d0 + z_t(i+1))
       else
          z_t(i+1) = z_t(i) + (z_0 - z_end_rec) / (n2 - 1.0)
          x_t(i+1) = -log(1.d0 + z_t(i+1))
       end if
       a_t(i+1) = 1.0/(1 + z_t(i+1))

    end do


    ! ==========================================================
    ! Compute the conformal time:
    allocate(x_eta(n_eta))
    allocate(eta(n_eta))
    allocate(eta2(n_eta))

    ! Make grid of x_eta:
    x_eta(1) = x_eta1
    do i=1, n_eta

       x_eta(i+1) = x_eta(i) + (x_eta2 - x_eta1)/(n_eta)

    end do
    x_eta(n_eta) = x_eta2

    ! Integration of eta, using odeint from ode_solver.f90
    dx = abs((x_eta(2) - x_eta(1)))
    eta(1) = a_init*c / (H_0*sqrt(omega_r))

    do i=2, n_eta

       eta(i) =  eta(i-1) ! + dx*c/get_H_p(x_eta(i-1))
       call odeint(eta(i:i), x_eta(i-1), x_eta(i), 1.0d-8*dx, dx, dx*0, d_eta_dx, bsstep, output)

    end do

    call spline(x_eta, eta, 1.d30, 1.d30, eta2)


    ! ==========================================================
    ! Find the density parameters as function of x and z?
    allocate(z_eta(n_eta))
    allocate(Omega_Radiation(n_eta))
    allocate(Omega_Mass(n_eta))
    allocate(Omega_Baryons(n_eta))
    allocate(Omega_L(n_eta))

    do i=1, n_eta

       z_eta(i) = exp(-x_eta(i)) - 1.0
       Omega_Radiation(i) = Omega_r*(1.0+z_eta(i))**4
       Omega_Mass(i) = Omega_m*(1.0+z_eta(i))**3
       Omega_Baryons(i) = Omega_b*(1.0+z_eta(i))**3
       Omega_L(i) = Omega_lambda

    end do

  end subroutine initialize_time_mod


  ! Task: Write a function that computes H at given x
  function get_H(x)
    implicit none

    real(dp), intent(in) :: x
    real(dp)             :: get_H

    get_H = H_0*sqrt((Omega_b + Omega_m)*exp(-3*x) + Omega_r*exp(-4*x) + Omega_lambda)

  end function get_H

  ! Task: Write a function that computes H' = a*H  at given x
  function get_H_p(x)
    implicit none

    real(dp), intent(in) :: x
    real(dp)             :: get_H_p

    get_H_p = exp(x) * get_H(x)

  end function get_H_p

  ! Task: Write a function that computes dH'/dx at given x
  function get_dH_p(x)
    implicit none

    real(dp), intent(in) :: x
    real(dp)             :: get_dH_p, dH

    dH = (-3*(Omega_m + Omega_b)*exp(-3*x) - 4*Omega_r*exp(-4*x))/(2*get_H(x))
    get_dH_p = exp(x)*(dH + get_H(x))

  end function get_dH_p

  ! Task: Write a function that computes eta(x), using the previously precomputed splined function
  function get_eta(x_in)
    implicit none

    real(dp), intent(in) :: x_in
    real(dp)             :: get_eta

    get_eta = splint(x_eta, eta, eta2, x_in) ! c / get_H_p(x_in)

  end function get_eta

  subroutine output(x, y)
    ! Input in the odeint solver
    use healpix_types
    implicit none
    real(dp),               intent(in) :: x
    real(dp), dimension(:), intent(in) :: y
  end subroutine output

  subroutine d_eta_dx(x, eta, dydx)
    ! Input parameter of the derivative of eta in odeint solver
    use healpix_types
    implicit none
    real(dp),               intent(in) :: x
    real(dp), dimension(:), intent(in) :: eta
    real(dp), dimension(:), intent(out) :: dydx

    dydx = c/get_H_p(x)

  end subroutine d_eta_dx


end module time_mod
