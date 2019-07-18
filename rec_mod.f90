module rec_mod
  use healpix_types
  use params
  use time_mod
  use ode_solver
  use spline_1D_mod
  implicit none

  integer(i4b),                        private :: n                           ! Number of grid points
  real(dp), allocatable, dimension(:), private :: x_rec, a, z_rec             ! Grid
  real(dp), allocatable, dimension(:), private :: tau, tau2, tau22            ! tau and second derivatives
  real(dp), allocatable, dimension(:), private :: n_e, n_e2                   ! Electron density
  real(dp), allocatable, dimension(:), private :: log_ne, log_ne2             ! Splined (log of) electron density, n_e
  real(dp), allocatable, dimension(:), private :: g, g2, g22                  ! visibility function
  real(dp), allocatable, dimension(:), private :: s_g, s_g2, s_g22            ! Splined of visibility function
  real(dp), allocatable, dimension(:), private :: s_tau, s_tau2, s_tau22      ! Splined optical depth
  real(dp), allocatable, dimension(:), private :: logtau, logtau2, logtau22   ! log of optical depth
  
  

contains

  subroutine initialize_rec_mod
    implicit none
    
    real(dp), allocatable, dimension(:) :: X_e ! Fractional electron density, n_e / n_H
    integer(i4b) :: i, j, k
    real(dp)     :: saha_limit, y, T_b, n_b, dydx, xmin, xmax, dx, f, n_e0, X_e0 
    real(dp)     :: xstart, xstop, xe_const, step, step_min, eps, yp1, yp2
    logical(lgt) :: use_saha
    
    
    saha_limit = 0.99d0       ! Switch from Saha to Peebles when X_e < 0.99
    xstart     = log(1.d-10)  ! Start grids at a = 10^-10
    xstop      = 0.d0         ! Stop  grids at a = 1
    n          = 1000         ! Number of grid points between xstart and xstopo

    ! Constants used in functions:
    eps        = 1.d-8
    step_min   = 1.d-5
    yp1        = 1.d30
    yp2        = 1.d30
    
    ! Allocate arrays
    allocate(x_rec(n))
    allocate(X_e(n))
    allocate(tau(n))
    allocate(tau2(n))
    allocate(tau22(n))
    allocate(s_tau(n))
    allocate(s_tau2(n))
    allocate(s_tau22(n))
    allocate(logtau(n))
    allocate(logtau2(n))
    allocate(logtau22(n))
    allocate(n_e(n))
    allocate(n_e2(n))
    allocate(log_ne(n))
    allocate(log_ne2(n))
    allocate(g(n))
    allocate(g2(n))
    allocate(g22(n))
    allocate(s_g(n))
    allocate(s_g2(n))
    allocate(s_g22(n))

    
    ! Task: Fill in x (rec) grid and scale factor and red shift
    ! ==============================================
    ! Filling in the x grid, and for a and z as well
    ! ==============================================

    x_rec(1) = xstart
    dx = (xstop - xstart)/(n - 1.0)
    do i=1, n-1
       
       x_rec(i+1) = x_rec(i) + dx
       
    end do
    x_rec(n)   = xstop
    
    a          = exp(x_rec) 
    z_rec      = 1.d0/a - 1.d0 


    ! ==================================
    ! Computing the electron density
    ! ==================================
    ! Task: Compute X_e and n_e at all grid times
    step       = abs(x_rec(1) - x_rec(2))
    use_saha = .true.
    do i = 1, n
       
       n_b = Omega_b*rho_c/(m_H*a(i)**3.d0)
       if (use_saha) then 
          ! Use the Saha equation
          
          T_b = T_0/a(i)
          xe_const = (1.d0/(n_b))*(m_e*k_b*T_b/(2.d0*pi*hbar*hbar))**(1.5d0)*exp(-epsilon_0/(k_b*T_b))
     
          X_e(i) = 2.d0/(sqrt(1.d0 + 4.d0/xe_const) + 1.d0)                    ! compphys variant
          !X_e(i) = 0.5d0*(-xe_const + sqrt(xe_const**2.d0 + 4.d0*xe_const))    ! abc formula
          
          if (X_e(i) < saha_limit) use_saha = .false.
       else
          ! Use the Peebles equation
          X_e(i) = X_e(i-1)
   
          call odeint(X_e(i:i), x_rec(i-1), x_rec(i), eps, step, step_min, dXe_dx, rkqs, output)
          
       end if
       n_e(i) = n_b*X_e(i)

    end do

    ! Task: Compute splined (log of) electron density function
    log_ne = log(n_e)
    call spline(x_rec, log_ne, yp1, yp2, log_ne2)
    write(*,*) 'Done with X_e, calculating tau and g:'
    


    ! =============================================
    ! Compute optical depth and visibility function
    ! =============================================

    ! Task: Compute optical depth at all grid points
    tau(n) = 0.0d0                     ! optical depth today (i=n) locally

    ! Integrate tau: 
    do i=n, 1, -1

       tau(i) = tau(i+1)
       call odeint(tau(i:i), x_rec(i+1), x_rec(i), eps, step, step_min, dtaudx, rkqs, output)

    end do
    tau(n) = tau(n-1)
    ! Task: Compute splined (log of) optical depth
    ! Task: Compute splined second derivative of (log of) optical depth
    logtau = log(tau)

    call spline(x_rec, tau, yp1, yp2, tau2)             ! spline of tau
    call spline(x_rec, tau2, yp1, yp2, tau22)           ! Spline of tau' 

    call spline(x_rec, logtau, yp1, yp2, logtau2)       ! Spline of log of tau
    call spline(x_rec, logtau2, yp1, yp2, logtau22)     ! Spline of log of tau'

    ! Get the tau values
    do i = 1, n
       
       s_tau(i)    = get_tau(x_rec(i))
       s_tau2(i)   = get_dtau(x_rec(i))
       s_tau22(i)  = get_ddtau(x_rec(i))

    end do
    
    ! Task: Compute splined visibility function
    ! Task: Compute splined second derivative of visibility function
    g = -s_tau2*exp(-s_tau)
    
    call spline(x_rec, g, yp1, yp2, g2)
    call spline(x_rec, g2, yp1, yp2, g22)

    ! Get the g values
    do i = 1, n
       
       s_g(i)    = get_g(x_rec(i))
       s_g2(i)   = get_dg(x_rec(i))
       s_g22(i)  = get_ddg(x_rec(i))
       
    end do   


    ! =========================
    ! Write to files:
    ! =========================

    write(*,*) 'Writing X_e, tau and g to file:'
    
    open(unit=1, file='electron_density.dat', status='replace')
    do i=1, n
      
       write(unit=1, fmt='(10(ES15.7E3, 1x))') z_rec(i), x_rec(i), X_e(i), n_e(i), s_tau(i), s_tau2(i), s_tau22(i), s_g(i), s_g2(i), s_g22(i)

    end do   
    close(unit=1) 
  end subroutine initialize_rec_mod



  
! ================================
! Extra needed subroutines
! ================================


  subroutine dXe_dx(x, X_e, dXe)
    ! Subroutine to compute dX_e/dx:
    
    use healpix_types
    implicit none
    real(dp), intent(in)               :: x
    real(dp), dimension(:), intent(in) :: X_e
    real(dp)                           :: n1s, L_alpha, C_r, a, H     
    real(dp)                           :: L_2to1, beta2, beta, alpha2, phi2, T_b, n_H
    real(dp), dimension(:), intent(out):: dXe
    
    !
    H    = get_H(x)
    a    = exp(x)
    n_H  = Omega_b*rho_c/(m_H*a**3.d0)           ! number density of hydrogen
    T_b  = T_0/a                             ! baryon temperature
    
    ! Atomic physics parameters
    phi2    = 0.448d0 * log(epsilon_0/(k_b*T_b))        
    alpha2  = (64.d0*pi/sqrt(27.d0*pi))*(alpha/m_e)**2.d0*sqrt(epsilon_0/(k_b*T_b))*phi2 *(hbar**2.d0/c) ! [m³/s]
    beta    = alpha2 * (m_e*k_b*T_b/(2.d0*pi*hbar**2.d0))**(1.5d0) * exp(-epsilon_0/(k_b*T_b))           ! [s⁻¹]
    !beta2   = beta*exp(0.75d0*epsilon_0/(k_b*T_b))                                                      ! [s⁻¹]
    beta2   = alpha2*(m_e*k_b*T_b/(2.d0*pi*hbar**2.d0))**(1.5d0) * exp(-0.25d0*epsilon_0/(k_b*T_b))
    n1s     = (1.d0 - X_e(1))*n_H                                                                        ! [s⁻¹]
    L_alpha = get_H(x)*(3.0d0*epsilon_0)**3.d0 / ((n1s*(8.d0*pi)**2.d0) * (hbar*c)**3.d0)                ! [s⁻¹]
    L_2to1  = 8.227d0                                                                                    ! [s⁻¹]
    C_r     = (L_2to1 + L_alpha) / (L_2to1 + L_alpha + beta2)
    
    ! Peebles eqaution
    dXe = C_r*(beta*(1.d0 - X_e(1)) - n_H*alpha2*X_e(1)**2.d0) / (get_H(x))
    
  end subroutine dXe_dx  
 

  subroutine dtaudx(x, tau, dtau)
    ! Compute the optical depth

    use healpix_types
    implicit none
    real(dp), intent(in)                :: x
    real(dp), dimension(:), intent(in)  :: tau
    real(dp), dimension(:), intent(out) :: dtau

    dtau = -(sigma_T*c*get_n_e(x))/(get_H(x))

  end subroutine dtaudx



! ==========================================
! Functions used in calculation of tau and g
! ==========================================

  ! Task: Complete routine for computing n_e at arbitrary x, using precomputed information
  ! Hint: Remember to exponentiate...
  function get_n_e(x)
    implicit none

    real(dp), intent(in) :: x
    real(dp)             :: get_n_e
    
    ! Need to use slint since using log_ne2, which is splined
    get_n_e = exp(splint(x_rec, log_ne, log_ne2, x))

  end function get_n_e


  ! Task: Complete routine for computing tau at arbitrary x, using precomputed information
  function get_tau(x)
    implicit none

    real(dp), intent(in) :: x
    real(dp)             :: get_tau

    get_tau = exp(splint(x_rec, logtau, logtau2, x))

  end function get_tau

  ! Task: Complete routine for computing the derivative of tau at arbitrary x, 
  ! using precomputed information
  function get_dtau(x)
    implicit none

    real(dp), intent(in) :: x
    real(dp)             :: get_dtau

    get_dtau = splint_deriv(x_rec, logtau, logtau2, x)*get_tau(x)

  end function get_dtau

  ! Task: Complete routine for computing the second derivative of tau at arbitrary x, 
  ! using precomputed information
  function get_ddtau(x)
    implicit none

    real(dp), intent(in) :: x
    real(dp)             :: get_ddtau

    get_ddtau = splint(x_rec, logtau2, logtau22, x)*get_tau(x) + (get_dtau(x)**2.d0)/get_tau(x)

  end function get_ddtau


  ! Task: Complete routine for computing the visibility function, g, at arbitray x
  function get_g(x)
    implicit none

    real(dp), intent(in) :: x
    real(dp)             :: get_g

    get_g = splint(x_rec, g, g2, x)

  end function get_g

  ! Task: Complete routine for computing the derivative of the visibility function, g, 
  ! at arbitray x
  function get_dg(x)
    implicit none

    real(dp), intent(in) :: x
    real(dp)             :: get_dg

    get_dg = splint_deriv(x_rec, g, g2, x)

  end function get_dg

  ! Task: Complete routine for computing the second derivative of the visibility function, g, 
  ! at arbitray x
  function get_ddg(x)
    implicit none

    real(dp), intent(in) :: x
    real(dp)             :: get_ddg

    get_ddg = splint(x_rec, g2, g22, x)

  end function get_ddg

 
end module rec_mod
