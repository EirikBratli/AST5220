module evolution_mod
  use healpix_types
  use params
  use time_mod
  use ode_solver
  use rec_mod
  implicit none

  ! Accuracy parameters
  real(dp),     parameter, private :: a_init   = 1.d-8
  real(dp),                private :: x_init   = log(a_init)
  real(dp),     parameter, private :: k_min    = 0.1d0 * H_0 / c
  real(dp),     parameter, private :: k_max    = 1.d3  * H_0 / c
  integer(i4b), parameter          :: n_k      = 100
  integer(i4b), parameter, private :: lmax_int = 6

  ! Perturbation quantities
  real(dp), allocatable, dimension(:,:,:) :: Theta
  real(dp), allocatable, dimension(:,:)   :: delta
  real(dp), allocatable, dimension(:,:)   :: delta_b
  real(dp), allocatable, dimension(:,:)   :: Phi
  real(dp), allocatable, dimension(:,:)   :: Psi
  real(dp), allocatable, dimension(:,:)   :: v
  real(dp), allocatable, dimension(:,:)   :: v_b
  real(dp), allocatable, dimension(:,:)   :: dPhi
  real(dp), allocatable, dimension(:,:)   :: dPsi
  real(dp), allocatable, dimension(:,:)   :: dv_b
  real(dp), allocatable, dimension(:,:,:) :: dTheta

  real(dp), allocatable, dimension(:)     :: x_t2
  ! Fourier mode list
  real(dp), allocatable, dimension(:) :: ks

  ! Book-keeping variables
  real(dp),     private :: k_current
  integer(i4b), private :: npar = 6+lmax_int

contains


  ! NB!!! New routine for 4th milestone only; disregard until then!!!
  subroutine get_hires_source_function(k, x, S)
    implicit none

    real(dp), pointer, dimension(:),   intent(out) :: k, x
    real(dp), pointer, dimension(:,:), intent(out) :: S

    integer(i4b) :: i, j
    real(dp)     :: g, dg, ddg, tau, dt, ddt, H_p, dH_p, ddHH_p, Pi, dPi, ddPi
    real(dp), allocatable, dimension(:,:) :: S_lores

    ! Task: Output a pre-computed 2D array (over k and x) for the 
    !       source function, S(k,x). Remember to set up (and allocate) output 
    !       k and x arrays too. 
    !
    ! Substeps:
    !   1) First compute the source function over the existing k and x
    !      grids
    !   2) Then spline this function with a 2D spline
    !   3) Finally, resample the source function on a high-resolution uniform
    !      5000 x 5000 grid and return this, together with corresponding
    !      high-resolution k and x arrays

  end subroutine get_hires_source_function




  ! Routine for initializing and solving the Boltzmann and Einstein equations
  subroutine initialize_perturbation_eqns
    implicit none

    integer(i4b) :: l, i
    real(dp)     :: k_step, ckHp0, dtau0
    real(dp)     :: x_init

    ! Allocate arrays for perturbation quantities
    allocate(Theta(0:n_t, 0:lmax_int, n_k))
    allocate(delta(0:n_t, n_k))
    allocate(delta_b(0:n_t, n_k))
    allocate(v(0:n_t, n_k))
    allocate(v_b(0:n_t, n_k))
    allocate(Phi(0:n_t, n_k))
    allocate(Psi(0:n_t, n_k))
    allocate(dPhi(0:n_t, n_k))
    allocate(dPsi(0:n_t, n_k))
    allocate(dv_b(0:n_t, n_k))
    allocate(dTheta(0:n_t, 0:lmax_int, n_k))
    allocate(ks(n_k))

    x_init = log(a_init)
    ! ===============================================================
    ! Task: Initialize k-grid, ks; quadratic between k_min and k_max
    ! Quadratic grid of k values
    ! ===============================================================
    
    do i=1, n_k
       
      ks(i) = k_min + (k_max-k_min)*((i-1.d0)/(n_k-1.d0))**2 
    
    end do

    ! ========================================================================
    ! Task: Set up initial conditions for the Boltzmann and Einstein equations
    ! ========================================================================

    Phi(0,:)         = 1.d0
    delta(0,:)       = 1.5d0*Phi(0,:)
    delta_b(0,:)     = delta(0,:)
    dPhi(0,:)        = 0.d0
    dPsi(0,:)        = 0.d0
    dv_b(0,:)        = 0.d0
    dTheta(0,:,:)    = 0.d0

    dtau0            = get_dtau(x_init)   
    do i = 1, n_k
       ! Optimizing for efficiency
       ckHp0        = c*ks(i)/get_H_p(x_init)
            
       v(0,i)       = ckHp0*Phi(0,i) / (2.d0)
       v_b(0,i)     = v(0,i)
       Theta(0,0,i) = 0.5d0*Phi(0,i)
       Theta(0,1,i) = -ckHp0*Phi(0,i) / (6.d0)
       Theta(0,2,i) = -20.d0*ckHp0*Theta(0,1,i) / (45.d0*dtau0)

       do l = 3, lmax_int
          Theta(0,l,i) = -(l/(2.d0*l+1.d0)) * ckHp0*Theta(0, l-1, i) / (dtau0)
       end do
       Psi(0,i)     = -Phi(0,i) - 12.d0*H_0**2.d0 * Omega_r*Theta(0,2,i) / (c*ks(i)*a_t(1))**2.d0 

    end do

  end subroutine initialize_perturbation_eqns


  subroutine integrate_perturbation_eqns
    implicit none

    integer(i4b) :: i, j_tc, k, l
    real(dp)     :: x1, x2, x_init
    real(dp)     :: eps, hmin, h1, x_tc, Hp, dt, t1, t2, ck, ckHp, dtau, eta

    real(dp), allocatable, dimension(:) :: y, y_tight_coupling, dydx

    x_init = log(a_init)
    eps    = 1.d-8
    hmin   = 1.d-10
    h1     = 1.d-5
   
    allocate(y(npar))
    allocate(dydx(npar))
    allocate(y_tight_coupling(7))
    
    ! ===================================
    ! Propagate each k-mode independently
    ! ===================================
    do k = 1, n_k
      
       k_current = ks(k)            ! Store k_current as a global module variable
       ck        = c*k_current      ! Get fewer calculations

       ! Initialize equation set for tight coupling
       y_tight_coupling(1) = delta(0,k)
       y_tight_coupling(2) = delta_b(0,k)
       y_tight_coupling(3) = v(0,k)
       y_tight_coupling(4) = v_b(0,k)
       y_tight_coupling(5) = Phi(0,k)
       y_tight_coupling(6) = Theta(0,0,k)
       y_tight_coupling(7) = Theta(0,1,k)

       ! Find the time to which tight coupling is assumed, 
       ! and integrate equations to that time
       x_tc = get_tight_coupling_time(k_current)
       

       ! Task: Integrate from x_init until the end of tight coupling, using
       !       the tight coupling equations
       ! =====================================================
       ! Integrate from initial x to the end of Tight coupling
       ! =====================================================
       
       do j_tc=1, n_t
          x1 = x_t(j_tc) ! x_t array start at index 1, while y start at index 0
          x2 = x_t(j_tc+1)
          if (x_t(j_tc) < x_tc) then  
            
             ! Integrate the tight coupling equations
             call odeint(y_tight_coupling, x1, x2, eps, h1, hmin, dfdx, bsstep, output)
          
             ! Increase effectivness
             Hp   = get_H_p(x_t(j_tc))
             eta  = get_eta(x_t(j_tc))
             ckHp = ck/Hp
             dtau = get_dtau(x_t(j_tc)) 
         
             ! Save variables for a given time
             delta(j_tc, k)    = y_tight_coupling(1)
             delta_b(j_tc, k)  = y_tight_coupling(2)
             v(j_tc, k)        = y_tight_coupling(3)
             v_b(j_tc, k)      = y_tight_coupling(4)
             Phi(j_tc, k)      = y_tight_coupling(5)
             Theta(j_tc, 0, k) = y_tight_coupling(6)
             Theta(j_tc, 1, k) = y_tight_coupling(7)
             Theta(j_tc, 2, k) = -(20.d0*ckHp*Theta(j_tc, 1, k))/(45.d0*dtau)
             
             do l=3, lmax_int
                Theta(j_tc, l, k) = -(l/(2.d0*l+1.d0))*ckHp*Theta(j_tc, l-1, k)/dtau                
             end do

             Psi(j_tc, k)      = -Phi(j_tc, k) - 12.d0*(H_0/(ck*exp(x_t(j_tc))))**2 * Omega_r*Theta(j_tc, 2, k) 
             
             ! Store the derivatives:
             call dfdx(x_t(j_tc), y_tight_coupling, dydx)
             
             dv_b(j_tc, k)      = dydx(4)
             dPhi(j_tc, k)      = dydx(5)
             dTheta(j_tc, 0, k) = dydx(6)
             dTheta(j_tc, 1, k) = dydx(7)
             dTheta(j_tc, 2, k) = 0.4d0*ckHp*Theta(j_tc, 1, k) - 0.6d0*ckHp*Theta(j_tc, 3, k) &
                                + 0.9d0*dtau*Theta(j_tc, 2, k)             
             ! loop for higher order theta
             do l=3, lmax_int-1
                dTheta(j_tc, l, k) = (l/(2.d0*l+1.d0))*ckHp*Theta(j_tc, l-1, k) &
                                   - ((l+1.d0)/(2.d0*l+1.d0))*ckHp*Theta(j_tc,l+1,k) + dtau*Theta(j_tc,l,k) 
             end do
             dTheta(j_tc, lmax_int, k) = ckHp*Theta(j_tc, lmax_int-1, k)&
                                       - c*(lmax_int+1.d0)*Theta(j_tc, lmax_int, k)/(Hp*eta) &
                                       + dtau*Theta(j_tc, lmax_int, k)

             dPsi(j_tc, k)      = -dPhi(j_tc, k) - (12.d0*H_0**2.d0*Omega_r/(ck*exp(x_t(j_tc)))**2.d0) &
                                * (dTheta(j_tc, 2, k) - 2.d0*Theta(j_tc, 2, k))
             
          else   
             i = j_tc +1
             exit
          end if  ! x condition
           
       end do   
       
       ! Set up variables for integration from the end of tight coupling until today
       
       y(1:7) = y_tight_coupling(1:7)
       y(8)   = Theta(j_tc-1, 2, k)
       do l = 3, lmax_int
          y(6+l) = Theta(j_tc-1, l, k) 
       end do

       ! =====================================
       ! Integrate from after tight coupling
       ! =====================================
       do i = j_tc, n_t
            
          ! Integrate equations from tight coupling to today
          call odeint(y, x_t(i-1), x_t(i),eps, h1, hmin, derivs, bsstep, output)

          ! Store variables at time step i in global variables
          delta(i,k)   = y(1)
          delta_b(i,k) = y(2)
          v(i,k)       = y(3)
          v_b(i,k)     = y(4)
          Phi(i,k)     = y(5)

          do l = 0, lmax_int
            Theta(i,l,k) = y(6+l) 
          end do
          Psi(i,k)     = -Phi(i,k) - 12.d0*(H_0/(ck*exp(x_t(i))))**2 * Omega_r*Theta(i, 2, k)
          
          ! Store derivatives that are required for C_l estimation
          call derivs(x_t(i), y, dydx)
          dv_b(i,k)     = dydx(4)
          dPhi(i,k)     = dydx(5)
          
          do l=0, lmax_int            
             dTheta(i,l,k) = dydx(l+6)
          end do 
          dPsi(i,k)     = -dPhi(i,k) - 12.d0*(H_0/(ck*exp(x_t(i))))**2 &
                        * Omega_r * (dTheta(i,2,k) - 2.d0*Theta(i,2,k))  
         
       end do      ! ends i loop
    end do         ! ends k loop

    
    deallocate(y_tight_coupling)
    deallocate(y)
    deallocate(dydx)
    
  end subroutine integrate_perturbation_eqns


  subroutine dfdx(x, y_tc, dydx)
    ! Solve tight coupling equations, using derivatives for intergration of y_tight_coupling
    use healpix_types
    implicit none

    real(dp), intent(in)                :: x
    real(dp), dimension(:), intent(in)  :: y_tc
    real(dp), dimension(:), intent(out) :: dydx
    real(dp)                            :: delta, delta_b, v, vb, Phi, Theta0, Theta1 ! Variables for derivation
    real(dp)                            :: dtau, ddtau, a, Hp, dHp, ckHp, k           ! Efficiency
    real(dp)                            :: Theta2, R, Psi, dPhi, dTheta0, d_delta     ! Calculation derivatives 
    real(dp)                            :: d_delta_b, dv, q, dv_b, dTheta1            ! Tight coupling eqns.

   
    ! Variables for derivation:
    delta    = y_tc(1)
    delta_b  = y_tc(2)
    v        = y_tc(3)
    vb       = y_tc(4)
    Phi      = y_tc(5)
    Theta0   = y_tc(6)
    Theta1   = y_tc(7)
    
    ! Increase effectivness:
    k        = k_current 
    dtau     = get_dtau(x)
    ddtau    = get_ddtau(x)
    a        = exp(x)
    Hp       = get_H_p(x)
    dHp      = get_dH_p(x)
    ckHp     = c*k/Hp
    
    Theta2    = -20.d0*ckHp/(45.d0*dtau)*Theta1
    R         = (4.d0*Omega_r)/(3.d0*Omega_b*a)
    Psi       = -Phi - 12.d0*(H_0/(c*k*a))**2 * Omega_r*Theta2 
    
   
    ! Calculate derivatives:
    !dPhi
    dydx(5) = Psi - ckHp**2/3.d0*Phi + H_0**2/(2.d0*Hp**2) & 
            * (Omega_m/a*delta + Omega_b/a*delta_b + 4.d0*Omega_r/a**2*Theta0)
    ! d_delta
    dydx(1) = ckHp*v - 3.d0*dydx(5)
    ! d_delta_b
    dydx(2) = ckHp*vb - 3.d0*dydx(5)
    ! v
    dydx(3) = -v - ckHp*Psi
    ! Theta0
    dydx(6) = -ckHp*Theta1 - dydx(5)

    ! Tight coupling:
    q       = ( -((1.d0-2.d0*R)*dtau + (1.d0+R)*ddtau) * (3.d0*Theta1+vb) - ckHp*Psi &
            + (1.d0-dHp/Hp)*ckHp*(-Theta0+2.d0*Theta2) - ckHp*dydx(6) ) &
            / ( (1.d0+R)*dtau + dHp/Hp - 1.d0)
    
    ! v_b
    dydx(4) = 1.d0/(1.d0+R) * ( -vb - ckHp*Psi + R*(q + ckHp*(-Theta0 + 2.d0*Theta2) - ckHp*Psi) )
    ! Theta1
    dydx(7) = 1.d0/3.d0*(q - dydx(4))
    
  end subroutine dfdx

  subroutine derivs(x, y, dydx)
    use healpix_types
    implicit none

    real(dp),               intent(in)  :: x
    real(dp), dimension(:), intent(in)  :: y
    real(dp), dimension(:), intent(out) :: dydx

    integer(i4b) :: l
    real(dp)     :: delta, delta_b, v, vb, Phi, R, Psi
    real(dp)     :: Theta0, Theta1, Theta2, Theta3, Theta4, Theta5, Theta6
    real(dp)     :: dtheta0, dTheta1, dTheta2, dv, dv_b, dPhi, d_delta, d_delta_b
    real(dp)     :: Hp, ckHp, dtau, a

    
    ! Define variables to the input array
    delta    = y(1)
    delta_b  = y(2)
    v        = y(3)
    vb       = y(4)
    Phi      = y(5)
    Theta0   = y(6)
    Theta1   = y(7)
    Theta2   = y(8)
    Theta3   = y(9)
    Theta4   = y(10)
    Theta5   = y(11)
    Theta6   = y(12)
    
    ! Increase effectivness
    Hp   = get_H_p(x)
    ckHp = c*k_current/Hp
    dtau = get_dtau(x)
    a    = exp(x)
    !print *, ckHp, 'jada', k_current
    R         = (4.d0*Omega_r)/(3.d0*Omega_b*a)
    Psi       = -Phi - 12.d0*(H_0/(c*k_current*a))**2 * Omega_r*Theta2

    ! Derivation
    
    !dPhi       
    dydx(5) = Psi - ckHp**2/3.d0*Phi + H_0**2/(2.d0*Hp**2) &
            * ( Omega_m/a*delta + Omega_b/a*delta_b + 4.d0*Omega_r/a**2*Theta0 )
    !d_delta  
    dydx(1) = ckHp*v - 3.d0*dydx(5)
    !d_delta_b 
    dydx(2) = ckHp*vb - 3.d0*dydx(5)
    ! dv
    dydx(3) = -v - ckHp*Psi 
    ! dv_b
    dydx(4) = -vb - ckHp*Psi + dtau*R*(3.d0*Theta1 + vb)
    ! dTheta0
    dydx(6) = -ckHp*Theta1 - dydx(5)
    !dTheta1   
    dydx(7) = ckHp/3.d0*Theta0 - 2.d0*ckHp/3.d0*Theta2 + ckHp/3.d0*Psi + dtau*(Theta1 + 1.d0/3.d0*vb)

    !dTheta_l    
    do l=2, lmax_int-1       
       dydx(6+l) = l*ckHp/(2.d0*l+1.d0)*y(6+l-1) - (l+1.d0)*ckHp/(2.d0*l+1.d0)*y(6+l+1) &
                 + dtau*(y(6+l) - 1.d0/10.d0*y(6+l)*abs(l==2))
    end do   
    
    dydx(12) = ckHp*Theta5 - (7.d0*c/(Hp*get_eta(x)))*Theta6 + dtau*Theta6


  end subroutine derivs 



  ! Task: Complete the following routine, such that it returns the time at which
  !       tight coupling ends. In this project, we define this as either when
  !       dtau < 10 or c*k/(H_p*dt) > 0.1 or x > x(start of recombination)
  function get_tight_coupling_time(k)
    implicit none

    real(dp), intent(in)  :: k
    real(dp)              :: get_tight_coupling_time, x, x_start
    integer(i4b)          :: i,j, n
    
    n = 1000
    x_start = -log(1.d0+1630.4d0)
    x = x_init
    do i = 1, n
       x = x + (x_start - x_init)/(n-1.d0)
       
       if ((abs(c*k/(get_H_p(x))*get_dtau(x))) >= 0.1d0 .and. abs(get_dtau(x)) <= 10.d0 .or. (x >= x_start)) then
          get_tight_coupling_time = x
          exit
       end if
    end do    
  end function get_tight_coupling_time


  subroutine write_to_file_mk3
    ! Subroutine for write quantities to file
    use healpix_types
    implicit none
    integer(i4b)                            :: i,j
    integer(i4b), allocatable, dimension(:) :: kv
    character(*), parameter                 :: path = 'Documents/Master/AST5220/src/data_mk3/'      
    allocate(kv(n_t))
    kv(1:6) = (/ 1, 10, 30, 50, 80, 100 /)
    
    ! Open and write the quantities to files
    open(1, file='phi.dat', status='replace')
    open(2, file='psi.dat', status='replace')
    open(3, file='delta.dat', status='replace')
    open(4, file='delta_b.dat', status='replace')
    open(5, file='v.dat', status='replace')
    open(6, file='v_b.dat', status='replace')
    open(7, file='Theta0.dat', status='replace')
    open(8, file='Theta1.dat', status='replace')
    open(9, file='k_values.dat', status='replace')
    open(10, file='x_t.dat', status='replace')
    !open(10, file='derivatives.dat', status='replace')
    
    do i=1, n_t
       write(1, fmt='(6(ES15.5E3, 1x))') Phi(i,kv(1)), Phi(i,kv(2)),Phi(i,kv(3)),Phi(i,kv(4)),&
            Phi(i,kv(5)), Phi(i,kv(6))
       write(2, fmt='(6(ES15.5E3, 1x))') Psi(i,kv(1)), Psi(i,kv(2)),Psi(i,kv(3)),Psi(i,kv(4)),&
            Psi(i,kv(5)),Psi(i,kv(6))
       write(3, fmt='(6(ES15.5E3, 1x))') delta(i,kv(1)), delta(i,kv(2)),delta(i,kv(3)),&
            delta(i,kv(4)),delta(i,kv(5)),delta(i,kv(6))
       write(4, fmt='(6(ES15.5E3, 1x))') delta_b(i,kv(1)), delta_b(i,kv(2)),delta_b(i,kv(3)),&
            delta_b(i,kv(4)),delta_b(i,kv(5)),delta_b(i,kv(6))
       write(5, fmt='(6(ES15.5E3, 1x))') v(i,kv(1)), v(i,kv(2)),v(i,kv(3)),v(i,kv(4)),&
            v(i,kv(5)),v(i,kv(6))
       write(6, fmt='(6(ES15.5E3, 1x))') v_b(i,kv(1)), v_b(i,kv(2)),v_b(i,kv(3)),v_b(i,kv(4)),&
            v_b(i,kv(5)),v_b(i,kv(6))
       write(7, fmt='(6(ES15.5E3, 1x))') Theta(i,0,kv(1)), Theta(i,0,kv(2)),Theta(i,0,kv(3)),&
            Theta(i,0,kv(4)),Theta(i,0,kv(5)),Theta(i,0,kv(6)) 
       write(8, fmt='(6(ES15.5E3, 1x))') Theta(i,1,kv(1)), Theta(i,1,kv(2)),Theta(i,1,kv(3)),&
            Theta(i,1,kv(4)),Theta(i,1,kv(5)),Theta(i,1,kv(6))
       write(10, fmt='(3(ES15.5E3, 1x))') x_t(i), ks(i), ks(i)*c/H_0

       !do j=1,6
       !   write(9 , '(13(ES15.5E3,1x))') dPhi(i, kv(j)), dPsi(i, kv(j)), dv_b(i,kv(j)), &
        !       dTheta(i,0,kv(j)), dTheta(i,1,kv(j)), dTheta(i,2,kv(j)), &
         !      dTheta(i,3,kv(j)), dTheta(i,4,kv(j)), dTheta(i,5,kv(j)), dTheta(i,6,kv(j)) 
       !end do   
    end do   

    do i=1,6
       write(9,*) ks(kv(i)), ks(kv(i))*c/H_0, kv(i)
    end do

    ! close the files
    do i=1, 9
       close(i)
    end do   
    
    deallocate(kv)

  end subroutine write_to_file_mk3


end module evolution_mod
