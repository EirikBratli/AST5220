module cl_mod
  use healpix_types
  use evolution_mod
  use sphbess_mod
  implicit none

  real(dp),     pointer, dimension(:)       :: x_hires, k_hires, z_spline
  real(dp),     pointer, dimension(:,:)     :: S, S2, j_l, j_l2


contains

  ! Driver routine for (finally!) computing the CMB power spectrum
  subroutine compute_cls
    implicit none

    integer(i4b) :: i, j, l, l_num, x_num, n_spline, n_hires, n_l
    real(dp)     :: dx, dk, S_func, j_func, z, eta, eta0, x0, x_min, x_max, d, e, j_l_splint
    integer(i4b), allocatable, dimension(:)       :: ls
    real(dp),     allocatable, dimension(:,:)       :: integrand, integrand2
    !real(dp),     allocatable,     dimension(:,:)     :: j_l, j_l2 ! pointer
    real(dp),     pointer,     dimension(:)       :: x_arg, cls, cls2 ! pointer
    real(dp),     pointer,     dimension(:)       :: k, x, int_cl       ! pointer
    real(dp),     pointer,     dimension(:,:,:,:) :: S_coeff    ! pointer
    !real(dp),     allocatable,     dimension(:,:)     :: S, S2      ! pointer
    real(dp),     allocatable, dimension(:,:)     :: Theta, int_arg
    real(dp),     allocatable, dimension(:)       :: j_l_spline, j_l_spline2, besseltest
    real(dp),     allocatable, dimension(:)       :: ls_hires, ls_dp, cls_hires

    real(dp)           :: t1, t2, integral
    logical(lgt)       :: exist
    character(len=128) :: filename
    real(dp), allocatable, dimension(:) :: y, y2
    real(dp), allocatable, dimension(:) :: ckH0, c5, c10, c17,c25, c35, c44
    real(dp), allocatable, dimension(:) :: T5, T10, T17,T25, T35, T44

    l_num     = 44
    n_hires   = 5000
    n_spline  = 5400

    allocate(ls(l_num))
    ! Set up which l's to compute
    ls = (/ 2, 3, 4, 6, 8, 10, 12, 15, 20, 30, 40, 50, 60, 70, 80, 90, 100, &
         & 120, 140, 160, 180, 200, 225, 250, 275, 300, 350, 400, 450, 500, 550, &
         & 600, 650, 700, 750, 800, 850, 900, 950, 1000, 1050, 1100, 1150, 1200 /)

    n_l = ls(l_num)

    allocate(x_hires(n_hires))
    allocate(k_hires(n_hires))
    allocate(S(n_hires, n_hires))

    allocate(z_spline(n_spline))    ! Note: z is *not* redshift, but simply the dummy argument of j_l(z)
    allocate(j_l(n_spline, l_num))
    allocate(j_l2(n_spline, l_num))
    allocate(x_arg(n_hires))

    allocate(integrand(l_num, n_hires))
    allocate(integrand2(l_num, n_hires))
    allocate(Theta(l_num, n_hires))
    allocate(int_arg(l_num, n_hires))
    allocate(int_cl(l_num))
    allocate(cls(l_num))
    allocate(cls2(l_num))

    allocate(ls_dp(l_num))
    allocate(ls_hires(n_l))
    allocate(cls_hires(n_l))


    ! Create double precision and high resolution l arrays
    do l=1, l_num
       ls_dp(l) = ls(l)
    end do
    do i=1, n_l
       ls_hires(i) = i
    end do

    !open(unit=41, file='cl_data.dat', status='replace')   ! default (h = 0.7, Om=0.224, Ob=0.046, n=1)
    !open(unit=41, file='cl_h08.dat', status='replace')    ! h = 0.6, 0.8   xx done
    !open(unit=41, file='cl_Ob005.dat', status='replace')  ! Ob=0.04, 0.05   xx
    !open(unit=41, file='cl_n096.dat', status='replace')    ! n =0.9, (0.96), 1.1  xxx
    !open(unit=41, file='cl_Om025.dat', status='replace')   ! Om=0.2, 0.25  xx
    !open(unit=41, file='cl_Or11e4.txt', status='replace')    ! Or=5e-5, 1.1e-4 x
    !open(unit=41, file='cl_fit5.dat', status='replace')    ! h=0.4, n=0.96, Ob=0.09, Om=0.25 x good
    !open(unit=41, file='cl_fit8.dat', status='replace')    ! h=0.4, n=0.96, Ob=0.09, Om=0.225  x good
    open(unit=41, file='cl_fit10.dat', status='replace')    ! h=0.4, n=0.96, Ob=0.1, Om=0.275, Or=1.0e-4  x
    !open(unit=41, file='cl_fit6.dat', status='replace')    ! h=0.4, n=0.96, Ob=0.07, Om=0.25  x

    open(unit=42, file='Cl_file.dat', status='replace')
    open(unit=43, file='Transfer_func.dat', status='replace')
    open(unit=44, file='Source.dat', status='replace')
    open(unit=45, file='S2.dat', status='replace')
    ! Task: Get source function from evolution_mod
    call get_hires_source_function(x_hires, k_hires, S) ! have switched order of x and k!!
    do i=1, n_hires
       write(45, fmt='(6(ES15.5E3,1x))') x_hires(i), S(i,100), S(i,1000), S(i, 1700), S(i,2500), S(i,4000)
    end do
    close(45)

    print *, '   - Source function computed'
    !print *, S(3000,:)
    ! Task: Initialize spherical Bessel functions for each l; use 5400 sampled points between
    !       z = 0 and 3500. Each function must be properly splined
    ! Hint: It may be useful for speed to store the splined objects on disk in an unformatted
    !       Fortran (= binary) file, so that these only has to be computed once. Then, if your
    !       cache file exists, read from that; if not, generate the j_l's on the fly.


    ! Initialize z from 0 to 3500
    do i=1, n_spline
       z_spline(i) = 0.d0 + (i-1.d0)*3500.d0/(n_spline-2.d0)
    end do

    ! find the Bessel functions
    filename = 'j_l.unf'
    inquire(file=filename, exist = exist)
    !if (exist) then
    !   print *, '   - Use pre-calculated Bessel functions'
    !   open(100, form='unformatted', file=filename)
    !   read(100) j_l !, j_l2
    !   close(100)

    !else
       print *, '   - Compute Bessel functions'
       do i=1, n_spline
          do l=1, l_num
            j_l(1,l) = 0.d0
            if (z_spline(i) > 2.d0) then
               call sphbes(ls(l), z_spline(i), j_l(i,l))
            end if
          end do
       end do

       ! Write bessel function to binary file to save computation time
       !open(100, form='unformatted', file=filename)
       !write(100) j_l!, j_l2
       !close(100)

       ! Need to spline the bessel function for each l
       do l=1, l_num
          call spline(z_spline, j_l(:,l), 1.d30, 1.d30, j_l2(:,l))
       end do

    !end if

    allocate(besseltest(n_hires))
    do i=1, n_hires
       besseltest(i) = splint(z_spline,j_l(:,l),j_l2(:,l), k_hires(1700)*(get_eta(0.d0)-get_eta(x_hires(i))))
       !j_l_func(17,k_hires(1700), x_hires(i))
    end do
    open(unit=101, file='besseltest.dat', status='replace')
    do i=1, n_hires
       write(101, fmt='(3(ES14.6E3,1x))') besseltest(i), S(i,1700), S(i,1700)*besseltest(i)
    end do
    close(101)


    ! =================================================
    ! Overall task: Compute the C_l's for each given l
    ! =================================================
    !do i=1,n_hires
       !write(44, fmt='(2(ES14.6E3,1x))') x_hires(i), S(1700,i)*splint(z_spline, j_l(:,17),j_l2(:,17),&
       !     (k_hires(1700)*c/H_0)*(get_eta(0.d0)-get_eta(x_hires(i))))
    !end do
    dx = (x_hires(n_hires)-x_hires(1))/(n_hires-1.d0)
    dk = (k_hires(n_hires)-k_hires(1))/(n_hires-1.d0)
    do l = 1, l_num
       ! Task: Compute the transfer function, Theta_l(k)
       do j=1, n_hires  ! k loop

          do i=1, n_hires  ! x loop

             !integrand(l,i) = S(i,j) * j_l_func(l, k_hires(j), x_hires(i))
             integrand(l,i) = S(i,j)*splint(z_spline,j_l(:,l),j_l2(:,l), k_hires(j)*(get_eta(0.d0)-get_eta(x_hires(i)) ) )

             !Theta(l,j) = Theta(l,j) + integrand(l,i) ! Integrate Theta

          end do ! x loop
          if (ls(l) == 100) then

             if (j == 1700) then
                print *, 'Write S to file for l=100, k=340 H0', k_hires(j)*c/H_0

                do i=1,n_hires
                   write(44, fmt='(2(ES14.6E3, 1x))') x_hires(i), integrand(100, i)
                end do
             end if
          end if

          do i=1, n_hires
             Theta(l,j) = Theta(l,j) + integrand(l,i) ! Eulers method
          end do
          Theta(l,j) = dx*Theta(l,j)
          ! integrate with trapeziodal method
          !call trapez(k_hires(:), integrand(l,:), Theta(l,j))

          ! Task: Integrate P(k) * (Theta_l^2 / k) over k to find un-normalized C_l's
          int_arg(l,j) = (c*k_hires(j)/H_0)**(n_s-1.d0) * (Theta(l,j)**2.d0)/(c*k_hires(j))
          int_cl(l) = int_cl(l) + int_arg(l,j)  ! integrate C_l

       end do    ! k loop

       !call trapez(k_hires(:), int_arg(l,:), int_cl(l))


       ! Task: Store C_l in an array. Optionally output to file
       cls(l) = int_cl(l)*ls(l)*(ls(l) + 1.d0)/(2.d0*pi) *dk

       print *, 'l=', ls(l)

    end do ! l loop

    ! Task: Spline C_l's found above, and output smooth C_l curve for each integer l

    call spline(ls_dp, cls, 1.d30, 1.d30, cls2)

    do i=1, n_l
       cls_hires(i) = splint(ls_dp, cls, cls2, ls_hires(i))

       write(41,fmt='(2(ES15.5E3, 1x))') ls_hires(i), cls_hires(i)
    end do
    close(41)

    ! ==================================================
    ! Write C_l and Theta to file for 6 different ls,
    ! ls: 5, 10, 17, 25, 35, 44
    ! ==================================================

    allocate(ckH0(n_hires))
    allocate(c5(n_hires), T5(n_hires))
    allocate(c10(n_hires), T10(n_hires))
    allocate(c17(n_hires), T17(n_hires))
    allocate(c25(n_hires), T25(n_hires))
    allocate(c35(n_hires), T35(n_hires))
    allocate(c44(n_hires), T44(n_hires))

    ckH0 = c*k_hires/H_0
    c5   = ls(5)*(ls(5)+1.d0)*Theta(5,:)**2/(ckH0)
    c10  = ls(10)*(ls(10)+1.d0)*Theta(10,:)**2/(ckH0)
    c17  = ls(17)*(ls(17)+1.d0)*Theta(17,:)**2/(ckH0)
    c25  = ls(25)*(ls(25)+1.d0)*Theta(25,:)**2/(ckH0)
    c35  = ls(35)*(ls(35)+1.d0)*Theta(35,:)**2/(ckH0)
    c44  = ls(44)*(ls(44)+1.d0)*Theta(44,:)**2/(ckH0)

    T5   = ls(5)*(ls(5)+1.d0)*Theta(5,:)/ckH0
    T10  = ls(10)*(ls(10)+1.d0)*Theta(10,:)/ckH0
    T17  = ls(17)*(ls(17)+1.d0)*Theta(17,:)/ckH0
    T25  = ls(25)*(ls(25)+1.d0)*Theta(25,:)/ckH0
    T35  = ls(35)*(ls(35)+1.d0)*Theta(35,:)/ckH0
    T44  = ls(44)*(ls(44)+1.d0)*Theta(44,:)/ckH0


    print *, 'Write Theta and Theta^2 to file:'
    do i=1, n_hires
       !print *, c5(i), c10(i), c44(i)
       write(42, fmt='(7(ES14.6E3))') ckH0(i), c5(i), c10(i), c17(i), c25(i), c35(i), c44(i)
       write(43, fmt='(7(ES14.6E3, 1x))') ckH0(i), T5(i), T10(i), T17(i), T25(i), T35(i), T44(i)

    end do

    close(42)
    close(43)

  end subroutine compute_cls

  subroutine trapez(x, y, integral)
    implicit none
    real(dp), dimension(:), intent(in)  :: x, y
    real(dp),               intent(out) :: integral
    integer(i4b)                        :: N, i
    real(dp)                            :: dx

    if (size(x) .ne. size(y)) then
       print *, 'x and y does not have the same shape'
       return
    end if

    N = size(x)
    integral = 0.d0
    dx = (x(N) - x(1))/N
    do i=1, N-1
       !dx = x(i+1) - x(i)
       integral = integral + dx*(y(i+1) + y(i))/2.d0
    end do

  end subroutine trapez

  function j_l_func(l, x, k)
    implicit none

    integer(i4b), intent(in) :: l
    real(dp),     intent(in) :: x, k
    real(dp)                 :: j_l_func

    j_l_func = splint(z_spline,j_l(:,l),j_l2(:,l), k*(get_eta(0.d0)-get_eta(x)))

  end function j_l_func

end module cl_mod
