! ------------------------------------------------------------------
! -                                                                -
! -  Solution and particles information                            -
! -                                                                -
! ------------------------------------------------------------------
module ParticleData

  type Particle
     real(8):: XX(3)     ! particle position at time step t+1
     real(8):: Xp(3)     ! particle position at time step t
     real(8):: Disp_p(3) ! displacement
     real(8):: H         ! depth of the particle below surface
     real(8):: Hw        ! depth of the particle below the groundwater table
     real(8):: Cp_apic(3, 3)        ! coefficient of APIC
     real(8):: Volu      ! initial volume
     real(8):: Vol       ! current volume
     !for solid    
     real(8):: Mass      ! particle mass
     real(8):: FXp(3)    ! load
     real(8):: VXp(3)    ! velocity
     real(8):: sig_y     ! yield stress
     real(8):: Seqv      ! Von Mises Stress(Equivalent stress)
     real(8):: Press     ! Mean stress (Tension is positive)
     real(8):: SolidPress_con, SolidPress_coll
     real(8):: SDxx, SDyy, SDzz, SDyz, SDzx, SDxy           ! deviatoric stress
     real(8):: epeff     ! effective plastic strain
     real(8):: Sxx_p, Syy_p, Szz_p, Syz_p, Szx_p, Sxy_p     ! deviatoric plastic strain
     real(8):: edpeff                                       ! effective deviatoric plastic strain
     real(8):: dsp
     real(8):: fract
     real(8):: Jacobian
     real(8):: Fp(3,3)   ! the gradient deformation
     real(8):: Jacobi
     real(8):: Fp_F(3,3) ! the modified gradient deformation
     real(8):: Jacobi_F
     ! for fluid
     real(8):: mass_F
     real(8):: mass_FF
     real(8):: Den_F
     real(8):: VXp_F(3)
     real(8):: PXp_F(3)   ! load for fluid
     real(8):: PorePress  ! pore pressure
     real(8):: dPorePress
     real(8):: Nf         ! porosity
     real(8):: Kf         ! hydraulic conductivity

     real(8):: de_(6)
     real(8):: vort_(3)
     real(8):: de_f_(6)
     real(8):: grad_vel(3)

     real(8):: Qp
     real(8):: ap
     
     logical:: SkipThis   ! for plot less
     logical:: failure    ! failure
     integer:: icell(3)   ! cell number
     ! basenode number
     integer:: basenode1(3)
     integer:: boundaryflag
     real(8):: ie         ! internal energy
     real(8):: cp         ! sound speed
     real(8):: DMG        ! damage
     
  end type Particle

  type Body
     integer:: mat       ! material set number (1 ~ nb_mat)
     integer:: comID     ! component set number (1 ~ nb_component) 
     real(8):: Gravp(3)  ! gravity
     integer:: par_begin ! the beginning of particle of body
     integer:: par_end   ! the end of particle of body
  end type Body
  
  type(Particle), target, allocatable:: particle_list(:)    ! particles list
  type(Body), target, allocatable:: body_list(:)            ! bodies list
  
  integer:: nb_particle = 0    ! number of particles
  integer:: nb_component = 1   ! number of components for contact simulating
  integer:: nb_body = 0        ! number of bodies

  character(256):: Title
  integer:: iTitle(64)
  equivalence (Title, iTitle)

  logical:: GIMP = .false.      ! use GIMP method?
  logical:: CUBI = .false.      ! use Cubic B spline shape function?
  ! initial shape function setting for linear shape function
  integer:: InfluInit = 0
  integer:: InfluEnd = 1
  integer:: Offset = 1
  logical:: contact = .false.   ! use contact method?
  
  integer(8):: istep = 0        ! Current time step
  
  real(8):: DT = 0.0            ! time step interval
  real(8):: CurrentTime = 0.0   ! Time for current step
  real(8):: loadtime = 0.1
  real(8):: loadrate(3)
  real(8):: EndTime = 0.0       ! Time for end of solution
  real(8):: DTScale = 0.5       ! time step size factor (<= 1.0)

  real(8):: EngInternal = 0.0   ! Internal energy
  real(8):: EngKinetic = 0.0    ! Kinetic energy
  real(8):: Momentum(3) 
  real(8):: Mombody1(3)
  real(8):: Mombody2(3)

contains

  subroutine InitParticle()
!------------------------------------------------------------------
!    purpose: initialize particle_list                            -
!    used after particle_list space allocated                     -
!------------------------------------------------------------------
    implicit none
    integer :: i, j

    particle_list%SolidPress_con = 0.0d0
    particle_list%SolidPress_coll = 0.0d0
    particle_list%Press = 0.0d0     ! mean stress
    particle_list%seqv = 0.0d0      ! Mises stress
    particle_list%SDxx = 0.0d0      ! deviatoric stress
    particle_list%SDyy = 0.0d0
    particle_list%SDzz = 0.0d0
    particle_list%SDyz = 0.0d0
    particle_list%SDzx = 0.0d0
    particle_list%SDxy = 0.0d0
    
    particle_list%epeff = 0.0d0     ! effective plastic strain
    particle_list%Sxx_p = 0.0d0
    particle_list%Syy_p = 0.0d0
    particle_list%Szz_p = 0.0d0
    particle_list%Syz_p = 0.0d0
    particle_list%Szx_p = 0.0d0
    particle_list%Sxy_p = 0.0d0
    particle_list%edpeff = 0.0d0
    
    particle_list%dsp = 0.0d0
    particle_list%fract = 0.0d0
    particle_list%PorePress = 0.0d0
    particle_list%dPorePress = 0.0d0
    particle_list%Nf = 0.0d0
    particle_list%Kf = 0.0d0

    particle_list%Qp = 1.0d0
    particle_list%ap = 0.0d0
    
    !ËÙ¶È ºÉÔØ Î»ÒÆ
    do i = 1, 3
        particle_list%VXp(i) = 0.0d0
        particle_list%FXp(i) = 0.0d0
        particle_list%Disp_p(i) = 0.0d0
        particle_list%VXp_F(i) = 0.0d0
        particle_list%PXp_F(i) = 0.0d0
        particle_list%icell(i) = 0
        particle_list%basenode1(i) = 0
    end do

    do i = 1, 3
        do j = 1, 3
            particle_list%Cp_apic(i, j) = 0.0d0
            if (i == j) then
                particle_list%Fp(i, j) = 1.0d0
            else
                particle_list%Fp(i, j) = 0.0d0
            end if
        end do
    end do

    particle_list%boundaryflag = 0
    particle_list%SkipThis = .false.
    particle_list%failure = .false.
    particle_list%ie = 0.0d0
    particle_list%cp = 0.0d0
    particle_list%DMG = 0.0d0

  end subroutine InitParticle

  subroutine InitBody
!------------------------------------------------------------------
!    purpose: initialize body_list                                -
!    used after body_list space allocated                         -
!------------------------------------------------------------------
	  implicit none
	  integer :: i

	  do i = 1, 3
        body_list%Gravp(i) = 0.0d0
	  end do

  end subroutine InitBody

  subroutine calcEnergy()
!-------------------------------------------------------------------
!-   purpose: calculate kinematic energy  and internal energy      -
!----------------------------------------------------------------- -
    implicit none
    integer:: p    ! loop counter
    real(8):: delta_EngKinetic, iener, delta_Momentum(3)
    type(Particle), POINTER:: pt
    
    EngKinetic  = 0.0      ! initial Kinetic energy 
    EngInternal = 0.0      ! initial internal energy
    Momentum    = 0.0      ! set initial momentum equal to  zero
  
    !$OMP PARALLEL DO PRIVATE(pt,delta_EngKinetic,iener,delta_Momentum) &
    !$OMP REDUCTION(+:EngKinetic, EngInternal, Momentum) SCHEDULE(STATIC)
    do p = 1, nb_particle
        pt => particle_list(p)
        delta_EngKinetic = DOT_PRODUCT(pt%VXp, pt%VXp)*pt%mass*0.5d0
        iener = pt%ie
        delta_Momentum = pt%mass*pt%VXp
      
        EngKinetic = EngKinetic + delta_EngKinetic
        EngInternal = EngInternal + iener
        Momentum = Momentum + delta_Momentum
    end do
    !$OMP END PARALLEL DO

  end subroutine calcEnergy

end module ParticleData
