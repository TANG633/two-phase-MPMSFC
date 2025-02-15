module ConstitutiveModel
  use MaterialData
  use ParticleData

  integer:: mid       ! material set id
  integer:: mtype_    ! Type of material model
  integer:: flagboundary
  real(8):: de(6)
  real(8):: vort(3)
  real(8):: young_, poisson_, tangmod_
  ! the shear modulus, the bulk modulus, the plastic modulus
  real(8):: Gmod, Kmod, PlaMod
  
  real(8):: yield0_
  real(8):: sig_y_    ! Current yield stress
  real(8):: SolidPress_con_, SolidPress_coll_
  real(8):: sm        ! Mean stress
  real(8):: sd(6)     ! deviatoric stress
  real(8):: sig(6)    ! stress components
  real(8):: sold(6)   ! deviatoric stress of step n
  real(8):: seqv_     ! Equivalent stress
  real(8):: epeff_    ! Effective plastic strain
  real(8):: depeff    ! increment of equivalent plastic strain
  real(8):: ratio     ! for hardening caculation
  real(8):: strain_rate(6)
  real(8):: Norm_DeStrain_rate
  real(8):: mp_       ! particle mass
  real(8):: den0_     ! Initial density
  real(8):: den_      ! Current density
  real(8):: fract_
  real(8):: ds
  real(8):: Jacobi_
  !*****************************
  logical:: Strain_softening = .false.
  !deviatoric plastic strain
  real(8):: d_strain_p(6)
  real(8):: edpeff_   !effective deviatoric plastic strain
  !*****************************
  real(8):: vol0_     ! Initial volume
  real(8):: vold      ! volume of step n
  real(8):: vol_      ! Current volume
  real(8):: iener     ! internal energy
  real(8):: ieinc     ! internal energy increment
  real(8):: cp        ! sound speed

  real(8):: de_f(6)
  real(8):: PorePress_
  real(8):: dPorePress_
  real(8):: BulkMod_F_
  real(8):: mp_f
  real(8):: mp_ff
  real(8):: den0_f
  real(8):: den_f_
  real(8):: Nf0_
  real(8):: Nf_
  real(8):: Kf0_
  real(8):: Kf_

  real(8):: I_inertial, J_viscous, K_visciner
  real(8):: mu
  real(8):: mu_sta_, mu_dyn_, I0_ref_

contains

  subroutine Constitution(b, p)
!------------------------------------------------------------------
!-  Purpose                                                       -
!-      Update stresses by using a constitution model             -
!-  Inputs:                                                       -
!-      b     - body index                                        -
!-      p     - particle index                                    -
!-    Outputs:                                                    -
!-        stress component                                        -
!-    Note:                                                       -
!-        sd(i) and de(i) comply the Voigt rule                   -
!------------------------------------------------------------------
    implicit none
    
    integer, intent(in):: b, p
    
    logical:: failure
    ! Pick parameters    
    mid = body_list(b)%mat
    mtype_ = mat_list(mid)%MatType   
    young_ = mat_list(mid)%Young       ! Young's Modulus
    poisson_ = mat_list(mid)%Poisson   ! Poisson Ratio
    tangmod_ = mat_list(mid)%TangMod   ! Tangential modulus
    yield0_ = mat_list(mid)%Yield0     ! Initial yield strength
    BulkMod_F_ = mat_list(mid)%BulkMod_F
    den_ = mat_list(mid)%Density
    den0_f = mat_list(mid)%Density_F
    Nf0_ = mat_list(mid)%porosity
    Kf0_ = mat_list(mid)%HydrCon
    
    flagboundary = particle_list(p)%boundaryflag
    vol0_ = particle_list(p)%Volu      ! initial volume
    vold = particle_list(p)%Vol        ! Volume at time step t
    mp_ = particle_list(p)%mass
    mp_f = particle_list(p)%mass_F
    mp_ff = particle_list(p)%mass_FF
    de = particle_list(p)%de_
    vort = particle_list(p)%vort_
    de_f = particle_list(p)%de_f_
    iener = particle_list(p)%ie        ! internal energy
    den_f_ = particle_list(p)%Den_F
    Nf_ = particle_list(p)%Nf          ! current porosity
    Kf_ = particle_list(p)%Kf
    fract_ = particle_list(p)%fract
    ds = particle_list(p)%dsp
    Jacobi_ = particle_list(p)%Jacobian

    Gmod = young_ / (2.0*(1 + poisson_))
    Kmod = young_ / (3.0*(1 - 2.0*poisson_))
    PlaMod = young_ * tangmod_ / (young_ - tangmod_)

    sm = particle_list(p)%Press        ! Mean stress
    SolidPress_con_ = particle_list(p)%SolidPress_con
    SolidPress_coll_ = particle_list(p)%SolidPress_coll
    sd(1) = particle_list(p)%SDxx      ! deviatoric stress
    sd(2) = particle_list(p)%SDyy
    sd(3) = particle_list(p)%SDzz
    sd(4) = particle_list(p)%SDyz
    sd(5) = particle_list(p)%SDzx
    sd(6) = particle_list(p)%SDxy
    sig(1) = sd(1) + sm                ! cauchy stress
    sig(2) = sd(2) + sm
    sig(3) = sd(3) + sm
    sig(4) = sd(4)
    sig(5) = sd(5)
    sig(6) = sd(6)
    sold = sd                          ! cauchy stress at time step t
    sig_y_ = particle_list(p)%sig_y    ! current yield stress
    seqv_ = particle_list(p)%seqv      ! Von Mises Stress(Equivalent stress)
    
    epeff_ = particle_list(p)%epeff    ! Effective plastic strain
    ! deviatoric plastic strain
    d_strain_p(1) = particle_list(p)%Sxx_p
    d_strain_p(2) = particle_list(p)%Syy_p
    d_strain_p(3) = particle_list(p)%Szz_p
    d_strain_p(4) = particle_list(p)%Syz_p
    d_strain_p(5) = particle_list(p)%Szx_p
    d_strain_p(6) = particle_list(p)%Sxy_p
    edpeff_ = particle_list(p)%edpeff

    !PorePress_ = particle_list(p)%PorePress
    
    call isotropic_f()                  ! incompressible fluid
    
    ! volume, porosity update
    ! det(pt%Fp)
    call Vol_upd1(p)
    ! volumetric strain
    !call Vol_upd2(p)

    ! Select material model
    select case(mtype_)

    case(1) 
       ! elas: elastic model
       call sigrot(vort, sig, sm, sd)       ! Rotate stress
       call M3DM1()
       call lieupd()

    case(2) 
       ! pla1: elastic-perfectly plastic
       call sigrot(vort, sig, sm, sd)
       call M3DM2()
       call lieupd()

    case(3) 
       ! pla2: isotropic hardening
       call sigrot(vort, sig, sm, sd)
       call M3DM3()
       call lieupd()
       
    case(4) 
       ! Drucker-Prager elastic-perfectly plastic
       call sigrot(vort, sig, sm, sd)
       call M3DM4(mat_list(mid))
       call lieupd()

    case(5) 
       ! Drucker-Prager elastic-perfectly plastic with rheology
       call sigrot(vort, sig, sm, sd)
       call M3DM5(mat_list(mid))
       call lieupd()
       
    case(6) 
       ! rheological model
       call sigrot(vort, sig, sm, sd)
       call M3DM6()
       call lieupd()
       
    case default 
       write(*, 10) mtype_
10     format(1x,'*** Stop *** material type ', i2, &
              ' has not been implemented !')
       stop

    end select

    ! Write stress result
    particle_list(p)%Press = sm
    particle_list(p)%SolidPress_con = SolidPress_con_
    particle_list(p)%SolidPress_coll = SolidPress_coll_
    particle_list(p)%SDxx = sd(1)
    particle_list(p)%SDyy = sd(2)
    particle_list(p)%SDzz = sd(3)
    particle_list(p)%SDyz = sd(4)
    particle_list(p)%SDzx = sd(5)
    particle_list(p)%SDxy = sd(6)
    
    particle_list(p)%sig_y = sig_y_
    particle_list(p)%Seqv = seqv_
    particle_list(p)%ie = iener
    
    particle_list(p)%epeff = epeff_
    particle_list(p)%Sxx_p = d_strain_p(1)
    particle_list(p)%Syy_p = d_strain_p(2)
    particle_list(p)%Szz_p = d_strain_p(3)
    particle_list(p)%Syz_p = d_strain_p(4)
    particle_list(p)%Szx_p = d_strain_p(5)
    particle_list(p)%Sxy_p = d_strain_p(6)
    particle_list(p)%edpeff = edpeff_
 
    !particle_list(p)%PorePress = PorePress_
    particle_list(p)%dPorePress = dPorePress_

    particle_list(p)%Nf = Nf_
    particle_list(p)%Kf = Kf_
    particle_list(p)%fract = fract_
    particle_list(p)%VOL = vol_
    particle_list(p)%dsp = ds
    particle_list(p)%mass_F = mp_f
    particle_list(p)%mass_FF = mp_ff

  end subroutine Constitution

  subroutine sigrot(vort, sig, sm, sd)
!------------------------------------------------------------------
!-  Purpose                                                       -
!-      Rotate stresses, and then update the mean stress sm       -
!-      and deviate stresses sd                                   -
!-  Input                                                         -
!-      vort - Vorticity increments (W32, W13, W21)*DT            -
!-      sig  - Caucy stresses at time step t                      -
!-             (S11, S22, S33, S23, S13, S12)                     -
!-  Output                                                        -
!-      sig  - Rotated caucy stresses                             -
!-      sd   - Rotated deviate stress                             -
!-      sm   - Rotated mean stress                                -
!------------------------------------------------------------------
    implicit none
    
    real(8), intent(in):: vort(3)      ! Vorticity increment
    real(8), intent(inout):: sig(6)    ! Caucy stresses
    real(8), intent(out):: sm, sd(6)   ! Mean stress, deviatoric stress

    real(8):: rot(6), q(3)

    q(1) = 2.0*vort(1)*sig(4)
    q(2) = 2.0*vort(2)*sig(5)
    q(3) = 2.0*vort(3)*sig(6)

    rot(1) = - q(3) + q(2)
    rot(2) = + q(3) - q(1)
    rot(3) = - q(2) + q(1)
    rot(4) = vort(1)*(sig(2)-sig(3)) + vort(3)*sig(5) - vort(2)*sig(6)
    rot(5) = vort(2)*(sig(3)-sig(1)) + vort(1)*sig(6) - vort(3)*sig(4)
    rot(6) = vort(3)*(sig(1)-sig(2)) + vort(2)*sig(4) - vort(1)*sig(5)

    sig = sig + rot * DT
    sm = (sig(1)+sig(2)+sig(3))/3d0     !Rotated mean stress
    !Rotated deviate stress
    sd(1) = sig(1) - sm  
    sd(2) = sig(2) - sm
    sd(3) = sig(3) - sm
    sd(4) = sig(4)
    sd(5) = sig(5)
    sd(6) = sig(6)

  end subroutine sigrot
  
  subroutine elastic_devi()
!------------------------------------------------------------------
!-  Purpose                                                       -
!-      Update devitoric stress by elastic relation               -
!-  Inputs                                                        -
!-      de   - strain increment                                   -
!-             (D11, D22, D33, 2D23, 2D13, 2D12)*DT               -
!-  Outputs                                                       -
!-    sd     - devitoric stress component                         -
!------------------------------------------------------------------
    implicit none
    real(8):: dem

    dem = (de(1) + de(2) + de(3)) / 3.0 

    sd(1) = sd(1) + 2.0*Gmod*(de(1)-dem) * DT
    sd(2) = sd(2) + 2.0*Gmod*(de(2)-dem) * DT
    sd(3) = sd(3) + 2.0*Gmod*(de(3)-dem) * DT
    sd(4) = sd(4) + 2.0*Gmod*de(4) * DT
    sd(5) = sd(5) + 2.0*Gmod*de(5) * DT
    sd(6) = sd(6) + 2.0*Gmod*de(6) * DT

  end subroutine elastic_devi

  subroutine elastic_p()
!------------------------------------------------------------------
!-  Purpose                                                       -
!-      Update pressure by elastic relation                       -
!-  Inputs                                                        -
!       dinc - strain increment                                   -
!   Outputs                                                       -
!       sm   - mean stress (pressure) positive for tension        -
!------------------------------------------------------------------
    implicit none
    real(8):: dem, dsm

    ! the bulk strain rate
    dem = de(1) + de(2) + de(3)
    dsm = Kmod*dem * DT
    sm = sm + dsm

  end subroutine elastic_p

  subroutine stress_shear()
!------------------------------------------------------------------
!-  Purpose                                                       -
!-         Update solid shear stress                              -
!-  Inputs                                                        -
!          sm   - mean stress (solid pressure)                    -
!          friction                                               -
!   Outputs                                                       -
!-         sd   - devitoric stress component                      -
!------------------------------------------------------------------
    implicit none
    real(8):: dem
    
    dem = (de(1) + de(2) + de(3)) / 3.0 
    
    sd(1) = mu*sm*(de(1)-dem) / Norm_DeStrain_rate
    sd(2) = mu*sm*(de(2)-dem) / Norm_DeStrain_rate
    sd(3) = mu*sm*(de(3)-dem) / Norm_DeStrain_rate
    sd(4) = mu*sm*de(4) / Norm_DeStrain_rate
    sd(5) = mu*sm*de(5) / Norm_DeStrain_rate
    sd(6) = mu*sm*de(6) / Norm_DeStrain_rate
    
  end subroutine stress_shear

  subroutine isotropic_f()
!------------------------------------------------------------------
!-  Purpose                                                       -
!-                isotropic pore pressure                         -
!------------------------------------------------------------------
    implicit none

    dPorePress_ = (BulkMod_F_ / Nf_)*(fract_*(de(1) + de(2) + de(3)) + &
                  Nf_*(de_f(1) + de_f(2) + de_f(3))) * DT

    !PorePress_ = PorePress_ - dPorePress_

  end subroutine isotropic_f

  real function SecDeStreInva()
!------------------------------------------------------------------
!-  Purpose                                                       -
!-      Calculate the second invariant of deviatoric stress       -
!-  Input                                                         -
!-      sd    - the deviatoric stress components                  -
!-  Return                                                        -
!-      J2    - the equivalent stress                             -
!------------------------------------------------------------------
    implicit none

    real(8):: J2  ! the second deviatoric stress invariant

    J2 = 0.5d0*(sd(1)**2 + sd(2)**2 + sd(3)**2) + sd(4)**2 + &
         sd(5)**2 + sd(6)**2
    SecDeStreInva = J2
    
    return

  end function SecDeStreInva
  
  real function SecDeStrainInva()
!------------------------------------------------------------------
!-  Purpose                                                       -
!-      Calculate the second invariant of deviatoric strain       -
!-  Input                                                         -
!-      de     - the deviatoric strain components                 -
!-  Return                                                        -
!-      S2     - the second invariant of deviatoric stress        -
!------------------------------------------------------------------
    implicit none

    real(8):: dem, strain_rate_dev(6)
    real(8):: S2

    dem = (de(1) + de(2) + de(3)) / 3.0
    
    strain_rate_dev(1) = (de(1) - dem)
    strain_rate_dev(2) = (de(2) - dem)
    strain_rate_dev(3) = (de(3) - dem)
    strain_rate_dev(4) = de(4)
    strain_rate_dev(5) = de(5)
    strain_rate_dev(6) = de(6)

    S2 = 0.5*(strain_rate_dev(1)**2 + strain_rate_dev(2)**2 + strain_rate_dev(3)**2) + &
         strain_rate_dev(4)**2 + strain_rate_dev(5)**2 + strain_rate_dev(6)**2
    SecDeStrainInva = S2

    return
    
  end function SecDeStrainInva
  
  subroutine DimenNum_upd()
!------------------------------------------------------------------
!-  Purpose                                                       -
!-      Calculate the dimensionless number for different regime   -
!------------------------------------------------------------------
    implicit none

    Norm_DeStrain_rate = sqrt(SecDeStrainInva())
    
    ! for ineritial regime
    ! ineritial number
    I_inertial = Norm_DeStrain_rate*ds / sqrt(sm / den_)
    mu = mu_sta_ + (mu_dyn_ - mu_sta_)/((I0_ref_ / I_inertial) + 1.0)

  end subroutine DimenNum_upd

  subroutine Vol_upd1(p)
!------------------------------------------------------------------
!-  Purpose                                                       -
!-   Update particle volume using the Jacobian of                 -
!-   the deformation gradient                                     -
!-   Update porosity and solid volume fraction                    -
!------------------------------------------------------------------
    implicit none
    
    integer, intent(in):: p
    
    ! volume, porosity update
    ! det(pt%Fp)
    vol_ = Jacobi_ * vol0_             ! Current volume
    if (vol_ .lt. 0) then
        write(*,*) '=== warning: negative volume at particle', p, vol_
    end if    
    Nf_ = 1.0 - (1.0 - Nf0_) / Jacobi_
    if (Nf_<0.0 .OR. Nf_>1.0) then
        write(*,*) '=== warning: wrong porosity at solid particle', p, Nf_
    end if
    fract_ = 1.0 - Nf_

    ds = (fract_*vol_)**(1.0/3.0)
    mp_f = Nf_*den_f_*vol_
    mp_ff = den_f_*vol_
    !Kf_ = Kf0_*(((1-Nf0_)/(1-Nf_))**2)
    !Kf_ = 2.0*((Nf_**3)/((1-Nf_)**2))

  end subroutine Vol_upd1
  
  subroutine Vol_upd2(p)
!------------------------------------------------------------------
!-  Purpose                                                       -
!-   Update particle volume using volumetric strain               -
!-   Update porosity and solid volume fraction                    -
!------------------------------------------------------------------
    implicit none
    
    integer, intent(in):: p

    ! volumetric strain
    vol_ = vold * (1.0 + (de(1) + de(2) + de(3))*DT)    ! Current volume of solid
    if (vol_ .lt. 0) then
        write(*,*) '=== warning: negative volume at particle', p, vol_
    end if 
    Nf_ = Nf_ + fract_*(de(1) + de(2) + de(3))*DT
    if (Nf_ .lt. 0.0 .OR. Nf_ .gt. 1.0) then
        write(*,*) '=== warning: wrong porosity at solid particle', p, Nf_
    end if
    fract_ = 1.0 - Nf_
    
    ds = (fract_*vol_)**(1.0/3.0)
    mp_f = Nf_*den_f_*vol_
    mp_ff = den_f_*vol_
    !Kf_ = Kf0_*(((1-Nf0_)/(1-Nf_))**2)
    !Kf_ = 2.0*((Nf_**3)/((1-Nf_)**2))
    
  end subroutine Vol_upd2

  subroutine lieupd()
!------------------------------------------------------------------
!-  Purpose                                                       -
!-      Update energy for material models                         -
!------------------------------------------------------------------
    implicit none

    real:: vavg

    vavg = vol_ + vold
    iener = iener + 0.25d0 * ieinc * vavg

  end subroutine lieupd
  
  subroutine M3DM1()
!------------------------------------------------------------------
!-  Purpose                                                       -
!-      Elastic material model                                    -
!------------------------------------------------------------------
    implicit none

    ieinc = (sig(1)*de(1) + sig(2)*de(2) + sig(3)*de(3) + &
            sig(4)*de(4) + sig(5)*de(5) + sig(6)*de(6))*DT

    call elastic_devi()
    call elastic_p()

    seqv_ = sqrt(3.0*SecDeStreInva())

    ieinc = ieinc + ((sd(1)+sm)*de(1) + (sd(2)+sm)*de(2) + (sd(3)+sm)*de(3) + &
            sd(4)*de(4) + sd(5)*de(5) + sd(6)*de(6))*DT

  end subroutine M3DM1

  subroutine M3DM2()
!------------------------------------------------------------------
!-  Purpose                                                       -
!-      Elastic-perfectly plastic material model                  -
!------------------------------------------------------------------
    implicit none

    ieinc = (sig(1)*de(1) + sig(2)*de(2) + sig(3)*de(3) + &
            sig(4)*de(4) + sig(5)*de(5) + sig(6)*de(6))*DT

    call elastic_devi()
    call elastic_p()

    seqv_ = sqrt(3.0*SecDeStreInva())

    if (seqv_ .GT. yield0_) then
       depeff = (seqv_ - sig_y_) / (3.0*Gmod)       ! increment of equivalent plastic strain
       epeff_ = epeff_ + depeff                     ! Effective plastic strain

       ratio = yield0_/seqv_

       sd(1) = sd(1)*ratio
       sd(2) = sd(2)*ratio
       sd(3) = sd(3)*ratio
       sd(4) = sd(4)*ratio
       sd(5) = sd(5)*ratio
       sd(6) = sd(6)*ratio

       seqv_ = seqv_*ratio
    end if

    ieinc = ieinc + ((sd(1)+sm)*de(1) + (sd(2)+sm)*de(2) + (sd(3)+sm)*de(3) + &
            sd(4)*de(4) + sd(5)*de(5) + sd(6)*de(6))*DT

  end subroutine M3DM2

  subroutine M3DM3()
!------------------------------------------------------------------
!-  Purpose                                                       -
!-      Isotropic hardening plastic material model                -
!------------------------------------------------------------------
    implicit none

    ieinc = (sig(1)*de(1) + sig(2)*de(2) + sig(3)*de(3) + &
            sig(4)*de(4) + sig(5)*de(5) + sig(6)*de(6))*DT

    call elastic_devi()
    call elastic_p()
    
    seqv_ = sqrt(3.0*SecDeStreInva())

    if (seqv_ .GT. sig_y_) then
       depeff = (seqv_ - sig_y_) / (3.0*Gmod + PlaMod)
       epeff_ = epeff_ + depeff

       sig_y_ = sig_y_ + PlaMod*depeff
       ratio = sig_y_/seqv_

       sd(1) = sd(1)*ratio
       sd(2) = sd(2)*ratio
       sd(3) = sd(3)*ratio
       sd(4) = sd(4)*ratio
       sd(5) = sd(5)*ratio
       sd(6) = sd(6)*ratio

       seqv_ = seqv_*ratio                
    end if

    ieinc = ieinc + ((sd(1)+sm)*de(1) + (sd(2)+sm)*de(2) + (sd(3)+sm)*de(3) + &
            sd(4)*de(4) + sd(5)*de(5) + sd(6)*de(6))*DT

  end subroutine M3DM3

  subroutine M3DM4(mat)
!------------------------------------------------------------------
!-  Purpose                                                       -
!-      perfect elastic-plastic Drucker-Prager model for soil     -
!------------------------------------------------------------------
    implicit none
    
    type(material), intent(in):: mat
    real(8):: cohen_, fai_, psi_
    real(8):: qfai, kfai, qpsi, tenf_, tenf_max 
    real(8):: J2, Tau, dpFi, dpsig, dlamd, newTau
    real(8):: ratio
    integer:: iplas, inum

    if (Strain_softening == .true.) then
        fai_ = mat%fai_res + (mat%fai_peak - mat%fai_res)*exp(-mat%Para_yita*edpeff_)
        cohen_= mat%cohen_res + (mat%cohen_peak - mat%cohen_res)*exp(-mat%Para_yita*edpeff_)
        psi_= mat%psi_int + (mat%psi_crit - mat%psi_int)*exp(-mat%Para_yita*edpeff_)
        
        qfai = 6.0*sind(fai_)/(sqrt(3.0)*(3.0-sind(fai_)))
        kfai = 6.0*cohen_*cosd(fai_)/(sqrt(3.0)*(3.0-sind(fai_)))
        qpsi = 6.0*sind(psi_)/(sqrt(3.0)*(3.0-sind(psi_)))
    else
        fai_ = mat%fai
        cohen_ = mat%cohen
        psi_ = mat%psi
        
        qfai = 6.0*sind(fai_)/(sqrt(3.0)*(3.0-sind(fai_)))
        kfai = 6.0*cohen_*cosd(fai_)/(sqrt(3.0)*(3.0-sind(fai_)))
        qpsi = 6.0*sind(psi_)/(sqrt(3.0)*(3.0-sind(psi_)))
    end if
    ! internal energy incremental
    ieinc = (sig(1)*de(1) + sig(2)*de(2) + sig(3)*de(3) + &
            sig(4)*de(4) + sig(5)*de(5) + sig(6)*de(6))*DT
    ! Give the tension stress value
    if (qfai == 0.0) then
        tenf_ = 0.0
    else
        tenf_max = kfai / qfai
        tenf_ = min(mat%ten_f, tenf_max)
    end if
    ! --- trial elastic stresses ---
    iplas = 0               ! elastic calculation
    call elastic_devi()     ! Update devitoric stress by elastic relation
    call elastic_p()        ! Update pressure by elastic relation

    J2 = SecDeStreInva()
    Tau = sqrt(J2)
    seqv_ = sqrt(3.0*J2)
    
    dpFi = Tau + qfai*sm - kfai     ! shear yield surface
    dpsig = sm - tenf_              ! tensile yield surface
    if (dpsig < 0.0) then
        if (dpFi > 0.0) then
            iplas = 1                             ! shear plastic flow
            ! plastic flow coefficient
            dlamd = dpFi/(Gmod + Kmod*qfai*qpsi)
            ! correct spherical stress
            sm = sm - Kmod*qpsi*dlamd
            ! correct shear stress
            newTau = kfai - qfai*sm
            ratio = newTau / Tau
            ! correct deviatoric stress
            sd(1) = sd(1)*ratio
            sd(2) = sd(2)*ratio
            sd(3) = sd(3)*ratio
            sd(4) = sd(4)*ratio
            sd(5) = sd(5)*ratio
            sd(6) = sd(6)*ratio
            
            seqv_ = seqv_*ratio                  ! correct the mises stress
            
            ! calculate the plastic strain
            depeff = dlamd*sqrt(1.0/3.0 + (2.0/9.0)*(qpsi**2))
            epeff_ = epeff_ + depeff
        end if
    else   ! (dpsig >= 0.0)
        iplas = 2                                ! tension plastic flow
        ! plastic flow coefficient
        dlamd = dpsig / Kmod
        sm = tenf_
        ! calculate the plastic strain
        depeff = dlamd*sqrt(2.0)/3.0
        epeff_ = epeff_ + depeff
    end if

    if (iplas /= 0 .and. newTau > epsilon(1.0d-15)) then
        do inum = 1, 6
            d_strain_p(inum) = d_strain_p(inum) + (dlamd * sd(inum) / (2 * newTau))
        end do
        edpeff_ = sqrt(((d_strain_p(1)**2 + d_strain_p(2)**2 + d_strain_p(3)**2) + &
                  2 * d_strain_p(4)**2 + 2*d_strain_p(5)**2 + 2*d_strain_p(6)**2)*(2.0/3.0))
    end if
    
    ieinc = ieinc + ((sd(1)+sm)*de(1) + (sd(2)+sm)*de(2) + (sd(3)+sm)*de(3) + &
            sd(4)*de(4) + sd(5)*de(5) + sd(6)*de(6))*DT
  end subroutine M3DM4

  subroutine M3DM5(mat)
!---------------------------------------------------------------------
!-  Purpose                                                          -
!-   perfect elastic-plastic Drucker-Prager model with rtheology law -
!---------------------------------------------------------------------
    implicit none
    type(material), intent(in):: mat
    ! internal energy increment
    ieinc = sig(1)*de(1) + sig(2)*de(2) + sig(3)*de(3) + &
            sig(4)*de(4) + sig(5)*de(5) + sig(6)*de(6)

    ieinc = ieinc + ((sd(1)+sm)*de(1) + (sd(2)+sm)*de(2) + (sd(3)+sm)*de(3) + &
            sd(4)*de(4) + sd(5)*de(5) + sd(6)*de(6))*DT
  end subroutine M3DM5
  
  subroutine M3DM6()
!-----------------------------------------------------
!-                 Rheology model                    -
!-----------------------------------------------------
    implicit none
    ! internal energy increment
    ieinc = (sig(1)*de(1) + sig(2)*de(2) + sig(3)*de(3) + &
            sig(4)*de(4) + sig(5)*de(5) + sig(6)*de(6))*DT

    ieinc = ieinc + ((sd(1)+sm)*de(1) + (sd(2)+sm)*de(2) + (sd(3)+sm)*de(3) + &
            sd(4)*de(4) + sd(5)*de(5) + sd(6)*de(6))*DT
  end subroutine M3DM6

end module ConstitutiveModel