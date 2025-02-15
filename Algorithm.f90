! -----------------------------------------------------------------
! -                                                               -
! -  Calculation procedures                                       -
! -                                                               -
! -----------------------------------------------------------------
module update_step

contains

  subroutine GravLoad()
  !------------------------------------------------------------------
  !-  Purpose                                                       -
  !-    linear gravity loading                                      -
  !------------------------------------------------------------------
    use ParticleData
    implicit none
  
    integer:: b
  
    type(body), POINTER:: bd

    do b = 1, nb_body
        bd =>body_list(b)
        bd%Gravp = loadrate*CurrentTime
    end do

  end subroutine GravLoad

  subroutine GridInformationInitial()
  !-------------------------------------------------------------------
  !-  Purpose                                                        -
  !-      1. The variables of particle is mapped to the grid node    -
  !-------------------------------------------------------------------
    use ParticleData
    use GridData
    use FFI, only: iomsg
    implicit none

    integer:: i, b, p, parBegin, parEnd     ! loop counter
    integer:: ix, iy, iz, comID_ = 1
    real(8):: mp, mp_f, mp_ff, vp(3), vp_f(3)
    real(8):: pos(3), Cp(3, 3), dis_det(3)
    real(8), dimension(2, 3, 4):: shi
    real(8):: SHPi

    type(Particle), POINTER:: pt
    type(GridNode), POINTER:: gn
    ! Calculate the grid nodal masses, moemntum only
    ! Reset Grid data
    do i = 1, 3
        node_list%PXg(i) = 0.0d0
        node_list%PXg_F(i) = 0.0d0
        node_list%deltx(i) = 0.0d0
    end do
    node_list%Mg = 0.0d0
    node_list%Mg_F = 0.0d0
    node_list%Mg_FF = 0.0d0
    node_list%interp = .false.
    cell_list%pnum = 0
    cell_list%index = 0

    do b = 1, nb_body
        parBegin = body_list(b)%par_begin
        parEnd = body_list(b)%par_End
        if(contact) comID_ = body_list(b)%comID   ! Get comID from body
        !$OMP PARALLEL DO SCHEDULE(DYNAMIC,3) &
        !$OMP PRIVATE(pt,mp,mp_f,mp_ff,vp,vp_f,Cp,shi) &
        !$OMP PRIVATE(ix,iy,iz,SHPi,gn,dis_det)
        do p = parBegin, parEnd                   ! Loop over all particles of the body
            pt => particle_list(p)
            pt%icell = InWhichCell(pt%Xp, p)      ! Coordinate of the cell in which the particle located
            ! Particle p is out of the computational region
            if (any(pt%icell < 0)) cycle
            mp = pt%Mass; vp = pt%VXp
            mp_f = pt%Mass_F; mp_ff = pt%Mass_FF; vp_f = pt%VXp_F
            pos = pt%Xp; Cp = pt%Cp_apic
            ! find the base node
            pt%basenode1 = BaseNode(pt%Xp)
            ! Calculate the shape functions and their derivatives
            if (CUBI) then
                shi = NShape_Cubic(comID_, pt%basenode1, pt%Xp)
            else if (GIMP) then
                shi = NShape_GIMP(comID_, pt%basenode1, pt%Xp)
            else
                shi = NShape(comID_, pt%basenode1, pt%Xp)
            end if
            ! loop all the neighbor nodes
            do iz = InfluInit, InfluEnd
                if (pt%basenode1(3)+iz < 1 .OR. pt%basenode1(3)+iz > NGz) cycle
                do iy = InfluInit, InfluEnd
                    if (pt%basenode1(2)+iy < 1 .OR. pt%basenode1(2)+iy > NGy) cycle
                    do ix = InfluInit, InfluEnd
                        if (pt%basenode1(1)+ix < 1 .OR. pt%basenode1(1)+ix > NGx) cycle
                      
                        SHPi = shi(1, 1, ix+Offset)*shi(1, 2, iy+Offset)*shi(1, 3, iz+Offset)
                        if (SHPi <= 0.0) cycle
                        gn => node_list(comID_, pt%basenode1(1)+ix, pt%basenode1(2)+iy, pt%basenode1(3)+iz)
                        gn%interp = .true.
                        dis_det = gn%Xg - pt%Xp
                        !$OMP CRITICAL
                        gn%Mg = gn%Mg + SHPi*mp               ! the nodal mass
                        gn%Mg_F = gn%Mg_F + SHPi*mp_f
                        gn%Mg_FF = gn%Mg_FF + SHPi*mp_ff
                        gn%PXg = gn%PXg + SHPi*mp*vp          ! the nodal momentum
                        ! APIC
                        !gn%PXg = gn%PXg + SHPi*mp*(vp + MATMUL(Cp, dis_det))
                        gn%PXg_F = gn%PXg_F + SHPi*mp_f*vp_f
                        if (contact) then
                            gn%deltx = gn%deltx + SHPi*mp*pos ! the nodal position
                        end if
                        !$OMP END CRITICAL
                    end do    !ix
                end do    !iy
            end do    !iz
        end do    !p
        !$OMP END PARALLEL DO
    end do    !b

  end subroutine GridInformationInitial

  subroutine GridNodalForce()
  !------------------------------------------------------------------
  !-  Purpose                                                       -
  !-      1. calculate the background grid nodal force              -
  !------------------------------------------------------------------
    use ParticleData
    use GridData
    use MaterialData
    use FFI, only: iomsg
    implicit none

    integer:: i, b, p, parBegin, parEnd     ! loop counter
    integer:: ix, iy, iz, comID_ = 1
    real(8):: sxx, syy, szz, syz, szx, sxy
    real(8):: fx(3), fx_int(3), fx_ext(3)
    real(8):: PorePress_
    real(8):: fx_f(3), fxf_int(3), fxf_ext(3), fx_drag(3)
    real(8):: Drag, mp, vol_
    real(8):: vs(3), vf(3)
  
    real(8), dimension(2, 3, 4):: shi
    real(8):: SHPi, DNDXi, DNDYi, DNDZi
  
    type(Particle), POINTER:: pt
    type(GridNode), POINTER:: gn
    type(ContactGridNodeProperty), POINTER:: CP
    
    ! Calculate the grid nodal forces only
    ! Reset nodal forces
    do i = 1, 3
        node_list%FXg(i) = 0.0d0
        node_list%FXg_F(i) = 0.0d0
    end do

    if (contact) then
       CP_list%ndir(1) = 0.0d0    
       CP_list%ndir(2) = 0.0d0
       CP_list%ndir(3) = 0.0d0
    end if
  
    do b = 1, nb_body             ! Loop over all bodies
        parBegin = body_list(b)%par_begin
        parEnd = body_list(b)%par_End
      
        if(contact) comID_ = body_list(b)%comID   ! Get comID from body
        !$OMP PARALLEL DO SCHEDULE(DYNAMIC,3) &
        !$OMP PRIVATE(pt,mp,vol_,sxx,syy,szz,syz,szx,sxy,PorePress_,Drag,shi) &
        !$OMP PRIVATE(ix,iy,iz,SHPi,DNDXi,DNDYi,DNDZi,gn,CP) &
        !$OMP PRIVATE(vs,vf,fx,fx_int,fx_ext,fx_f,fxf_int,fxf_ext,fx_drag)
        do p = parBegin, parEnd   ! Loop over all particles
            pt => particle_list(p)
            ! old particle position
            ! Particle p is out of the computational region
            if (any(pt%icell < 0)) cycle
            mp = pt%Mass; vol_ = pt%Vol
            ! solid stress
            sxx = pt%SDxx + pt%Press
            syy = pt%SDyy + pt%Press
            szz = pt%SDzz + pt%Press
            syz = pt%SDyz
            szx = pt%SDzx
            sxy = pt%SDxy
            ! fluid stress
            PorePress_ = pt%PorePress
            ! inter-phase force
            Drag = pt%Mass_F*9.8 / pt%Kf

            ! External forces
            fx_f = - pt%PXp_F
            fx = pt%FXp - pt%PXp_F
          
            fx_f = fx_f + pt%Mass_FF * (body_list(b)%Gravp)
            fx = fx + (pt%Mass + pt%Mass_F) * (body_list(b)%Gravp)
            ! Calculate the shape functions and their derivatives
            if (CUBI) then
                shi = NShape_Cubic(comID_, pt%basenode1, pt%Xp)
            else if (GIMP) then
                shi = NShape_GIMP(comID_, pt%basenode1, pt%Xp)
            else         
                shi = NShape(comID_, pt%basenode1, pt%Xp)
            end if
            ! loop all the neighbor nodes
            do iz = InfluInit, InfluEnd
                if (pt%basenode1(3)+iz < 1 .OR. pt%basenode1(3)+iz > NGz) cycle
                do iy = InfluInit, InfluEnd
                    if (pt%basenode1(2)+iy < 1 .OR. pt%basenode1(2)+iy > NGy) cycle
                    do ix = InfluInit, InfluEnd
                        if (pt%basenode1(1)+ix < 1 .OR. pt%basenode1(1)+ix > NGx) cycle
                      
                        gn => node_list(comID_, pt%basenode1(1)+ix, pt%basenode1(2)+iy, pt%basenode1(3)+iz)
                        if (.NOT.gn%interp) cycle
                        SHPi = shi(1, 1, ix+Offset)*shi(1, 2, iy+Offset)*shi(1, 3, iz+Offset)
                        DNDXi = shi(2, 1, ix+Offset)*shi(1, 2, iy+Offset)*shi(1, 3, iz+Offset)
                        DNDYi = shi(1, 1, ix+Offset)*shi(2, 2, iy+Offset)*shi(1, 3, iz+Offset)
                        DNDZi = shi(1, 1, ix+Offset)*shi(1, 2, iy+Offset)*shi(2, 3, iz+Offset)
                        ! The nodal mass is not too small
                        if (gn%Mg_F > 0.0) then
                            ! particle  velocity
                            vs = gn%PXg/gn%Mg; vf = gn%PXg_F/gn%Mg_F
                            fx_drag = SHPi*Drag*(vs - vf)
                            ! for fluid
                            fxf_int(1) = (DNDXi*PorePress_)*vol_
                            fxf_int(2) = (DNDYi*PorePress_)*vol_
                            fxf_int(3) = (DNDZi*PorePress_)*vol_
                            fxf_ext = SHPi*fx_f
                            ! for solid
                            fx_int(1) = - (DNDXi*sxx + DNDYi*sxy + DNDZi*szx)*vol_
                            fx_int(2) = - (DNDXi*sxy + DNDYi*syy + DNDZi*syz)*vol_
                            fx_int(3) = - (DNDXi*szx + DNDYi*syz + DNDZi*szz)*vol_
                            fx_ext = SHPi*fx
                            !$OMP CRITICAL
                            gn%FXg_F = gn%FXg_F + fxf_int + fxf_ext + fx_drag
                            gn%FXg = gn%FXg + fx_int + fxf_int + fx_ext
                            !$OMP END CRITICAL
                            if(contact) then
                                CP => CP_list(comID_, pt%basenode1(1)+ix, pt%basenode1(2)+iy, pt%basenode1(3)+iz)
                                !$OMP ATOMIC
                                CP%ndir(1) = CP%ndir(1) + DNDXi*mp
                                !$OMP ATOMIC
                                CP%ndir(2) = CP%ndir(2) + DNDYi*mp
                                !$OMP ATOMIC
                                CP%ndir(3) = CP%ndir(3) + DNDZi*mp
                            end if
                          
                        end if

                    end do    !ix
                end do    !iy
            end do    !iz
        end do    !p
        !$OMP END PARALLEL DO
    end do    !b

  end subroutine GridNodalForce

  subroutine ParticleUpdate()
  !------------------------------------------------------------------
  !-  Purpose                                                       -
  !-            Update particle acceleration                        -
  !------------------------------------------------------------------
    use ParticleData
    use GridData
    implicit none

    integer:: i, b, p, parBegin, parEnd       ! loop counter
    integer:: ix, iy, iz, comID_ = 1
    real(8):: node_vel(3), node_vel_f(3), node_vel_(3), node_vel_f_(3)
    real(8):: node_acc(3), node_acc_f(3)
    real(8):: vx(3), vx_f(3)
    real(8):: ax(3), ax_f(3)
    real(8):: vg(3), vxg(3), vg_f(3), vxg_f(3)
    real(8):: para_FLIP = 1.0, para_FLIP_f = 1.0
    real(8):: para_PIC = 0.1, para_PIC_f = 0.1
    real(8):: Cp(3, 3), dis_det(3)
    real(8), dimension(2, 3, 4):: shi
    real(8):: SHPi
  
    type(Particle), POINTER:: pt
    type(GridNode), POINTER:: gn
  
    ! Update particle position and velocity
    do b = 1, nb_body
        parBegin = body_list(b)%par_begin
        parEnd = body_list(b)%par_End
        ! Get comID from body
        if(contact)  comID_ = body_list(b)%comID
        !$OMP PARALLEL DO SCHEDULE(DYNAMIC,3) &
        !$OMP PRIVATE(pt,vx,vx_f,ax,ax_f,vxg,vxg_f,Cp,shi) &
        !$OMP PRIVATE(ix,iy,iz,SHPi,gn,i,dis_det) &
        !$OMP PRIVATE(node_vel,node_vel_f,node_vel_,node_vel_f_,node_acc,node_acc_f)
        do p = parBegin, parEnd       ! Loop over all particles
            pt => particle_list(p)
            ! old particle position
            ! Particle p is out of the computational region
            if (any(pt%icell < 0)) cycle
            ax = 0.0d0; ax_f = 0.0d0
            vx = 0.0d0; vx_f = 0.0d0
            vxg = 0.0d0; vxg_f = 0.0d0
            Cp = 0.0d0
            ! Calculate the shape functions and their derivatives
            if (CUBI) then
                shi = NShape_Cubic(comID_, pt%basenode1, pt%Xp)
            else if (GIMP) then
                shi = NShape_GIMP(comID_, pt%basenode1, pt%Xp)
            else
                shi = NShape(comID_, pt%basenode1, pt%Xp)
            end if
            ! loop all the neighbor nodes
            do iz = InfluInit, InfluEnd
                if (pt%basenode1(3)+iz < 1 .OR. pt%basenode1(3)+iz > NGz) cycle
                do iy = InfluInit, InfluEnd
                    if (pt%basenode1(2)+iy < 1 .OR. pt%basenode1(2)+iy > NGy) cycle
                    do ix = InfluInit, InfluEnd
                        if (pt%basenode1(1)+ix < 1 .OR. pt%basenode1(1)+ix > NGx) cycle
                      
                        gn => node_list(comID_, pt%basenode1(1)+ix, pt%basenode1(2)+iy, pt%basenode1(3)+iz)
                        if (.NOT.gn%interp) cycle
                        SHPi = shi(1, 1, ix+Offset)*shi(1, 2, iy+Offset)*shi(1, 3, iz+Offset)

                        if (gn%Mg_F > 0.0) then
                            ! velocity and acceleration on the node
                            node_vel_f = gn%PXg_F/gn%Mg_F; node_vel = gn%PXg/gn%Mg
                            node_acc_f = gn%FXg_F/gn%Mg_FF; node_acc = (gn%FXg-(gn%FXg_F/gn%Mg_FF)*gn%Mg_F)/gn%Mg
                            ! velocity at t + 1
                            node_vel_f_ = node_vel_f + node_acc_f * DT
                            node_vel_ = node_vel + node_acc * DT
                            
                            ax_f = ax_f + SHPi*node_acc_f
                            ax = ax + SHPi*node_acc
                            vx_f = vx_f + SHPi*node_vel_f
                            vx = vx + SHPi*node_vel
                            vxg_f = vxg_f + SHPi*node_vel_f_
                            vxg = vxg + SHPi*node_vel_
                            ! APIC
                            !dis_det = gn%Xg - pt%Xp
                            !do i = 1, 3
                            !    Cp(1, i) = Cp(1, i) + SHPi*node_vel_(1)*dis_det(i)
                            !    Cp(2, i) = Cp(2, i) + SHPi*node_vel_(2)*dis_det(i)
                            !    Cp(3, i) = Cp(3, i) + SHPi*node_vel_(3)*dis_det(i)
                            !end do
                        end if
          
                    end do    !ix
                end do    !iy
            end do    !iz
            !$OMP CRITICAL
            ! FLIP and PIC
            ! Update Particle Velocity at time step t + 1
            pt%VXp = para_FLIP * (pt%VXp + ax * DT) + (1.0 - para_FLIP) * vxg
            pt%VXp_F = para_FLIP_f * (pt%VXp_F + ax_f * DT) + (1.0 - para_FLIP_f) * vxg_f
            ! the next particle position
            pt%XX = pt%Xp + vxg * DT
            pt%Disp_p = pt%Disp_p + vxg * DT
            ! FLIP with PIC damping
            !pt%XX = pt%Xp + vxg*DT - 0.5*para_PIC*(pt%VXp - vx)*DT - 0.5 * ax * DT * DT         ! the next particle position
            !pt%Disp_p = pt%Disp_p + vxg*DT - 0.5*para_PIC*(pt%VXp - vx)*DT - 0.5 * ax * DT * DT
            !pt%VXp = pt%VXp + ax*DT - para_PIC*(pt%VXp - vx)        ! Update Particle Velocity at time step t+1
            !pt%VXp_F = pt%VXp_F + ax_f*DT - para_PIC_f*(pt%VXp_F - vx_f)

            !pt%Cp_apic = Cp*(3.0/(DCell**2))
            !$OMP END CRITICAL
        end do    ! p
        !$OMP END PARALLEL DO
    end do    !b
  
  end subroutine  ParticleUpdate

  subroutine GridMomentumUpdate()
  !------------------------------------------------------------------
  !-  Purpose                                                       -
  !-      1. recalculate the grid node momentum by mapping          -
  !-         the updated particle information                       -
  !------------------------------------------------------------------
    use ParticleData
    use GridData
    use FFI, only: iomsg
    implicit none

    integer:: i, b, p, parBegin, parEnd    ! loop counter
    integer:: ix, iy, iz, comID_ = 1
    real(8):: mp, mp_f, vp(3), vp_f(3)
    real(8):: Cp(3, 3), dis_det(3)
    real(8), dimension(2, 3, 4):: shi
    real(8):: SHPi

    type(Particle), POINTER:: pt
    type(GridNode), POINTER:: gn
    ! Calculate the grid nodal masses, moemntum only 
    ! Reset Grid data
    do i = 1, 3
        node_list%PXg(i) = 0.0d0
        node_list%PXg_F(i) = 0.0d0
    end do
  
    ! Recalculate the grid node momentum
    do b = 1, nb_body     ! Loop over all bodies
        parBegin = body_list(b)%par_begin
        parEnd = body_list(b)%par_End
        if(contact) comID_ = body_list(b)%comID ! Get comID from body
        !$OMP PARALLEL DO SCHEDULE(DYNAMIC,3) &
        !$OMP PRIVATE(pt,mp,mp_f,vp,vp_f,Cp,shi) &
        !$OMP PRIVATE(ix,iy,iz,SHPi,gn,dis_det)
        do p = parBegin, parEnd    ! Loop over all particles of the body
            pt => particle_list(p)
            ! old particle position
            ! Particle p is out of the computational region
            if (any(pt%icell < 0)) cycle
            mp = pt%Mass; vp = pt%VXp
            mp_f = pt%Mass_F; vp_f = pt%VXp_F
            Cp = pt%Cp_apic
            ! Calculate the shape functions and their derivatives
            if (CUBI) then
                shi = NShape_Cubic(comID_, pt%basenode1, pt%Xp)
            else if (GIMP) then
                shi = NShape_GIMP(comID_, pt%basenode1, pt%Xp)
            else
                shi = NShape(comID_, pt%basenode1, pt%Xp)
            end if
            ! loop all the neighbor nodes
            do iz = InfluInit, InfluEnd
                if (pt%basenode1(3)+iz < 1 .OR. pt%basenode1(3)+iz > NGz) cycle
                do iy = InfluInit, InfluEnd
                    if (pt%basenode1(2)+iy < 1 .OR. pt%basenode1(2)+iy > NGy) cycle
                    do ix = InfluInit, InfluEnd
                        if (pt%basenode1(1)+ix < 1 .OR. pt%basenode1(1)+ix > NGx) cycle
                      
                        gn => node_list(comID_, pt%basenode1(1)+ix, pt%basenode1(2)+iy, pt%basenode1(3)+iz)
                        if (.NOT.gn%interp) cycle
                        SHPi = shi(1, 1, ix+Offset)*shi(1, 2, iy+Offset)*shi(1, 3, iz+Offset)
                        dis_det = gn%Xg - pt%Xp
                        !$OMP CRITICAL
                        gn%PXg = gn%PXg + SHPi*mp*vp
                        ! APIC
                        !gn%PXg = gn%PXg + SHPi*mp*(vp + MATMUL(Cp, dis_det))
                        gn%PXg_F = gn%PXg_F + SHPi*mp_f*vp_f
                        !$OMP END CRITICAL
                    end do    !ix
                end do    !iy
            end do    !iz
        end do    !p
        !$OMP END PARALLEL DO
    end do    !b
  
  end subroutine GridMomentumUpdate

  subroutine ParticleStrainUpdate()
  !------------------------------------------------------------------
  !-  Purpose                                                       -
  !-         Calculate the strain rate and spin tensor              -
  !------------------------------------------------------------------
    use ParticleData
    use MaterialData
    use GridData
    use FFI, only: iomsg
    implicit none
  
    integer:: b, p, parBegin, parEnd     ! loop counter
    integer:: ix, iy, iz, comID_ = 1
    real(8):: vs(3), vf(3)
    real(8):: de(6), vort(3), vel_grad(3,3), de_f(6)
    real(8):: I_iden(3,3), det_defgrad(3,3), def_grad(3,3), Jacobi_
    real(8), dimension(2, 3, 4):: shi
    real(8):: SHPi, DNDXi, DNDYi, DNDZi
  
    type(Particle), POINTER:: pt
    type(GridNode), POINTER:: gn
  
    data I_iden /1.0, 0.0, 0.0, &
                 0.0 ,1.0, 0.0, &
                 0.0, 0.0, 1.0/
  
    ! Calculate the increment strain and vorticity
    ! de(i) comply the Voigt rule (d11, d22, d33, d23, d13, d12)
    do b = 1, nb_body
        parBegin = body_list(b)%par_begin
        parEnd = body_list(b)%par_End
        
        if(contact) comID_ = body_list(b)%comID       ! Get comID from body
        !$OMP PARALLEL DO SCHEDULE(DYNAMIC,3) &
        !$OMP PRIVATE(pt,de,vort,vel_grad,de_f,shi) &
        !$OMP PRIVATE(ix,iy,iz,gn,SHPi,DNDXi,DNDYi,DNDZi) &
        !$OMP PRIVATE(vs,vf,det_defgrad,def_grad,Jacobi_)
        do p = parBegin, parEnd                       ! Loop over all particles
            pt => particle_list(p)
            ! old particle position
            ! Particle p is out of the computational region
            if (any(pt%icell < 0)) cycle
            de = 0d0      ! Incremental strain
            vort = 0d0    ! Incremental vorticity
            de_f = 0d0
            vel_grad = 0d0
            ! Calculate the shape functions and their derivatives
            if (CUBI) then
                shi = NShape_Cubic(comID_, pt%basenode1, pt%Xp)
            else if (GIMP) then
                shi = NShape_GIMP(comID_, pt%basenode1, pt%Xp)
            else
                shi = NShape(comID_, pt%basenode1, pt%Xp)
            end if
            ! Mapping from grid to particle
            ! loop all the neighbor nodes
            do iz = InfluInit, InfluEnd
                if (pt%basenode1(3)+iz < 1 .OR. pt%basenode1(3)+iz > NGz) cycle
                do iy = InfluInit, InfluEnd
                    if (pt%basenode1(2)+iy < 1 .OR. pt%basenode1(2)+iy > NGy) cycle
                    do ix = InfluInit, InfluEnd
                        if (pt%basenode1(1)+ix < 1 .OR. pt%basenode1(1)+ix > NGx) cycle
                      
                        gn => node_list(comID_, pt%basenode1(1)+ix, pt%basenode1(2)+iy, pt%basenode1(3)+iz)
                        if (.NOT.gn%interp) cycle
                        SHPi = shi(1, 1, ix+Offset)*shi(1, 2, iy+Offset)*shi(1, 3, iz+Offset)
                        DNDXi = shi(2, 1, ix+Offset)*shi(1, 2, iy+Offset)*shi(1, 3, iz+Offset)
                        DNDYi = shi(1, 1, ix+Offset)*shi(2, 2, iy+Offset)*shi(1, 3, iz+Offset)
                        DNDZi = shi(1, 1, ix+Offset)*shi(1, 2, iy+Offset)*shi(2, 3, iz+Offset)
                        if (gn%Mg_F > 0.0) then
                            ! Grid nodal velocity
                            vs = gn%PXg / gn%Mg; vf = gn%PXg_F / gn%Mg_F
                            ! The rate of deformation tensor
                            de(1) = de(1) + DNDXi*vs(1)   ! D11
                            de(2) = de(2) + DNDYi*vs(2)   ! D22
                            de(3) = de(3) + DNDZi*vs(3)   ! D33
                            de(4) = de(4) + (DNDZi*vs(2) + DNDYi*vs(3))*0.5   ! D23
                            de(5) = de(5) + (DNDXi*vs(3) + DNDZi*vs(1))*0.5   ! D31
                            de(6) = de(6) + (DNDYi*vs(1) + DNDXi*vs(2))*0.5   ! D12
                            ! The spin tensor
                            vort(1) = vort(1) + (DNDYi*vs(3) - DNDZi*vs(2))   ! W32
                            vort(2) = vort(2) + (DNDZi*vs(1) - DNDXi*vs(3))   ! W13
                            vort(3) = vort(3) + (DNDXi*vs(2) - DNDYi*vs(1))   ! W21
                            ! The velocity gradient tensor
                            vel_grad(1,1) = vel_grad(1,1) + DNDXi*vs(1)
                            vel_grad(1,2) = vel_grad(1,2) + DNDYi*vs(1)
                            vel_grad(1,3) = vel_grad(1,3) + DNDZi*vs(1)
                            vel_grad(2,1) = vel_grad(2,1) + DNDXi*vs(2)
                            vel_grad(2,2) = vel_grad(2,2) + DNDYi*vs(2)
                            vel_grad(2,3) = vel_grad(2,3) + DNDZi*vs(2)
                            vel_grad(3,1) = vel_grad(3,1) + DNDXi*vs(3)
                            vel_grad(3,2) = vel_grad(3,2) + DNDYi*vs(3)
                            vel_grad(3,3) = vel_grad(3,3) + DNDZi*vs(3)
                            ! For fluid
                            de_f(1) = de_f(1) + DNDXi*vf(1)
                            de_f(2) = de_f(2) + DNDYi*vf(2)
                            de_f(3) = de_f(3) + DNDZi*vf(3)
                            de_f(4) = de_f(4) + (DNDZi*vf(2) + DNDYi*vf(3))*0.5
                            de_f(5) = de_f(5) + (DNDXi*vf(3) + DNDZi*vf(1))*0.5
                            de_f(6) = de_f(6) + (DNDYi*vf(1) + DNDXi*vf(2))*0.5
                        end if
                    end do    !ix
                end do    !iy
            end do    !iz
            det_defgrad = I_iden + vel_grad * DT
            def_grad = matmul(det_defgrad, pt%Fp)
            Jacobi_ = def_grad(1, 1)*(def_grad(2, 2)*def_grad(3, 3) - def_grad(2, 3)*def_grad(3, 2)) - &
                      def_grad(1, 2)*(def_grad(2, 1)*def_grad(3, 3) - def_grad(2, 3)*def_grad(3, 1)) + &
                      def_grad(1, 3)*(def_grad(2, 1)*def_grad(3, 2) - def_grad(2, 2)*def_grad(3, 1))
          
            pt%de_ = de
            pt%vort_ = vort * 0.5
            pt%de_f_ = de_f
            pt%Fp = def_grad
            pt%Jacobi = Jacobi_
        end do    !p
        !$OMP END PARALLEL DO
    end do    !b
  
  end subroutine ParticleStrainUpdate
   
  subroutine ParticleStressUpdate()
  !------------------------------------------------------------------
  !-  Purpose                                                       -
  !-         Update stresses by appropriate constitution law        -
  !------------------------------------------------------------------
    use ParticleData
    use ConstitutiveModel, only: Constitution
    implicit none

    integer:: b, p, parBegin, parEnd      ! loop counter
    real(8):: Fp_F_(3, 3), Jacobian_
    type(Particle), POINTER:: pt

    do b = 1, nb_body
        parBegin = body_list(b)%par_begin
        parEnd = body_list(b)%par_End
        !$OMP PARALLEL DO PRIVATE(Fp_F_,Jacobian_) SCHEDULE(DYNAMIC,3)
        do p = parBegin, parEnd        ! Loop over all particles
            pt => particle_list(p)
            ! use old position
            ! Particle p is out of the computational region
            if (any(pt%icell < 0)) cycle
            Fp_F_ = ((pt%Jacobi_F/pt%Jacobi)**(1.0/3.0))*pt%Fp
            Jacobian_ = Fp_F_(1, 1)*(Fp_F_(2, 2)*Fp_F_(3, 3) - Fp_F_(2, 3)*Fp_F_(3, 2)) - &
                        Fp_F_(1, 2)*(Fp_F_(2, 1)*Fp_F_(3, 3) - Fp_F_(2, 3)*Fp_F_(3, 1)) + &
                        Fp_F_(1, 3)*(Fp_F_(2, 1)*Fp_F_(3, 2) - Fp_F_(2, 2)*Fp_F_(3, 1))
            
            pt%Fp_F = Fp_F_
            pt%Jacobian = Jacobian_
            ! Update stress by constitution law
            !$OMP CRITICAL
            call Constitution(b, p)
            !$OMP END CRITICAL
        end do !p
        !$OMP END PARALLEL DO
    end do    !b
  
  end subroutine ParticleStressUpdate
  
  subroutine GridMVPInformation()
  !--------------------------------------------------------
  !-  Purpose                                             -
  !-      1. Calculate the MVP parameters on grid node    -
  !--------------------------------------------------------
    use ParticleData
    use GridData
    implicit none

    integer:: b, p, parBegin, parEnd     ! loop counter
    integer:: ix, iy, iz, comID_ = 1
    real(8), dimension(2, 3, 4):: shi
    real(8):: SHPi
    
    type(Particle), POINTER:: pt
    type(GridNode), POINTER:: gn
  
    ! Reset Grid data
    node_list%Qi = 0.0d0
    node_list%Hi = 0.0d0

    do b = 1, nb_body
        parBegin = body_list(b)%par_begin
        parEnd = body_list(b)%par_End

        if(contact) comID_ = body_list(b)%comID   ! Get comID from body
        !$OMP PARALLEL DO SCHEDULE(DYNAMIC,3) &
        !$OMP PRIVATE(pt,shi) &
        !$OMP PRIVATE(ix,iy,iz,gn,SHPi)
        do p = parBegin, parEnd                   ! Loop over all particles of the body
            pt => particle_list(p)
            ! Particle p is out of the computational region
            if (any(pt%icell < 0)) cycle
            ! Calculate the shape functions and their derivatives
            if (CUBI) then
                shi = NShape_Cubic(comID_, pt%basenode1, pt%Xp)
            else if (GIMP) then
                shi = NShape_GIMP(comID_, pt%basenode1, pt%Xp)
            else
                shi = NShape(comID_, pt%basenode1, pt%Xp)
            end if
            ! loop all the neighbor nodes
            do iz = InfluInit, InfluEnd
                if (pt%basenode1(3)+iz < 1 .OR. pt%basenode1(3)+iz > NGz) cycle
                do iy = InfluInit, InfluEnd
                    if (pt%basenode1(2)+iy < 1 .OR. pt%basenode1(2)+iy > NGy) cycle
                    do ix = InfluInit, InfluEnd
                        if (pt%basenode1(1)+ix < 1 .OR. pt%basenode1(1)+ix > NGx) cycle
                      
                        SHPi = shi(1, 1, ix+Offset)*shi(1, 2, iy+Offset)*shi(1, 3, iz+Offset)
                        gn => node_list(comID_, pt%basenode1(1)+ix, pt%basenode1(2)+iy, pt%basenode1(3)+iz)
                        if (.NOT.gn%interp) cycle
                      
                        gn%Qi = gn%Qi + SHPi*pt%Qp*pt%dPorePress*pt%Vol
                        gn%Hi = gn%Hi + SHPi*pt%Qp*pt%Qp*pt%Vol
                      
                    end do    !ix
                end do    !iy
            end do    !iz
            
        end do    !p
        !$OMP END PARALLEL DO
    end do    !b

  end subroutine GridMVPInformation

  subroutine GridCoefficient()
  !-------------------------------------------------------------------
  !-  Purpose                                                        -
  !-      1. calculate the MSL coefficient on the node               -
  !-------------------------------------------------------------------
    use GridData
    implicit none

    integer:: i, ix, iy, iz, comID_ = 1  ! loop counter

    type(GridNode), POINTER:: gn
  
    node_list%ai = 0.0d0
    !$OMP PARALLEL DO COLLAPSE(3) PRIVATE(gn) SCHEDULE(DYNAMIC,3)
    do iz = 1, NGz
        do iy = 1, NGy
            do ix = 1, NGx
                gn => node_list(comID_, ix, iy, iz)
                if (gn%Hi > 0) then
                    gn%ai = - gn%Qi / gn%Hi
                end if
            end do
        end do
    end do
    !$OMP END PARALLEL DO
    
  end subroutine GridCoefficient
  
  subroutine PorePressUpdate()
  !-----------------------------------------------------
  !-  Purpose                                          -
  !-         Calculate the pore pressure               -
  !-----------------------------------------------------
    use ParticleData
    use GridData
    implicit none

    integer:: b, p, parBegin, parEnd     ! loop counter
    integer:: ix, iy, iz, comID_ = 1
    real(8):: ap_
    real(8), dimension(2, 3, 4):: shi
    real(8):: SHPi
  
    type(Particle), POINTER:: pt
    type(GridNode), POINTER:: gn

    do b = 1, nb_body
        parBegin = body_list(b)%par_begin
        parEnd = body_list(b)%par_End
      
        if(contact) comID_ = body_list(b)%comID       ! Get comID from body
        !$OMP PARALLEL DO SCHEDULE(DYNAMIC,3) &
        !$OMP PRIVATE(pt,ap_,shi) &
        !$OMP PRIVATE(ix,iy,iz,gn,SHPi)
        do p = parBegin, parEnd                       ! Loop over all particles
            pt => particle_list(p)
            ! old particle position
            ! Particle p is out of the computational region
            if (any(pt%icell < 0)) cycle
            ap_ = 0.0d0
            ! Calculate the shape functions and their derivatives
            if (CUBI) then
                shi = NShape_Cubic(comID_, pt%basenode1, pt%Xp)
            else if (GIMP) then
                shi = NShape_GIMP(comID_, pt%basenode1, pt%Xp)
            else
                shi = NShape(comID_, pt%basenode1, pt%Xp)
            end if
            ! loop all the neighbor nodes
            do iz = InfluInit, InfluEnd
                if (pt%basenode1(3)+iz < 1 .OR. pt%basenode1(3)+iz > NGz) cycle
                do iy = InfluInit, InfluEnd
                    if (pt%basenode1(2)+iy < 1 .OR. pt%basenode1(2)+iy > NGy) cycle
                    do ix = InfluInit, InfluEnd
                        if (pt%basenode1(1)+ix < 1 .OR. pt%basenode1(1)+ix > NGx) cycle
                      
                        SHPi = shi(1, 1, ix+Offset)*shi(1, 2, iy+Offset)*shi(1, 3, iz+Offset)
                        gn => node_list(comID_, pt%basenode1(1)+ix, pt%basenode1(2)+iy, pt%basenode1(3)+iz)
                        if (.NOT.gn%interp) cycle
                      
                        ap_ = ap_ + SHPi*gn%ai
                      
                    end do    !ix
                end do    !iy
            end do    !iz
          
            pt%ap = ap_
            pt%PorePress = pt%PorePress + pt%Qp*pt%ap

            ! the next particle position
            pt%Xp = pt%XX
   
        end do    !p
        !$OMP END PARALLEL DO
    end do    !b

  end subroutine PorePressUpdate
  
  subroutine BoundaryCellSearch()
  !------------------------------------------------------------------
  !-  Purpose                                                       -
  !-     1. serach for boundary cell                                -
  !------------------------------------------------------------------
    use ParticleData
    use GridData
    implicit none
    
    integer:: ix, iy, iz      ! loop counter
    integer:: id1, id2, id3, id4
    integer:: p, id
    
    type(CellElement), POINTER:: ce
    type(Particle), POINTER:: pt

    !$OMP PARALLEL DO COLLAPSE(3) SCHEDULE(DYNAMIC,3) &
    !$OMP PRIVATE(ce,id1,id2,id3,id4,p,id,pt)
    do iz = 1, NumCellz
        do iy = 1, NumCelly
            do ix = 1, NumCellx
                ce => cell_list(ix, iy, iz)
                if (ce%index == 1) then
                    id1 = ix - 1
                    if (id1 < 1) then
                        ce%index  = 2
                    else if (cell_list(id1, iy, iz)%index == 0) then
                        ce%index  = 2
                    end if
                    id2 = iz - 1
                    if (id2 < 1) then
                        ce%index  = 3
                    else if (cell_list(ix, iy, id2)%index == 0) then
                        ce%index  = 3
                    end if
                  
                    id3 = ix + 1
                    if (id3 > NumCellx) then
                        ce%index  = 4
                    else if (cell_list(id3, iy, iz)%index == 0) then
                        ce%index  = 4
                    end if
                    id4 = iz + 1
                    if (id4 > NumCellz) then
                        ce%index  = 5
                    else if (cell_list(ix, iy, id4)%index == 0) then
                        ce%index  = 5
                    end if

                    ! Loop over all particles in cell
                    do p = 1, ce%pnum
                        id = ce%Cell_par(p)
                        pt => particle_list(id)
                        pt%boundaryflag = ce%index
                    end do
                end if
            end do
        end do
    end do
    !$OMP END PARALLEL DO
  end subroutine BoundaryCellSearch
    
  subroutine CellAverage1()
  !------------------------------------------
  !-  Purpose                               -
  !-  1. volumetric strain rate smoothing   -
  !-  2. F-bar method                       -
  !------------------------------------------
    use ParticleData
    use GridData
    implicit none
 
    integer:: ix, iy, iz, p, id   ! loop counter
    real(8):: strain_vol, strain_vol_f
    real(8):: sum_strain_vol, sum_strain_vol_F
    real(8):: sum_mass, sum_mass_F
    real(8):: sum_Jacobi
    real(8):: sum_vol
  
    type(Particle), POINTER:: pt
    type(CellElement), POINTER:: ce
    !$OMP PARALLEL DO COLLAPSE(3) SCHEDULE(DYNAMIC,3) &
    !$OMP PRIVATE(ce,p,id,pt) &
    !$OMP PRIVATE(strain_vol,strain_vol_f,sum_strain_vol,sum_strain_vol_F) &
    !$OMP PRIVATE(sum_mass,sum_mass_F,sum_Jacobi,sum_vol)
    do iz = 1, NumCellz
        do iy = 1, NumCelly
            do ix = 1, NumCellx
                ce => cell_list(ix, iy, iz)
                if (ce%pnum .GT. 0) then
                    sum_strain_vol = 0.0
                    sum_strain_vol_F = 0.0
                    sum_Jacobi = 0.0
                    sum_mass = 0.0
                    sum_mass_F = 0.0
                    sum_vol = 0.0
                    ! sum up
                    do p = 1, ce%pnum    ! Loop over all particles in cell
                        id = ce%Cell_par(p)
                        pt => particle_list(id)
                        ! volumetric strain rate and Jacobi
                        sum_strain_vol = sum_strain_vol + (pt%de_(1) + pt%de_(2) + pt%de_(3))*pt%mass
                        sum_strain_vol_F = sum_strain_vol_F + (pt%de_f_(1)+pt%de_f_(2)+pt%de_f_(3))*pt%mass_F
                        sum_Jacobi = sum_Jacobi + pt%Jacobi*pt%Volu
                        ! mass
                        sum_mass = sum_mass + pt%mass
                        sum_mass_F = sum_mass_F + pt%mass_F
                        ! volume
                        sum_vol = sum_vol + pt%Volu
                    end do    ! p
                    ! cell average
                    do p = 1, ce%pnum
                        id = ce%Cell_par(p)
                        pt => particle_list(id)
                        strain_vol = (pt%de_(1) + pt%de_(2) + pt%de_(3)) / 3.0
                        strain_vol_f = (pt%de_f_(1) + pt%de_f_(2) + pt%de_f_(3)) / 3.0
                        
                        !$OMP ATOMIC
                        pt%de_(1) = pt%de_(1) - strain_vol + (sum_strain_vol/(sum_mass*3.0))
                        !$OMP ATOMIC
                        pt%de_(2) = pt%de_(2) - strain_vol + (sum_strain_vol/(sum_mass*3.0))
                        !$OMP ATOMIC
                        pt%de_(3) = pt%de_(3) - strain_vol + (sum_strain_vol/(sum_mass*3.0))
                        !$OMP ATOMIC
                        pt%de_f_(1) = pt%de_f_(1) - strain_vol_f + (sum_strain_vol_F/(sum_mass_F*3.0))
                        !$OMP ATOMIC
                        pt%de_f_(2) = pt%de_f_(2) - strain_vol_f + (sum_strain_vol_F/(sum_mass_F*3.0))
                        !$OMP ATOMIC
                        pt%de_f_(3) = pt%de_f_(3) - strain_vol_f + (sum_strain_vol_F/(sum_mass_F*3.0))
                        
                        pt%Jacobi_F = sum_Jacobi / sum_vol
                    end do    ! p
                end if

            end do    !ix
        end do    !iy
    end do    !iz
    !$OMP END PARALLEL DO
  end subroutine CellAverage1

  subroutine ApplyBoundaryConditions()
  !------------------------------------------------------------------
  !-  Purpose                                                       -
  !-      1. apply boundary condition                               -
  !------------------------------------------------------------------
    use GridData
    use ParticleData, only: nb_component
    implicit none

    integer:: ix, iy, iz, comID_ = 1  ! loop counter
    type(GridNode), POINTER:: gn
    !$OMP PARALLEL DO COLLAPSE(3) PRIVATE(gn) SCHEDULE(DYNAMIC,3)
    do iz = 1, NGz
        do iy = 1, NGy
            do ix = 1, NGx
                if (nb_component == 1) then
                    gn => node_list(comID_, ix, iy, iz)
                    if (.NOT.gn%isboundary) cycle
                    if (.NOT.gn%interp) cycle
                  
                    if (gn%Fix_x) then    ! Grid node n is fixed in x direction
                        gn%PXg(1) = 0.0
                        gn%FXg(1) = 0.0
                    end if
                    if (gn%Fix_y) then    ! Grid node n is fixed in y direction
                        gn%PXg(2) = 0.0
                        gn%FXg(2) = 0.0
                    end if
                    if (gn%Fix_z) then    ! Grid node n is fixed in z direction
                        gn%PXg(3) = 0.0
                        gn%FXg(3) = 0.0
                    end if
                  
                    if (gn%Fix_xf) then    ! Grid node n is fixed in x direction
                        gn%PXg_F(1) = 0.0
                        gn%FXg_F(1) = 0.0
                    end if
                    if (gn%Fix_yf) then    ! Grid node n is fixed in y direction
                        gn%PXg_F(2) = 0.0
                        gn%FXg_F(2) = 0.0
                    end if
                    if (gn%Fix_zf) then    ! Grid node n is fixed in z direction
                        gn%PXg_F(3) = 0.0
                        gn%FXg_F(3) = 0.0
                    end if
                else
                    do comID_ = 1, nb_component
                        gn => node_list(comID_, ix, iy, iz)
                        if (.NOT.gn%isboundary) cycle
                        if (.NOT.gn%interp) cycle
                      
                        if (gn%Fix_x) then    ! Grid node n is fixed in x direction
                            gn%PXg(1) = 0.0
                            gn%FXg(1) = 0.0
                        end if
                        if (gn%Fix_y) then    ! Grid node n is fixed in y direction
                            gn%PXg(2) = 0.0
                            gn%FXg(2) = 0.0
                        end if
                        if (gn%Fix_z) then    ! Grid node n is fixed in z direction
                            gn%PXg(3) = 0.0
                            gn%FXg(3) = 0.0
                        end if
                      
                        if (gn%Fix_xf) then    ! Grid node n is fixed in x direction
                            gn%PXg_F(1) = 0.0
                            gn%FXg_F(1) = 0.0
                        end if
                        if (gn%Fix_yf) then    ! Grid node n is fixed in y direction
                            gn%PXg_F(2) = 0.0
                            gn%FXg_F(2) = 0.0
                        end if
                        if (gn%Fix_zf) then    ! Grid node n is fixed in z direction
                            gn%PXg_F(3) = 0.0
                            gn%FXg_F(3) = 0.0
                        end if
                    end do
                end if
              
            end do    !ix
        end do    !iy
    end do    !iz
    !$OMP END PARALLEL DO
    
  end subroutine ApplyBoundaryConditions

  subroutine Lagr_NodContact()
  !------------------------------------------------------------------
  !-  Purpose                                                       -
  !-      1. Establishing the nodal contact criteria,          and  -
  !-      2. correct  the normal vectors                       and  -
  !-      3. Apply the contact force and adjust nodal velocities    -
  !------------------------------------------------------------------
    use  ParticleData
    use  GridData
    implicit none

    integer:: ix, iy, iz        ! loop counter
    real(8):: nx, ny, nz, tt,crit
    real(8):: delt_dis
    real(8):: val_fslip, val_tanforce, val_nomforce
    real(8):: nomforce(3), tanforce(3), cforce(3)

    type(GridNode), POINTER:: gn1
    type(GridNode), POINTER:: gn2
    type(ContactGridNodeProperty), POINTER:: CP1
    type(ContactGridNodeProperty), POINTER:: CP2

    tot_cont_for = 0.0 ! the total contact force between of bodies
  
    ! calculate contact force and adjust the nodal force and momentum
    !$OMP PARALLEL DO COLLAPSE(3) SCHEDULE(DYNAMIC,3) &
    !$OMP PRIVATE(CP1,CP2,gn1,gn2,nx,ny,nz,tt,crit,delt_dis) &
    !$OMP PRIVATE(val_tanforce,tanforce,val_nomforce,nomforce,cforce,val_fslip) &
    !$OMP REDUCTION(+:tot_cont_for)
    do iz = 1, NGz
        do iy = 1, NGy
            do ix = 1, NGx
                CP1 =>  CP_list(1, ix, iy, iz)
                CP2 =>  CP_list(2, ix, iy, iz)
                gn1 => node_list(1, ix, iy, iz)
                gn2 => node_list(2, ix, iy, iz)

                ! recalculate the nodal normal direction    
                ! if normbody 0 then using average method 
                ! if 1,using abody; if 2,using bbody
                if (normbody == 0) then
                    nx = CP1%ndir(1)  - CP2%ndir(1)
                    ny = CP1%ndir(2)  - CP2%ndir(2)
                    nz = CP1%ndir(3)  - CP2%ndir(3)
                else if (normbody == 1) then         
                    nx = CP1%ndir(1)
                    ny = CP1%ndir(2)
                    nz = CP1%ndir(3)
                else if (normbody == 2) then
                    nx = - CP2%ndir(1)
                    ny = - CP2%ndir(2)
                    nz = - CP2%ndir(3)
                end if

                ! unitize normal vector        
                tt = sqrt(nx*nx + ny*ny + nz*nz)        
                if(tt > epsilon(tt)) then
                    nx = nx / tt
                    ny = ny / tt
                    nz = nz / tt
                end if

                CP1%ndir(1) = nx; ! Nodal direction for contact 保证接触点的法向量共线
                CP1%ndir(2) = ny;    
                CP1%ndir(3) = nz;        
                CP2%ndir = - CP1%ndir        

                crit = 0.0

                ! contact criteria using the unit normal vectors

                if (gn1%Mg > CutOff .AND. gn2%Mg > CutOff) then
                    crit = (gn1%Pxg(1)*gn2%Mg - gn2%Pxg(1)*gn1%Mg)*nx +&
                           (gn1%Pxg(2)*gn2%Mg - gn2%Pxg(2)*gn1%Mg)*ny +&  
                           (gn1%Pxg(3)*gn2%Mg - gn2%Pxg(3)*gn1%Mg)*nz

                    gn1%deltx =  gn1%deltx / gn1%Mg;
                    gn2%deltx =  gn2%deltx / gn2%Mg;

                    delt_dis = (gn1%deltx(1)-gn2%deltx(1))*nx + &
                               (gn1%deltx(2)-gn2%deltx(2))*ny + &
                               (gn1%deltx(3)-gn2%deltx(3))*nz       
                end if
                !calculate contact force for body 1, for body 2 contact force is negative
                if (crit > epsilon(crit) .and. abs(delt_dis) < 0.8*Dcell) then
                    tt = (gn1%Mg + gn2%Mg)*Dt
                    ! calculate the normal contact force
                    val_nomforce = crit/tt            !法向接触力大小
                    nomforce = val_nomforce*CP1%ndir  !法向接触力
                    ! for friction contact
                    if (fricfa > epsilon(1.0d-15)) then
                        ! calculate the contact force
                        cforce = (gn1%Pxg*gn2%Mg - gn2%Pxg*gn1%Mg)/tt         !接触力
                        ! calculate the tangential force
                        tanforce = cforce - nomforce       !切向接触力
                        val_tanforce = sqrt(DOT_PRODUCT(tanforce, tanforce))  !切向接触力大小
                        CP1%sdir = tanforce /val_tanforce  !切向量

                        val_fslip = fricfa * val_nomforce

                        if (val_fslip < val_tanforce) then
                            cforce = nomforce + val_fslip*CP1%sdir  !计算滑移接触力
                        end if
                    end if
                    ! for contact without friction 
                    if (fricfa == 0) then
                        cforce = nomforce
                    end if
                    ! add contact force to nodal force
                    ! adjust the nodal component by contact force
                    if (.not. gn1%Fix_x) then
                        !$OMP ATOMIC
                        gn1%Fxg(1) = gn1%Fxg(1) - cforce(1)
                        !$OMP ATOMIC
                        gn2%Fxg(1) = gn2%Fxg(1) + cforce(1)
                    end if
                  
                    if (.not. gn1%Fix_y) then
                        !$OMP ATOMIC
                        gn1%Fxg(2) = gn1%Fxg(2) - cforce(2)
                        !$OMP ATOMIC
                        gn2%Fxg(2) = gn2%Fxg(2) + cforce(2)
                    end if
              
                    if (.not. gn1%Fix_z) then
                        !$OMP ATOMIC
                        gn1%Fxg(3) = gn1%Fxg(3) - cforce(3)
                        !$OMP ATOMIC
                        gn2%Fxg(3) = gn2%Fxg(3) + cforce(3)
                    end if

                    tot_cont_for = tot_cont_for + cforce
          
                end if

            end do    !ix
        end do    !iy
    end do    !iz
    !$OMP END PARALLEL DO
  end subroutine Lagr_NodContact

end module update_step