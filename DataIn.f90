!------------------------------------------------------------------
!-   Data input procedures                                        -
!------------------------------------------------------------------
module DataIn

  use ParticleData
  use GridData
  use DataOut
  use MaterialData

  integer:: bodyCounter = 0 ! body counter
  integer:: comCounter = 0  ! component counter
  integer:: parCounter = 0  ! particle counter

contains

  subroutine InputPara()
!------------------------------------------------------------------
!-    purpose: input data both from file and by program           -
!------------------------------------------------------------------  
    use FFI
    implicit none

    integer:: idx

    call GETARG(1,jobname)

    idx = index(jobname,".")
    if (idx > 0) then
        fName = jobname(1:index(jobname,".")-1)
    else
        fName = jobname
    endif

    FileInp = trim(fName) // "example.mpm"          ! Input data
    FileOut = trim(fName) // "example.out"          ! Message file
    FileCurv = trim(fName) // "example_curv.dat"    ! Curve data (TecPlot)

    open(iow2, file = FileCurv, status = 'unknown')

    call FFIOpen()      ! Open the input file and initialization
    call GetData()
    close(iord)         ! Close input file

    ! Write results in TecPlot format
    if (WriteTecPlot) then
        FileAnim = trim(fName) // "example_anim.dat"    ! Animation data (TecPlot)
        open(iow1, file = FileAnim, status = 'unknown')
    endif

    open(iow03, file = 'EnergyPlot.dat', status = 'unknown')
    open(iow04, file = 'MomentumPlot.dat', status = 'unknown')
    open(iow05, file = 'ContforcPlot.dat', status = 'unknown')

    call Initial()

    call SetGridData()

    ! Create mutiple grid for contact
    if(contact) call SetContact_GridNodeData()

    call SetDT()

    call statinfo()

  end subroutine InputPara

  subroutine GetData()
!------------------------------------------------------------------
!-  purpose: Get input data using FFI module                      -
!-           data are stored in ParticleData module               -
!------------------------------------------------------------------
    use FFI
    implicit none
    
    integer:: key, i, tempI

    integer, parameter:: nbkw = 29
    character(4), parameter:: kw(nbkw) = (/&
        'endi','mpmc','nbmp','endt','grid','spx ', &
        'spy ','spz ','dcel','dtsc','outt','rptt', &
        'fixs','nmat','mate','part','load','velo', &
        'outr','curv','pt2d','curx','gimp','cont', &
        'nbco','nbbo','fixf','cubi','curt'/)

    do while(.true.)
       key = keyword(kw,nbkw)
       select case(key)
       
       case(1)    ! end input
          exit    ! Terminates execution of Do loop

       case(2)    ! mpm3 (title)
          call GetString(Title)

       case(3)    ! nbmp - number of material particles
          nb_particle = GetInt()
          write(*,"(a,i12)") 'Number of particles = ', nb_particle
          write(iomsg,"(a,i12)") 'Number of particles = ', nb_particle
          allocate(particle_list(nb_particle))
          call InitParticle()

       case(4)    ! endtime - Simulation end time
          EndTime = GetReal()
          OutTime = EndTime

       case(5)    ! grid - Define computational grid
          SpanX(1) = GetReal()
          SpanX(2) = GetReal()

          SpanY(1) = GetReal()
          SpanY(2) = GetReal()

          SpanZ(1) = GetReal()
          SpanZ(2) = GetReal()

       case(6)    ! spx - Define computational grid in x direction
          SpanX(1) = GetReal()
          SpanX(2) = GetReal()

       case(7)    ! spy - Define computational grid in y direction
          SpanY(1) = GetReal()
          SpanY(2) = GetReal()

       case(8)    ! spz - Define computational grid in z direction
          SpanZ(1) = GetReal()
          SpanZ(2) = GetReal()

       case(9)    ! dcel - Define computational grid space
          DCell = GetReal()
          if (abs(DCell) .le. 1.0e-15) &
              write(*,"(a50)") '***Warning*** Computational grid is too small'

       case(10)   ! DTScale - Define scale for time step size
          DTScale = GetReal()
          write(*,"(a,f10.3)") 'DTScale = ', DTScale
          write(iomsg,"(a,f10.3)") 'DTScale = ', DTScale
          if (DTScale .gt. 2.5) then
              call ErrorMsg()
              stop '*** Error *** DTScale too large'
          end if

       case(11)   ! outTime - Define time interval for writting plot file
          OutTime = GetReal()

       case(12)   ! ReportTime - Define time interval for reporting status
          ReportTime = GetReal()

       case(13)   ! fixs - Define fixed rigid plane for solid
          do i = 1, 6
              FixS(i) = GetInt()
          end do

       case(14)   ! nmat - Number of material sets
          nb_mat = GetInt()
          write(*,"(a,i12)") 'Number of material sets  = ', nb_mat
          write(iomsg,"(a,i12)") 'Number of material sets  = ', nb_mat
          allocate(mat_list(nb_mat))
          call InitMaterial()

       case(15)   ! material - Read in material properties
          call SetMaterial()

       case(16)   ! particle - Read in particles
          call SetParticle()

       case(17)   ! load - Define external load
          call SetLoad()

       case(18)   ! velo - Read in initial velocity
          call SetVelocity()

       case(19)   ! outres - Define variables which will be saved to plot file        
          nAnimate = nAnimate + 1
          AnimOption(nAnimate) = SetResOption()

       case(20)   ! curve - Define variables which will be saved to plot file
          nCurves = nCurves + 1
          CurveOption(nCurves) = SetResOption()
          if (nb_word .gt. 0) then
              tempI = GetInt()
              if (tempI.le.0 .OR. tempI.gt.nb_particle) then
                  call ErrorMsg()
                  stop '*** Error *** Invalid curved particle ID!'
              else
                  CurvePoint(nCurves) = tempI
              end if
          else    ! No curve point specified
              CurvePoint(nCurves) = 1
          end if
          
       case(21)   ! pt2d
          plot2dTrue = .true.
          plot2d(1) = GetReal()
          plot2d(2) = GetReal()
          plot2d(3) = GetReal()
          plot2d(4) = GetReal()
          plot2d(5) = GetReal()
          plot2d(6) = GetReal()

       case(22)   ! curx - curve x y z
          call SetCurX()
          
       case(23)   ! GIMP
          GIMP = .true.
          InfluInit = -1
          InfluEnd = 2
          Offset = 2
          write(*,"(a)") 'GIMP shape function is used'
          write(iomsg,"(a)") 'GIMP shape function is used'

       case(24)   ! contact ( = 0 by default )
          contact = .true.      ! using contact method
          if(contact)then
             call Setcontact()
          end if

       case(25)   ! nbco - number of components
          nb_component = GetInt() 
          write(*,"(a,i12)") 'Number of components = ', nb_component
          write(iomsg,"(a,i12)") 'Number of components = ', nb_component

       case(26)   ! nbbody - number of bodies
          nb_body = GetInt()
          write(*,"(a,i12)") 'Number of bodies = ', nb_body
          write(iomsg,"(a,i12)") 'Number of bodies = ', nb_body
          allocate(body_list(nb_body))
          call InitBody()

       case(27)   ! fixf - Define fixed rigid plane for fluid
          do i = 1, 6
             FixF(i) = GetInt()
          end do
          
       case(28)   ! Cubic B spline shape function
          CUBI = .true.
          InfluInit = -1
          InfluEnd = 2
          Offset = 2
          write(*,"(a)") 'Cubic B spline shape function is used'
          write(iomsg,"(a)") 'Cubic B spline shape function is used'
          
       case(29)   ! CurveTime - Define time interval for curve status
          CurvTime = GetReal()
          
       case default ! error
          stop 'STOP - Error encountered in reading data'

       end select
    end do

  end subroutine GetData

  subroutine Initial()
!------------------------------------------------------------------
!-  purpose: Initialize specific variable used after reading all  -
!-           particle data                                        -
!-    including: VOL, sig_y,  ie                               -
!------------------------------------------------------------------
    use FFI
    implicit none

    integer:: b, p, i, m
    integer:: parBegin, parEnd
    real(8), dimension(3):: Grav = (/0.0, 0.0, -9.8/)
    real(8):: sig1, sig2, sig3
    real(8):: volume

    type(particle), pointer:: pt
    type(particle), pointer:: pt1
    type(body), pointer:: bd

    loadrate = Grav / loadtime

    do b = 1, nb_body
       bd => body_list(b)
       m = bd%mat
       parBegin = bd%par_begin
       parEnd = bd%par_end
       !$OMP PARALLEL DO PRIVATE(pt,pt1,i,volume) SCHEDULE(STATIC)
       do p = parBegin, parEnd
          pt => particle_list(p)
          pt%sig_y = mat_list(m)%Yield0
          pt%Nf = mat_list(m)%porosity
          pt%fract = 1.0 - mat_list(m)%porosity
          pt%Kf = mat_list(m)%HydrCon
          pt%Den_F = mat_list(m)%Density_F
          ! particle depth
          !do i = 1, nb_particle  
          !    pt1 => particle_list(i)
          !    if ((pt%Xp(1)-pt1%Xp(1)) < epsilon(1.0_8)) then
          !        if (pt%Xp(3) <= pt1%Xp(3)) then  
          !            pt%H = pt1%Xp(3)
          !        end if
          !    end if
          !end do    ! p
          ! PorePressure initialization
          !if((pt%Xp(3)-pt%H) < epsilon(1.0_8)) pt%PorePress = 0.0
          pt%PorePress = 1.0d4
          ! stress initialization
          !sig3 = - (mat_list(m)%Density - mat_list(m)%Density_F)*9.8*(pt%H - pt%Xp(3))
          !sig1 = sig3*0.3 / (1 - 0.3)
          !sig2 = sig3*0.3 / (1 - 0.3)
          !pt%Press = (sig1 + sig2 + sig3) / 3
          
          volume = DCell**3 / 8       ! particle volume 4 material points per grid cell
          pt%Volu = volume
          pt%Vol = volume
          pt%mass_F = pt%Nf * volume * mat_list(m)%Density_F
          pt%mass_FF = volume * mat_list(m)%Density_F
          pt%mass = pt%fract * volume * mat_list(m)%Density
          pt%dsp = (pt%fract * volume)**(1.0/3.0)

       end do
       !$OMP END PARALLEL DO
    end do
 
  end subroutine Initial

  subroutine SetDT()
!------------------------------------------------------------------
!-  purpose:  Computational time step size                        -
!------------------------------------------------------------------
    implicit none
    
    integer:: m
    real(8):: cp, nu, E, ro

    DT = 1.0e-1

    do m = 1, nb_mat      
       E = mat_list(m)%Young
       nu = mat_list(m)%Poisson
       ro = mat_list(m)%Density
       cp = sqrt(E*(1-nu)/(1+nu)/(1-2*nu)/ro)   ! sound speed
       DT = min(DT, DCell/cp)
    end do
    DT = 0.00001   ! DT * DTScale

  end subroutine SetDT

  subroutine statinfo()
!------------------------------------------------------------------
!-   purpose: Output statistical data to message file             -
!------------------------------------------------------------------
    use FFI
    implicit none
    
    integer:: i, ibody, parBegin, parEnd
    real(8):: cmass, cke, cie, tmass=0, tke=0, tie=0
    type(Particle), pointer:: pt

    write(iomsg,*)

    do ibody = 1, nb_body
        write(iomsg,"(a13,i3,a4)") '--- body', ibody, ' ---'

        ! write mass and energy information
        cmass = 0    ! mass
        cke = 0      ! kinetic energy
        cie = 0      ! internal energy
        
        parBegin = body_list(ibody)%par_begin
        parEnd = body_list(ibody)%par_end
        !$OMP PARALLEL DO PRIVATE(pt) &
        !$OMP REDUCTION(+:cmass, cke, cie) SCHEDULE(STATIC)
        do i = parBegin, parEnd
            pt => particle_list(i)
            cmass = cmass + pt%mass
            cke = cke + dot_product(pt%VXp, pt%VXp)*pt%mass*0.5d0
            cie = cie + pt%ie
        end do
        !$OMP END PARALLEL DO
        tmass = tmass + cmass
        tke = tke + cke
        tie = tie + cie
        write(iomsg,"(a25,e12.4)"), ' mass:', cmass
        write(iomsg,"(a25,e12.4)") '  kinetic energy:', cke
        write(iomsg,"(a25,e12.4)") ' internal energy:', cie
    end do

    ! write total infomation
    write(iomsg,"(a20)") '--- total ---'
    write(iomsg,"(a25,e12.4)") 'total mass:', tmass
    write(iomsg,"(a25,e12.4)") ' total kinetic energy:', tke
    write(iomsg,"(a25,e12.4)") 'total internal energy:', tie

  end subroutine statinfo
  
  subroutine SetMaterial()
!------------------------------------------------------------------
!-    purpose: Get material data using FFI module                 -
!-             data are stored in mat_list                        -
!------------------------------------------------------------------
    use FFI
    implicit none

    integer:: i, t, matcount
    integer, parameter:: nbkw = 6
    character(4), parameter:: kw(nbkw) = (/&
        'elas','pla1','pla2','dpm ','dpr ','rheo'/)

    if (nb_mat.eq.0 .OR. nb_particle.eq.0) then
        stop '*** Error *** nb_mat/nb_particle must be defined in advance!'
    end if

    write(iomsg,"(a)") 'Material:'
    i = 0
    do while(i.lt.nb_mat)
        i = GetInt()
        write(*,*) trim(sss)
        write(iomsg,"(a)") trim(sss)

        if (i.gt.nb_mat) then
            call ErrorMsg()
            print *, '*** Error *** Too many material sets'
            print *, 'required : ',nb_mat
            print *, 'inputed  : ',i
            stop
        end if
        
        t = KeyWord(kw,nbkw)
        select case(t)
        
        case(1) ! elastic
            mat_list(i)%MatType = 1
            mat_list(i)%Density = GetReal()
            mat_list(i)%Young = GetReal()
            mat_list(i)%Poisson = GetReal()
        
        case(2) ! pla1: elastic-perfectly plastic
            mat_list(i)%MatType = 2
            mat_list(i)%Density = GetReal()
            mat_list(i)%Young = GetReal()
            mat_list(i)%Poisson = GetReal()
            mat_list(i)%Yield0 = GetReal()
        
        case(3) ! pla2: isotropic hardening
            mat_list(i)%MatType = 3
            mat_list(i)%Density = GetReal()
            mat_list(i)%Young = GetReal()
            mat_list(i)%Poisson = GetReal()
            mat_list(i)%Yield0 = GetReal()
            mat_list(i)%TangMod = GetReal()
        
        case(4) ! drucker-prager model
            mat_list(i)%MatType = 4
            mat_list(i)%Density = GetReal()
            mat_list(i)%Young = GetReal()
            mat_list(i)%Poisson = GetReal()
            mat_list(i)%fai = GetReal()
            mat_list(i)%cohen = GetReal()
            mat_list(i)%psi = GetReal()
            mat_list(i)%ten_f = GetReal()
        
        case(5) ! drucker-prager model with rheology
            mat_list(i)%MatType = 4
            mat_list(i)%Density = GetReal()
            mat_list(i)%Young = GetReal()
            mat_list(i)%Poisson = GetReal()
            mat_list(i)%fai = GetReal()
            mat_list(i)%cohen = GetReal()
            mat_list(i)%psi = GetReal()
            mat_list(i)%ten_f = GetReal()
        
        case(6) ! rheological model
            mat_list(i)%MatType = 4
            mat_list(i)%Density = GetReal()
            mat_list(i)%Young = GetReal()
            mat_list(i)%Poisson = GetReal()
            mat_list(i)%fai = GetReal()
            mat_list(i)%cohen = GetReal()
            mat_list(i)%psi = GetReal()
            mat_list(i)%ten_f = GetReal()
        
        case default
            call ErrorMsg()
            stop '*** Error *** Invalid Material Type!'

        end select
    end do
    ! set CutOff value
    CutOff = 0.0
    matcount = 0
    do i = 1, nb_mat
        CutOff = mat_list(i)%density_F + CutOff
        matcount = matcount + 1
    end do
    CutOff = CutOff * DCell**3 * 1e-6 / matcount

    ! print *, 'Grid mass CutOff value: ', CutOff
    write(*,"(a,i12)") 'the effective number of material', matcount
    write(iomsg,"(a,i12)") 'the effective number of material', matcount
    write(*,"(a,e12.4)") 'Grid mass CutOff value: ', CutOff
    write(iomsg,"(a,e12.4)") 'Grid mass CutOff value: ', CutOff

  end subroutine SetMaterial

  subroutine SetParticle()
!------------------------------------------------------------------
!-  purpose: Get particle data using FFI module                   -
!-           data are stored in particle_list                     -
!------------------------------------------------------------------
    use FFI
    implicit none

    integer:: p, newi, pnum
    integer:: matID

    if (nb_particle .eq. 0) then
        stop '*** Error *** nbmp must be defined in advance!'
    end if

    ! particle (comID, pnum / pID, matID, x, y, z)
    bodyCounter = bodyCounter + 1
    body_list(bodyCounter)%comID = GetInt()

    pnum = GetInt()     ! number of particle in the body
    body_list(bodyCounter)%par_begin = parCounter + 1
    body_list(bodyCounter)%par_end = parCounter + pnum
    do p = 1, pnum
        newi = GetInt() ! particle ID
        parCounter = parCounter + 1     ! count the particle number
        if (newi .lt. parCounter) then
            call ErrorMsg()
            stop 'particle ID conflict'
        end if

        if (parCounter .gt. nb_particle) then
            call ErrorMsg()
            print *, '*** Error *** Too many particles'
            print *, 'required : ',nb_particle
            print *, 'inputed  : ',parCounter
            stop
        end if

        matID = GetInt()  ! material ID of particle

        if (matID .gt. nb_mat .or. matID .le. 0) then
            call ErrorMsg()
            stop '*** Error *** Improper material ID'
        end if

        if (matID.gt.nb_mat) then
            call ErrorMsg()
            stop '*** Error *** material ID greater than nb_mat'
        end if

        particle_list(parCounter)%Xp(1) = GetReal()
        particle_list(parCounter)%Xp(2) = GetReal()
        particle_list(parCounter)%Xp(3) = GetReal()

        if (plot2dTrue) then
            call setskip(parCounter)
        end if
    end do
    
    body_list(bodyCounter)%mat = matID

    if (bodyCounter == nb_body) then
        write(*,"(a,i12)") 'particles read:', parCounter
    end if 

  end subroutine SetParticle

  subroutine setskip(p)
!------------------------------------------------------------------
!-             purpose: set display range                         -
!------------------------------------------------------------------
    implicit none
    
    integer:: p

    if (particle_list(p)%Xp(1)<plot2d(1) .or. &
        particle_list(p)%Xp(1)>plot2d(2)) then
        particle_list(p)%SkipThis = .true.
    else if (particle_list(p)%Xp(2)<plot2d(3) .or. &
             particle_list(p)%Xp(2)>plot2d(4)) then
        particle_list(p)%SkipThis = .true.
    else if (particle_list(p)%Xp(3)<plot2d(5) .or. &
             particle_list(p)%Xp(3)>plot2d(6)) then 
        particle_list(p)%SkipThis = .true.
    end if

  end subroutine setskip

  subroutine SetLoad()
!------------------------------------------------------------------
!-   purpose: Get Load by component or by node                    -
!------------------------------------------------------------------
    use FFI
    implicit none

    integer:: k, inode, i, ibody
    real(8):: fxp, fyp, fzp
    real(8):: pxp_F, pyp_F, pzp_F
    integer, parameter:: nbkw = 4
    character(4), parameter:: kw(nbkw) = &
        (/'endl','node','body','grav'/)

    if (nb_body.eq.0) then
        stop '*** Error *** nbby must be defined !'
    end if

    do while(.true.)
        k = keyword(kw,nbkw)
        select case(k)
        case(1)    ! end load
            exit

        case(2)    ! by node 
            inode = GetInt()
            particle_list(inode)%FXp(1) = GetReal()
            particle_list(inode)%FXp(2) = GetReal()
            particle_list(inode)%FXp(3) = GetReal()
            particle_list(inode)%PXp_F(1) = GetReal()
            particle_list(inode)%PXp_F(2) = GetReal()
            particle_list(inode)%PXp_F(3) = GetReal()
        
        case(3)    ! by body
            ibody = GetInt()    ! body number
            fxp = GetReal()
            fyp = GetReal()
            fzp = GetReal()
            pxp_F = GetReal()
            pyp_F = GetReal()
            pzp_F = GetReal()
            do i = body_list(ibody)%par_begin, body_list(ibody)%par_end
                particle_list(i)%FXp(1) = fxp
                particle_list(i)%FXp(2) = fyp
                particle_list(i)%FXp(3) = fzp
                particle_list(i)%PXp_F(1) = pxp_F
                particle_list(i)%PXp_F(2) = pyp_F
                particle_list(i)%PXp_F(3) = pzp_F
            end do
        
        case(4) ! gravity
            ibody = GetInt()    ! body number
            body_list(ibody)%Gravp(1) = GetReal()
            body_list(ibody)%Gravp(2) = GetReal()
            body_list(ibody)%Gravp(3) = GetReal()
            
        case default    ! error
            call ErrorMsg()
            stop 'error GetLoad'

       end select
    end do

  end subroutine SetLoad

  subroutine SetVelocity()
!------------------------------------------------------------------
!-  purpose: Get Initial Velocity by component or by node         -
!------------------------------------------------------------------
    use FFI
    implicit none

    integer:: k, inode, ibody, cpl, i
    real(8):: vxp, vyp, vzp
    real(8):: vxp_F, vyp_F, vzp_F
    integer, parameter:: nbkw = 3
    character(4), parameter:: kw(nbkw) = (/'endv','node','body'/)

    if (nb_body.eq.0) then
        stop '*** Error *** nbby must be defined !'
    end if

    do while(.true.)
        k = keyword(kw,nbkw)
       
        select case(k)
        case(1)    ! endload
            exit

        case(2)    ! by node
            inode = GetInt()
            particle_list(inode)%VXp(1) = GetReal()
            particle_list(inode)%VXp(2) = GetReal()
            particle_list(inode)%VXp(3) = GetReal()
            particle_list(inode)%VXp_F(1) = GetReal()
            particle_list(inode)%VXp_F(2) = GetReal()
            particle_list(inode)%VXp_F(3) = GetReal()
            
        case(3)    ! by body
            ibody = GetInt()      ! body number
            vxp = GetReal()
            vyp = GetReal()
            vzp = GetReal()
            vxp_F = GetReal()
            vyp_F = GetReal()
            vzp_F = GetReal()
            do i = body_list(ibody)%par_begin, body_list(ibody)%par_end
                particle_list(i)%VXp(1) = vxp
                particle_list(i)%VXp(2) = vyp
                particle_list(i)%VXp(3) = vzp
                particle_list(i)%VXp_F(1) = vxp_F
                particle_list(i)%VXp_F(2) = vyp_F
                particle_list(i)%VXp_F(3) = vzp_F
            end do
            
        case default    ! error
            call ErrorMsg()
            stop 'error GetVelocity'
            
        end select
    end do

  end subroutine SetVelocity

  logical function SetOnOff()
!------------------------------------------------------------------
!-   purpose: Get On/Off switch                                   -
!------------------------------------------------------------------
    use FFI
    implicit none

    integer:: k
    integer, parameter:: nbkw = 2
    character(4), parameter:: kw(nbkw) = (/'on', 'off'/)

    k = keyword(kw,nbkw)

    select case(k)
    case(1)    ! on
        SetOnOff = .true.

    case(2)    ! off
        SetOnOff = .false.

    case default    ! error
        call ErrorMsg()
        stop 'switch should be ON/OFF'

    end select

  end function SetOnOff

  integer function SetResOption()
!------------------------------------------------------------------
!-    purpose: Get result option                                  -
!------------------------------------------------------------------
    use FFI
    implicit none
    
    integer:: k
    integer, parameter:: nbkw = 30
    character(4), parameter:: kw(nbkw) = (/&
        'engk','engi','mat ','seqv','epef',&
        'vexS','veyS','vezS','velS','disx',&
        'disy','disz','disp','Pres','sigx',&
        'sigy','sigz','SDyz','SDzx','SDxy',&
        'vexF','veyF','vezF','velF','Porp',&
        'volu','volf','poro','cont','coll'/)

    SetResOption = keyword(kw,nbkw)

    if (SetResOption.le.0 .OR. SetResOption.gt.nbkw) then
        call ErrorMsg()
        stop 'Result option error'
    end if

  end function SetResOption
  
  subroutine SetCurX()
!-------------------------------------------------------------------
!-  purpose:                                                       -
!-     find particle ID according to coordinates                   -
!-------------------------------------------------------------------
    use FFI
    implicit none

    integer:: pid, bid
    real(8):: tx, ty, tz, maxx

    if (parCounter .ne. nb_particle) then
        call ErrorMsg()
        stop 'curx should be used after particle input'
    end if

    nCurves = nCurves + 1
    CurveOption(nCurves) = SetResOption()

    tx = GetReal()
    ty = GetReal()
    tz = GetReal()

    do bid = 1, nb_body
        do pid = body_list(bid)%par_begin, body_list(bid)%par_end
            maxx = MAX(ABS(particle_list(pid)%Xp(1)-tx), &
                       ABS(particle_list(pid)%Xp(2)-ty), &
                       ABS(particle_list(pid)%Xp(3)-tz))

            if (maxx .le. (0.25*DCell)) then
                CurvePoint(nCurves) = pid
                exit
            end if
        end do
    end do
    
  end subroutine SetCurX

  subroutine ErrorPt()
!-------------------------------------------------------------------
!-                                                                 -
!-  purpose:                                                       -
!-     output error message and line number                        -
!-                                                                 -
!-------------------------------------------------------------------
    implicit none
    
    print *, '*** Error *** Too many particles'
    print *, 'required : ', nb_particle
    print *, 'inputed  : ', parCounter
    stop
    
  end subroutine ErrorPt

  subroutine Setcontact()
!------------------------------------------------------------------
!-    purpose: Get contact data using FFI module                  -
!------------------------------------------------------------------
    use FFI
    implicit none

    integer:: k, t
    integer,parameter:: nbkw = 2
    character(4), parameter:: kw(nbkw) = (/'lagr', 'others'/)

    k = keyword(kw,nbkw)
    
    select case(k)
    case(1)
       contact_type= 1     ! using Bardenhagens method
       fricfa = GetReal()  ! the frictional coefficient
       ! the computational method of contact normal
       normbody = GetInt()    
       write(*,"(a)") 'Notice: Bardenhagens contact method is used'
       write(iomsg,"(a)") 'Notice: Bardenhagens contact method is used'
       write(*,"(a,e10.3)") 'frictional factor = ', fricfa
       write(iomsg,"(a,e10.3)") 'frictional factor = ', fricfa

    case(2)            
       !using another contact method by users 
       contact_type= 2    
       fricfa = GetReal() ! the frictional coefficient
       ! the computational method of contact normal
       normbody = GetInt()    
       write(*,"(a)") 'Notice: another contact method is used'
       write(iomsg,"(a)") 'Notice: another contact method is used '
       write(*,"(a,e10.3)") 'frictional factor = ', fricfa
       write(iomsg,"(a,e10.3)") 'frictional factor = ', fricfa

    case default
       call ErrorMsg()
       stop 'Error set contact parameter'

    end select

  end subroutine Setcontact
  
end module DataIn