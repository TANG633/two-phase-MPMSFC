!-------------------------------------------------------------------
!-                                                                 -
!-   Result output procedures                                      -
!-                                                                 -
!-------------------------------------------------------------------
module DataOut

  use ParticleData
  use GridData

  integer:: nCvs = 0        ! 时间历程曲线的数据总点数
  integer:: nAnm = 0        ! 后处理数据文件的数据集总数

  integer, parameter:: nVariables = 30
  character(4), parameter:: OutputName(nVariables) = (/&
      'engk','engi','mat ','seqv','epef',&
      'vexS','veyS','vezS','velS','disx',&
      'disy','disz','disp','Pres','sigx',&
      'sigy','sigz','SDyz','SDzx','SDxy',&
      'vexF','veyF','vezF','velF','Porp',&
      'volu','volf','poro','cont','coll'/)

  logical:: WriteTecPlot = .true.

  real(8):: OutTime = 0.0     ! time interval for plot data output
  real(8):: CurvTime = 0.0    ! time interval for plot curves
  real(8):: ReportTime = 0.0  ! time interval for status report

  ! Number of curves defined, default 2 curves
  integer:: nCurves  = 2
  ! Number of animation variables defined
  integer:: nAnimate = 1

  ! Limited 20 curves
  integer, parameter:: MaxCurves = 20
  ! Limited 20 animation variables
  integer, parameter:: MaxAnim = 20

  ! variable of each curve 时间历程曲线的变量在outputname中的编号
  integer CurveOption(MaxCurves) /1, 2, 18*0/       ! 赋值
  ! animation variables后处理文件中 输出变量在outputname中的编号
  integer AnimOption(MaxAnim) /3, 19*0/         
  integer:: CurvePoint(MaxCurves) = 1

  real(8):: plot2d(6) = 0.0      ! plot less particles
  logical:: plot2dTrue = .false.

contains

  subroutine OutCurve()
!------------------------------------------------------------------
!-  purpose: output result data for plotting time-history curve   -
!------------------------------------------------------------------
    use FFI, only: iow2, iow03, iow04, iow05
    implicit none
    
    integer:: i    
    character(20):: istr
    
    type(Particle), pointer:: pt

    nCvs = nCvs + 1
    ! write variable name
    if (nCvs .eq. 1) then
        write(iow2,"(a5)",advance='no') 'Time '
        print *, 'ncurves=', nCurves
        do i = 1, nCurves
            write(istr,*) CurvePoint(i)
            istr = trim(adjustl(istr))
            write(iow2,"(a4,a8)",advance='no') OutputName(CurveOption(i)), istr
        end do
        write(iow2,*)
    end if

    write(iow2,"(e12.4)", advance='no') CurrentTime

    !$OMP PARALLEL DO PRIVATE(pt) SCHEDULE(STATIC)
    do i = 1, nCurves
       pt => particle_list(CurvePoint(i))

       select case(CurveOption(i))
       case(1)
           write(iow2,"(e12.4)", advance='no') EngKinetic
       case(2)
           write(iow2,"(e12.4)", advance='no') EngInternal
       case(4)
           write(iow2,"(e12.4)", advance='no') pt%seqv
       case(5)
           write(iow2,"(e12.4)", advance='no') pt%epeff
       case(6)
           write(iow2,"(e12.4)", advance='no') pt%VXp(1)
       case(7)
           write(iow2,"(e12.4)", advance='no') pt%VXp(2)
       case(8)
           write(iow2,"(e12.4)", advance='no') pt%VXp(3)
       case(9)
           write(iow2,"(e12.4)", advance='no') sqrt(DOT_PRODUCT(pt%VXp, pt%VXp))
       case(10)
           write(iow2,"(e12.4)", advance='no') pt%Disp_P(1)
       case(11)
           write(iow2,"(e12.4)", advance='no') pt%Disp_P(2)
       case(12)
           write(iow2,"(e12.4)", advance='no') pt%Disp_P(3)
       case(13)
           write(iow2,"(e12.4)", advance='no') sqrt(DOT_PRODUCT(pt%Disp_P, pt%Disp_P))
       case(14)
           write(iow2,"(e12.4)", advance='no') pt%Press
       case(15)
           write(iow2,"(e12.4)", advance='no') pt%SDxx + pt%Press
       case(16)
           write(iow2,"(e12.4)", advance='no') pt%SDyy + pt%Press
       case(17)
           write(iow2,"(e12.4)", advance='no') pt%SDzz + pt%Press
       case(18)
           write(iow2,"(e12.4)", advance='no') pt%SDyz   
       case(19)
           write(iow2,"(e12.4)", advance='no') pt%SDzx
       case(20)
           write(iow2,"(e12.4)", advance='no') pt%SDxy
       case(21)
           write(iow2,"(e12.4)", advance='no') pt%VXp_F(1) 
       case(22)
           write(iow2,"(e12.4)", advance='no') pt%VXp_F(2)
       case(23)
           write(iow2,"(e12.4)", advance='no') pt%VXp_F(3)
       case(24)
           write(iow2,"(e12.4)", advance='no') sqrt(DOT_PRODUCT(pt%VXp_F, pt%VXp_F))
       case(25)
           write(iow2,"(e12.4)", advance='no') pt%PorePress
       case(26)
           write(iow2,"(e12.4)", advance='no') pt%vol
       case(27)
           write(iow2,"(e12.4)", advance='no') pt%Nf*pt%vol
       case(28)
           write(iow2,"(e12.4)", advance='no') pt%Nf
       case(29)
           write(iow2,"(e12.4)", advance='no') pt%SolidPress_con
       case(30)
           write(iow2,"(e12.4)", advance='no') pt%SolidPress_coll
       end select
       
    end do
    !$OMP END PARALLEL DO
    write(iow2,*)

    ! add code for plot energy and momentum
    write(iow03,"(4e12.4)") CurrentTime, EngKinetic+EngInternal, &
                            EngKinetic, EngInternal
    write(iow04,"(4e12.4)") CurrentTime, Mombody1(1), Mombody1(2), &
                            Mombody1(3)

    if (contact) then
        ! output contact force
        write(iow05,"(4e12.4)") CurrentTime, tot_cont_for(1), &
                                tot_cont_for(2), tot_cont_for(3) 
    end if

  end subroutine OutCurve

  subroutine OutAnim()
!------------------------------------------------------------------
!-  purpose: output result data for animation                     -
!------------------------------------------------------------------
    use FFI, only: iow1, iomsg
    implicit none

    integer:: p, i, b, parBegin, parEnd ! loop counter
    integer:: fCounter                  ! failure particle counter
    integer:: matID                     ! material number
    
    type(Particle), POINTER :: pt

    fCounter = 0

    nAnm = nAnm + 1
    if (nAnm .eq. 1) then
       write(iow1,10) Title
       write(iow1,20) ('"', OutputName(AnimOption(i)), '"', i=1,nAnimate)
    end if
    
10  format('TITLE = "', a60, '"')
20  format('VARIABLES= "X"   "Y"   "Z"  ', 50(A2, A4, A1))

    write(iow1,"(a12,e10.3,a1)") 'ZONE T="time', CurrentTime, '"'
    do b = 1, nb_body
       parBegin = body_list(b)%par_begin
       parEnd = body_list(b)%par_end
       matID = body_list(b)%mat
       !$OMP PARALLEL DO PRIVATE(pt) SCHEDULE(STATIC)
       do p = parBegin, parEnd
          pt => particle_list(p)
          ! omit particles out of the computational region
          if (pt%icell(1)<0 .or. pt%icell(2)<0 .or.pt%icell(3)<0) cycle
          if (pt%SkipThis) cycle

          write(iow1, "(3e12.4)", advance='no') pt%Xp

          do i = 1, nAnimate
              select case(AnimOption(i))
                
             case(3) !mat
                 write(iow1,"(i3)", advance='no') matID
             case(4) !seqv
                 write(iow1,"(e12.4)", advance='no') pt%seqv
             case(5) !epeff
                 write(iow1,"(e12.4)", advance='no') pt%epeff
             case(6) !vx
                 write(iow1,"(e12.4)", advance='no') pt%VXp(1)
             case(7) !vy
                 write(iow1,"(e12.4)", advance='no') pt%VXp(2)
             case(8) !vz
                 write(iow1,"(e12.4)", advance='no') pt%VXp(3)
             case(9) !Pvel
                 write(iow1,"(e12.4)", advance='no') sqrt(DOT_PRODUCT(pt%VXp, pt%VXp))
             case(10) !disx
                 write(iow1,"(e12.4)", advance='no') pt%Disp_P(1)
             case(11) !disy
                 write(iow1,"(e12.4)", advance='no') pt%Disp_P(2)
             case(12) !disz
                 write(iow1,"(e12.4)", advance='no') pt%Disp_P(3)
             case(13) !disp
                 write(iow1,"(e12.4)", advance='no') sqrt(DOT_PRODUCT(pt%Disp_P, pt%Disp_P)) 
             case(14) !soild pressure
                 write(iow1,"(e12.4)", advance='no') pt%Press
             case(15) !sigx
                 write(iow1,"(e12.4)", advance='no') pt%SDxx + pt%Press
             case(16) !sigy
                 write(iow1,"(e12.4)", advance='no') pt%SDyy + pt%Press
             case(17) !sigz
                 write(iow1,"(e12.4)", advance='no') pt%SDzz + pt%Press
             case(18) !SDyz
                 write(iow1,"(e12.4)", advance='no') pt%SDyz
             case(19) !SDzx
                 write(iow1,"(e12.4)", advance='no') pt%SDzx
             case(20) !SDxy
                 write(iow1,"(e12.4)", advance='no') pt%SDxy
             case(21)
                 write(iow1,"(e12.4)", advance='no') pt%VXp_F(1)
             case(22)
                 write(iow1,"(e12.4)", advance='no') pt%VXp_F(2)
             case(23)
                 write(iow1,"(e12.4)", advance='no') pt%VXp_F(3)
             case(24)
                 write(iow1,"(e12.4)", advance='no') sqrt(DOT_PRODUCT(pt%VXp_F, pt%VXp_F))
             case(25) !pore pressure
                 write(iow1,"(e12.4)", advance='no') pt%PorePress
             case(26) !vol
                 write(iow1,"(e12.4)", advance='no') pt%vol
             case(27)
                 write(iow1,"(e12.4)", advance='no') pt%Nf*pt%vol
             case(28)
                 write(iow1,"(e12.4)", advance='no') pt%Nf
             case(29)
                 write(iow1,"(e12.4)", advance='no') pt%SolidPress_con
             case(30)
                 write(iow1,"(e12.4)", advance='no') pt%SolidPress_coll
             end select
             
          end do ! i
          write(iow1,*)
       end do ! p
       !$OMP END PARALLEL DO
    end do    ! b

  end subroutine OutAnim
  
end module DataOut