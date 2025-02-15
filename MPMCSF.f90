! ------------------------------------------------------------------
! -                                                                -
! -  Main procedures                                               -
! -                                                                -
! ------------------------------------------------------------------
program MPMCSF
  use omp_lib
  use ParticleData
  use FFI, only: iomsg, iow1, iow2
  use DataIn
  use DataOut

  implicit none

  real(8):: t_begin, t_end, t_cpu=0.0, TPZC
  real(8):: t_bg, t_ed, t_elapsed, t0, t1
  real(8):: t_sec, t_s
  integer:: t_min, t_m
  integer:: iplotstep = 1       ! Write results in TecPlot format

  real(8):: plt = 0.0           ! time interval for plot data output
  real(8):: pct = 0.0           ! time interval for plot Curves
  real(8):: prt = 0.0           ! time interval for status report

  call cpu_time( t_bg )         ! the initial time of preprocessing

  call InputPara()              ! Input data

  call calcEnergy()             ! Calculate kinetic energy

  write(iomsg,*)
  write(iomsg,"(a,e10.3)") 'DT = ', DT

  plt = plt + OutTime           ! time interval for plot data output time interval for plot data output
  pct = pct + CurvTime
  ! Write results in TecPlot format: initial step
  if (WriteTecPlot) call OutAnim()
  call OutCurve()
  iplotstep = iplotstep + 1
  
  call cpu_time( t_ed )         ! the end time of preprocessing
  print *, '** Time for preprocessing is', t_ed - t_bg, ' seconds'

  write(*,"(a,e10.3)") 'DT = ', DT
  write(*,*) 'solving...'

  call cpu_time( t_bg )
  t0 = secnds(0.0)              ! reference time

  ! Solving
  do while(CurrentTime .le. EndTime)
     call cpu_time( t_begin )
     istep = istep + 1          ! Current time step
     CurrentTime = CurrentTime + DT
     EngInternal = 0.0

     !if (CurrentTime <= loadtime) call GravLoad()
     
     call GridInformationInitial()
     !call BoundaryCellSearch()
     call GridNodalForce()
     call ApplyBoundaryConditions()

     call ParticleUpdate()
     call GridMomentumUpdate()
     call ApplyBoundaryConditions()

     call ParticleStrainUpdate()
     call CellAverage1()
     call ParticleStressUpdate()
     
     call GridMSLInformation()
     call GridCoefficient()
     call PorePressUpdate()

     call calcEnergy()          ! Calculate kinetic energy

     call cpu_time( t_end )
     ! time for one step calculation
     t_cpu = t_cpu + t_end - t_begin
     
     ! out put curve
     if (CurrentTime .ge. pct) then
         pct = pct + CurvTime
         call OutCurve()
     end if
     ! out put animation data
     if (CurrentTime .ge. plt) then
         plt = plt + OutTime
         write(*,*) 'Write output data'
         write(iomsg,*) 'Write output data'

         ! Write results in TecPlot format: current step
         if (WriteTecPlot) call OutAnim()
         iplotstep = iplotstep + 1
     end if

     ! report current computational progress
     if (CurrentTime .ge. prt) then
        prt = prt + ReportTime
  
        write(*, 50) t_cpu
        write(iomsg, 50) t_cpu
50      format("** The elapsed time : ", f12.4, " seconds")

        write(*,100) istep, CurrentTime, EngKinetic, EngKinetic + EngInternal
        write(iomsg,100) istep, CurrentTime, EngKinetic, EngKinetic + EngInternal
100     format(1x, "Step=", i10, 1p, "  T=", e10.3, "  K.E.=", e10.3, "  T.E.=", e10.3)

        write(*,300) Momentum(1), Momentum(2), Momentum(3)
        write(iomsg,300) Momentum(1), Momentum(2), Momentum(3)
300     format(18x,"Mx=", e10.3, "   My=", e10.3, "   Mz=", e10.3)
     end if

  end do

  call cpu_time( t_ed )
  t1 = secnds(t0)   ! the elapsed time

  t_elapsed = t_ed - t_bg

  t_min = floor(t_cpu/60)
  t_sec = t_cpu - t_min*60
  t_m = floor(t_elapsed/60)
  t_s = t_elapsed - t_m*60
  
  TPZC = t_cpu / nb_particle / istep

  write(*,200) t_cpu, t_min, t_sec
  write(*,201) t_elapsed, t_m, t_s
  write(iomsg,*)
  write(iomsg,200) t_cpu, t_min, t_sec
  write(iomsg,201) t_elapsed, t_m, t_s
  write(*,202) TPZC *1.0e9
  write(iomsg,202) TPZC*1.0e9
  write(*,203) t1
  write(iomsg,203) t1

  ! Close OutPut Files
  if (WriteTecPlot) close(iow1)  
  close(iow2)
  close(iomsg)

200 format("** Total CPU Time is ", f12.4, " seconds (", i5, " minutes", f8.4, " seconds)")
201 format("** Elapsed time is   ", f12.4, " seconds (", i5, " minutes", f8.4, " seconds)")
202 format("** Time per particle cycle: ", f12.4, " nanoseconds")
203 format("** Total elasped time : ", f12.4, " seconds")

end program MPMCSF