! ------------------------------------------------------------------
! -                                                                -
! -  Background grid procedures                                    -
! -                                                                -
! -  NOTE: only 8 node cube element is available                   -
! -        only cuboid computational region is available           -
! -  computational region is:                                      -
! -   [SpanX(1),SpanX(2)]*[SpanY(1),SpanY(2)]*[SpanZ(1),SpanZ(2)]  -
! -  Grid,Cell(element) infomation will be computed by program     -
! -  GridNode Numbering is keypoint of searching algorithm         -
! ------------------------------------------------------------------
module GridData

  type CellElement
     real(8):: Cg(3)     ! Cell element centre coordinate
     integer:: index     ! 0 for empty / 1 for with particle
     integer:: pnum      ! Number of particle in cell element
     integer:: Cell_par(100)
  end type CellElement

  type GridNode
     real(8):: Xg(3)     ! grid node coordinate
     integer:: NumNode(3)
     
     real(8):: Mg        ! mass on grid node
     real(8):: PXg(3)    ! momentum on grid node
     real(8):: FXg(3)    ! internal/external force on gride node

     real(8):: Mg_F      ! mass of liquid on grid node
     real(8):: Mg_FF
     real(8):: PXg_F(3)  ! momentum of fluid on grid node
     real(8):: FXg_F(3)  ! internal/external force of fluid on gride node

     real(8):: Qi
     real(8):: Hi
     real(8):: ai
     
     real(8):: deltx(3)
     logical:: interp
     logical:: isboundary
     logical:: Fix_x, Fix_y, Fix_z      ! BC
     logical:: Fix_xf, Fix_yf, Fix_zf   ! BC
  end type GridNode

  type ContactGridNodeProperty
     ! the normal direction of contact grid node
     real(8):: ndir(3)
     ! the tangential unit vetors of contact grid node
     real(8):: sdir(3)  
  end type ContactGridNodeProperty

  type(CellElement), target, allocatable:: cell_list(:,:,:)
  type(GridNode), target, allocatable:: node_list(:,:,:,:)
  type(ContactGridNodeProperty), target, allocatable:: CP_list(:,:,:,:)
  
  real(8):: fricfa = 0.0    ! the frictional coefficient
  integer:: normbody = 0    ! the flag of computaional normal
  ! the flag of contact type: 0-no, 1-langrange, 2-penalty
  integer:: contact_type = 0
  ! the total contact force between of bodies
  real(8):: tot_cont_for(3)
  ! computational region
  real(8):: SpanX(2), SpanY(2), SpanZ(2)    

  real(8):: DCell = 0.0       ! Grid node interval
  ! Number of cells
  integer:: NumCell=0, &
            NumCellx=0, NumCelly=0, NumCellz=0
  ! Number of nodes
  integer:: NGx, NGy, NGz
  integer:: nb_gridnode = 0    ! number of gridnodes
  
  real(8):: iJacobi, iDCell
  real(8):: CutOff = 0.0       ! Grid mass cutoff value
  ! boundary condition
  integer:: FixS(6)
  integer:: FixF(6)

contains

  subroutine SetGridData()
! -----------------------------------------------------------------
! - Purpose                                                       -
! -    Create computational grid                                  -
! -----------------------------------------------------------------
    use ParticleData
    use FFI
    implicit none

    integer:: i, ix, iy, iz        ! loop counter
    real(8):: spx, spy, spz

    ! backgroud grid setting
    spx = SpanX(2) - SpanX(1)
    spy = SpanY(2) - SpanY(1)
    spz = SpanZ(2) - SpanZ(1)

    if (DCell.eq.0) then
       stop '*** Error *** DCell must be defined !'
    end if
    
    iDCell = 1.0 / DCell
    
    if (spx.le.0 .or. spy.le.0 .or. spz.le.0) then
       stop '*** Error *** SPX/SPY/SPZ must be defined !'
    end if
    ! Number of Cell
    NumCellx = floor(spx/DCell + 0.5)
    NumCelly = floor(spy/DCell + 0.5) + 1
    NumCellz = floor(spz/DCell + 0.5)

    SpanX(2) = SpanX(1) + NumCellx*DCell
    SpanY(2) = SpanY(1) + NumCelly*DCell
    SpanZ(2) = SpanZ(1) + NumCellz*DCell

    NumCell = NumCellx * NumCelly * NumCellz
    ! Number of node
    NGx = NumCellx + 1
    NGy = NumCelly + 1
    NGz = NumCellz + 1
    
    nb_gridnode = NGx*NGy*NGz
    
    node_list%NumNode(1) = NGx
    node_list%NumNode(2) = NGy
    node_list%NumNode(3) = NGz
    
    print *, 'Number of grid nodes in x axis = ', NGx
    print *, 'Number of grid nodes in y axis = ', NGy
    print *, 'Number of grid nodes in z axis = ', NGz
    print *, 'Number of grid nodes = ', nb_gridnode
    
    print *, 'Number of cell elements in x axis = ', NumCellx
    print *, 'Number of cell elements in y axis = ', NumCelly
    print *, 'Number of cell elements in z axis = ', NumCellz
    print *, 'Number of cell elements = ', NumCell
    write(iomsg,*)
    write(iomsg,"(a14,i10)"), 'Number of grid nodes = ', nb_gridnode
    write(iomsg,"(a14,i10)"), 'Number of cell elements = ', NumCell
    write(iomsg,"(a14,3i4)"), 'Nodes (x,y,z) ', NGx, NGy, NGz
    write(iomsg,"(a14,3i4)"), 'Cells (x,y,z) ', NumCellx, NumCelly, NumCellz

    allocate(node_list(nb_component, NGx, NGy,NGz))
    allocate(cell_list(NumCellx, NumCelly, NumCellz))
    
    if (contact) then
        CP_list%ndir(1) = 0.0d0    
        CP_list%ndir(2) = 0.0d0
        CP_list%ndir(3) = 0.0d0
    end if

    node_list%isboundary = .false.
    
    node_list%Fix_x = .false.
    node_list%Fix_y = .false.
    node_list%Fix_z = .false.
    
    node_list%Fix_xf = .false.
    node_list%Fix_yf = .false.
    node_list%Fix_zf = .false.

    ! create grid node info
    ! loop over every node
    !$OMP PARALLEL DO SIMD COLLAPSE(3) SCHEDULE(STATIC)
    do iz = 1, NGz
        do iy = 1, NGy
            do ix = 1, NGx
                ! Position of node
                node_list(:, ix, iy, iz)%Xg(1) = (ix - 1)*DCell + SpanX(1)
                node_list(:, ix, iy, iz)%Xg(2) = (iy - 1)*DCell + SpanY(1)
                node_list(:, ix, iy, iz)%Xg(3) = (iz - 1)*DCell + SpanZ(1)

                ! Solid boundary condition
                if ((ix==1 .and. FixS(1)==1).or.(ix==NGx .and. FixS(2)==1).or.&
                    (iy==1 .and. FixS(3)==1).or.(iy==NGy .and. FixS(4)==1).or.&
                    (iy==NGy-1 .and. FixS(4)==1).or.&
                    (iz==1 .and. FixS(5)==1).or.(iz==NGz .and. FixS(6)==1)) then
                    node_list(:, ix, iy, iz)%Fix_x = .true.
                    node_list(:, ix, iy, iz)%Fix_y = .true.
                    node_list(:, ix, iy, iz)%Fix_z = .true.
                end if

                if ((ix==1 .and. FixS(1)==2).or.(ix==NGx .and. FixS(2)==2)) then
                    node_list(:, ix, iy, iz)%Fix_x = .true.
                end if

                if ((iy==1 .and. FixS(3)==2).or.(iy==NGy .and. FixS(4)==2).or.&
                    (iy==NGy-1 .and. FixS(4)==2)) then
                    node_list(:, ix, iy, iz)%Fix_y = .true.
                end if

                if ((iz==1 .and. FixS(5)==2).or.(iz==NGz .and. FixS(6)==2)) then
                    node_list(:, ix, iy, iz)%Fix_z = .true.
                end if

                ! Fluid boundary condition
                if ((ix==1 .and. FixF(1)==1).or.(ix==NGx .and. FixF(2)==1).or.&
                    (iy==1 .and. FixF(3)==1).or.(iy==NGy .and. FixF(4)==1).or.&
                    (iy==NGy-1 .and. FixS(4)==1).or.&
                    (iz==1 .and. FixF(5)==1).or.(iz==NGz .and. FixF(6)==1)) then
                    node_list(:, ix, iy, iz)%Fix_xf = .true.
                    node_list(:, ix, iy, iz)%Fix_yf = .true.
                    node_list(:, ix, iy, iz)%Fix_zf = .true.
                end if

                if ((ix==1 .and. FixF(1)==2).or.(ix==NGx .and. FixF(2)==2)) then
                    node_list(:, ix, iy, iz)%Fix_xf = .true.
                end if

                if ((iy==1 .and. FixF(3)==2).or.(iy==NGy .and. FixF(4)==2).or.&
                    (iy==NGy-1 .and. FixF(4)==2)) then
                    node_list(:, ix, iy, iz)%Fix_yf = .true.
                end if

                if ((iz==1 .and. FixF(5)==2).or.(iz==NGz .and. FixF(6)==2)) then
                    node_list(:, ix, iy, iz)%Fix_zf = .true.
                end if
                
                ! Mark as boundary if necessary
                node_list(:, ix, iy, iz)%isboundary = node_list(:, ix, iy, iz)%Fix_x .or. &
                                                      node_list(:, ix, iy, iz)%Fix_y .or. &
                                                      node_list(:, ix, iy, iz)%Fix_z .or. &
                                                      node_list(:, ix, iy, iz)%Fix_xf .or. &
                                                      node_list(:, ix, iy, iz)%Fix_yf .or. &
                                                      node_list(:, ix, iy, iz)%Fix_zf
            end do
        end do
    end do
    !$OMP END PARALLEL DO SIMD

    ! create cell node info
    ! loop over every cell element
    !$OMP PARALLEL DO SIMD COLLAPSE(3) SCHEDULE(STATIC)
    do iz = 1, NumCellz
        do iy = 1, NumCelly
            do ix = 1, NumCellx
                ! Position of cell element
                cell_list(ix, iy, iz)%Cg(1) = (ix - 1)*DCell + SpanX(1) + 0.5*DCell
                cell_list(ix, iy, iz)%Cg(2) = (iy - 1)*DCell + SpanY(1) + 0.5*DCell
                cell_list(ix, iy, iz)%Cg(3) = (iz - 1)*DCell + SpanZ(1) + 0.5*DCell
            end do
        end do
    end do
    !$OMP END PARALLEL DO SIMD

  end subroutine SetGridData

  subroutine SetContact_GridNodeData()
! -----------------------------------------------------------------
! - Purpose                                                       -
! -    Create computational grid for simulating contact problem   -
! -----------------------------------------------------------------
    use ParticleData, only:nb_component
    use FFI
    implicit none

    allocate(CP_list(nb_component, NGx, NGy, NGz))

  end subroutine SetContact_GridNodeData

  function InWhichCell(xx, p) result(cell)
! -----------------------------------------------------------------
! - Purpose                                                       -
! -    Determine which cell the point (xx,yy,zz) is located in    -
! -                                                               -
! - Input                                                         -
! -    xx(3) - Coordinates of a point                             -
! -                                                               -
! - Out values                                                    -
! -    cell(ix, iy,iz) : Cell coordinate in which the point       -
!-                       is located (MP)                          -
! -    cell(-1, -1, -1) : point out off the cell                  -
! -----------------------------------------------------------------
    implicit none
    
    real(8), intent(in):: xx(3)
    integer, intent(in):: p
    real(8):: dx, dy, dz
    integer:: idx, idy, idz
    integer:: cell(3)   ! result
    real(8), parameter:: EPSILON = 1.0d-12
    
    type(CellElement), POINTER:: ce

    if (xx(1) < SpanX(1) - EPSILON .or. xx(1) > SpanX(2) + EPSILON .or. &
        xx(2) < SpanY(1) - EPSILON .or. xx(2) > SpanY(2) + EPSILON .or. &
        xx(3) < SpanZ(1) - EPSILON .or. xx(3) > SpanZ(2) + EPSILON) then
        cell = (/-1, -1, -1/)
        return  
    end if
    ! Precompute the distances to avoid recalculating multiple times
    dx = xx(1) - SpanX(1)
    dy = xx(2) - SpanY(1)
    dz = xx(3) - SpanZ(1)
    ! Calculate indices
    idx = floor(dx / DCell) + 1
    idy = floor(dy / DCell) + 1
    idz = floor(dz / DCell) + 1

    cell = (/idx, idy, idz/)
    
    if (idx > 0 .and. idx <= NumCellx .and. &
        idy > 0 .and. idy <= NumCelly .and. &
        idz > 0 .and. idz <= NumCellz) then
        ce => cell_list(idx, idy, idz)
        ce%index = 1
        ce%pnum = ce%pnum + 1
        ce%Cell_par(ce%pnum) = p
    end if
    
  end function InWhichCell
  
  function BaseNode(xx) result(base)
! -----------------------------------------------------------------
! - Purpose                                                       -
! -    Determine which base node that point (xx,yy,zz) influences -
! - Input                                                         -
! -    xx(3) - Coordinates of a point                             -
! -                                                               -
! - Out values                                                    -
! -    node(ix, iy,iz): base node coordinate in which the 
! -                     point located                             -
! -----------------------------------------------------------------
    implicit none
    
    real(8), intent(in):: xx(3)
    integer:: idx, idy, idz
    real(8):: dx, dy, dz
    integer:: base(3)   ! result
    
    ! Precompute the distances to avoid recalculating multiple times
    dx = xx(1) - SpanX(1)
    dy = xx(2) - SpanY(1)
    dz = xx(3) - SpanZ(1)

    ! Compute the base node indices based on cell size
    idx = floor(dx / DCell) + 1
    idy = floor(dy / DCell) + 1
    idz = floor(dz / DCell) + 1

    base(1) = idx
    base(2) = idy
    base(3) = idz

  end function BaseNode

  function NShape(IDcom, node1, xx) result(res)
!------------------------------------------------------------------
!-  Purpose                                                       -
!-      Evaluate the shape functions and their derivatives        -
!-      at particle p associated with nodes of the cell in        -
!-      which the particle p is located                           -
!-  Inputs                                                        -
!-      node1 - the base node  that the particle can influence    -
!-      xx    - particle position                                 -
!------------------------------------------------------------------
    implicit none

    integer, intent(in):: IDcom, node1(3)
    real(8), intent(in):: xx(3)

    integer:: i
    real(8):: x1(3), x2(3)
    real(8), dimension(3, 4):: sh, dh  ! Shape functions and their derivatives
    real(8), dimension(2, 3, 4):: res  ! Return array: res(1,:,:) -> sh, res(2,:,:) -> dh
    
    type(GridNode), POINTER:: gn

    gn => node_list(IDcom, node1(1), node1(2), node1(3))
    ! Initialize
    sh = 0.0
    dh = 0.0
    res = 0.0
    
    x1 = (xx - gn%Xg) * iDCell
    x2 = x1 - 1.0

    ! 1>___2
    !$OMP PARALLEL DO SIMD SCHEDULE(STATIC)
    do i = 1, 3
        ! node 1
        sh(i, 1) = 1.0d0 - abs(x1(i))
        dh(i, 1) = - sign(real(1.0, 8), x1(i)) * iDCell
        ! node 2
        sh(i, 2) = 1.0d0 - abs(x2(i))
        dh(i, 2) = - sign(real(1.0, 8), x2(i)) * iDCell
    end do
    !$OMP END PARALLEL DO SIMD
    ! Combine results into the return array
    res(1, :, :) = sh
    res(2, :, :) = dh

  end function NShape

  function NShape_GIMP(IDcom, node1, xx) result(res)
!------------------------------------------------------------------
!-  Purpose                                                       -
!-      Evaluate the shape functions and their derivatives        -
!-      at particle p associated with nodes of the cell in        -
!-      which the particle p is located                           -
!-  Inputs                                                        -
!-      node1 - the base node  that the particle can influence    -
!-      xx    - particle position                                 -
!------------------------------------------------------------------
    implicit none

    integer, intent(in):: IDcom, node1(3)
    real(8), intent(in):: xx(3)
    
    integer:: i
    real(8):: x1(3), x2(3), x3(3), x4(3), lp
    
    real(8), dimension(3, 4):: sh, dh  ! Shape functions and their derivatives
    real(8), dimension(2, 3, 4):: res  ! Return array: res(1,:,:) -> sh, res(2,:,:) -> dh
    
    type(GridNode), POINTER:: gn
    ! Initialize
    sh = 0.0
    dh = 0.0
    res = 0.0
    
    gn => node_list(IDcom, node1(1), node1(2), node1(3))

    x2 = xx - gn%Xg
    lp = 0.5 * DCell
    x1 = x2 + DCell
    x3 = x2 - DCell
    x4 = x2 - 2.0*DCell
    ! 1___2>___3___4
    !$OMP PARALLEL DO SIMD SCHEDULE(STATIC)
    do i = 1, 3
        ! node 1
        if (x1(i) <= (DCell + 0.5*lp)) then
            sh(i, 1) = ((DCell + 0.5 * lp - abs(x1(i)))**2) / (2.0 * DCell * lp)
            dh(i, 1) = - sign(real(1.0, 8), x1(i)) * (DCell + 0.5 * lp - abs(x1(i))) / (DCell * lp)
        else
            sh(i, 1) = 0.0d0
            dh(i, 1) = 0.0d0
        end if
        ! node 2
        if (x2(i) >= (DCell - 0.5*lp)) then
            sh(i, 2) = ((DCell + 0.5 * lp - abs(x2(i)))**2) / (2.0 * DCell * lp)
            dh(i, 2) = - sign(real(1.0, 8), x2(i)) * (DCell + 0.5 * lp - abs(x2(i))) / (DCell * lp)
        else if (x2(i) >= 0.5*lp .and. x2(i) < (DCell - 0.5*lp)) then
            sh(i, 2) = 1.0 - (abs(x2(i)) / DCell)
            dh(i, 2) = - sign(real(1.0, 8), x2(i)) / DCell
        else
            sh(i, 2) = 1.0 - ((4*(x2(i)**2) + lp**2) / (4.0 * DCell * lp))
            dh(i, 2) = - 2.0*x2(i) / (DCell * lp)
        end if
        ! node 3
        if (abs(x3(i)) >= (DCell - 0.5*lp)) then
            sh(i, 3) = ((DCell + 0.5 * lp - abs(x3(i)))**2) / (2.0 * DCell * lp)
            dh(i, 3) = - sign(real(1.0, 8), x3(i)) * (DCell + 0.5 * lp - abs(x3(i))) / (DCell * lp)
        else if (abs(x3(i)) >= 0.5*lp .and. abs(x3(i)) < (DCell - 0.5*lp)) then
            sh(i, 3) = 1.0 - (abs(x3(i)) / DCell)
            dh(i, 3) = - sign(real(1.0, 8), x3(i)) / DCell
        else
            sh(i, 3) = 1.0 - ((4*(x3(i)**2) + lp**2) / (4.0 * DCell * lp))
            dh(i, 3) = - 2.0*x3(i) / (DCell * lp)
        end if
        ! node 4
        if (abs(x4(i)) <= (DCell + 0.5*lp)) then
            sh(i, 4) = ((DCell + 0.5 * lp - abs(x4(i)))**2) / (2.0 * DCell * lp)
            dh(i, 4) = - sign(real(1.0, 8), x4(i)) * (DCell + 0.5 * lp - abs(x4(i))) / (DCell * lp)
        else
            sh(i, 4) = 0.0d0
            dh(i, 4) = 0.0d0
        end if
    end do
    !$OMP END PARALLEL DO SIMD
    ! Combine results into the return array
    res(1, :, :) = sh
    res(2, :, :) = dh
    
  end function NShape_GIMP

  function NShape_Cubic(IDcom, node1, xx) result(res)
!------------------------------------------------------------------
!-  Purpose                                                       -
!-      Evaluate the shape functions and their derivatives        -
!-      at particle p associated with nodes of the cell in        -
!-      which the particle p is located                           -
!-  Inputs                                                        -
!-      node1 - the base node that the particle can influence     -
!-      xx    - particle position                                 -
!------------------------------------------------------------------
    implicit none
    
    integer, intent(in):: IDcom, node1(3)
    real(8), intent(in):: xx(3)
    
    integer:: i
    real(8):: x1(3), x2(3), x3(3), x4(3)
    
    real(8), dimension(3, 4):: sh, dh  ! Shape functions and their derivatives
    real(8), dimension(2, 3, 4):: res  ! Return array: res(1,:,:) -> sh, res(2,:,:) -> dh
    
    type(GridNode), POINTER:: gn
    ! Initialize
    sh = 0.0
    dh = 0.0
    res = 0.0
    
    gn => node_list(IDcom, node1(1), node1(2), node1(3))
    
    x2 = (xx - gn%Xg) * iDCell
    x1 = x2 + 1.0
    x3 = x2 - 1.0
    x4 = x2 - 2.0
    !$OMP PARALLEL DO SIMD SCHEDULE(STATIC)
    do i = 1, 3
        if (node1(i) == 1) then
            ! node 1 (out of the boundary)
            sh(i, 1) = 0.0
            dh(i, 1) = 0.0
            ! node 2 (boundary)
            sh(i, 2) = ((x2(i)**3)/6.0) - x2(i) + 1.0
            dh(i, 2) = (0.5*(x2(i)**2) - 1.0) * iDCell
            ! node 2 (right side of closest boundary)
            sh(i, 3) = - ((x3(i)**3)/3.0) - (x3(i)**2) + (2.0/3.0)
            dh(i, 3) = (- (x3(i)**2) - 2.0*x3(i)) * iDCell
            
            sh(i, 4) = ((x4(i)**3)/6.0) + (x4(i)**2) + 2.0*x4(i) + (4.0/3.0)
            dh(i, 4) = (0.5*(x4(i)**2) + 2.0*x4(i) + 2.0) * iDCell
        else if (node1(i) == (gn%NumNode(i)-1)) then
            sh(i, 1) = - ((x1(i)**3)/6.0) + (x1(i)**2) - 2.0*x1(i) + (4.0/3.0)
            dh(i, 1) = (- 0.5*(x1(i)**2) + 2.0*x1(i) - 2.0) * iDCell
            
            ! node 2 (left side of closest boundary)
            sh(i, 2) = ((x2(i)**3)/3.0) - (x2(i)**2) + (2.0/3.0)
            dh(i, 2) = ((x2(i)**2) - 2.0*x2(i)) * iDCell
            ! node 3 (boundary)
            sh(i, 3) = - ((x3(i)**3)/6.0) + x3(i) + 1.0
            dh(i, 3) = (- 0.5*(x3(i)**2) + 1.0) * iDCell
            ! node 4 (out of boundary)
            sh(i, 4) = 0.0
            dh(i, 4) = 0.0
        else
            ! node 1
            sh(i, 1) =  - ((x1(i)**3)/6.0) + (x1(i)**2) - 2.0*x1(i) + (4.0/3.0)
            dh(i, 1) = (- 0.5*(x1(i)**2) + 2.0*x1(i) - 2.0) * iDCell
            ! node 2
            sh(i, 2) = 0.5*(x2(i)**3) - (x2(i)**2) + (2.0/3.0)
            dh(i, 2) = (1.5*(x2(i)**2) - 2.0*x2(i)) * iDCell
            ! node 3
            sh(i, 3) = - 0.5*(x3(i)**3) - (x3(i)**2) + (2.0/3.0)
            dh(i, 3) = (- 1.5*(x3(i)**2) - 2.0*x3(i)) * iDCell
            ! node 4
            sh(i, 4) = ((x4(i)**3)/6.0) + (x4(i)**2) + 2.0*x4(i) + (4.0/3.0)
            dh(i, 4) = (0.5*(x4(i)**2) + 2.0*x4(i) + 2.0) * iDCell
        end if
    end do  ! i
    !$OMP END PARALLEL DO SIMD
    ! Combine results into the return array
    res(1, :, :) = sh
    res(2, :, :) = dh
    
  end function NShape_Cubic
  
end module GridData