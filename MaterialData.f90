! ------------------------------------------------------------------
! -                                                                -
! -  Material parameters                                           -
! -                                                                -
! ------------------------------------------------------------------
module MaterialData

  type Material
     integer:: MatType    ! Material type
     real(8):: Density    ! Initial density
     real(8):: Young      ! Young's modulus
     real(8):: Poisson    ! Poisson's ratio
     real(8):: Yield0     ! Initial yield stress
     real(8):: TangMod    ! Tangential modulus
     ! Parameters for Drucker-Prager soil material
     real(8):: fai, cohen, psi, ten_f
     ! 应变软化参数
     real(8):: cohen_peak, cohen_res
     real(8):: fai_peak, fai_res
     real(8):: psi_int, psi_crit
     real(8):: Para_yita
     
     real(8):: Density_F
     real(8):: BulkMod_F
     real(8):: porosity
     real(8):: HydrCon
  end type Material
  
  integer:: nb_mat = 0    ! number of material sets
  type(Material), allocatable:: mat_list(:)

  logical:: Jaum = .true.    ! use Jaumann rate ?
  
contains
  
  subroutine InitMaterial()
!------------------------------------------------------------------
!    purpose: initialize particle_list                            -
!    used after particle_list space allocated                     -
!------------------------------------------------------------------
    implicit none
    
    integer :: i
    
    mat_list%Density = 0.0d0
    mat_list%Young = 0.0d0
    mat_list%Poisson = 0.0d0
    mat_list%Yield0 = 0.0d0
    mat_list%TangMod = 0.0d0
    mat_list%fai = 0.0d0
    mat_list%cohen = 0.0d0
    mat_list%psi = 0.0d0
    mat_list%ten_f = 0.0d0
    !应变软化初始化     
    mat_list%cohen_peak = 0.0d0
    mat_list%cohen_res = 0.0d0
    mat_list%fai_peak= 30.0
    mat_list%fai_res= 25.0
    mat_list%psi_int = 0.0d0
    mat_list%psi_crit = 0.0d0
    mat_list%Para_yita = 500

    mat_list%Density_F = 1.0d3
    mat_list%BulkMod_F = 2.2d9
    mat_list%porosity = 0.3d0
    mat_list%HydrCon = 0.001d0
    
  end subroutine InitMaterial  
  
end module MaterialData