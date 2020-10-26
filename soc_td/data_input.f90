module data_input

  use prec_mod
  implicit none
  
  !private

  !public  :: contrl, soc_ao, ci_coeff, mo_coeff, mo_ene

  type  contrl
    integer               :: n_singl, n_tripl
    integer               :: num_bov(3), tbcart_bov(3)
    integer,allocatable   :: root_s(:,:), root_t(:,:)
    real(dpr),allocatable :: ene_s(:), ene_t(:)
    character(len=10)     :: do_s0
    character(len=10)     :: do_dip
    character(len=20)     :: qm_flag
    real(dpr)             :: coeff_thresh
  end type
  
  type ao
    integer                   :: nc, nb
    integer, allocatable      :: nsh(:)
  end type

  type soc_ao
    real(dpr),allocatable     :: gto(:)    
    real(dpr),allocatable     :: tbgto(:)  
  end type

  type dip_ao
    real(dpr),allocatable     :: gto(:)    
    real(dpr),allocatable     :: tbgto(:)  
  end type

  type overlp_ao
    real(dpr),allocatable     :: gto(:)    !from third party QM code
    real(dpr),allocatable     :: tbgto(:)  !base matching output
    real(dpr),allocatable     :: molsoc(:) !from MolSOC 
    real(dpr),allocatable     :: molout(:) !base matching output
  end type

  type ci_coeff
    real(dpr),allocatable     :: mix_all(:) 
    real(dpr),allocatable     :: xpy(:)   
    real(dpr),allocatable     :: xmy(:)   
  end type

  real(dpr),allocatable       :: mo_coeff(:) 
  real(dpr),allocatable       :: mo_ene(:)   
  
end module data_input

