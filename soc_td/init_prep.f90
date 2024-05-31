module init_prep
  
  use prec_mod
  use data_input

  implicit none

  private
  public read_input, read_soc_ao, read_dip_ao, read_overlp_ao
  public read_ci_coeff, read_mo, read_ao_orb
  public input, aoint, ao_orb, dipint, overlpint, td_coeff
  public t_dim, base_dim, ndim1
  
  integer    :: t_dim, base_dim, ndim1
  type(contrl)    :: input
  type(ao)        :: ao_orb
  type(overlp_ao) :: overlpint
  type(soc_ao)    :: aoint
  type(dip_ao)    :: dipint
  type(ci_coeff)  :: td_coeff
  
  contains

  subroutine read_input()
    
    open(unit=171, file='soc_td_input.dat', status='old')
    read(171, *) input%qm_flag
    read(171, *) input%do_dip
    read(171, *) input%do_s0
    read(171, *) input%coeff_thresh
    read(171, *) input%n_singl
    read(171, *) input%n_tripl

    allocate(input%root_s(input%n_singl,2))
    allocate(input%ene_s(input%n_singl))
    allocate(input%root_t(input%n_tripl,2))
    allocate(input%ene_t(input%n_tripl))

    read(171, *) input%root_s(:,1)
    read(171, *) input%root_s(:,2)
    read(171, *) input%ene_s(:)
    read(171, *) input%root_t(:,1)
    read(171, *) input%root_t(:,2)
    read(171, *) input%ene_t(:)
    read(171, *) input%num_bov(:)

    if (input%qm_flag == 'tddftb') then
      read(171, *) input%tbcart_bov(:) !for tddftb
      print *, "tbcart_bov", input%tbcart_bov(:)
    endif
    
    open(unit=16, file='ene_out.dat', status='replace') 
    write(16, 1705) input%ene_s(:)
    write(16, 1706) input%ene_t(:)
    1705 format('Singlet(eV):', 20f12.4)
    1706 format('Triplet(eV):', 20f12.4)
    1707 format(20f12.4)

    close(16)
    close(171)

  end subroutine read_input
  
  subroutine read_ao_orb()
    
    open(unit=178, file='ao_basis.dat', status='old')
    read(178, *) ao_orb%nc, ao_orb%nb 
    allocate(ao_orb%nsh(ao_orb%nc))
    read(178, *) ao_orb%nsh 
    
    close(178)
   
  end subroutine read_ao_orb
  
  subroutine read_soc_ao()
    
    if (input%qm_flag == 'tddftb') then
      allocate(aoint%tbgto(input%tbcart_bov(1)**2*3)) 
      allocate(aoint%gto(input%num_bov(1)**2*3)) 
      open(unit=172, file='soc_ao.dat', status='old')
      read(172, *) aoint%tbgto(:)
      aoint%gto = 0.0
    else
      allocate(aoint%gto(input%num_bov(1)**2*3)) 
      open(unit=172, file='soc_ao.dat', status='old')
      read(172, *) aoint%gto(:)

      print *, "soc", aoint%gto(4:8)
      print *, "soc_size", size(aoint%gto)
      print *, "soc_ao_last", aoint%gto(input%num_bov(1)**2*3-1)
    endif

    close(172)

  end subroutine read_soc_ao
  
  subroutine read_dip_ao()
    ! read dipole integral from MolSOC 
    if (input%qm_flag == 'tddftb') then
      allocate(dipint%tbgto(input%tbcart_bov(1)**2*3)) 
      allocate(dipint%gto(input%num_bov(1)**2*3)) 
      open(unit=176, file='dip_ao.dat', status='old')
      read(176, *) dipint%tbgto(:)
      dipint%gto = 0.0
    else
      allocate(dipint%gto(input%num_bov(1)**2*3)) 
      open(unit=176, file='dip_ao.dat', status='old')
      read(176, *) dipint%gto(:)

      print *, "dip_x(1,1)", dipint%gto(1)
      print *, "dip_y(1,1)", dipint%gto(input%num_bov(1)**2+1)
      print *, "dip_z(1,1)", dipint%gto(input%num_bov(1)**2*2+1)
      print *, "dip_size", size(dipint%gto)
    endif

    close(176)

  end subroutine read_dip_ao

  subroutine read_overlp_ao()
    ! read ovlap integral from the third party QM code

    integer    :: dim_half ! down triangle including diag 
    integer    :: dim_full ! in tddftb

    if (input%qm_flag == 'tddftb') then
      !allocate(overlpint%tbgto(dim_full))
      dim_full = input%num_bov(1)*input%num_bov(1)
      allocate(overlpint%gto(dim_full))
      allocate(overlpint%molout(dim_full))
      open(unit=177, file='ao_overlap.dat', status='old')
      read(177, *) overlpint%gto(:)

      dim_full = input%tbcart_bov(1)*input%tbcart_bov(1)
      allocate(overlpint%molsoc(dim_full))
      open(unit=1771, file='s_matr.dat', status='old')
      read(1771, *) overlpint%molsoc(:)
      overlpint%molout = 0.0
    else
      dim_half = input%num_bov(1)*(input%num_bov(1)+1)/2
      dim_half = input%num_bov(1)**2
      allocate(overlpint%gto(dim_half))
      !open(unit=177, file='ao_overlap.dat', status='old')
      open(unit=177, file='s_matr.dat', status='old')
      read(177, *) overlpint%gto(:)
    endif
   
    close(177)
   
  end subroutine read_overlp_ao
  
  subroutine read_mo()
    
    integer    :: nroot, ndim0, ni, nf
    real(dpr)  :: norm_coeff, tolrence=0.5, error=1.0E-5
    
    ndim1 = input%num_bov(2)+input%num_bov(3)
    allocate(mo_ene(ndim1))
    ndim0 = ndim1*input%num_bov(1)
    allocate(mo_coeff(ndim0))

    open(unit=173, file='mo_ene.dat', status='old')
    read(173, *) mo_ene(:)
    if (input%qm_flag /= 'tddftb') then 
      mo_ene = mo_ene*au2ev
    endif
    open(unit=174, file='mo_coeff.dat', status='old')
    read(174, *) mo_coeff(:)
    
    do nroot = 1, ndim1

      ni = (nroot-1)*input%num_bov(1)
      nf = nroot*input%num_bov(1)
      norm_coeff = sum(mo_coeff(ni+1:nf)**2)
      !print *, "mo_ene", mo_ene(nroot)
      !print *, "norm_mo_coeff", nroot, norm_coeff
     
    enddo

    close(173)
    close(174)

  end subroutine read_mo

  subroutine read_ci_coeff()
    
    integer    :: nroot, total_root, ni, nf
    real(dpr)  :: norm_coeff, norm_xpy
    real(dpr)  :: tolrence=1.0, error=1.0E-5
   
    if (input%qm_flag /= 'tddftb') then 
      !alpha and beta
      base_dim = 2*input%num_bov(2)*input%num_bov(3) 
      total_root = input%n_singl+input%n_tripl
      t_dim = base_dim*total_root*2
      allocate(td_coeff%mix_all(t_dim))
      allocate(td_coeff%xpy(t_dim/2))
      allocate(td_coeff%xmy(t_dim/2))

      open(unit=175, file='ci_coeff.dat', status='old')
      read(175, *) td_coeff%mix_all(:)
      td_coeff%xpy = td_coeff%mix_all(1:t_dim/2) 
      td_coeff%xmy = td_coeff%mix_all(t_dim/2+1:t_dim)
      print *, "x+y(a),x+y(b)", td_coeff%xpy(1), td_coeff%xpy(base_dim/2+1)
      print *, "xpy", td_coeff%mix_all(1), td_coeff%mix_all(t_dim/2)
      print *, "xmy", td_coeff%mix_all(t_dim/2+1), td_coeff%mix_all(t_dim)
      !check normalization
    
      do nroot = 1, total_root 
        ni = (nroot-1)*base_dim
        nf = nroot*base_dim
        norm_coeff = sum(td_coeff%xpy(ni+1:nf)*td_coeff%xmy(ni+1:nf))
        norm_xpy = sum(td_coeff%xpy(ni+1:nf)*td_coeff%xpy(ni+1:nf))
        
        print *, "xpy*xpy from output", sum(td_coeff%xpy(ni+1:nf)*td_coeff%xpy(ni+1:nf))
        td_coeff%xpy(ni+1:nf) = td_coeff%xpy(ni+1:nf)/sqrt(norm_xpy)
        print *, "xpy*xpy normlized", sum(td_coeff%xpy(ni+1:nf)*td_coeff%xpy(ni+1:nf))
        !print *, "xmy*xmy", sum(td_coeff%xmy(ni+1:nf)*td_coeff%xmy(ni+1:nf))
        print *, "nroot=", nroot, "norm_coeff=", norm_coeff
        print *, "x+y(a),x+y(b)", td_coeff%xpy(ni+1), td_coeff%xpy(ni+base_dim/2+1)

        if (abs(norm_coeff - tolrence) > error) then
          print *, "ci coeff are not normalized,,,"
          exit
        endif
      enddo
    
      if(.false.) then
      !obtain X and Y
      td_coeff%xpy = td_coeff%mix_all(1:t_dim/2) + td_coeff%mix_all(t_dim/2+1:t_dim)
      td_coeff%xmy = td_coeff%mix_all(1:t_dim/2) - td_coeff%mix_all(t_dim/2+1:t_dim)
      td_coeff%xpy = td_coeff%xpy/2.0_dpr
      td_coeff%xmy = td_coeff%xmy/2.0_dpr
      endif

    else
      base_dim = input%num_bov(2)*input%num_bov(3) 
      total_root = input%n_singl+input%n_tripl
      t_dim = base_dim*total_root
      allocate(td_coeff%mix_all(t_dim))
      allocate(td_coeff%xpy(t_dim))
      !allocate(td_coeff%xmy(t_dim/2))

      open(unit=175, file='ci_coeff.dat', status='old')
      read(175, *) td_coeff%mix_all(:)
      td_coeff%xpy = td_coeff%mix_all(1:t_dim) 
      do nroot=1, total_root
        ni = (nroot-1)*base_dim
        nf = nroot*base_dim
        norm_xpy = sum(td_coeff%xpy(ni+1:nf)*td_coeff%xpy(ni+1:nf))
        print *, "xpy*xpy from output", sum(td_coeff%xpy(ni+1:nf)*td_coeff%xpy(ni+1:nf))
        td_coeff%xpy(ni+1:nf) = td_coeff%xpy(ni+1:nf)/sqrt(norm_xpy)
        print *, "xpy*xpy normlized", sum(td_coeff%xpy(ni+1:nf)*td_coeff%xpy(ni+1:nf))
      enddo
      td_coeff%xpy = td_coeff%xpy/sqrt(2.0_dpr)
      !td_coeff%xmy = td_coeff%mix_all(t_dim/2+1:t_dim)
      print *, "x+y0", td_coeff%xpy(1)
      print *, "x+y-1", td_coeff%xpy(t_dim)
    endif
       
    close(175)
  end subroutine read_ci_coeff
  
end module init_prep
