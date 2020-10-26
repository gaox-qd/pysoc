program soc_td
  
  use init_prep
  use soc_mo
  use dip_mo
  use norm_mo
  use basis_match
  use soc_state
  use trans_dipole

  implicit none

  !read soc_td_input.dat
  call read_input()
  print *, "input", input%ene_t(:)
  print *, "root_s", input%root_s
  print *, "root_t", input%root_t
  !read basis information
  call read_ao_orb()
  !soc_int.dat
  call read_soc_ao()
  !mo_coeff.dat
  !mo_ene.dat
  call read_mo()
  call read_ci_coeff()
  !ci_coeff.dat
  print *, "ci_size", size(td_coeff%mix_all)
  print *, "mo_coeff_i", mo_coeff(1)
  print *, "mo_coeff_i", mo_coeff(ndim1*input%num_bov(1))

  call read_overlp_ao() !from QM code
  print *, "hello0"
  !call check_norm_mo(overlpint%gto, mo_coeff)
  print *, "hello1"
  if (input%qm_flag /= 'tddftb') then
    !match AO integral(from molsoc) orders with MO coeffs(from QM code)
    !dxx,dyy,dzz,dxy,dxz,dyz in Gaussian
    !while dxx,dxy,dxz,dyy,dyz,dzz in MolSOC
    call basmatch_matr(overlpint%gto, input%num_bov(1), 1)
    call check_norm_mo(overlpint%gto, mo_coeff)
    call basmatch_matr(aoint%gto, input%num_bov(1), 3)
  else
    !Cartesian to spherical GTOs for tddftb
    call basmatch_tb(aoint%tbgto, input%tbcart_bov(1), 3,&
         aoint%gto, input%num_bov(1))
    !convert overlap integral from molsoc, then compare it with that from DFTB+ output(oversqr.dat)
    call basmatch_tb(overlpint%molsoc, input%tbcart_bov(1), 1,&
         overlpint%molout, input%num_bov(1))
    
  endif
  !stop
  !MO integral for SOC
  call cal_soc_moint(aoint%gto, mo_coeff) 
  
  call statint()
  !stop

  !transition dipole
  if(input%do_dip == "True") then
    print *, "calculation of transition dipole"
    call read_dip_ao() !from molsoc
    if (input%qm_flag /= 'tddftb') then
      call basmatch_matr(dipint%gto, ao_orb%nb, 3) !3 means x,y,z
    else
      call basmatch_tb(dipint%tbgto, input%tbcart_bov(1), 3,&
           dipint%gto, input%num_bov(1))
      !checked with Gaussian output
    endif
    call cal_dip_moint(dipint%gto, mo_coeff)
    call dip_state()
    deallocate(dip_moint, dip_singl, dip_tripl)
  endif

  deallocate(input%root_s,input%ene_s)
  deallocate(input%root_t,input%ene_t)
  deallocate(soc_statint)
end program
