module dip_mo
  
  use prec_mod
  use data_input
  use init_prep

  implicit none
    
  real(dpr), allocatable    :: dip_moint(:,:,:)
  
  contains

  subroutine cal_dip_moint(aoint_temp, mo_coeff_temp0)
  
    integer      :: i, k
    real(dpr)    :: aoint_temp(input%num_bov(1),input%num_bov(1),3)     
    real(dpr)    :: mo_coeff_temp0(input%num_bov(1),ndim1)    
    real(dpr)    :: mo_coeff_temp(input%num_bov(1),ndim1)    
    real(dpr), allocatable    :: tp_mo_coeff_temp(:,:) 
    
    allocate(dip_moint(ndim1, ndim1, 3))
    allocate(tp_mo_coeff_temp(ndim1, input%num_bov(1)))
    
    tp_mo_coeff_temp = transpose(mo_coeff_temp0)
    
    k = input%num_bov(1)
    do i=1, 3
      mo_coeff_temp = matmul(aoint_temp(:,:,i), mo_coeff_temp0)
      dip_moint(:,:,i) = matmul(tp_mo_coeff_temp, mo_coeff_temp)
      print *, "aoint_temp_diag", aoint_temp(k,k,i)
      print *, "aoint_temp_symm", aoint_temp(k,k-1,i), aoint_temp(k-1,k,i)
      print *, "aoint_temp_size", size(aoint_temp(:,1,i)), size(aoint_temp(:,:,i))
      print *, "dip_moint_diag",  dip_moint(1,1,i), dip_moint(ndim1,ndim1,i)
      print *, "dip_moint_size",  size(dip_moint(:,1,i)), size(dip_moint(:,:,i))
      print *, "dip_moint_symm",  dip_moint(ndim1,ndim1-1,i), dip_moint(ndim1-1,ndim1,i)
    enddo
    dip_moint = au2debye*dip_moint
     
    deallocate(tp_mo_coeff_temp)
  end subroutine cal_dip_moint

end module dip_mo
