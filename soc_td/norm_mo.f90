module norm_mo
  
  use prec_mod 
  use data_input
  use init_prep
  
  implicit none

  real(dpr), allocatable    :: overlap_mo(:,:)

  contains
  
  subroutine check_norm_mo(overlap_temp, mo_coeff_temp)
    
    integer      :: i, j, k, np
    real(dpr)    :: overlap_temp(*)
    real(dpr)    :: mo_coeff_temp(input%num_bov(1),ndim1)
    !real(dpr)    :: mo_coeff_temp0(input%num_bov(1),ndim1)
    real(dpr)    :: overlap_direct(ndim1,ndim1)
    real(dpr)    :: diag_ele

    real(dpr), allocatable    :: overlap_matr(:,:)
    real(dpr), allocatable    :: tp_mo_coeff_temp(:,:)
    
    print *, "input%num_bov(1),ndim1", input%num_bov(1),ndim1
    allocate(overlap_mo(ndim1,ndim1))
    allocate(overlap_matr(input%num_bov(1),input%num_bov(1)))
    allocate(tp_mo_coeff_temp(ndim1, input%num_bov(1)))
    print *, "hello1, norm_mo.f90" 
    overlap_mo = 0.0
    overlap_matr = 0.0
    if (input%qm_flag /= 'tddftb') then
      k = input%num_bov(1)**2
      overlap_matr = reshape(overlap_temp(1:k), (/input%num_bov(1),input%num_bov(1)/))
      np = 0
      do i=1, input%num_bov(1)
        do j=1, i
          np = np + 1
          !print *, "np,j,i", np, j, i
          !overlap_matr(i,j) = overlap_temp(np)
          !overlap_matr(j,i) = overlap_matr(i,j)
        enddo
      enddo
    else
      k = input%num_bov(1)**2
      overlap_matr = reshape(overlap_temp(1:k), (/input%num_bov(1),input%num_bov(1)/))
    endif
    print *, "hello2, norm_mo.f90" 

    !print *, "mo_coeff_temp1", mo_coeff_temp(1,1), mo_coeff_temp(2,1)
    !print *, "mo_coeff_temp1", mo_coeff_temp(input%num_bov(1)-1,ndim1), mo_coeff_temp(input%num_bov(1),ndim1)
    tp_mo_coeff_temp = transpose(mo_coeff_temp)
    !print *, "mo_coeff_temp2", mo_coeff_temp(1,1), mo_coeff_temp(2,1)
    !print *, "mo_coeff_temp2", mo_coeff_temp(input%num_bov(1)-1,ndim1), mo_coeff_temp(input%num_bov(1),ndim1)
    !print *, "tp_mo_coeff_temp", tp_mo_coeff_temp(1,1), tp_mo_coeff_temp(1,2)
    !print *, "tp_mo_coeff_temp", tp_mo_coeff_temp(ndim1,input%num_bov(1)-1), tp_mo_coeff_temp(ndim1, input%num_bov(1))
    !mo_coeff_temp = matmul(overlap_matr, mo_coeff_temp)
    tp_mo_coeff_temp = matmul(tp_mo_coeff_temp, overlap_matr)
    overlap_mo = matmul(tp_mo_coeff_temp, mo_coeff_temp)
    print *, "hello3, norm_mo.f90" 
    print *, "overlap_ao_diag" , overlap_matr(1,1), overlap_matr(input%num_bov(1),input%num_bov(1))
    print *, "overlap_ao_offdiag" , overlap_matr(1,2), overlap_matr(2,1)
    print *, "overlap_mo_diag", overlap_mo(1,1), overlap_mo(ndim1,ndim1)
    print *, "overlap_mo_offdiag", overlap_mo(1,2), overlap_mo(2,1)
    
    do i=1, ndim1
      diag_ele = overlap_mo(i,i) - 1.0
      if( abs(diag_ele) > 1.0e-3) then
        print *, "mo normalization problem, check the MO reading", diag_ele
        stop
      endif
    enddo
    
    deallocate(overlap_matr, tp_mo_coeff_temp)

  end subroutine check_norm_mo

end module norm_mo
