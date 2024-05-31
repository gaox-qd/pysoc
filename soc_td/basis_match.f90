module basis_match

  use prec_mod
  use data_input
  use init_prep
  
  implicit none
   
!  interface basmath
!    module procedure basmatch_matr
!    module procedure basmatch_tb
!  end interface

  contains
  
  subroutine shell_matrix(nsh, matrix_n, flag)
    
    integer      :: nsh, flag
    real(dpr)    :: matrix_n(nsh,nsh)

    matrix_n = 0.0
    
    select case(nsh)
    case(6)
      flag = 1
      matrix_n(1,1) = 1.0
      matrix_n(2,4) = 1.0
      matrix_n(3,6) = 1.0
      matrix_n(4,2) = 1.0
      matrix_n(5,3) = 1.0
      matrix_n(6,5) = 1.0
    case(10)
      flag = 1
      matrix_n(1,1)   = 1.0
      matrix_n(2,7)   = 1.0
      matrix_n(3,10)  = 1.0
      matrix_n(4,4)   = 1.0
      matrix_n(5,2)   = 1.0
      matrix_n(6,3)   = 1.0
      matrix_n(7,6)   = 1.0
      matrix_n(8,9)   = 1.0
      matrix_n(9,8)   = 1.0
      matrix_n(10,5)  = 1.0
          
    case default
      flag = -1

    end select
    
  end subroutine shell_matrix

  subroutine basmatch_matr(matrix, ndim, compont) 
  
    integer    :: flag, compont, ndim
    integer    :: k, i, j, dim_i, dim_j, nshi, nshj
    real(dpr)  :: matrix(ndim,ndim,compont)
    real(dpr), allocatable    :: matrix_temp(:,:)
    real(dpr), allocatable    :: matrix_i(:,:), matrix_j(:,:)
   
    do k=1, compont 
      !print *,"DIM_old=", k
      !print *, matrix(10:20,10,k)
      dim_i = 0
      do i=1, ao_orb%nc
        nshi = ao_orb%nsh(i)
        dim_i = dim_i + nshi
        dim_j = 0
        do j=1, ao_orb%nc
          nshj = ao_orb%nsh(j)
          dim_j = dim_j + nshj
          if(nshi == 6 .or. nshj == 6 .or. nshi == 10 .or. nshj == 10) then
            !print *, "nshi,nshj", nshi, nshj
            allocate(matrix_temp(nshj, nshi))
            allocate(matrix_i(nshi, nshi))
            allocate(matrix_j(nshj, nshj))
            matrix_temp(:,:) = matrix(dim_j-nshj+1:dim_j,dim_i-nshi+1:dim_i,k)
            call shell_matrix(nshi, matrix_i, flag)
            if(flag == 1) then
              matrix_i = transpose(matrix_i)
              matrix_temp = matmul(matrix_temp, matrix_i)
            endif
            call shell_matrix(nshj, matrix_j, flag)
            if(flag == 1) then
              matrix_temp = matmul(matrix_j, matrix_temp)
            endif
            matrix(dim_j-nshj+1:dim_j,dim_i-nshi+1:dim_i,k) = matrix_temp
            deallocate(matrix_temp, matrix_i, matrix_j)
          endif       
        enddo
      enddo
      !print *,"DIM_new=", k
      !print *, matrix(10:20,10,k)
    enddo 
    
  end subroutine basmatch_matr

  subroutine shell_tb(nsh1, nsh2, matrix_n, flag)
    
    integer      :: nsh1, nsh2, flag
    real(dpr)    :: matrix_n(nsh1,nsh2)
    real(dpr)    :: cs, cp, cd1, cd2, cf1, cf2, cf3, cf4

    cs = r_pi
    cp = r_pi * sqrt(3.0) 
    cd1 = r_pi * sqrt(15.0) 
    cd2 = 0.5 * r_pi * sqrt(5.0)
    cf1 = r_pi * sqrt(35.0/8.0)
    cf2 = r_pi * sqrt(105.0/4.0)
    cf3 = r_pi * sqrt(21.0/8.0)
    cf4 = r_pi * sqrt(7.0/4.0)
    matrix_n = 0.0
    
    select case(nsh1)
    case(1)
      flag = 1
      matrix_n(1,1) = cs
    case(3)
      flag = 1
      matrix_n(1,2) = cp
      matrix_n(2,3) = cp
      matrix_n(3,1) = cp
    case(5)
      flag = 1
      matrix_n(1,2) = cd1
      matrix_n(2,5) = cd1
      matrix_n(3,1) = -cd2
      matrix_n(3,4) = -cd2
      matrix_n(3,6) = 2.0 * cd2
      matrix_n(4,3) = cd1
      matrix_n(5,1) = 0.5 * cd1
      matrix_n(5,4) = -0.5 * cd1
    case(7)
      flag = 1
      matrix_n(1,2)  = 3.0 * cf1 
      matrix_n(1,7)  = -cf1 
      matrix_n(2,5)  = 2.0 * cf2
      matrix_n(3,2)  = -cf3
      matrix_n(3,7)  = -cf3
      matrix_n(3,9)  = 4.0 * cf3
      matrix_n(4,3)  = -3.0 * cf4
      matrix_n(4,8)  = -3.0 * cf4
      matrix_n(4,10) = 2.0 * cf4
      matrix_n(5,1)  = -cf3
      matrix_n(5,4)  = -cf3
      matrix_n(5,6)  = 4.0 * cf3
      matrix_n(6,3)  = cf2
      matrix_n(6,8)  = -cf2
      matrix_n(7,1)  = cf1
      matrix_n(7,4)  = -3.0 * cf1
    case default
      flag = -1

    end select
    
  end subroutine shell_tb

  subroutine basmatch_tb(matrix, ndim, compont, matr_out, dim_out) 
    !from Cartesian to spherical AOs 
    integer    :: flag, compont, ndim, dim_out
    integer    :: k, i, j, dim_i, dim_j, nshi, nshj
    integer    :: dim_i1, dim_i2, dim_j1, dim_j2
    integer    :: nshi1, nshi2, nshj1, nshj2
    real(dpr)  :: matrix(ndim,ndim,compont)
    real(dpr)  :: matr_out(dim_out,dim_out,compont)
    real(dpr), allocatable    :: matrix_temp1(:,:), matrix_temp2(:,:)
    real(dpr), allocatable    :: matrix_temp3(:,:)
    real(dpr), allocatable    :: matrix_i(:,:), matrix_j(:,:)
    real(dpr), allocatable    :: tp_matrix_i(:,:)

    !print *, "mo_coeff_i", mo_coeff(1)
    !print *, "mo_coeff_i", mo_coeff(ndim1*input%num_bov(1))
   
    do k=1, compont 
      !print *,"DIM_old=", k
      !print *, matrix(52:68,68,k)
      dim_i1 = 0
      dim_i2 = 0
      do i=1, ao_orb%nc
        nshi = ao_orb%nsh(i)
        if(nshi == 6) then
          nshi1 = 5
          nshi2 = 6
        else if(nshi == 10) then
          nshi1 = 7
          nshi2 = 10
        else 
          nshi1 = nshi
          nshi2 = nshi
        endif
 
        dim_i1 = dim_i1 + nshi1
        dim_i2 = dim_i2 + nshi2
        dim_j1 = 0
        dim_j2 = 0
        do j=1, ao_orb%nc
          nshj = ao_orb%nsh(j)
          !print *, "nshi, nshj", nshi, nshj
          if(nshj == 6 ) then
            !print *, "nshi, nshj", nshi, nshj
            nshj1 = 5
            nshj2 = 6
          else if(nshj == 10) then
            nshj1 = 7
            nshj2 = 10
          else
            nshj1 = nshj
            nshj2 = nshj
          endif

          dim_j1 = dim_j1 + nshj1
          dim_j2 = dim_j2 + nshj2

          allocate(matrix_temp1(nshj2, nshi2))
          allocate(matrix_temp2(nshj2, nshi1))
          allocate(matrix_temp3(nshj1, nshi1))
          allocate(matrix_i(nshi1, nshi2))
          allocate(tp_matrix_i(nshi2, nshi1))
          allocate(matrix_j(nshj1, nshj2))
          ![J]_nshj1*nshj2 * [Martrix]_nshj2*nshi2 * [transpose(I)]_nshi2*nshi1
          matrix_temp1(:,:) = matrix(dim_j2-nshj2+1:dim_j2,dim_i2-nshi2+1:dim_i2,k)
          call shell_tb(nshi1, nshi2, matrix_i, flag)
          tp_matrix_i = transpose(matrix_i) !nshi2*nshi1
          matrix_temp2 = matmul(matrix_temp1, tp_matrix_i)
          call shell_tb(nshj1, nshj2, matrix_j, flag)
          matrix_temp3 = matmul(matrix_j, matrix_temp2)
          matr_out(dim_j1-nshj1+1:dim_j1,dim_i1-nshi1+1:dim_i1,k) = matrix_temp3
          !print *, "label1", dim_j2-nshj2+1,dim_j2,dim_i2-nshi2+1,dim_i2
          !print *, "label2", dim_j1-nshj1+1,dim_j1,dim_i1-nshi1+1,dim_i1
          !print *, "matrix_i", matrix_i
          !print *, "matrix_j", matrix_j
          !print *, "matrix_temp1", matrix_temp1
          !print *, "matrix_temp3", matrix_temp3
          deallocate( matrix_i, matrix_j, tp_matrix_i)
          deallocate(matrix_temp1, matrix_temp2, matrix_temp3)
          if(j == 4) then
          !  stop
          endif
        enddo
      enddo
      print *, "k=", k
      print *, "matr_in diag", matrix(1,1,k)
      print *, "matr_in off1", matrix(1,2,k), matrix(2,1,k)
      print *, "matr_in off2", matrix(ndim,ndim-1,k), matrix(ndim-1,ndim,k)
      do i=1, dim_out
        !write(*, "(500f12.7)") matr_out(:,i,k)
      enddo
      !print *, "matr_out diag2", matr_out(dim_out,dim_out,k)
      !print *, "matr_out off2", matr_out(dim_out,dim_out-1,k), matr_out(dim_out-1,dim_out,k)
    enddo 
    
  end subroutine basmatch_tb

end module
