module soc_state
  
  use prec_mod
  use data_input
  use init_prep
  use soc_mo

  implicit none 

  complex(dpr), allocatable    :: soc_statint(:,:,:)
  
  interface cal_socci
    module procedure cal_socci_ex
    module procedure cal_socci_eg
  end interface 

  contains

  subroutine statint()
    
    integer       :: k, ki
    integer       :: si, ti
    real          :: sum_soc, soc_sub(3)
    complex(dpr)  :: soc_statint_temp(3)
        
    soc_statint_temp = (0.0, 0.0)
    open(unit=17, file='soc_out.dat', status='replace')

    if(input%do_s0 == 'True') then
      allocate(soc_statint(input%n_singl+1,input%n_tripl,3))
      do si=0, input%n_singl
        do ti=1, input%n_tripl
          if(si == 0) then
            call cal_socci(ti, soc_statint_temp)
            soc_statint(si+1, ti, :) = soc_statint_temp(:)*au2wavnum
          else
            call cal_socci(si, ti, soc_statint_temp)
            !print *, "soc_statint_temp0", soc_statint_temp(1:2)
            soc_statint(si+1, ti, :) = soc_statint_temp(:)*au2wavnum
            !print *, "soc_statint"
            !write(*, 1717) soc_statint(si+1, ti, :)
          endif
        enddo
        do k=1, input%n_tripl
          soc_sub(:) = sqrt(soc_statint(si+1, k, :)*conjg(soc_statint(si+1, k, :)))
          sum_soc=sum(soc_statint(si+1, k, :)*conjg(soc_statint(si+1, k, :)))
          sum_soc=sqrt(sum_soc)
          write(17, 1715) si, k, sum_soc, (soc_sub(ki), ki=1,3)
          write(*, 1716) si, k, sum_soc, (soc_statint(si+1, k, ki), ki=1,3)
        enddo
        !write(*, 1717) (soc_statint(si+1, k, 1), k=1,input%n_tripl)
        !write(*, 1717) (soc_statint(si+1, k, 2), k=1,input%n_tripl)
        !write(*, 1717) (soc_statint(si+1, k, 3), k=1,input%n_tripl)
      enddo
    else
      allocate(soc_statint(input%n_singl,input%n_tripl,3))
      do si=1, input%n_singl
        do ti=1, input%n_tripl
          call cal_socci(si, ti, soc_statint_temp)
          soc_statint(si, ti, :) = soc_statint_temp(:)*au2wavnum
        enddo
        do k=1, input%n_tripl
          sum_soc=sum(soc_statint(si+1, k, :)*conjg(soc_statint(si+1, k, :)))
          sum_soc=sqrt(sum_soc)
          write(*, 1716) sum_soc, si, k, (soc_statint(si+1, k, ki), ki=1,3)
        enddo
        write(*, 1716) si, k, sum_soc, (soc_statint(si, k, ki), ki=1,3)
        !write(*, 1717) (soc_statint(si, k, 1), k=1,input%n_tripl)
        !write(*, 1717) (soc_statint(si, k, 2), k=1,input%n_tripl)
        !write(*, 1717) (soc_statint(si, k, 3), k=1,input%n_tripl)
      enddo
    endif
    1715 format('sum_soc, ','<S',I1,'|Hso|T',I1,',1,0,-1> (cm-1): ',4f15.5)
    1716 format('sum_soc, ','<S',I1,'|Hso|T',I1,',1,0,-1> (cm-1): ',7f15.5)  
    1717 format(50f15.5)  
    close(17)   
     
  end subroutine statint

  subroutine cal_socci_eg(ti, soc_statint_temp)

    integer        :: ti
    integer        :: i, a
    integer        :: nt
    real(dpr)      :: pre_all, pre_t, coeff_t
    complex(dpr)   :: sum_st_temp(2), soc_st_temp, soc_statint_temp(3)
    
    soc_statint_temp = (0.0, 0.0)
    nt = base_dim*(input%n_singl+ti-1)
    sum_st_temp = (0.0, 0.0)
    do a=1, input%num_bov(3)
      do i=1, input%num_bov(2)
        coeff_t = td_coeff%xpy(nt+(i-1)*input%num_bov(3)+a)  
        if(abs(coeff_t) > input%coeff_thresh) then
          pre_t = sqrt((mo_ene(input%num_bov(2)+a)-mo_ene(i))/input%ene_t(ti))
          pre_all = pre_t*coeff_t
          soc_st_temp = pre_all*cmplx(soc_moint(i,input%num_bov(2)+a,1), soc_moint(i,input%num_bov(2)+a,2))
          sum_st_temp(1) = sum_st_temp(1) + soc_st_temp
          soc_st_temp = pre_all*cmplx(soc_moint(i,input%num_bov(2)+a,3), 0.0)
          sum_st_temp(2) = sum_st_temp(2) + soc_st_temp
        endif
      enddo
    enddo
    soc_statint_temp(1:2) = soc_statint_temp(1:2) + sum_st_temp(:)
    soc_statint_temp(1) = -1.0/sqrt(2.0)*soc_statint_temp(1)
    soc_statint_temp(3) = -1.0*soc_statint_temp(1)
          
  end subroutine cal_socci_eg

  subroutine cal_socci_ex(si, ti, soc_statint_temp)
    
    integer        :: si, ti
    integer        :: ns, nt
    integer        :: i, j, a, b
    real(dpr)      :: pre_all, pre_s, pre_t, coeff_s, coeff_t
    complex(dpr)   :: sum_st_temp(2), soc_st_temp, soc_statint_temp(3)
    
    soc_statint_temp = (0.0, 0.0)
    ns = base_dim*(si-1)
    nt = base_dim*(input%n_singl+ti-1)
    print *, "s, t", si, ti
    sum_st_temp = (0.0, 0.0)
    do i=1, input%num_bov(2)
      do j=1, input%num_bov(2)
        do a=1, input%num_bov(3)
          coeff_s = td_coeff%xpy(ns+(i-1)*input%num_bov(3)+a)  
          coeff_t = td_coeff%xpy(nt+(j-1)*input%num_bov(3)+a)  
          if(abs(coeff_s) > input%coeff_thresh .or. abs(coeff_t) > input%coeff_thresh &
            & .and. i /= j) then
            pre_s = sqrt((mo_ene(input%num_bov(2)+a)-mo_ene(i))/input%ene_s(si))
            pre_t = sqrt((mo_ene(input%num_bov(2)+a)-mo_ene(j))/input%ene_t(ti))
            pre_all = pre_s*pre_t*coeff_s*coeff_t
            soc_st_temp = pre_all*cmplx(soc_moint(j,i,1), soc_moint(j,i,2))
            sum_st_temp(1) = sum_st_temp(1) + soc_st_temp
            soc_st_temp = pre_all*cmplx(soc_moint(j,i,3), 0.0)
            sum_st_temp(2) = sum_st_temp(2) + soc_st_temp
            if(sqrt(real(soc_st_temp)**2 + aimag(soc_st_temp)**2) > 1.0e-5) then
              print *, "soc_st_temp1", soc_st_temp
              print *, "coeff", coeff_s, coeff_t
              print *, "ene", input%ene_s(si), input%ene_t(ti)
              print *, "ene mo", mo_ene(input%num_bov(2)+a), mo_ene(i), mo_ene(j)
              print *, "soc_moint", soc_moint(j,i,1:2)
            endif
          endif
        enddo
      enddo
    enddo
    soc_statint_temp(1:2) = soc_statint_temp(1:2) + sum_st_temp(:)
    soc_statint_temp(2) = -1.0*soc_statint_temp(2)
    print *, "sum_st_temp1", sum_st_temp(:)*au2wavnum
    
    sum_st_temp = (0.0, 0.0)
    do i=1, input%num_bov(2)
      do b=1, input%num_bov(3)
        do a=1, input%num_bov(3)
          coeff_s = td_coeff%xpy(ns+(i-1)*input%num_bov(3)+a)  
          coeff_t = td_coeff%xpy(nt+(i-1)*input%num_bov(3)+b)  
          if(abs(coeff_s) > input%coeff_thresh .or. abs(coeff_t) > input%coeff_thresh &
             & .and. a /= b) then
            pre_s = sqrt((mo_ene(input%num_bov(2)+a)-mo_ene(i))/input%ene_s(si))
            pre_t = sqrt((mo_ene(input%num_bov(2)+b)-mo_ene(i))/input%ene_t(ti))
            pre_all = pre_s*pre_t*coeff_s*coeff_t
            soc_st_temp = pre_all*cmplx(soc_moint(input%num_bov(2)+a,input%num_bov(2)+b,1), soc_moint(input%num_bov(2)+a,input%num_bov(2)+b,2))
            sum_st_temp(1) = sum_st_temp(1) + soc_st_temp
            soc_st_temp = pre_all*cmplx(soc_moint(input%num_bov(2)+a,input%num_bov(2)+b,3), 0.0)
            sum_st_temp(2) = sum_st_temp(2) + soc_st_temp
            if(sqrt(real(soc_st_temp)**2 + aimag(soc_st_temp)**2) > 1.0e-3) then
              print *, "soc_st_temp2", soc_st_temp
              print *, "coeff", coeff_s, coeff_t
              print *, "ene", input%ene_s(si), input%ene_t(ti)
              print *, "ene mo", mo_ene(input%num_bov(2)+a), mo_ene(input%num_bov(2)+b), mo_ene(i)
              print *, "soc_moint", soc_moint(input%num_bov(2)+a,input%num_bov(2)+b,1:2)
            endif
          endif
        enddo
      enddo
    enddo
    soc_statint_temp(1) = soc_statint_temp(1) - sum_st_temp(1)
    soc_statint_temp(1) = 1.0/sqrt(2.0)*soc_statint_temp(1)
    soc_statint_temp(2) = soc_statint_temp(2) + sum_st_temp(2)
    soc_statint_temp(3) = -1.0*soc_statint_temp(1)
    print *, "sum_st_temp2", sum_st_temp(:)*au2wavnum

  end subroutine cal_socci_ex

end module soc_state
