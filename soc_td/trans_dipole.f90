module trans_dipole
  
  use prec_mod
  use data_input
  use init_prep
  use dip_mo

  implicit none 

  real(dpr), allocatable    :: dip_singl(:,:,:)
  real(dpr), allocatable    :: dip_tripl(:,:,:)
  
  interface cal_dipci
    module procedure cal_dipci_ex
    module procedure cal_dipci_eg
  end interface 

  contains

  subroutine dip_state()
    
    integer       :: k
    integer       :: si, sj, ti, tj, ki
    character(10) :: multip
    real(dpr)     :: dip_statint_temp(3)
    
    dip_singl = 0.0
    dip_tripl = 0.0    
    dip_statint_temp = 0.0
    print *, "input%ene_s", input%ene_s
    print *, "input%ene_t", input%ene_t
    open(unit=17, file='trans_dip_out.dat', status='replace')
    if(input%do_s0 == 'True') then
      multip = "singlet"
      allocate(dip_singl(input%n_singl,input%n_singl,3))
      !singlets
      do si=0, input%n_singl-1
        do sj=si+1, input%n_singl
          if(si == 0) then
            call cal_dipci(sj, dip_statint_temp)
            dip_singl(si+1, sj, :) = dip_statint_temp(:)
          else
            call cal_dipci(si, sj, multip, input%ene_s, dip_statint_temp)
            dip_singl(si+1, sj, :) = dip_statint_temp(:)
          endif
          write(17, 1715) si, sj, (dip_singl(si+1, sj, ki), ki=1,3)
        enddo
      enddo
      !triplets
      multip = "triplet"
      allocate(dip_tripl(input%n_tripl,input%n_tripl,3))
      do ti=1, input%n_tripl
        do tj=ti+1, input%n_tripl
          call cal_dipci(ti, tj, multip, input%ene_t, dip_statint_temp)
          dip_tripl(ti, tj, :) = dip_statint_temp(:)
          write(17, 1717) ti, tj, (dip_tripl(ti, tj, ki), ki=1,3)
        enddo
      enddo
      
      !write(*, 1716) (dip_singl(1, k, 1), k=2,input%n_singl+1)
      !write(*, 1716) (dip_singl(1, k, 2), k=2,input%n_singl+1)
      !write(*, 1716) (dip_singl(1, k, 3), k=2,input%n_singl+1)

    else
      multip = "singlet"
      allocate(dip_singl(input%n_singl,input%n_singl,3))
      !singlets
      do si=1, input%n_singl
        do sj=si+1, input%n_singl
          call cal_dipci(si, sj, multip, input%ene_s, dip_statint_temp)
          dip_singl(si, sj, :) = dip_statint_temp(:)
          write(17, 1715) si, sj, (dip_singl(si, sj, ki), ki=1,3)
        enddo
      enddo
      !triplets
      multip = "triplet"
      allocate(dip_tripl(input%n_tripl,input%n_tripl,3))
      do ti=1, input%n_tripl
        do tj=ti+1, input%n_tripl
          call cal_dipci(ti, tj, multip, input%ene_t, dip_statint_temp)
          dip_tripl(ti, tj, :) = dip_statint_temp(:)
          write(17, 1717) ti, tj, (dip_tripl(ti, tj, ki), ki=1,3)
        enddo
      enddo
    endif
    1715 format('<S',I1,'|u|S',I1,',x,y,z> (Debye): ',3f15.5)
    1717 format('<T',I1,'|u|T',I1,',x,y,z> (Debye): ',3f15.5)
    1716 format(50f15.5)  
    close(17)

  end subroutine dip_state

  subroutine cal_dipci_eg(si, dip_statint_temp)

    integer     :: si
    integer     :: ni
    integer     :: i, a
    real(dpr)   :: pre_all, pre_e, coeff_e
    real(dpr)   :: sum_st_temp(3), dip_st_temp(3), dip_statint_temp(3)
    
    print *, "si", si
    dip_statint_temp = 0.0
    ni = base_dim*(si-1)
    sum_st_temp = 0.0
    do a=1, input%num_bov(3)
      do i=1, input%num_bov(2)
        coeff_e = td_coeff%xpy(ni+(i-1)*input%num_bov(3)+a)  
        !coeff_e = td_coeff%xpy(ni+(a-1)*input%num_bov(3)+i)  
        if(abs(coeff_e) > input%coeff_thresh) then
          !pre_e = sqrt((mo_ene(input%num_bov(2)+a)-mo_ene(i))/input%ene_s(si))
          !pre_all = pre_e*coeff_e
          pre_all = coeff_e
          dip_st_temp = pre_all*dip_moint(i,a+input%num_bov(2),:)
          !print *, "pre_all, dip_moint", pre_all, dip_moint(i,a,:)
          sum_st_temp = sum_st_temp + dip_st_temp
          !print *, "sum_st_temp, dip_st_temp", sum_st_temp, dip_st_temp
          !print *, "i->a", i, input%num_bov(2)+a 
          !print *, "coeff_e", coeff_e
          !print *, "ene_s", input%ene_s(si), si
          !print *, "mo_ene", mo_ene(input%num_bov(2)+a), mo_ene(i)
        endif
      enddo
    enddo
    
    dip_statint_temp = dip_statint_temp + 2.0*sum_st_temp
    write(*, "(A20,3f15.7)" ) "dip_statint_temp_eg", dip_statint_temp

  end subroutine cal_dipci_eg

  subroutine cal_dipci_ex(ei, ej, multip, ene_temp, dip_statint_temp)
    
    integer        :: ei, ej
    integer        :: ni, nf
    integer        :: i, j, a, b
    character(10)  :: multip
    real(dpr)      :: ene_temp(*)
    real(dpr)      :: pre_all, pre_e1, pre_e2, coeff_e1, coeff_e2
    real(dpr)      :: sum_st_temp(3), dip_st_temp(3), dip_statint_temp(3)
    real(dpr)      :: sum_dip_mo(3)
    
    print *, "multip:  ", multip
    print *, "ei,ej", ei, ej
    print *, "ene_temp", ene_temp(ei), ene_temp(ej)
    dip_statint_temp = 0.0
    if(multip == 'singlet') then
      ni = base_dim*(ei-1)
      nf = base_dim*(ej-1)
    else if(multip == 'triplet') then
      ni = base_dim*(input%n_singl+ei-1)
      nf = base_dim*(input%n_singl+ej-1)
    end if
    sum_st_temp = 0.0
    do i=1, input%num_bov(2)
      do j=1, input%num_bov(2)
        do a=1, input%num_bov(3)
          coeff_e1 = td_coeff%xpy(ni+(i-1)*input%num_bov(3)+a)  
          coeff_e2 = td_coeff%xpy(nf+(j-1)*input%num_bov(3)+a)  
          if(abs(coeff_e1) > input%coeff_thresh .or. abs(coeff_e2) &
             & > input%coeff_thresh) then
            !pre_e1 = sqrt((mo_ene(input%num_bov(2)+a)-mo_ene(i))/ene_temp(ei))
            !pre_e2 = sqrt((mo_ene(input%num_bov(2)+a)-mo_ene(j))/ene_temp(ej))
            !pre_all = pre_e1*pre_e2*coeff_e1*coeff_e2
            pre_all = coeff_e1*coeff_e2
            dip_st_temp = pre_all*dip_moint(j,i,:)
            sum_st_temp = sum_st_temp + dip_st_temp
          endif
        enddo
      enddo
    enddo
    dip_statint_temp = dip_statint_temp - 2.0*sum_st_temp
    !print *, "dip_statint_temp1", 2.0*sum_st_temp

    sum_st_temp = 0.0
    do i=1, input%num_bov(2)
      do b=1, input%num_bov(3)
        do a=1, input%num_bov(3)
          coeff_e1 = td_coeff%xpy(ni+(i-1)*input%num_bov(3)+a)  
          coeff_e2 = td_coeff%xpy(nf+(i-1)*input%num_bov(3)+b)  
          if(abs(coeff_e1) > input%coeff_thresh .or. abs(coeff_e2) &
             & > input%coeff_thresh) then
            !pre_e1 = sqrt((mo_ene(input%num_bov(2)+a)-mo_ene(i))/ene_temp(ei))
            !pre_e2 = sqrt((mo_ene(input%num_bov(2)+b)-mo_ene(i))/ene_temp(ej))
            !pre_all = pre_e1*pre_e2*coeff_e1*coeff_e2
            pre_all = coeff_e1*coeff_e2
            dip_st_temp = pre_all*dip_moint(a+input%num_bov(2),b+input%num_bov(2),:)
            sum_st_temp = sum_st_temp + dip_st_temp
          endif
        enddo
      enddo
    enddo
    dip_statint_temp = dip_statint_temp + 2.0*sum_st_temp
    !print *, "dip_statint_temp2", 2.0*sum_st_temp

    sum_st_temp = 0.0
    do i=1, input%num_bov(2)
      do a=1, input%num_bov(3)
        coeff_e1 = td_coeff%xpy(ni+(i-1)*input%num_bov(3)+a)  
        coeff_e2 = td_coeff%xpy(nf+(i-1)*input%num_bov(3)+a)  
        if(abs(coeff_e1) > input%coeff_thresh .or. abs(coeff_e2) &
           & > input%coeff_thresh) then
          !pre_e1 = sqrt((mo_ene(input%num_bov(2)+a)-mo_ene(i))/ene_temp(ei))
          !pre_e2 = sqrt((mo_ene(input%num_bov(2)+a)-mo_ene(i))/ene_temp(ej))
          !pre_all = pre_e1*pre_e2*coeff_e1*coeff_e2
          pre_all = coeff_e1*coeff_e2
          dip_st_temp = pre_all
          sum_st_temp = sum_st_temp + dip_st_temp
        endif
      enddo
    enddo

    sum_dip_mo = 0.0
    do j=1, input%num_bov(2)
      sum_dip_mo = sum_dip_mo + dip_moint(j,j,:)
    enddo

    dip_statint_temp = dip_statint_temp + 4.0*sum_st_temp*sum_dip_mo
    !print *, "4.0*sum_st_temp*sum_dip_mo", 4.0*sum_st_temp*sum_dip_mo
    write(*, "(A20,3f15.7)" ) "dip_statint_temp_ex", dip_statint_temp

  end subroutine cal_dipci_ex

end module trans_dipole
