subroutine RK4()
  use variables
  implicit none
  integer i,j
  !guardo valores de las variables en auxiliares
  aux1 = g_rr
  aux2 = g_T
  aux3 = K_rr
  aux4 = K_T
  aux5 = f_rrr
  aux6 = f_rt
     
  do i=1,4
     if(i==1) then
        call derivadas()
        call rhs(i)
        call fronteras(i)
     else if (i==2) then
        do j=0,n
           g_rr(j) = aux1(j) + ht*kg_rr(i-1,j)*0.5d0
           g_T(j)= aux2(j) + ht*kg_T(i-1,j)*0.5d0
           K_rr(j)= aux3(j) + ht*kK_rr(i-1,j)*0.5d0
           K_T(j)= aux4(j) + ht*kK_T(i-1,j)*0.5d0
           f_rrr(j)= aux5(j) + ht*kf_rrr(i-1,j)*0.5d0
           f_rT(j)= aux6(j) + ht*kf_rT(i-1,j)*0.5d0
        end do
        grr=1.0d0/g_rr
        ing_T=1.0d0/g_T
        call derivadas()
        call rhs(i)
        call fronteras(i)
     else if (i==3) then
        do j=0,n
           g_rr(j) = aux1(j) + ht*kg_rr(i-1,j)*0.5d0
           g_T(j)= aux2(j) + ht*kg_T(i-1,j)*0.5d0
           K_rr(j)= aux3(j) + ht*kK_rr(i-1,j)*0.5d0
           K_T(j)= aux4(j) + ht*kK_T(i-1,j)*0.5d0
           f_rrr(j)= aux5(j) + ht*kf_rrr(i-1,j)*0.5d0
           f_rT(j)= aux6(j) + ht*kf_rT(i-1,j)*0.5d0
        end do
        grr=1.0d0/g_rr
        ing_T=1.0d0/g_T
        call derivadas()
        call rhs(i)
        call fronteras(i)
     else if(i==4) then
        do j=0,n
           g_rr(j) = aux1(j) + ht*kg_rr(i-1,j)
           g_T(j)= aux2(j) + ht*kg_T(i-1,j)
           K_rr(j)= aux3(j) + ht*kK_rr(i-1,j)
           K_T(j)= aux4(j) + ht*kK_T(i-1,j)
           f_rrr(j)= aux5(j) + ht*kf_rrr(i-1,j)
           f_rT(j)= aux6(j) + ht*kf_rT(i-1,j)
        end do
        grr=1.0d0/g_rr
        ing_T=1.0d0/g_T
        call derivadas()
        call rhs(i)
        call fronteras(i)
     end if

     do j=0,n
        g_rr(j) = aux1(j) + ht*in6*(kg_rr(1,j) + 2.0d0*kg_rr(2,j) + 2.0d0*kg_rr(3,j) + kg_rr(4,j))
        
        g_T(j) = aux2(j) + ht*in6*(kg_T(1,j) + 2.0d0*kg_T(2,j) + 2.0d0*kg_T(3,j) + kg_T(4,j))
        
        K_rr(j) = aux3(j) + ht*in6*(kK_rr(1,j) + 2.0d0*kK_rr(2,j) + 2.0d0*kK_rr(3,j) + kK_rr(4,j))
        
        K_T(j) = aux4(j) + ht*in6*(kk_T(1,j) + 2.0d0*kk_T(2,j) + 2.0d0*kk_T(3,j) + kk_T(4,j))
        
        f_rrr(j) = aux5(j) + ht*in6*(kf_rrr(1,j) + 2.0d0*kf_rrr(2,j) + 2.0d0*kf_rrr(3,j) + kf_rrr(4,j))
        
        f_rT(j) = aux6(j) + ht*in6*(kf_rT(1,j) + 2.0d0*kf_rT(2,j) + 2.0d0*kf_rT(3,j) + kf_rT(4,j))
     end do
     
     ene=alpha*g_T*sqrt(g_rr)
     grr=1.0d0/g_rr
     ing_T=1.0d0/g_T
  end do
end subroutine RK4
