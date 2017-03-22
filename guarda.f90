subroutine guarda(i)
  use variables
  implicit none
  integer i,j
  
  if(mod(i,us)==0) then
     print*,i*ht
     call derivadas()
     !constricciones
     C_rrr=dg_rr + 8.0d0*g_rr*f_rT*ing_T - 2.0d0*f_rrr 
     C_rT=dg_T + 2.0d0*g_T*inr - 2.0d0*f_rT
     C_r=dK_t/g_T + 2.0d0*K_T/(r*g_T) - f_rT*(grr*k_rr + K_T/g_T)/g_T
     CH=df_rT/(g_rr*g_T) - 1.0d0/(2.0d0*g_T*r**2) + f_rT*(2.0d0/r + &
          3.5d0*f_rT/g_T - grr*f_rrr)/(g_rr*g_T) - K_T*(grr*K_rr +&
          K_T/2.0d0*g_T)/g_T
     
     write(22,*)'#tiempo=',i*ht
     do j=0,n
        write(22,*)r(j),g_rr(j),g_T(j),K_rr(j),K_T(j),C_r(j),CH(j)
     end do
     write(22,*)
     write(22,*)
     norma=0.0d0
     do j=0,n
        norma=norma + CH(j)**2
     end do
     norma=sqrt(norma)
     write(11,*)i*ht,norma
     call evalua(C_r,0.05d0,r_min,r_max,i*ht,21)
     call evalua(CH,0.05d0,r_min,r_max,i*ht,20)
  end if
  !====================================================================
end subroutine guarda

