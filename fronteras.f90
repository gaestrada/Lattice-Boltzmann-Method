subroutine fronteras(i)
  use variables
  implicit none
  integer i
  dk_T(0)=df_rT(0) / sqrt(g_rr(0)) - sqrt(g_rr(0)) / (2.0d0*r(0)**2)  +&
       f_rT(0) / sqrt(g_rr(0)) * (2.0d0/r(0) + (3.5d0*f_rT(0)/g_T(0)) -&
       f_rrr(0)/g_rr(0)) - sqrt(g_rr(0))*K_T(0) * (K_rr(0)/g_rr(0) + &
       K_T(0)/(2.0d0*g_T(0))) - 2.0d0*K_T(0) / r(0)  + f_rT(0) *&
       (K_rr(0)/g_rr(0) + K_T(0)/g_T(0))
  kg_rr(i,0)=2.0d0*beta_r(0)*f_rrr(0)- 8.0d0*g_rr(0)*f_rT(0)*beta_r(0)*ing_T(0) -&
       2.0d0*ene(0)*K_rr(0) + 2.0d0*g_rr(0)*dbeta_r(0)
  kg_T(i,0)=2.0d0*beta_r(0)*f_rT(0) - 2.0d0*ene(0)*K_T(0)
  kK_T(i,0)=beta_r(0)*dk_T(0) - ene(0)*grr(0)*df_rT(0) + ene(0)*(K_T(0)*grr(0)*K_rr(0) +&
       inr(0)**2 - 2.0d0*grr(0)*ing_T(0)*&
       f_rT(0)**2 - f_rT(0)*grr(0)*dlog_alpha(0)) + 2.0d0*beta_r(0)*inr(0)*K_T(0)
  
  kf_rT(i,0) = beta_r(0)*df_rT(0) - ene(0)*dk_T(0) + ene(0)*(K_T(0)*(2.0d0*f_rT(0)*&
       ing_T(0) - grr(0)*f_rrr(0) - &
       dlog_alpha(0))) + (dbeta_r(0) + 2.0d0*beta_r(0)*inr(0))*f_rT(0)
  
  !kg_rr(i,0)=0.0d0
  !kg_T(i,0)=0.0d0
  kK_rr(i,0)=kf_rrr(i,0)/sqrt(g_rr(0))
  !kK_T(i,0)=k4f_rt(0)/sqrt(g_rr(0))
end subroutine fronteras
