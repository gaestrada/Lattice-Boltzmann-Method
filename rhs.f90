subroutine rhs(j)
  use variables
  implicit none
  integer j,i
  do i=0,n
     kg_rr(j,i) = -2.0d0*alpha(i)*g_T(i)*sqrt(g_rr(i))*K_rr(i) +&
          beta_r(i)*dg_rr(i) + 2.0d0*g_rr(i)*dbeta_r(i)
     
     kg_T(j,i) = -2.0d0*alpha(i)*g_T(i)*sqrt(g_rr(i))*K_T(i) + &
          2.0d0 * beta_r(i) *g_T(i) / r(i) + beta_r(i)*dg_T(i)
     
     kK_rr(j,i) = alpha(i)*g_T(i)*sqrt(g_rr(i)) * (2.0d0 * (f_rrr(i)/g_rr(i))*&
          ((f_rrr(i)/g_rr(i)) +  1.0d0 /r(i) -	4.0d0*f_rT(i)/g_T(i))&
          - 6.0d0/(r(i)**2)					&
          + K_rr(i)*(						&
          2.0d0*K_T(i)/g_T(i)			&
          - K_rr(i)/g_rr(i))			&
          - 6.0d0*(f_rT(i)/g_T(i))**2			&
          - ddlog_alpha(i)					&
          - (dlog_alpha(i))**2					&
          +((4.0d0/r(i)) - (f_rrr(i)/g_rr(i))) * dlog_alpha(i))	&
          + 2.0d0 * K_rr(i) * dbeta_r(i) 					&
          + beta_r(i) * dK_rr(i)						&
          - (alpha(i)*g_T(i)/sqrt(g_rr(i))) * df_rrr(i)
     
     kK_T(j,i) = alpha(i)*g_T(i)*sqrt(g_rr(i)) * (						&
          K_T(i) * K_rr(i)/g_rr(i)			&
          + 1.0d0 / r(i)**2				&
          - 2.0d0 * (f_rT(i)**2)/(g_rr(i)*g_T(i)) 	&
          - f_rT(i)* dlog_alpha(i)/g_rr(i))		&
          + 2.0d0 * beta_r(i) * K_T(i) / r(i)					&
          + beta_r(i) * dK_T(i)						&
          - alpha(i)*g_T(i)*df_rT(i)/sqrt(g_rr(i))
     
     kf_rrr(j,i) = alpha(i)*g_T(i)*sqrt(g_rr(i)) * (   &
          4.0d0 * g_rr(i) * K_T(i) / g_T(i) * (&
          3.0d0*f_rT(i)/g_T(i)	&
          - (f_rrr(i)/g_rr(i))	&
          + 2.0d0 / r(i)	&
          - dlog_alpha(i))	&
          - K_rr(i) * (				&
          10.0d0 * f_rT(i)/g_T(i)  				&
          + f_rrr(i)/g_rr(i)				&
          - 2.0d0/r(i)					&
          + dlog_alpha(i)))				&
          + 3.0d0*f_rrr(i)*dbeta_r(i)				&
          + g_rr(i) * ddbeta_r(i)				&
          + beta_r(i) * df_rrr(i)				&
          - alpha(i)*g_T(i)*sqrt(g_rr(i)) * dK_rr(i)
     
     kf_rT(j,i) = alpha(i)*g_T(i)*sqrt(g_rr(i)) *(K_T(i) * (	&
          2.0d0 * f_rT(i) /g_T(i)		&
          - f_rrr(i)/g_rr(i)			&
          - dlog_alpha(i)))			&
          + f_rT(i) * (dbeta_r(i) + 2.0d0*beta_r(i)/r(i))		&
          + beta_r(i)* df_rT(i)					&
          - alpha(i)*g_T(i)*sqrt(g_rr(i))*dK_T(i)
  end do  
end subroutine rhs
