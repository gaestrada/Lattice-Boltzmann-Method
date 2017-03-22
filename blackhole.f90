module variables
  implicit none
  real(8),allocatable,dimension(:)::x,r
  real(8),allocatable,dimension(:,:)::Dn
  real(8) coef_escala,r_min,r_max
  integer n
end module variables

program e3
  use variables
  implicit none
  real(8),allocatable,dimension(:)::g_rr,g_T,K_rr,K_T,f_rrr,f_rT
  real(8),allocatable,dimension(:)::dg_rr,dg_T,dK_rr,dK_T,df_rrr,df_rT
  real(8),allocatable,dimension(:)::alpha,log_alpha,beta_r
  real(8),allocatable,dimension(:)::dbeta_r,ddbeta_r,dlog_alpha,ddlog_alpha
  real(8),allocatable,dimension(:)::aux1,aux2,aux3,aux4,aux5,aux6
  real(8),allocatable,dimension(:)::ene,grr,ing_T,inr
  real(8),allocatable,dimension(:)::CH,C_r,C_rrr,C_rT
  real(8),allocatable,dimension(:)::k1g_rr,k2g_rr,k3g_rr,k4g_rr
  real(8),allocatable,dimension(:)::k1g_T,k2g_T,k3g_T,k4g_T
  real(8),allocatable,dimension(:)::k1K_rr,k2K_rr,k3K_rr,k4K_rr
  real(8),allocatable,dimension(:)::k1K_T,k2K_T,k3K_T,k4K_T
  real(8),allocatable,dimension(:)::k1f_rrr,k2f_rrr,k3f_rrr,k4f_rrr
  real(8),allocatable,dimension(:)::k1f_rT,k2f_rT,k3f_rT,k4f_rT
  real(8) M,error,norma
  real(8) pi,ht,in6,ttiempo,unseg
  real(8) T,U,DT,L,start,finish,A,r0,sigma
  integer tiempo,us
  integer i,j,coordenadas,condicionfrontera

  !definiendo parametros
  n=10
  pi=acos(-1.0d0)
  r_min=1.9d0
  r_max=11.9d0
  ht= 1.0d0/real(n)**2
  M=1.0d0
  coef_escala=2.0d0/(r_max-r_min)
  in6=1.0d0/6.0d0
  ttiempo=1000.0d0/ht
  tiempo=int(ttiempo)
  unseg= 1.0d0/ht
  us=int(unseg)
  A=0.1d0
  r0=5.0d0
  sigma=0.05d0
  coordenadas= 1
  condicionfrontera= 1
  

  !abriendo archivos de datos
  !open(10,file='32g_rr.dat')
  !open(11,file='12norma.dat')
  open(20,file='10Hnorma2.dat')
  !open(21,file='24C.dat')
  open(22,file='10variables.dat')

  !asignando memora
  allocate(x(0:n),r(0:n))
  allocate(Dn(0:n,0:n))
  allocate(g_rr(0:n),g_T(0:n),K_rr(0:n),K_T(0:n),f_rrr(0:n),f_rT(0:n))
  allocate(dg_rr(0:n),dg_T(0:n),dK_rr(0:n),dK_T(0:n),df_rrr(0:n),df_rT(0:n))
  allocate(alpha(0:n),log_alpha(0:n),beta_r(0:n))
  allocate(dbeta_r(0:n),ddbeta_r(0:n),dlog_alpha(0:n),ddlog_alpha(0:n))
  allocate(aux1(0:n),aux2(0:n),aux3(0:n),aux4(0:n),aux5(0:n),aux6(0:n))
  allocate(ene(0:n),grr(0:n),ing_T(0:n),inr(0:n))
  allocate(CH(0:n),C_r(0:n),C_rrr(0:n),C_rT(0:n))
  allocate(k1g_rr(0:n),k2g_rr(0:n),k3g_rr(0:n),k4g_rr(0:n))
  allocate(k1g_T(0:n),k2g_T(0:n),k3g_T(0:n),k4g_T(0:n))
  allocate(k1K_rr(0:n),k2K_rr(0:n),k3K_rr(0:n),k4K_rr(0:n))
  allocate(k1K_T(0:n),k2K_T(0:n),k3K_T(0:n),k4K_T(0:n))
  allocate(k1f_rrr(0:n),k2f_rrr(0:n),k3f_rrr(0:n),k4f_rrr(0:n))
  allocate(k1f_rT(0:n),k2f_rT(0:n),k3f_rT(0:n),k4f_rT(0:n))

  !======definiendo los puntos de colocacion en las mallas==============
  !definiendo la malla x
  do i=0,n
     !if (i==(n/2)) then
      !  x(i)=0.0d0
     !else     
        x(i)=cos(pi*dble(i)/dble(n))
     !end if
  end do
  !definiendo la malla de r
  do i=0,n
     r(i)=(x(i)-1.0d0)/coef_escala + r_max
  end do
  
  !=====================================================================

  !==================calculo la matriz de derivada======================
  call matrizderivada()
  !=====================================================================

  !============================preguntando cosas varias=================
  !print*,'harmonicas 1; kerr 2'
  !read*,coordenadas
  !print*,'FC 1; CPBC 2'
  !read*,condicionfrontera
  !============================preguntando cosas varias=================
  
  !==========dando valores iniciales a las funciones====================
  !============================
  if (coordenadas==1) then
     g_rr= (1.0d0 + 2.0d0*M/r)*(1.0d0 + 4.0d0*M**2/r**2)
     g_T= 1.0d0
     alpha= (1.0d0/(1.0d0 + 2.0d0*M/r))*(1.0d0/(1.0d0 + 4.0d0*M**2/r**2))
     log_alpha= log(alpha)
     beta_r= 4.0d0*alpha*M**2/r**2
     K_rr=-(4.0d0*M**2/r**3)*sqrt(alpha)*(2.0d0 + 3.0d0*M/r + 4.0d0*M**2/r**2 +&
          4.0d0*M**3/r**3)
     K_T=(4.0d0*M**2/r**3)*sqrt(alpha)
     f_rrr=4.0d0/r + 7.0d0*M/r**2 + 12.0d0*M**2/r**3 + 20.0d0*M**3/r**4
     f_rT=1.0d0/r
  else if (coordenadas==2)then
     g_rr= 1.0d0 + 2.0d0*M/r
     g_T= 1.0d0
     alpha= 1.0d0! /(1.0d0 + 2.0d0*M/r)
     log_alpha= log(alpha)
     beta_r= (2.0d0*M/r)/(1.0d0 + 2.0d0*M/r)
     K_rr=-(2.0d0*M/r**2)*(1.0d0 + M/r)/sqrt((1.0d0 + 2.0d0*M/r))
     K_T=(2.0d0*M/r**2)/sqrt((1.0d0 + 2.0d0*M/r))
     f_rrr=(1.0d0/r)*(4.0d0 + 7.0d0*M/r)
     f_rT=1.0d0/r
  end if
  
  inr=1.0d0/r
  ene=alpha*g_T*sqrt(g_rr)
  grr=1.0d0/g_rr
  ing_T=1.0d0/g_T
  in6=1.0d0/6.0d0
  
  
  
  !=====================================================================

  !=========Calculando derivadas para los valores iniciales=============
  !derivadas de la norma
  call deriva(log_alpha,dlog_alpha)
  call deriva(dlog_alpha,ddlog_alpha)
  call deriva(beta_r,dbeta_r)
  call deriva(dbeta_r,ddbeta_r)
  
   
  !empezando a cotar el tiempo
  !call cpu_time(start)
  
  !=====================Evolucionando ecuaciones con RK4==================
  do i=0, tiempo

     !guardo valores de las variables en auxiliares
     aux1 = g_rr
     aux2 = g_T
     aux3 = K_rr
     aux4 = K_T
     aux5 = f_rrr
     aux6 = f_rt
     
     !derivadas de las variables fundamentales
     call deriva(g_rr,dg_rr)
     call deriva(g_T,dg_T)
     call deriva(k_rr,dK_rr)
     call deriva(K_T,dK_T)
     call deriva(f_rrr,df_rrr)
     call deriva(f_rT,df_rT)
     
     !=================Guarda los valores de la constriccion==============
     if(mod(i,us)==0) then
        print*,i*ht
        
        !constricciones
        !C_rrr=dg_rr + 8.0d0*g_rr*f_rT*ing_T - 2.0d0*f_rrr 
        !C_rT=dg_T + 2.0d0*g_T*inr - 2.0d0*f_rT
        !C_r=dK_t/g_T + 2.0d0*K_T/(r*g_T) - f_rT*(grr*k_rr + K_T/g_T)/g_T
        CH=df_rT/(g_rr*g_T) - 1.0d0/(2.0d0*g_T*r**2) + f_rT*(2.0d0/r + &
             3.5d0*f_rT/g_T - grr*f_rrr)/(g_rr*g_T) - K_T*(grr*K_rr +&
             K_T/2.0d0*g_T)/g_T
        
        !write(22,*)'#tiempo=',i*ht
        !do j=0,n
        !write(22,*)r(j),g_rr(j),alpha(j)
        !end do
        !write(22,*)
        !write(22,*)
        !norma=0.0d0
        !do j=0,n
         !  norma=norma + CH(j)**2
        !end do
        !norma=norma/10.0d0
        !norma=sqrt(norma)
        !write(11,*)i*ht,norma
        !call evaluaG_rr(CH,0.1d0,1.9d0,11.9d0,i*ht)
        call evaluaH(CH,0.1d0,1.9d0,11.9d0,i*ht)
     end if
     !====================================================================
     
     !calculando los K1's
     k1g_rr = beta_r*dg_rr - 2.0d0*ene*K_rr + 2.0d0*g_rr*dbeta_r
     !k1g_rr(0)=0.0d0

     k1g_T = beta_r*dg_T - 2.0d0*ene*K_T + 2.0d0*beta_r*inr*g_T
     !k1g_T(0)=0.0d0
     
     k1K_rr = beta_r*dK_rr - ene*grr*df_rrr + ene*(2.0d0*grr*f_rrr*(grr*f_rrr + inr -&
          4.0d0*f_rT*ing_T) - 6.0d0*inr*inr + K_rr*(2.0d0*K_T*ing_T - grr*K_rr) - &
          6.0d0*(f_rT*ing_T)**2 - ddlog_alpha - dlog_alpha**2 + (4.0d0*inr - grr*f_rrr)*&
          dlog_alpha) + 2.0d0*K_rr*dbeta_r
     
     K1K_T = beta_r*dK_T - ene*grr*df_rT + ene*(K_T*grr*K_rr + inr**2 - 2.0d0*grr*ing_T*&
          f_rT**2 - f_rT*grr*dlog_alpha) + 2.0d0*beta_r*inr*K_T
     

     k1f_rrr = beta_r*df_rrr - ene*dK_rr + ene*(4.0d0*g_rr*K_T*ing_T*(3.0d0*f_rT*ing_T - &
          grr*f_rrr + 2.0d0*inr - dlog_alpha) - K_rr*(10.0d0*f_rT*ing_T + grr*f_rrr - &
          2.0d0*inr + dlog_alpha)) + 3.0d0*f_rrr*dbeta_r + g_rr*ddbeta_r
     
     
     k1f_rt = beta_r*df_rT - ene*dK_T + ene*(K_T*(2.0d0*f_rT*ing_T - grr*f_rrr - &
          dlog_alpha)) + (dbeta_r + 2.0d0*beta_r*inr)*f_rT
     

     if(condicionfrontera==1) then
        !FC
        k1g_T(0)=0.0d0
        k1g_rr(0)=0.0d0
        k1K_rr(0)=k2f_rrr(0)/sqrt(g_rr(0))
        k1K_T(0)=k2f_rt(0)/sqrt(g_rr(0))
     else if (condicionfrontera==2)then
        !CPBC
        dk_T(0)=df_rT(0) / sqrt(g_rr(0)) - sqrt(g_rr(0)) / (2.0d0*r(0)**2)  +&
             f_rT(0) / sqrt(g_rr(0)) * (2.0d0/r(0) + (3.5d0*f_rT(0)/g_T(0)) -&
             f_rrr(0)/g_rr(0)) - sqrt(g_rr(0))*K_T(0) * (K_rr(0)/g_rr(0) + &
             K_T(0)/(2.0d0*g_T(0))) - 2.0d0*K_T(0) / r(0)  + f_rT(0) *&
             (K_rr(0)/g_rr(0) + K_T(0)/g_T(0))
        k1g_rr(0)=2.0d0*beta_r(0)*f_rrr(0)- 8.0d0*g_rr(0)*f_rT(0)*beta_r(0)*ing_T(0) -&
             2.0d0*ene(0)*K_rr(0) + 2.0d0*g_rr(0)*dbeta_r(0)
        k1g_T(0)=2.0d0*beta_r(0)*f_rT(0) - 2.0d0*ene(0)*K_T(0)
        k1K_T(0)=beta_r(0)*dk_T(0) - ene(0)*grr(0)*df_rT(0) + ene(0)*(K_T(0)*grr(0)*K_rr(0) +&
             inr(0)**2 - 2.0d0*grr(0)*ing_T(0)*&
             f_rT(0)**2 - f_rT(0)*grr(0)*dlog_alpha(0)) + 2.0d0*beta_r(0)*inr(0)*K_T(0)
        k1f_rT(0) = beta_r(0)*df_rT(0) - ene(0)*dk_T(0) + ene(0)*(K_T(0)*(2.0d0*f_rT(0)*&
             ing_T(0) - grr(0)*f_rrr(0) - &
             dlog_alpha(0))) + (dbeta_r(0) + 2.0d0*beta_r(0)*inr(0))*f_rT(0)
        k1K_rr(0)=k1f_rrr(0)/sqrt(g_rr(0))
     end if
     
     
     !paso intermedio
     g_rr = aux1 + ht*k1g_rr*0.5d0
     g_T= aux2 + ht*k1g_T*0.5d0
     K_rr= aux3 + ht*k1K_rr*0.5d0
     K_T= aux4 + ht*k1K_T*0.5d0
     f_rrr= aux5 + ht*k1f_rrr*0.5d0
     f_rT= aux6 + ht*k1f_rT*0.5d0
     ene=alpha*g_T*sqrt(g_rr)
     grr=1.0d0/g_rr
     ing_T=1.0d0/g_T

     !derivadas de las variables fundamentales
     call deriva(g_rr,dg_rr)
     call deriva(g_T,dg_T)
     call deriva(k_rr,dK_rr)
     call deriva(K_T,dK_T)
     call deriva(f_rrr,df_rrr)
     call deriva(f_rT,df_rT)

     !calculando los K2's
     k2g_rr = beta_r*dg_rr - 2.0d0*ene*K_rr + 2.0d0*g_rr*dbeta_r
     

     k2g_T = beta_r*dg_T - 2.0d0*ene*K_T + 2.0d0*beta_r*inr*g_T
     
     
     k2K_rr = beta_r*dK_rr - ene*grr*df_rrr + ene*(2.0d0*grr*f_rrr*(grr*f_rrr + inr -&
          4.0d0*f_rT*ing_T) - 6.0d0*inr*inr + K_rr*(2.0d0*K_T*ing_T - grr*K_rr) - &
          6.0d0*(f_rT*ing_T)**2 - ddlog_alpha - dlog_alpha**2 + (4.0d0*inr - grr*f_rrr)*&
          dlog_alpha) + 2.0d0*K_rr*dbeta_r
     
     K2K_T = beta_r*dK_T - ene*grr*df_rT + ene*(K_T*grr*K_rr + inr**2 - 2.0d0*grr*ing_T*&
          f_rT**2 - f_rT*grr*dlog_alpha) + 2.0d0*beta_r*inr*K_T
     

     k2f_rrr = beta_r*df_rrr - ene*dK_rr + ene*(4.0d0*g_rr*K_T*ing_T*(3.0d0*f_rT*ing_T - &
          grr*f_rrr + 2.0d0*inr - dlog_alpha) - K_rr*(10.0d0*f_rT*ing_T + grr*f_rrr - &
          2.0d0*inr + dlog_alpha)) + 3.0d0*f_rrr*dbeta_r + g_rr*ddbeta_r
     
     
     k2f_rt = beta_r*df_rT - ene*dK_T + ene*(K_T*(2.0d0*f_rT*ing_T - grr*f_rrr - &
          dlog_alpha)) + (dbeta_r + 2.0d0*beta_r*inr)*f_rT
     

     if(condicionfrontera==1) then
        !FC
        k2g_T(0)=0.0d0
        k2g_rr(0)=0.0d0
        k2K_rr(0)=k2f_rrr(0)/sqrt(g_rr(0))
        k2K_T(0)=k2f_rt(0)/sqrt(g_rr(0))
     else if (condicionfrontera==2)then
        !CPBC
        dk_T(0)=df_rT(0) / sqrt(g_rr(0)) - sqrt(g_rr(0)) / (2.0d0*r(0)**2)  +&
             f_rT(0) / sqrt(g_rr(0)) * (2.0d0/r(0) + (3.5d0*f_rT(0)/g_T(0)) -&
             f_rrr(0)/g_rr(0)) - sqrt(g_rr(0))*K_T(0) * (K_rr(0)/g_rr(0) + &
             K_T(0)/(2.0d0*g_T(0))) - 2.0d0*K_T(0) / r(0)  + f_rT(0) *&
             (K_rr(0)/g_rr(0) + K_T(0)/g_T(0))
        k2g_rr(0)=2.0d0*beta_r(0)*f_rrr(0)- 8.0d0*g_rr(0)*f_rT(0)*beta_r(0)*ing_T(0) -&
             2.0d0*ene(0)*K_rr(0) + 2.0d0*g_rr(0)*dbeta_r(0)
        k2g_T(0)=2.0d0*beta_r(0)*f_rT(0) - 2.0d0*ene(0)*K_T(0)
        k2K_T(0)=beta_r(0)*dk_T(0) - ene(0)*grr(0)*df_rT(0) + ene(0)*(K_T(0)*grr(0)*K_rr(0) +&
             inr(0)**2 - 2.0d0*grr(0)*ing_T(0)*&
             f_rT(0)**2 - f_rT(0)*grr(0)*dlog_alpha(0)) + 2.0d0*beta_r(0)*inr(0)*K_T(0)
        k2f_rT(0) = beta_r(0)*df_rT(0) - ene(0)*dk_T(0) + ene(0)*(K_T(0)*(2.0d0*f_rT(0)*&
             ing_T(0) - grr(0)*f_rrr(0) - &
             dlog_alpha(0))) + (dbeta_r(0) + 2.0d0*beta_r(0)*inr(0))*f_rT(0)
        k2K_rr(0)=k1f_rrr(0)/sqrt(g_rr(0))
     end if
     
     
     !paso intermedio
     g_rr = aux1 + ht*k2g_rr*0.5d0
     g_T= aux2 + ht*k2g_T*0.5d0
     K_rr= aux3 + ht*k2K_rr*0.5d0
     K_T= aux4 + ht*k2K_T*0.5d0
     f_rrr= aux5 + ht*k2f_rrr*0.5d0
     f_rT= aux6 + ht*k2f_rT*0.5d0
     ene=alpha*g_T*sqrt(g_rr)
     grr=1.0d0/g_rr
     ing_T=1.0d0/g_T

     !derivadas de las variables fundamentales
     call deriva(g_rr,dg_rr)
     call deriva(g_T,dg_T)
     call deriva(k_rr,dK_rr)
     call deriva(K_T,dK_T)
     call deriva(f_rrr,df_rrr)
     call deriva(f_rT,df_rT)

     !calculando los K3's
     k3g_rr = beta_r*dg_rr - 2.0d0*ene*K_rr + 2.0d0*g_rr*dbeta_r
     !k3g_rr(0)=0.0d0

     k3g_T = beta_r*dg_T - 2.0d0*ene*K_T + 2.0d0*beta_r*inr*g_T
     !k3g_T(0)=0.0d0
     
     k3K_rr = beta_r*dK_rr - ene*grr*df_rrr + ene*(2.0d0*grr*f_rrr*(grr*f_rrr + inr -&
          4.0d0*f_rT*ing_T) - 6.0d0*inr*inr + K_rr*(2.0d0*K_T*ing_T - grr*K_rr) - &
          6.0d0*(f_rT*ing_T)**2 - ddlog_alpha - dlog_alpha**2 + (4.0d0*inr - grr*f_rrr)*&
          dlog_alpha) + 2.0d0*K_rr*dbeta_r
     
     K3K_T = beta_r*dK_T - ene*grr*df_rT + ene*(K_T*grr*K_rr + inr**2 - 2.0d0*grr*ing_T*&
          f_rT**2 - f_rT*grr*dlog_alpha) + 2.0d0*beta_r*inr*K_T
     

     k3f_rrr = beta_r*df_rrr - ene*dK_rr + ene*(4.0d0*g_rr*K_T*ing_T*(3.0d0*f_rT*ing_T - &
          grr*f_rrr + 2.0d0*inr - dlog_alpha) - K_rr*(10.0d0*f_rT*ing_T + grr*f_rrr - &
          2.0d0*inr + dlog_alpha)) + 3.0d0*f_rrr*dbeta_r + g_rr*ddbeta_r
     
     
     k3f_rt = beta_r*df_rT - ene*dK_T + ene*(K_T*(2.0d0*f_rT*ing_T - grr*f_rrr - &
          dlog_alpha)) + (dbeta_r + 2.0d0*beta_r*inr)*f_rT
     

     if(condicionfrontera==1) then
        !FC
        k3g_T(0)=0.0d0
        k3g_rr(0)=0.0d0
        k3K_rr(0)=k2f_rrr(0)/sqrt(g_rr(0))
        k3K_T(0)=k2f_rt(0)/sqrt(g_rr(0))
     else if (condicionfrontera==2)then
        !CPBC
        dk_T(0)=df_rT(0) / sqrt(g_rr(0)) - sqrt(g_rr(0)) / (2.0d0*r(0)**2)  +&
             f_rT(0) / sqrt(g_rr(0)) * (2.0d0/r(0) + (3.5d0*f_rT(0)/g_T(0)) -&
             f_rrr(0)/g_rr(0)) - sqrt(g_rr(0))*K_T(0) * (K_rr(0)/g_rr(0) + &
             K_T(0)/(2.0d0*g_T(0))) - 2.0d0*K_T(0) / r(0)  + f_rT(0) *&
             (K_rr(0)/g_rr(0) + K_T(0)/g_T(0))
        k3g_rr(0)=2.0d0*beta_r(0)*f_rrr(0)- 8.0d0*g_rr(0)*f_rT(0)*beta_r(0)*ing_T(0) -&
             2.0d0*ene(0)*K_rr(0) + 2.0d0*g_rr(0)*dbeta_r(0)
        k3g_T(0)=2.0d0*beta_r(0)*f_rT(0) - 2.0d0*ene(0)*K_T(0)
        k3K_T(0)=beta_r(0)*dk_T(0) - ene(0)*grr(0)*df_rT(0) + ene(0)*(K_T(0)*grr(0)*K_rr(0) +&
             inr(0)**2 - 2.0d0*grr(0)*ing_T(0)*&
             f_rT(0)**2 - f_rT(0)*grr(0)*dlog_alpha(0)) + 2.0d0*beta_r(0)*inr(0)*K_T(0)
        k3f_rT(0) = beta_r(0)*df_rT(0) - ene(0)*dk_T(0) + ene(0)*(K_T(0)*(2.0d0*f_rT(0)*&
             ing_T(0) - grr(0)*f_rrr(0) - &
             dlog_alpha(0))) + (dbeta_r(0) + 2.0d0*beta_r(0)*inr(0))*f_rT(0)
        k3K_rr(0)=k1f_rrr(0)/sqrt(g_rr(0))
     end if
     
     
     !paso intermedio
     g_rr = aux1 + ht*k3g_rr
     g_T= aux2 + ht*k3g_T
     K_rr= aux3 + ht*k3K_rr
     K_T= aux4 + ht*k3K_T
     f_rrr= aux5 + ht*k3f_rrr
     f_rT= aux6 + ht*k3f_rT
     ene=alpha*g_T*sqrt(g_rr)
     grr=1.0d0/g_rr
     ing_T=1.0d0/g_T

     !derivadas de las variables fundamentales
     call deriva(g_rr,dg_rr)
     call deriva(g_T,dg_T)
     call deriva(k_rr,dK_rr)
     call deriva(K_T,dK_T)
     call deriva(f_rrr,df_rrr)
     call deriva(f_rT,df_rT)

     !calculando los K4's
     k4g_rr = beta_r*dg_rr - 2.0d0*ene*K_rr + 2.0d0*g_rr*dbeta_r
     !k4g_rr(0)=0.0d0

     k4g_T = beta_r*dg_T - 2.0d0*ene*K_T + 2.0d0*beta_r*inr*g_T
     !k4g_T(0)=0.0d0
     
     k4K_rr = beta_r*dK_rr - ene*grr*df_rrr + ene*(2.0d0*grr*f_rrr*(grr*f_rrr + inr -&
          4.0d0*f_rT*ing_T) - 6.0d0*inr*inr + K_rr*(2.0d0*K_T*ing_T - grr*K_rr) - &
          6.0d0*(f_rT*ing_T)**2 - ddlog_alpha - dlog_alpha**2 + (4.0d0*inr - grr*f_rrr)*&
          dlog_alpha) + 2.0d0*K_rr*dbeta_r
     
     K4K_T = beta_r*dK_T - ene*grr*df_rT + ene*(K_T*grr*K_rr + inr**2 - 2.0d0*grr*ing_T*&
          f_rT**2 - f_rT*grr*dlog_alpha) + 2.0d0*beta_r*inr*K_T
     

     k4f_rrr = beta_r*df_rrr - ene*dK_rr + ene*(4.0d0*g_rr*K_T*ing_T*(3.0d0*f_rT*ing_T - &
          grr*f_rrr + 2.0d0*inr - dlog_alpha) - K_rr*(10.0d0*f_rT*ing_T + grr*f_rrr - &
          2.0d0*inr + dlog_alpha)) + 3.0d0*f_rrr*dbeta_r + g_rr*ddbeta_r
     
     
     k4f_rt = beta_r*df_rT - ene*dK_T + ene*(K_T*(2.0d0*f_rT*ing_T - grr*f_rrr - &
          dlog_alpha)) + (dbeta_r + 2.0d0*beta_r*inr)*f_rT
     

     if(condicionfrontera==1) then
        !FC
        k4g_T(0)=0.0d0
        k4g_rr(0)=0.0d0
        k4K_rr(0)=k2f_rrr(0)/sqrt(g_rr(0))
        k4K_T(0)=k2f_rt(0)/sqrt(g_rr(0))
     else if (condicionfrontera==2)then
        !CPBC
        dk_T(0)=df_rT(0) / sqrt(g_rr(0)) - sqrt(g_rr(0)) / (2.0d0*r(0)**2)  +&
             f_rT(0) / sqrt(g_rr(0)) * (2.0d0/r(0) + (3.5d0*f_rT(0)/g_T(0)) -&
             f_rrr(0)/g_rr(0)) - sqrt(g_rr(0))*K_T(0) * (K_rr(0)/g_rr(0) + &
             K_T(0)/(2.0d0*g_T(0))) - 2.0d0*K_T(0) / r(0)  + f_rT(0) *&
             (K_rr(0)/g_rr(0) + K_T(0)/g_T(0))
        k4g_rr(0)=2.0d0*beta_r(0)*f_rrr(0)- 8.0d0*g_rr(0)*f_rT(0)*beta_r(0)*ing_T(0) -&
             2.0d0*ene(0)*K_rr(0) + 2.0d0*g_rr(0)*dbeta_r(0)
        k4g_T(0)=2.0d0*beta_r(0)*f_rT(0) - 2.0d0*ene(0)*K_T(0)
        k4K_T(0)=beta_r(0)*dk_T(0) - ene(0)*grr(0)*df_rT(0) + ene(0)*(K_T(0)*grr(0)*K_rr(0) +&
             inr(0)**2 - 2.0d0*grr(0)*ing_T(0)*&
             f_rT(0)**2 - f_rT(0)*grr(0)*dlog_alpha(0)) + 2.0d0*beta_r(0)*inr(0)*K_T(0)
        k4f_rT(0) = beta_r(0)*df_rT(0) - ene(0)*dk_T(0) + ene(0)*(K_T(0)*(2.0d0*f_rT(0)*&
             ing_T(0) - grr(0)*f_rrr(0) - &
             dlog_alpha(0))) + (dbeta_r(0) + 2.0d0*beta_r(0)*inr(0))*f_rT(0)
        k4K_rr(0)=k1f_rrr(0)/sqrt(g_rr(0))
     end if
     
     
     !calculando las variables al nuevo tiempo t=t+ht
     g_rr = aux1 + ht*in6*(k1g_rr + 2.0d0*k2g_rr + 2.0d0*k3g_rr + k4g_rr)
     
     g_T = aux2 + ht*in6*(k1g_T + 2.0d0*k2g_T + 2.0d0*k3g_T + k4g_T)
     
     K_rr = aux3 + ht*in6*(k1K_rr + 2.0d0*k2K_rr + 2.0d0*k3K_rr + k4K_rr)
     
     K_T = aux4 + ht*in6*(k1k_T + 2.0d0*k2k_T + 2.0d0*k3k_T + k4k_T)
     
     f_rrr = aux5 + ht*in6*(k1f_rrr + 2.0d0*k2f_rrr + 2.0d0*k3f_rrr + k4f_rrr)
     
     f_rT = aux6 + ht*in6*(k1f_rT + 2.0d0*k2f_rT + 2.0d0*k3f_rT + k4f_rT)

     ene=alpha*g_T*sqrt(g_rr)
     grr=1.0d0/g_rr
     ing_T=1.0d0/g_T
     
  end do
  !termina de contar
  !call cpu_time(finish)
  !open(1,file='40tiempo.dat')
  !write(1,*)n, finish-start
  !close(1)
  print*,'simulacion terminada'
  
  
  
  
  !close(10)
  !close(11)
  close(20)
  close(21)
  !close(22)
end program e3

subroutine matrizderivada()
  use variables
  implicit none
  integer k,j
  real(8) C
  do k=0,n
     do j=0,n
        if ((k/=j)) then
           Dn(k,j)= (C(k)* (-1)**(j+k)) / (C(j)*(x(k)-x(j)))
        else if ((k==j).and.(k/=0).and.(k/=n))then
           Dn(k,j)=- .5d0*x(k) / (1.0d0 - x(k)**2.0d0)
        end if
     end do
  end do
  Dn(0,0)=((2.0d0 * dble(n)**2)+1.0d0)/6.0d0
  Dn(n,n)=-Dn(0,0)
end subroutine matrizderivada

function C(i)
  use variables
  implicit none
  integer i
  real(8) C
  if ((i==0).or.(i==n)) then
     C=2.0d0
  else
     C=1.0d0
  end if
  return
end function C  

subroutine deriva(funcion,dfuncion)
  use variables
  implicit none
  real(8),dimension(0:n)::funcion,dfuncion
  integer i,j
  do i=0,n
     dfuncion(i)=0.0d0
     do j=0,n
        dfuncion(i)=dfuncion(i)+ Dn(i,j)*funcion(j)
     end do
  end do
  dfuncion=coef_escala*dfuncion
end subroutine deriva

subroutine evaluaH(funcion,h,rmin,rmax,time) 
  use variables
  implicit none
  real(8),dimension(0:n)::funcion
  real(8) z,L,h,rmin,rt,rmax,funcionaprox,normafuncion
  real(8)time
  integer i,j,k
  k=0
  rt=rmin
  normafuncion=0.0d0
  do while(rt<=rmax)
     z=(rt-r_max)*coef_escala + 1.0d0
     funcionaprox=0.0d0
     do i=0,n
        funcionaprox=funcionaprox + funcion(i)*L(i,z)
     end do
     !write(10,*)rt,funcionaprox
     k=k+1
     rt=rt+h
     normafuncion=normafuncion + h*funcionaprox**2
  end do
  normafuncion=sqrt(normafuncion)
  write(20,*)time,normafuncion
end subroutine evaluaH

subroutine evaluaG_rr(funcion,h,rmin,rmax,time) 
  use variables
  implicit none
  real(8),dimension(0:n)::funcion
  real(8) z,L,h,rmin,rt,rmax,funcionaprox,normafuncion
  real(8)time
  integer i,j,k
  k=0
  rt=rmin
  normafuncion=0.0d0
  write(21,*)'#time=',time
  do while(rt<=rmax)
     z=(rt-r_max)*coef_escala + 1.0d0
     funcionaprox=0.0d0
     do i=0,n
        funcionaprox=funcionaprox + funcion(i)*L(i,z)
     end do
     write(21,*)rt,funcionaprox
     k=k+1
     rt=rt+h
  end do
  write(21,*)
  write(21,*)
end subroutine evaluaG_rr

function L(i,z)
  use variables
  implicit none
  real(8) z,L,DT,C
  integer i,k,j
  L=((-1.0d0)**(i+1) *(1.0d0 - z**2)*DT(n,z))/(C(i)* dble(n)**2 *(z-x(i)))
  
  do k=0,n
     if((abs(z-x(k))<=10d-10).and.(k==i)) then
        L=1.0d0
     else if((abs(z-x(k))<=10d-10).and.(k/=i)) then
        L=0.0d0
     end if
  end do
  return
end function L

function T(i,z)  
  implicit none
  integer i
  real(8) z,T
  if (i==0) then
     T=1.0d0
  else if (i==1) then
     T=z
  else
     T=cos(dble(i)*acos(z))
  end if
  return
end function T

function DT(i,z)
  implicit none
  integer i
  real(8) j,z,DT,U
  j=dble(i)
  if (i==0) then
     DT=0.0d0
  else if(i==1) then
     DT=1.0d0
  else
     DT= j*U(i-1,z)
  end if
  return
end function DT

function U(i,z)
  implicit none
  real(8) U,z,j
  integer i
  j=dble(i)
  if(z==1.0d0) then
     U=j+1.0d0
  else if (z==-1.0d0) then
     U=(j+1.0d0)* (-1.0d0)**j
  else
     U=sin((j+1.0)*acos(z))/sin(acos(z))
  end if
end function U
