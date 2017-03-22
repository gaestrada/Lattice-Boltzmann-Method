module variables
  implicit none
  real(8),allocatable,dimension(:)::x,r
  real(8),allocatable,dimension(:,:)::Dn
  real(8),allocatable,dimension(:)::g_rr,g_T,K_rr,K_T,f_rrr,f_rT
  real(8),allocatable,dimension(:)::dg_rr,dg_T,dK_rr,dK_T,df_rrr,df_rT
  real(8),allocatable,dimension(:)::alpha,log_alpha,beta_r
  real(8),allocatable,dimension(:)::dbeta_r,ddbeta_r,dlog_alpha,ddlog_alpha
  real(8),allocatable,dimension(:)::aux1,aux2,aux3,aux4,aux5,aux6
  real(8),allocatable,dimension(:)::ene,grr,ing_T,inr
  real(8),allocatable,dimension(:)::CH,C_r,C_rrr,C_rT
  real(8),allocatable,dimension(:,:)::kg_rr,kg_T,kK_rr,kK_T,kf_rrr,kf_rT
  real(8) coef_escala,r_min,r_max
  integer n
  real(8) M,error,norma
  real(8) pi,ht,in6,ttiempo,unseg,tiempo2,tiempodeguardado
  integer tiempo,us
end module variables

program e3
  use variables
  implicit none
  real(8) T,U,DT,L
  integer i,j
  character(10) nfile,nfile1,nfile2,nfile3,nfile4,bla

  !aqui van los parametros
  read*,n
  read*,r_min
  read*,r_max
  read*,M
  read*,tiempo2
  read*,tiempodeguardado
  read*,nfile
  
  !definiendo parametros
  ht=1.0d0/real(n)**2
  pi=acos(-1.0d0)
  coef_escala=2.0d0/(r_max-r_min)
  in6=1.0d0/6.0d0
  ttiempo=tiempo2/ht
  tiempo=int(ttiempo)
  unseg=tiempodeguardado/ht
  us=int(unseg)
  bla='variables'
  nfile1=nfile//bla
  print*,nfile1,nfile//bla
  stop
  nfile2=nfile//'H.dat'
  nfile3=nfile // 'M.dat'
  
  
  !abriendo archivos de datos
  !open(10,file='60g_rr.dat')
  open(11,file=nfile1)
  open(20,file=nfile2)
  open(21,file='12Mnorma2.dat')
  open(22,file='12otras.dat')

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
  allocate(kg_rr(1:4,0:n),kg_T(1:4,0:n),kK_rr(1:4,0:n),kK_T(1:4,0:n))
  allocate(kf_rrr(1:4,0:n),kf_rT(1:4,0:n))

  !======definiendo los puntos de colocacion en las mallas==============
  !definiendo la malla x
  do i=0,n
     if (i==(n/2)) then
        x(i)=0.0d0
     else     
        x(i)=cos(pi*dble(i)/dble(n))
     end if
  end do
  !definiendo la malla de r
  do i=0,n
     r(i)=(x(i)-1.0d0)/coef_escala + r_max
  end do
  !=====================================================================

  !==================calculo la matriz de derivada======================
  call matrizderivada()
  !=====================================================================

  !==========dando valores iniciales a las funciones====================
  g_rr= 1.0d0 + 2.0d0*M/r
  g_T= 1.0d0
  alpha= 1.0d0/(1.0d0 + 2.0d0*M/r)
  log_alpha= log(alpha)
  beta_r= (2.0d0*M/r)/(1.0d0 + 2.0d0*M/r)
  K_rr=-(2.0d0*M/r**2)*(1.0d0 + M/r)/sqrt((1.0d0 + 2.0d0*M/r))
  K_T=(2.0d0*M/r**2)/sqrt((1.0d0 + 2.0d0*M/r))
  f_rrr=(1.0d0/r)*(4.0d0 + 7.0d0*M/r)
  f_rT=1.0d0/r
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
  
  
   
       
  
  !=====================Evolucionando ==================================
  do i=0, tiempo

     call guarda(i)
     
     call RK4()
     
  end do
  
  
  
  
  
  !close(10)
  close(11)
  close(20)
  close(21)
  close(22)
end program e3
