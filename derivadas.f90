subroutine derivadas()
  use variables
  implicit none
  call deriva(g_rr,dg_rr)
  call deriva(g_T,dg_T)
  call deriva(k_rr,dK_rr)
  call deriva(K_T,dK_T)
  call deriva(f_rrr,df_rrr)
  call deriva(f_rT,df_rT)
end subroutine derivadas

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
