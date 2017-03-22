subroutine evalua(funcion,h,rmin,rmax,time,k) 
  use variables
  implicit none
  real(8),dimension(0:n)::funcion
  real(8) z,L,h,rmin,rt,rmax,funcionaprox,normafuncion
  real(8)time
  integer i,j,k
  rt=rmin
  normafuncion=0.0d0
  do while(rt<=rmax)
     z=(rt-r_max)*coef_escala + 1.0d0
     funcionaprox=0.0d0
     do i=0,n
        funcionaprox=funcionaprox + funcion(i)*L(i,z)
     end do
     !write(10,*)rt,funcionaprox
     rt=rt+h
     normafuncion=normafuncion + h*funcionaprox**2
  end do
  normafuncion=sqrt(normafuncion)
  write(k,*)time,normafuncion
end subroutine evalua

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

function L(i,z)
  use variables
  implicit none
  real(8) z,L,DT,C
  integer i,k,j
  L=((-1.0d0)**(i+1) *(1.0d0 - z**2)*DT(n,z))/(C(i)* dble(n)**2 *(z-x(i)))
  
  do k=0,n
     if((abs(z-x(k))<=10d-16).and.(k==i)) then
        L=1.0d0
     else if((abs(z-x(k))<=10d-16).and.(k/=i)) then
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
