program main
  implicit none

  ! real, parameter :: RUNTIME = 1.0

  integer, parameter :: nParticles = 100
  integer, parameter :: tNumSteps = 10000
  ! real :: dt = RUNTIME / tNumSteps
  real :: dt = .00005
  integer :: i,time
  real, dimension (tNumSteps,nParticles,2) :: x
  real, dimension (nParticles) :: r
  real, dimension (nParticles,2) :: v
  real, dimension (2) :: V_AVE

  V_AVE = [0,0]
  ! character (len = 60) :: filepath
  ! filepath = "/Users/Jeff/Desktop/python/Helpful Code/Particle Box/Fortran"
  ! initialize particles
  call srand(1723)
  do i = 1,nParticles,1
    r(i) = .007 !(.01)*(2.**(1./6.))!
    x(1,i,:) = [rand(),rand()]
    v(i,:) = V_AVE + [100*rand()*((-1) ** aint(10*rand())),100*rand()*((-1) ** aint(10*rand()))]
  end do

  ! run through time steps
  do time = 1,tNumSteps,1
    do i = 1,nParticles,1
      ! print *,'time=',time,'i=',i
      call timeStep(time,i,x,v,r,nParticles,tNumSteps,dt)
      ! call sumEnergy(v,nParticles)
    end do
    print *, (real(time)/real(tNumSteps)*100),'% Complete'
  end do

  ! Write x data to file
  open(1, file = '/Users/Jeff/Desktop/python/Helpful Code/Particle Box/Fortran/xData.txt', status = 'old')
  open(2, file = '/Users/Jeff/Desktop/python/Helpful Code/Particle Box/Fortran/yData.txt', status = 'old')
    do time = 1,tNumSteps,1
      ! do i = 1,nParticles,1
        write(1,*) x(time,:,1)
        write(2,*) x(time,:,2)
      ! end do
    end do

   close(1)
   close(2)

end program main



!!!!!!!!!!!!!!!!!
subroutine timeStep(time,i,x,v,r,nParticles,tNumSteps,dt) !(integer(which timestep),which particle,positions,velocities)
  implicit none

  integer :: nParticles,tNumSteps
  integer :: i, time
  real :: dt
  real, dimension (nParticles) :: r
  real, dimension (tNumSteps,nParticles,2) :: x
  real, dimension (nParticles,2) :: v

  call hitparticle(time,i,x,v,r,nParticles,tNumSteps)
  ! call periodicBoundary(time,i,x,v,r,nParticles,tNumSteps)
  call wall(time,i,x,v,r,nParticles,tNumSteps)
  ! call wallhole(time,i,x,v,r,nParticles,tNumSteps)
  x(time+1,i,:) = x(time,i,:) + (v(i,:)*dt)

end subroutine timeStep


!!!!!!!!!!!!!!
subroutine hitParticle(time,i,x,v,r,nParticles,tNumSteps)
  implicit none

  integer :: i,j,time,nParticles,tNumSteps
  real :: vDist, dist
  real, dimension (nParticles) :: r
  real, dimension (tNumSteps,nParticles,2) :: x
  real, dimension (nParticles,2) :: v
  real, dimension (2) :: dx, dv, V_COM, g, gPrime

  do j = 1,nParticles,1
      if(i/=j) then
        dx = x(time,j,:) - x(time,i,:)
        dv = v(j,:) - v(i,:)
        dist = sqrt(dot_product(dx,dx))
        vDist = dot_product(dx,dv)
        if(dist<(r(i)+r(j)))then !if closer together than  radii
          if(vDist<0)then ! if derivitative of dist is negative
            ! compute new velocities
            V_COM = .5*(v(i,:)+v(j,:))
            g = .5*(v(i,:)-v(j,:))
            gPrime = g - (2*(dot_product(g,dx)*dx)/dot_product(dx,dx))
            ! update velocities
            v(j,:) = V_COM-gprime
            v(i,:) = V_COM+gprime
          end if
        end if
      end if
  end do
end subroutine hitParticle


!!!!!!!!!!!!!!!!!
subroutine periodicBoundary(time,i,x,v,r,nParticles,tNumSteps)
  implicit none

  integer :: i,time,nParticles,tNumSteps
  real, dimension (nParticles) :: r
  real, dimension (tNumSteps,nParticles,2) :: x
  real, dimension (nParticles,2) :: v

  ! Top
  if (x(time,i,2)+r(i)>=1 .and. v(i,2)>0)then
    x(time,i,2) = x(time,i,2)-1
  end if
  !Bottom
  if (x(time,i,2)-r(i)<=0 .and. v(i,2)<0)then
    x(time,i,2) = x(time,i,2)+1
  end if
  !Right
  if (x(time,i,1)+r(i)>=1 .and. v(i,1)>0)then
    x(time,i,1) = (x(time,i,1)-1)
    ! print *,(x(time,i,1)-1)
  end if
  !Left
  if (x(time,i,1)-r(i)<=0 .and. v(i,1)<0)then
    x(time,i,1) = x(time,i,1)+1
  end if
end subroutine periodicBoundary

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine wall(time,i,x,v,r,nParticles,tNumSteps)
  implicit none

  integer :: i,time,nParticles,tNumSteps
  real, dimension (nParticles) :: r
  real, dimension (tNumSteps,nParticles,2) :: x
  real, dimension (nParticles,2) :: v

  ! Top
  if (x(time,i,2)+r(i)>=1 .and. v(i,2)>0)then
    v(i,2) = -v(i,2)
  end if
  !Bottom
  if (x(time,i,2)-r(i)<=0 .and. v(i,2)<0)then
    v(i,2) = -v(i,2)
  end if
  !Right
  if (x(time,i,1)+r(i)>=1 .and. v(i,1)>0)then
    v(i,1) = -v(i,1)
    ! print *,(x(time,i,1)-1)
  end if
  !Left
  if (x(time,i,1)-r(i)<=0 .and.  v(i,1)<0)then
    v(i,1) = -v(i,1)
  end if
end subroutine wall

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine wallhole(time,i,x,v,r,nParticles,tNumSteps)
  implicit none

  integer :: i,time,nParticles,tNumSteps
  real, dimension (nParticles) :: r
  real, dimension (tNumSteps,nParticles,2) :: x
  real, dimension (nParticles,2) :: v

  ! Top
  if (x(time,i,2)+r(i)>=1 .and. v(i,2)>0)then
    v(i,2) = -v(i,2)
  end if
  !Bottom
  if (x(time,i,2)-r(i)<=0 .and. v(i,2)<0)then
    v(i,2) = -v(i,2)
  end if
  !Right
  if (x(time,i,1)+r(i)>=.5 .and. x(time,i,1)+r(i)<.55 .and. v(i,1)>0)then
    if (x(time,i,2)+r(i)>.54 .or. x(time,i,2)+r(i)<.46) then
      v(i,1) = -v(i,1)
    end if
  end if
  !Left
  if (x(time,i,1)-r(i)<=0 .and. v(i,1)<0)then
    v(i,1) = -v(i,1)
  end if
end subroutine wallhole


! sum Energy
subroutine sumEnergy(v,nParticles)
  implicit none

  integer :: nParticles,i
  real :: E,Msum
  real,dimension(2) :: M
  real, dimension(nParticles,2) :: v
  E = 0
  M = [0,0]
  do i = 1,nParticles,1
    E = E + (.5 * dot_product(v(i,:),v(i,:))) ! m is arbitrary
    M = M + [v(i,1),v(i,2)]
  end do
  Msum = sqrt(dot_product(M,M))

  print *,'E = ',E,'    M = ',M

end subroutine sumEnergy
