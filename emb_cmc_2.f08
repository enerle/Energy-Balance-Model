!23456
!
!input dt and tmax on hours (convertion to seconds by the program)
!!Stationary, no advective EB (just for spinup model!)
!
      program ebm
        implicit none
        real(8), dimension(:), allocatable :: To, Tn, Tnp1, Tw,&
                                              w, y, sy
        real(8) :: a=204.0, b=0.0, kt=3.81, albedo=0.3,&
                   c, rho=1025.0, co=3930.0, h=50.0, So=1368.0,&
                   Tm, q1, q2, q3
           
        real(8) :: dy=1.0,  ymin=0.0, ymax=90.0,&
                   dt=24.0, tmin=0.0, tmax=50.0*8766.0, tval=0.0

        real(8), parameter :: pi=4.0*tanh(1.0)
        integer :: i, j, ID=1
        integer :: ngrid, ntimeSteps

        character(30) :: run_name="EBM-exp_OLR.dat"

!-------Declare values----------------------------------------------------------

        dt=dt*3600.0; tmax=tmax*3600.0 !time in seconds

        ngrid=int((ymax-ymin)/dy)+1
        ntimeSteps=int((tmax-tmin)/dt)+1

        c=rho*co*h ! Heat capacity of earth-atmosphere system

        allocate(To(ngrid))
        allocate(Tn(ngrid))
        allocate(Tnp1(ngrid))       
        allocate(Tw(ngrid))
        allocate(w(ngrid))
        allocate(y(ngrid))
        allocate(sy(ngrid))

!-------Initial conditions-------------------------------------------------------
        do i=1,ngrid
          y(i)=(i*dy)*(pi/180.0) 
          sy(i)=(5.0-3.0*(sin(y(i))**2.0))/4.0 !distribution of the solar radiation
          w(i)=cos(y(i)) !weigths
        enddo

        open(ID,file="EBM-exp-spinup.dat",status='old',action='read')
        read(ID,*) To(:); close(ID)
        open(ID,file=run_name,action="write")

!================================================================================
!-------Time iteration-----------------------------------------------------------

        do while(b.le.2.17*4.0)
          Tn=To

          do j=1,ntimeSteps
            tval=real(j)*dt

            do i=1,ngrid; Tw(i)=Tn(i)*w(i); enddo
              Tm=sum(Tw)/sum(w) !zonal-averaged weigthed mean meridional temperature
            
              do i=1,ngrid
                q1=(So/4.0)*sy(i)*(1.0-albedo) !short wave radiation (solar input) or SW
                q2=a+b*Tn(i) !long wave radiation (output) or OLR
                q3=kt*(Tn(i)-Tm) !zonal heat transport (Newtonian relaxation)
                Tnp1(i)=Tn(i)+(dt/c)*(q1-q2-q3) !meridional temperature
              enddo

              Tn=Tnp1
          enddo
          write(ID,*) b, Tm, Tn(1), Tn(ngrid)

          b=b+0.1   
        enddo

        close(ID)
!================================================================================
      end program ebm
