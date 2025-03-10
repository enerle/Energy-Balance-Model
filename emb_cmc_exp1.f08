!23456
!
!input dt and tmax on hours (convertion to seconds by the program)
!!Stationary, no advective EB (just for spinup model!)
!
      program ebm
        implicit none
        real(8), dimension(:), allocatable :: To, Tn, Tnp1, Tw,&
                                              w, y, sy
        real(8) :: a=204.0, b=2.17, kt=3.81, albedo=0.3,&
                   c, rho=1025.0, co=3930.0, h=50.0,&
                   So=1368.0, r=0.0, eps=3.0/4.0, sigma=5.67E-8,&
                   eeps, lambda1, lambda2, Ts, Tm, q1, q2, q3,&
                   TcritVal=-10.15
           
        real(8) :: dy=1.0,  ymin=0.0, ymax=90.0,&
                   dt=24.0, tmin=0.0, tmax=50.0*8766.0, tval=0.0,&
                   dr=2.0,  rmin,     rmax 

        real(8), parameter :: pi=4.0*tanh(1.0)
        integer :: i, j, ID=1
        integer :: ngrid, ntimeSteps

        character(30) :: run_name

        logical :: spinup_flag=.false. !After spinup run use (.false.)      
        logical :: p_albedo_flag=.true. !Constant planetary albedo
        logical :: ice_albedo_flag=.false. !Temperture dependence albedo
        logical :: snowice_albedo_flag=.false. !Cuadratic change of albedo between thresholds
        logical :: TcritVal0_flag=.false.
        logical :: TcritVal15_flag=.false.
        logical :: snowball_flag=.false. !Ice covered earth
        logical :: kt0_flag=.false. !extreme values for transport coefficient
        logical :: kt2_flag=.false.
     
!-------Type of experiment---------------------------------------------------------
        if(spinup_flag) then
          rmin=0.0; rmax=0.0
        else
          rmin=-200.0; rmax=4.0
        endif

        if(ice_albedo_flag) then !each simulation includes toap analysis
          run_name="EBM-exp_ice-albedo.dat"
        endif

        if(snowice_albedo_flag) then
          run_name="EBM-exp_snowice-albedo.dat"
        endif
  
        if(p_albedo_flag) then
          run_name="EBM-exp_p-albedo.dat"
        endif

        if(snowball_flag) then
          run_name="EBM-exp_snowball.dat"
          snowice_albedo_flag=.true.
          rmax=50.0
        endif

        if(TcritVal0_flag) then
          run_name="EBM-exp_TcritVal_0.0.dat"; TcritVal=0.0
          ice_albedo_flag=.true.
        elseif(TcritVal15_flag) then
          run_name="EBM-exp_TcritVal_15.15.dat"; TcritVal=-15.15
          ice_albedo_flag=.true.
        endif

        if(kt0_flag) then
          run_name="EBM-exp_kt_0xkt.dat"; kt=kt*0.0
          ice_albedo_flag=.true.          
        elseif(kt2_flag) then
          run_name="EBM-exp_kt_2xkt.dat"; kt=kt*2.0          
          ice_albedo_flag=.true.
        endif

!-------Declare values----------------------------------------------------------

        dt=dt*3600.0; tmax=tmax*3600.0 !time in seconds

        ngrid=int((ymax-ymin)/dy)+1
        ntimeSteps=int((tmax-tmin)/dt)+1

        c=rho*co*h ! Heat capacity of earth-atmosphere system
        r=r+rmin

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

        if(spinup_flag) then
          do i=1,ngrid; To(i)=(-a+(So*sy(i)*(1.0-albedo))/4.0)/b; enddo
          open(ID,file="EBM-exp-spinup.dat",action="write")
        else
          open(ID,file="EBM-exp-spinup.dat",status='old',action='read')
          read(ID,*) To(:); close(ID)
          open(ID,file=run_name,action="write")
        endif

!================================================================================
!-------Time iteration--------------------------------------

        if(snowball_flag) then;  To(:)=-13.15; endif

        do while(r.le.rmax)
          Tn=To

          do j=1,ntimeSteps
            tval=real(j)*dt

            do i=1,ngrid; Tw(i)=Tn(i)*w(i); enddo
            Tm=sum(Tw)/sum(w) !zonal-averaged weigthed mean meridional temperature
            
            do i=1,ngrid
              if(ice_albedo_flag) then
                if(Tn(i).le.TcritVal) then
                  albedo=0.62
                else 
                  albedo=0.3
                endif
              endif

              if(snowice_albedo_flag) then
                if(Tn(i).le.-13.15) then
                  albedo=0.7 !ice-covered
                elseif(Tn(i).gt.-13.15 .and. Tn(i).lt.19.85) then
                  albedo=0.289+(0.411)*(((Tn(i)-19.85)**2.0)/(-33**2))
                elseif(Tn(i).ge.19.85) then
                  albedo=0.289 !ice-free
                endif
              endif
              q1=(So/4.0)*sy(i)*(1.0-albedo) !short wave radiation (solar input) or SW
              q2=a+b*Tn(i) !long wave radiation (output) or OLR
              q3=kt*(Tn(i)-Tm) !zonal heat transport (Newtonian relaxation)
              Tnp1(i)=Tn(i)+(dt/c)*(q1-q2-q3+r) !meridional temperature
            enddo

            Tn=Tnp1
          enddo

          if(spinup_flag) then
            write(ID,*) Tn
          else
            eeps=(2.0-eps)/2.0 !effective emissivity
            Ts=sqrt(sqrt((r+So*0.7)/(2.0*sigma*(2.0-eps)))) !Surface temperature
            lambda1=-4.0*eeps*sigma*Ts**3.0 !Planck feedback
            lambda2=-4.0*eeps*sigma*(Tm+273.15)**3.0 !Planck feedback
            write(ID,*) r, lambda1, lambda2, Ts-273.15, Tm,  Tn
          endif

          r=r+dr   
        enddo

        close(ID)
!================================================================================
      end program ebm
