!#########################################################################################
!#                    PHYCS426 - ADVANCED COMPUTATIONAL PHYSICS                          #
!#                                SPRING 2022                                            #
!#                                                                                       #
!#					From the Ising model to nano-granular systems                        #
!#				                                                                         #
!#		Project completed by:                                                            #
!#			Asif Bin Ayub                                                                #
!#			Kumail Abdulaziz Radhi                                                       #
!#				 		                                                             	 #
!#########################################################################################

!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!<><>							PROGRAM FOR 2D & 3D ISING MODELS					       <><>
!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

! program Ising_2D_3D
	! implicit none
	! !Variables Used:======================
	! ! imc 				--> Montecarlo steps
	! ! kB, JJ, mu 			--> Constants, which we take to be 1
	! ! B,dB, Bmax 			--> Magnetic Field variables
	! ! Eavg, E2avg 		--> Energy variables
	! ! M,Mth,Mavg,M2avg	--> Magnetization variables
	! ! start,finish		--> To record processor time
	! ! S(:,:), S(:,:,:)	--> Spin Grid in 2D,3D
	
	! ! Instructions
	
	! !Declaration:=========================
	! integer :: imc,iB,i,j,k,L
	! double precision :: kB,JJ,mu
	! double precision :: E,Eold,dE,P
	! double precision :: B,dB,Bmax,T,dT,Tmax
	! double precision :: Eavg,E2avg,C
	! double precision :: M,Mth,Mavg,M2avg,X
	! real :: start,finish
	! ! integer, allocatable :: S(:,:)
	! integer, allocatable :: S(:,:,:) ! FOR 3D ONLY
	
	! call cpu_time(start) !Starts measuring the processor time
	
	! call randomize()
	
	! !Open Data Files:======================
	! ! open(unit=10,file="M.txt")
	! ! open(unit=11,file="M-T(L=10)-3D.txt")
	! open(unit=12,file="M-B(L=8)-3D.txt")
	! ! open(unit=13,file="C-T(L=10)-3D.txt")
	! ! open(unit=14,file="X-T(L=10)-3D.txt")
	
	! !Parameters (2D):===========
	! ! ACTIVATE FOR 2D
	! ! JJ=0.5d0 !Interaction Strength
	! ! kB=1.0d0 !Boltzmann Constant
	! ! mu=1.0d0 !Magnetic Moment
	! ! L=50 !Grid Size
	! ! T=0.1d0; dT=0.01d0; Tmax=5.0d0; !J/kB
	! ! B=-0.0d0; dB=0.1d0; Bmax=5.0d0; !J/mu
	
	! !Parameters (3D):===========
	! ! ACTIVATE FOR 3D
	! JJ=1.0d0 	!Interaction Strength
	! kB=1.0d0 	!Boltzmann Constant
	! mu=1.0d0 	!Magnetic Moment
	! L=8			!Grid Size
	! T=0.5d0; dT=0.05d0; Tmax=8.0d0;	!J/kB
	! B=0.0d0; dB=0.1d0; Bmax=5.0d0; 	!J/mu 
	
	! ! ACTIVATE ONE AT A GIVEN TIME
	! ! allocate(S(0:L+1,0:L+1))	!FOR 2D ONLY
	! allocate(S(0:L+1,0:L+1,0:L+1)) ! FOR 3D ONLY
	
	! ! S(:,:)=1 !Setting All Spins To +1
	! S(:,:,:)=1 	!Setting All Spins To +1 [FOR 3D ONLY]
	
	! !Loop Over Temperature & Magnetic Field:=====================
	! do iB=-nint(2*Bmax/dB),nint(2*Bmax/dB)
	! ! For 2D, -5<B<5 
	! ! For 3D, -5<B<5
	
		! B=-abs(dB*iB)+Bmax
	! ! do while(B.le.Bmax) ! For Susceptibility vs B field, ignore
	
	! ! do while(T.le.Tmax)
	! ! For 2D, Tc ~ 2.27	, 0<T<5 
	! ! For 3D, Tc ~ 4.5	, 0<T<8
		
		! ! S(:,:)=1 	!Setting All Spins To +1
		! ! S(:,:,:)=1 	!Setting All Spins To +1 [FOR 3D ONLY --> Switch off for Hysteresis]
		! Mavg=0.0d0 	!Reset <M> "Thermal Average Magnetization"
		! Mth=0.0d0 	!Reset <M> "Thermal Average Magnetization" (Only last 100)
		! Eavg=0.0d0 	!Reset <E> "Thermal Average Energy"
		! M2avg=0.0d0 !Reset <M**2> "Thermal Average Square Magnetization"
		! E2avg=0.0d0 !Reset <E**2> "Thermal Average Square Energy"
		
		! !Loop Over MonteCarlo Steps:--------
		! ! For 2D choose imc ~ 1000-10000
		! ! For 3D choose imc ~ 500-1000
		! do imc=1,500
			
			! !!!!!!! ACTIVATE FOR 2D !!!!!!!!
			
			! !Loop Over Spins:-------------		
			! ! call Spins(L,S,B,T,JJ)
			
			! !Calculate Magnetization Per Spin
			! ! M=sum(S(1:L,1:L))/dble(L**2)
			
			! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			
			
			! !!!!!!! ACTIVATE FOR 3D !!!!!!!!
			
			! !Loop Over Spins (3D):-------------		
			! call Spins3D(L,S,B,T,JJ)
			
			! !Calculate Magnetization Per Spin (3D)
			! M=sum(S(1:L,1:L,1:L))/dble(L**3)
			
			! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			
			
			! !!!!!!! ACTIVATE FOR 2D !!!!!!!!
			
			! ! Heat Capacity (2D):------------
			! ! E=0.0d0 !Reset Energy
			! ! !Spin Loop
			! ! do i=1,L
				! ! do j=1,L
					! ! ! E=E-JJ*(S(i,j+1)+S(i,j-1)+S(i+1,j)+S(i-1,j))*S(i,j)
					
					! ! ! ACTIVATE THIS FOR 3D
					! ! do k=1,L
						! ! E=E-JJ*(S(i,j+1,k)+S(i,j-1,k)+S(i+1,j,k)+S(i-1,j,k)+S(i,j,k+1)+S(i,j,k-1))*S(i,j,k)-B*mu*S(i,j,k)
					! ! enddo
				! ! enddo
			! ! enddo
			! ! E=E/2.0d0 !Because of Overcounting
			! ! Eavg=Eavg+E/500
			! ! E2avg=E2avg+E**2/500
			
			! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

			! ! Magnetic Susceptibility:-------
			! ! Mavg=Mavg+M/500
			! ! M2avg=M2avg+M**2/500
			
			! !Writing Results
			! ! write(10,*) mc,M
			
			! !Calculate Thermal Average Magnetization
			! if(imc.gt.400) Mth=Mth+M/100
		! enddo
		
		! !Calculate Heat Capacity
		! ! C=(E2avg-Eavg**2)**2/(kB*T**2)
		! ! C=C/dble(L**3)
		
		! !Calculate Magnetic Susceptibility
		! ! X=(M2avg-Mavg**2)**2/(kB*T)
		
		! !Writing Results
		! write(11,*) T,Mth,1/Mth
		! write(12,*) B,Mth
		! ! write(13,*) T,C,1/C
		! ! write(14,*) B,X
		! ! write(14,*) T,X,1/X
		
		! !Update T,B
		! ! T=T+dT
		! B=B+dB
	! enddo
	
	
	
	! !Deallocation & Closing:========================
	! deallocate(S)
	! close(10)
	! close(11)
	! close(12)
	! close(13)
	! close(14)
	
	! call cpu_time(finish)
	! print *,"Execution time:",finish-start,"seconds"
	! print *,"Start/Finish",start,"/",finish
	
	! contains
	
	! !Spins Loop:==================
	! subroutine Spins(L,S,B,T,JJ)
		! implicit none
		! integer :: i,j,L,S(0:L+1,0:L+1)
		! double precision, intent(in) :: B,T,JJ
		! double precision :: dE,Eold,P,kB=1.0d0,mu=1.0d0,rnd
		
		! do i=1,L
			! do j=1,L
				! !Periodic Boundary Condition
				! S(0,:)=S(L,:)
				! S(:,0)=S(:,L)
				! S(L+1,:)=S(1,:)
				! S(:,L+1)=S(:,1)
				
				! !Calculate dE
				! Eold=-JJ*(S(i,j+1)+S(i,j-1)+S(i+1,j)+S(i-1,j))*S(i,j)-B*mu*S(i,j)
				! dE=-2*Eold !dE=Enew-Eold, Enew=-Eold
	
				! !Probability
				! P=exp(-dE/(kB*T))
				! call random_number(rnd)
				
				! !Checking Flip:
				! if(dE.lt.0)then
					! S(i,j)=-S(i,j) !Accept Flip 
				! else 
					! if(P.gt.rnd) S(i,j)=-S(i,j) !Accept Flip
				! endif
			! enddo
		! enddo
	! end subroutine Spins
	
	! !Spins Loop (3D):==================
	! subroutine Spins3D(L,S,B,T,JJ)
		! implicit none
		! integer :: i,j,k,L,S(0:L+1,0:L+1,0:L+1)
		! double precision, intent(in) :: B,T,JJ
		! double precision :: dE,Eold,P,kB=1.0d0,mu=1.0d0,rnd
		
		! do i=1,L
			! do j=1,L
				! do k=1,L
					! !Periodic Boundary Condition
					! S(0,:,:)=S(L,:,:)
					! S(:,0,:)=S(:,L,:)
					! S(:,:,0)=S(:,:,L)
					! S(L+1,:,:)=S(1,:,:)
					! S(:,L+1,:)=S(:,1,:)
					! S(:,:,L+1)=S(:,:,1)
						
					! !Calculate dE
					! Eold=-JJ*(S(i,j+1,k)+S(i,j-1,k)+S(i+1,j,k)+S(i-1,j,k)+S(i,j,k+1)+S(i,j,k-1))*S(i,j,k)-B*mu*S(i,j,k)
					! dE=-2*Eold !dE=Enew-Eold, Enew=-Eold
		
					! !Probability
					! P=exp(-dE/(kB*T))
					! call random_number(rnd)
					
					! !Checking Flip:
					! if(dE.lt.0)then
						! S(i,j,k)=-S(i,j,k) !Accept Flip 
					! else 
						! if(P.gt.rnd) S(i,j,k)=-S(i,j,k) !Accept Flip
					! endif
				! enddo
			! enddo
		! enddo
	! end subroutine Spins3D
	
	! !Randomize Subroutine:===================
	! subroutine randomize()
	! implicit none
	! integer :: i,n,clock
	! integer(kind=8) :: cc
	! integer, allocatable :: seed(:)
	
	! call random_seed(size=n)
	! allocate(seed(n))
	! call system_clock(count=clock,count_rate=cc)
	! ! print *, clock
	! seed=clock+37*(/(i-1,i=1,n)/)
	! call random_seed(put=seed)
	! deallocate(seed)
	
	! end subroutine randomize
	
! end program Ising_2D_3D



!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!<><>							PROGRAM FOR GENERAL SPIN DIRECTION					       	<><>
!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

program GeneralSpinDirection
	implicit none
	
	!Declaration:=============================================================
	integer :: imc,i,j,L,iB,counter=0
	double precision :: kB,JJ,mu,pi
	double precision :: E,Eold,Enew,dE,P,rnd
	double precision :: B,dB,Bmax,T,dT,Tmax
	real :: start,finish
	double precision, allocatable :: S(:,:)
	
	call cpu_time(start) ! Starts measuring the processor time
	
	call randomize()
	
	!Parameters:=============================================================
	pi=acos(-1.0d0)
	JJ=1.0d0 	!Interaction Strength
	kB=1.0d0 	!Boltzmann Constant
	mu=1.0d0 	!Magnetic Moment
	L=20 		!Grid Size
	T=0.5d0; dT=0.05d0; Tmax=4.0d0; !J/kB
	B=-3.0d0; dB=0.2d0; Bmax=3.0d0; !J/mu 
	
	allocate(S(0:L+1,0:L+1))
	
	!Assignment General Direction for Spins on 2D
	do i=1,L
		do j=1,L
			call random_number(rnd)
			S(i,j)=rnd*2*pi
		enddo
	enddo

	!Loop Over Temperature & Magnetic Field:=================================
	do iB=-nint(2*Bmax/dB),nint(2*Bmax/dB)
		B=-abs(dB*iB)+Bmax
	! do while(B.le.Bmax)
	! do while(T.le.Tmax)

		!Assignment General Direction for Spins on 2D:--------
		do i=1,L
			do j=1,L
				call random_number(rnd)
				S(i,j)=rnd*2*pi
			enddo
		enddo

		!Loop Over MonteCarlo Steps:--------------------------
		do imc=1,1000
			
			!Loop Over Spins:		
			call Spins(L,S,B,T,JJ)
			
		enddo
		
		!Writing Results:-------------------------------------
		counter=counter+1
		do i=1,L
			do j=1,L
				write(counter+10,*) i,j,S(i,j)
			enddo
		enddo

		!Update T,B:------------------------------------------
		! T=T+dT
		! B=B+dB
	enddo
	
	!Deallocation & Closing:==================================================
	deallocate(S)
	
	call cpu_time(finish)
	print *,"Execution time:",finish-start,"seconds"
	print *,"Start/Finish",start,"/",finish
	
	contains
	
	!Spins Loop:==================
	subroutine Spins(L,S,B,T,JJ)
	implicit none
	integer :: i,j,L
	double precision :: S(0:L+1,0:L+1),Snew
	double precision, intent(in) :: B,T,JJ
	double precision :: dE,Eold,P,kB=1.0d0,mu=1.0d0,rnd,rnd1
		
		do i=1,L
			do j=1,L
				!Periodic Boundary Condition
				S(0,:)=S(L,:)
				S(:,0)=S(:,L)
				S(L+1,:)=S(1,:)
				S(:,L+1)=S(:,1)
				
				call random_number(rnd1)

				!Calculate dE
				Eold=-JJ*(cos(S(i,j)-S(i+1,j))+cos(S(i,j)-S(i-1,j))+cos(S(i,j)-S(i,j+1))+cos(S(i,j)-S(i,j-1)))-B*mu*sin(S(i,j))
				Snew=rnd1*2.0d0*pi !The new direction of S(i,j)
				Enew=-JJ*(cos(Snew-S(i+1,j))+cos(Snew-S(i-1,j))+cos(Snew-S(i,j+1))+cos(Snew-S(i,j-1)))-B*mu*sin(Snew)
				dE=Enew-Eold
	
				!Probability
				P=exp(-dE/(kB*T))
				call random_number(rnd)
				
				!Checking Flip:
				if(dE.lt.0)then
					S(i,j)=Snew !Accept Flip 
				else 
					if(P.gt.rnd) S(i,j)=Snew !Accept Flip
				endif
			enddo
		enddo
	end subroutine
	
	!Randomize Subroutine:===================
	subroutine randomize()
	implicit none
	integer :: i,n,clock
	integer(kind=8) :: cc
	integer, allocatable :: seed(:)
	
	call random_seed(size=n)
	allocate(seed(n))
	call system_clock(count=clock,count_rate=cc)
	! print *, clock
	seed=clock+37*(/(i-1,i=1,n)/)
	call random_seed(put=seed)
	deallocate(seed)
	
	end subroutine randomize
	
end program GeneralSpinDirection

! SOME COMMENTS ABOUT THE CODE...
!
!	As we write this, there is less than 20 minutes to the presentation.
!	The code presented here can be arranged in a much better fashon.
!	Ideally, we would have a subroutine for each property, and all inputs
!	can be changed easily.
!
!	But due to lack of time, we ended up with this unorganized mess.
!	And we regret not being able to study the other modifications to Ising
!	model. 
