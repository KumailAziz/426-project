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
	JJ=1.0d0 !Interaction Strength
	kB=1.0d0 !Boltzmann Constant
	mu=1.0d0 !Magnetic Moment
	L=100 !Grid Size
	T=0.5d0; dT=0.05d0; Tmax=4.0d0; !J/kB
	B=-5.0d0; dB=0.1d0; Bmax=5.0d0; !J/mu 
	
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

