program Skeleton_2D
implicit none

!Declaration:=========================
integer :: imc,i,j,L
double precision :: kB,JJ,mu
double precision :: E,Eold,dE,M,Mavg,P
double precision :: B,dB,Bmax,T,dT,Tmax
double precision :: Eavg,E2avg,C
integer, allocatable :: S(:,:)

!Open Data File:======================
! open(unit=10,file="M.txt")
! open(unit=11,file="M-T(L=10).txt")
! open(unit=12,file="M-B(L=10).txt")
open(unit=13,file="C-T(L=10).txt")

!Parameters:===========
JJ=1.0d0 !Interaction Strength
kB=1.0d0 !Boltzmann Constant
mu=1.0d0 !Magnetic Moment
L=10 !Grid Size
T=0.5d0; dT=0.1d0; Tmax=5.0d0; !J/kB 
B=0.0d0; dB=0.1d0; Bmax=5.0d0; !J/mu 

allocate(S(0:L+1,0:L+1))

! S(:,:)=1 !Setting All Spins To +1

!Loop Over Temperature:=====================
! do while(B.le.Bmax)
do while(T.le.Tmax)
	S(:,:)=1 !Setting All Spins To +1
	Mavg=0.0d0 !Reset <M> "Thermal Average Magnetization"
	Eavg=0.0d0 !Reset <E> "Thermal Average Energy"
	E2avg=0.0d0 !Reset <E**2> "Thermal Average Square Energy"
	
	!Loop Over MonteCarlo Steps:--------
	do imc=1,1000
		!Loop Over Spins:-------------
		do i=1,L
			do j=1,L
				!Periodic Boundary Condition
				S(0,:)=S(L,:)
				S(:,0)=S(:,L)
				S(L+1,:)=S(1,:)
				S(:,L+1)=S(:,1)
				
				!Calculate dE
				Eold=-JJ*(S(i,j+1)+S(i,j-1)+S(i+1,j)+S(i-1,j))*S(i,j)-B*mu*S(i,j)
				dE=-2*Eold !dE=Enew-Eold, Enew=-Eold

				!Probability
				P=exp(-dE/(kB*T))
				
				!Checking Flip:
				if(dE.lt.0)then
					S(i,j)=-S(i,j) !Accept Flip 
				else 
					if(P.gt.rand()) S(i,j)=-S(i,j) !Accept Flip
				endif
			enddo
		enddo
		
		!Calculate Magnetization Per Spin
		M=sum(S(1:L,1:L))/dble(L**2)
		
		!Heat Capacity:-----
		E=0.0d0 !Reset Energy
		!Spin Loop
		do i=1,L
			do j=1,L
				E=E-JJ*(S(i,j+1)+S(i,j-1)+S(i+1,j)+S(i-1,j))*S(i,j)
			enddo
		enddo
		E=E/2.0d0 !Because of Overcounting
		Eavg=Eavg+E/1000
		E2avg=E2avg+E**2/1000
		
		!Writing Results
		! write(10,*) mc,M
		
		!Calculate Thermal Average Magnetization
		if(imc.gt.900) Mavg=Mavg+M/100
	enddo
	
	!Calculate Heat Capacity
	C=(E2avg-Eavg**2)**2/(kB*T**2)
	
	!Writing Results
	! write(11,*) T,Mavg
	! write(12,*) B,Mavg
	write(13,*) T,C
	
	!Update T,B
	T=T+dT
	! B=B+dB
enddo

!Deallocation & Closing:========================
deallocate(S)
! close(10)
! close(11)
! close(12)
close(13)

end program Skeleton_2D

