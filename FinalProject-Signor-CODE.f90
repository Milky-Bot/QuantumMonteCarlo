!###### DOCUMENTATION ######
!
! #REQUIREMENTS
!   To execute this program no particular packages or lybraries are needed.
!  
! #EXECUTION
!   gfortran FinalProject-Signor-CODE.f90
! 
! #Work report available https://www.overleaf.com/read/kxkysdhbwnjr
!
! ####MODULE DIFFUSION MONTE CARLO####
!   Module containing subroutines and functions useful to execute a 
!   diffusion monte carlo simulation to compute the ground state of 
!   many-body quantum systems, based on 
!   "Introduction to the Diffusion Monte Carlo Method", Kosztin et al.
!   (https://arxiv.org/abs/physics/9702023).
!   It contains also some potentials used to test the code:
!   1d harmonic oscillator,1d Morse oscillator, Hydrogen atom, H- ion,
!   H2+ ion, H3+ ion, H2 molecule
!
!   #SUBROUTINE RANDOM_NORMAL
!       SCOPE
!           Generate random normal numbers from a uniform distribution
!           using gasdev. See http://gi.ics.nara-wu.ac.jp/~takasu/lecture/old/archives/L2_random_variables-2.pdf
!       ARGUMENTS
!           u: a double precision real number
!       OUTPUT
!           assigns to the input variable a gaussian random number
!           with mean 0 and variance 1
!
!   #FUNCTION POTENTIAL
!       SCOPE
!           Compute the potential energy (in dimensionless units) of some test systems
!       ARGUMENTS
!           position: DOUBLE PRECISION(:). Array of coordinates
!           potential choice: INTEGER. Number of the potential we want to study
!               1:1D Harmonic Oscillator
!               2: Morse potential
!               3: Hidrogen atom
!               4: H- ion
!               5: H2+ ion
!               6: H3+ ion
!               7: H2 molecule 
!       OUTPUT
!           POTENTIAL: REAL*8. potential in dimensionless units
!      
!   #SUBROUTINE BRANCH
!       SCOPE 
!           Perform the branching process required by the DMC method
!
!       ARGUMENTS
!           m: INTEGER. computed within the DMC SUBROUTINE. given by m=min{int[W(x) +u],3} is the number of
!                       copies of the replica to keep. m=0,1,2.  
!           branched_psi: DOUBLE PRECISION(:,:) first column: existence flag 
!                                               second column: spatial distribution of replicas
!           j: INTEGER. replica index
!            
!       OUTPUT
!           branched_psi with new born or killed replicas
!
!   #SUBROUTINE DMC
!       SCOPE
!           It is the where the DMC algorithm is performed.
!           1)initialize the matrix (psi) containing in the first column the
!             existence flag and in the others the spatial coordinates
!           2)loop over time steps to get tau->inf
!           3)compute the potential calling the FUNCTION POTENTIAL
!           4)compute the weight
!           5)call the branch subroutine for every replica
!           6)after tau iterations(stationarity reached) write on a file
!             the coordinates for every alive replica
!           7)update the reference energy
!           8) if d=1 computes directly the spatial distribution of replicas, i.e. the wavefunction
!              if d>1 computes the radial wavefunction

!       ARGUMENTS
!           NO:INTEGER. Initial number of alive replicas
!           Nmax: INTEGER. Maximum number of replicas; i.e. the length of
!                          vector psi
!           tau: INTEGER. Number of time steps
!           delta_tau: REAL*8. Size of the time steps
!           D: INTEGER. Dimension of the system (e.g. for the H atom D=3)
!           initial_position: REAL*8. initial position of every replica
!           E_R: REAL*8. initial reference energy. Suggested to set it as the value
!                        of the potential energy at the initial position
!           xmin: INTEGER. Lower limit for coordinates for the spatial distribution histogramming
!           xmax: INTEGER. upper limit for coordinates for the spatial distribution histogramming
!           nb: INTEGER. Number of bins for sorting replicas
!           potential_choice: INTEGER. Number of potential chosen. See above 
!       
!       OUTPUT
!           Different outputs are written on different files
!               'N.txt' number of alive replicas after every time iteration
!               'ref_energy.txt' Reference energy after every time iteration
!               -'wavefunction.txt' the wavefunction of the ground state if d=1,
!                                   if d>1 the radial wavefunction
!                   
!        
!
!   ####PROGRAM QUANTUM_MONTE_CARLO####
!   
!   MAIN to run the diffusion monte carlo.
!   Here the inputs for the different DMC subroutines and function are set.
!   The potential is chosen by terminal by the user, and according to this choice,
!   the D, initial_position and E_R are computed authomatically
!   Furthermore the potential choice is written on "potential_choice.txt' and
!   tau on 'time_steps.txt'.
!   All the resulting files are useful to easily compute the ground state, being:
!   E_0 the mean of the reference energy after tau time iterations, its error
!   the standard deviation;
!   the groundstate wavefunction the spatial distribution of alive replicas 
!   
!   
!
!###AUTHOR###
!  Theosamuele Signor, Unipd
!  theosamuele.signor@studenti.unipd.it

!#######################################################################


MODULE DIFFUSION_MC
    USE ieee_arithmetic
    IMPLICIT NONE
    CONTAINS 
    
    !Module containing some useful subroutines and functions
    !1-generate a gaussian random variable
    !2-sample quantum system potentials
    !subroutines for achieveng a difffusion monte carlo 
    
    !subroutine for generating random normal numbers from a uniform distribution
    !using gasdev method
    SUBROUTINE RANDOM_NORMAL(u)
        REAL*8 xx,yy,u
        CALL RANDOM_NUMBER(xx)
        CALL RANDOM_NUMBER(yy)
        u=sqrt(-2d0*log(xx))*cos(2*4*atan(1d0)*yy)
    END SUBROUTINE RANDOM_NORMAL
    
    !examples of potentials used to test the algorithm
    FUNCTION POTENTIAL(position,potential_choice)
        REAL*8 HARMONIC_OSCILLATOR_1D_POTENTIAL,MORSE_POTENTIAL,HYDROGEN_ATOM_POTENTIAL,H2_molecule,H2_ion,H3_ion, H_ion
        REAL*8 POTENTIAL
        integer potential_choice
        REAL*8 ,DIMENSION(:):: position
        REAL*8,DIMENSION(:), ALLOCATABLE::R,R1,R2,R3

        !1D Harmonic Oscillator ohmega=hbar=1
        IF (POTENTIAL_CHOICE .EQ. 1) THEN
            HARMONIC_OSCILLATOR_1D_POTENTIAL=0.5d0*position(1)**2
            POTENTIAL=HARMONIC_OSCILLATOR_1D_POTENTIAL

        !morse potential D=1/2,A=1
        ELSEIF  (POTENTIAL_CHOICE .EQ. 2) THEN
            MORSE_POTENTIAL= 0.5*(exp(-2*position(1))-2*exp(-position(1)))
            POTENTIAL=MORSE_POTENTIAL

        !Hydrogen atom 3d potential 1/r
        ELSEIF (POTENTIAL_CHOICE .EQ. 3) THEN
            HYDROGEN_ATOM_POTENTIAL=-1d0/norm2(position)
            POTENTIAL=HYDROGEN_ATOM_POTENTIAL
        
        !H- ion 
        ELSEIF (potential_choice .EQ. 4) THEN
            H_ion=-1./norm2(position(:3))-1./norm2(position(4:))
            H_ion=H_ion+1./norm2(position(4:)-position(:3))
            POTENTIAL=H_ion

        !H2+ ion 
        ELSEIF (POTENTIAL_CHOICE .EQ. 5) THEN
            ALLOCATE(R1(3),R2(3))
            R1=([0.,0.,1.]) !proton 1 position
            R2=([0.,0.,-1.])!proton 2 position

            H2_ion=-1./norm2(position-R1) - 1./norm2(position-R2) + 1./2.
            POTENTIAL=H2_ion
            DEALLOCATE(R1,R2)

        !H3+ ion
        !using Bohr approximation: the protons position is fixed
        !protons are equally distanced by 1.66 bohr radius 
        !(https://www.huy.dev/Anderson-JChemPhys-1975.pdf)
        ELSEIF (POTENTIAL_CHOICE .EQ. 6) THEN
            ALLOCATE(R1(3))
            ALLOCATE(R2(3))
            ALLOCATE(R3(3))

            !three protons equally distanced
            R1=([0.,0.,1.])
            R2=([0.,0.,2.66])
            R3=([0.,1.4376,1.83])

            !first electron protons interaction
            H3_ion=-1./norm2(position(:3)-(R1))-1./norm2(position(:3)-(R2))-1./norm2(position(:3)-(R3))
            !second electron protons interaction
            H3_ion=H3_ion-1./norm2(position(4:)-(R1))-1./norm2(position(4:)-(R2))-1./norm2(position(4:)-(R3))
            !electron-electron repulsion
            H3_ion=H3_ion+1./norm2(position(4:)-position(:3))
            !proton-proton repulsion
            H3_ion = H3_ion + 1./norm2(R1-R2) + 1./norm2(R1-R3) + 1./norm2(R2-R3)
            POTENTIAL=H3_ion
            DEALLOCATE(R1,R2,R3)
        
        !H2 molecule
        !using Bohr approximation: the protons position is fixed
        !using known proton-proton equilibrium distance
        ELSEIF (POTENTIAL_CHOICE .EQ. 7) THEN
            ALLOCATE(R2(3),R1(3))
            R2=([0.,0.,0.699]) !proton 1 position
            R1=([0.,0.,-0.699]) !proton 2 position

            H2_molecule=-1./norm2(position(:3)-R1) - 1./norm2(position(:3)-R2) !electron 1/protons interaction term
            H2_molecule= H2_molecule -1./norm2(position(4:)-R1) - 1./norm2(position(4:)-R2) !electron 2/protons interaction term
            H2_molecule=H2_molecule + 1./norm2(position(:3)-position(4:)) !electron/electron repulsion term
            H2_molecule=H2_molecule + 1./norm2(R2-R1) !proton/proton repulsion term

            POTENTIAL=H2_molecule
            DEALLOCATE(R2,R1)
        END IF
        RETURN
    END  FUNCTION POTENTIAL
    
!#######actual diffusion monte carlo subroutines########

    SUBROUTINE BRANCH(m,branched_psi,j)
        integer m,j,ix
        DOUBLE PRECISION branched_psi(:,:)
        !finding where we have not born/dead replicas
        integer, allocatable:: indeces(:)
        indeces=pack([(ix,ix=1,size(branched_psi(:,1)))],abs(branched_psi(:,1))<=1e-15)
        !kill replica
        IF (M .EQ. 0) THEN
             branched_psi(j,1)=0.
            
        !if m=1 do nothing
        !if m=2 create a copy of the particle in the
        !first position we have a dead one
        ELSEIF (M .EQ. 2) THEN
             branched_psi(indeces(1),1)=1.
             branched_psi(indeces(1),2:)=branched_psi(j,2:)

        !if m=3 create 2 copies
        ELSEIF (M .EQ. 3) THEN
            !first copy
            branched_psi(indeces(1),1)=1.
            branched_psi(indeces(1),2:)=branched_psi(j,2:)

            !second copy
            branched_psi(indeces(2),1)=1.
            branched_psi(indeces(2),2:)=branched_psi(j,2:)
        END IF
        DEALLOCATE(indeces)
    
    END SUBROUTINE BRANCH
    
    SUBROUTINE DMC(N0,Nmax,tau,delta_tau,D,initial_position,E_R,xmin,xmax,nb,potential_choice)
        INTEGER N0,Nmax,tau,potential_choice,D,m
        integer*4 w, nb
        REAL*8 delta_tau,u,E_R,v,vmean,ww,xmin,xmax,spatial_bin_step,Radius,Z
        INTEGER ii,jj,DD,N,BB
        DOUBLE PRECISION, DIMENSION(:,:),ALLOCATABLE::psi,temp_psi
        DOUBLE PRECISION, DIMENSION(:)::initial_position
        DOUBLE PRECISION, DIMENSION(:),ALLOCATABLE::BUCKET,bin_centers !for 1D histograms
        ALLOCATE(psi(Nmax,D+1))
        ALLOCATE(temp_psi(Nmax,D+1))
        

        !spatial bin step & bucket: will be used for the 1d spatial distribution
        ALLOCATE(BUCKET(nb)) 
        BUCKET(:)=0
        
        !
        IF (D .EQ. 1) THEN
            spatial_bin_step=(xmax-xmin)/(nb*1d0)
        ELSE 
            xmin=0d0
            spatial_bin_step=xmax/(nb*1d0)
        END IF

        allocate(bin_centers(nb))
        DO BB=1,nb
            bin_centers(BB)=(xmin+BB*spatial_bin_step+xmin+(BB-1)*spatial_bin_step)/2d0
        END DO

        !!!!STARTING DMC!!!!    
        psi(:,:)=0
        !Setting first N0 replicas alive
        psi(:N0,1)=1
 
        !initial position
        DO DD=1,D
            psi(:N0,DD+1)=initial_position(DD)    
        END DO

        !loop over time
        DO ii=1,2*tau
            !number of alive replicas
            N=0  
            vmean=0d0
            temp_psi=psi
            DO jj=1,Nmax
                IF (ABS(psi(jj,1)-1d0) <= 1E-15) THEN
                    N=N+1  !enlarge the number of alive replicas
                    !walk
                    DO DD=1,D
                        CALL RANDOM_NORMAL(u)
                        u=u*sqrt(delta_tau)
                        temp_psi(jj,DD+1)=temp_psi(jj,DD+1)+u
                    END DO 
                    
                    !weight calculation

                    V=POTENTIAL(temp_psi(jj,2:),potential_choice)
                    !Handling with infinites
                    IF (v .lt. -2134483647) THEN
                        v=-2134483646.
                    END IF

                    vmean=((N-1)*vmean+v)/(N*1d0) !method to compute iteratively the mean of a sequence
                    CALL RANDOM_NUMBER(u)
                    ww=(exp((E_R-v)*delta_tau)+u)
                    !handling with infinites
                    IF (ww .gt. 2134483647) THEN
                        ww=2134483647
                    END IF
                    w=int(ww)
                    m=min(w,3)
                    !branching
                    CALL BRANCH(m,temp_psi,jj)
                    
                    !after tau iterations, we start to evaluate the spatial replica distribution
                    !if we are evaluating ions or molecules, save also the positions of the replicas on another file
                    IF (ii .gt. tau) THEN    
                        IF (potential_choice .gt. 3) THEN
                            WRITE(88,*) temp_psi(jj,2:)
                        ENDIF 

                        IF (D .EQ. 1) THEN
                            DO BB=1,nb  !loop over buckets
                                IF ((temp_psi(jj,2) .ge. xmin+(BB-1)*spatial_bin_step) .and. &
                                (temp_psi(jj,2) .lt. xmin+(BB)*spatial_bin_step))  THEN
                                    Bucket(BB) = Bucket(BB) + 1
                                END IF
                            END DO
                        ELSE
                            Radius=norm2(temp_psi(jj,2:))
                            DO BB=1,nb  !loop over buckets
                                IF ((Radius .ge. xmin+(BB-1)*spatial_bin_step) .and. &
                                (Radius .lt. xmin+(BB)*spatial_bin_step))  THEN
                                    Bucket(BB) = Bucket(BB) + 1
                                END IF
                            END DO
                    
                        END IF !way to make histogram depending on dimension
                            
                    END IF !if time iteration > tau, we start to evaluate psi
                END IF !if replics is alive
            END DO !loop over replicas
            WRITE(44,*) E_R
            
            !reference energy update
            E_R = vmean+(1-(N*1d0/N0))
            psi=temp_psi
            WRITE (2,*) N
        END DO !loop over time steps
        
        IF (D .EQ. 1) THEN
            BUCKET=BUCKET/(sum(BUCKET)*spatial_bin_step)
            DO BB=1,nb
                write(1,*) bin_centers(BB),bucket(BB)
            END DO  
        ELSE 
            Z=sum(BUCKET/(bin_centers**2))*spatial_bin_step!normalization
            DO BB=1,nb
                write(1,*) bin_centers(BB),bucket(BB)*1d0/bin_centers(BB)**2/Z
            END DO
        END IF
        DEALLOCATE(psi,temp_psi,bin_centers,bucket)
    END SUBROUTINE
END MODULE
    
PROGRAM QUANTUM_MONTE_CARLO
    USE DIFFUSION_MC
    IMPLICIT NONE
    
    !inputs
    INTEGER N0,tau,nb,hbar,potential_choice,D,Nmax
    REAL*8 delta_tau,xmin,xmax,E_R
    REAL T1,T2,T
    DOUBLE PRECISION, DIMENSION(:),ALLOCATABLE::initial_position
    LOGICAL sound
    N0 = 1000   			!inital number of replicas
    Nmax = 2000			    !maximum number of replicas
    tau = 4000 			    !number of time steps
    delta_tau = 0.05d0	    !time step size
    xmin=-10.               !minimum value for histogram bins
    xmax=10.                !maximum value for histogram bins
    nb=200                  !number of bins in the histogram
    
    !potential choice 
    PRINT*,'chose potential: '
    PRINT*,'[1:1D Harmonic Oscillator'
    PRINT*,'2: Morse potential'
    PRINT*,'3: Hidrogen atom'
    PRINT*,'4: H- ion'
    PRINT*,'5: H2+ ion'
    PRINT*,'6: H3+ ion'
    PRINT*,'7: H2 molecule]'
    READ*,potential_choice
    sound=.true. !set to false if you don't want the sound when completed
    
    !setting the dimension and potential according to the user's choice
    !if the number of potential is not comprehended in the module, the simulation
    !uses the Harmonic oscillator potential
    IF (potential_choice .EQ. 1) THEN
        print*,'Using 1D Harmonic Oscillator'
        D=1
        ALLOCATE(initial_position(D))
        initial_position=0.
            
    ELSEIF (potential_choice .EQ. 2) THEN
        print*,'Using Morse potential'
        D=1
        ALLOCATE(initial_position(D))
        initial_position=0.

    ELSEIF (potential_choice .EQ. 3) THEN
        print*,'Using Hidrogen atom'
        D=3
        ALLOCATE(initial_position(D))
        initial_position=([0.,0.,1.])
    
    ELSEIF (potential_choice .EQ. 4) THEN
        print*,'Using H- ion'
        D=3*2
        ALLOCATE(initial_position(D))
        initial_position=([0.,0.,1.,0.,0.,-1.])

    ELSEIF (potential_choice .EQ. 5) THEN
            print*,'Using H2+ ion'
            D=3
            ALLOCATE(initial_position(D))
            initial_position=([0.,0.,0.])
    
    ELSEIF (potential_choice .EQ. 6) THEN
            print*,'Using H3+ ion'
            D=3*2
            ALLOCATE(initial_position(D))
            initial_position=([0.,0.,0.,0.,0.,0.5])
    
    ELSEIF (potential_choice .EQ. 7) THEN
        print*,'Using H2 molecule'
        D=3*2
        ALLOCATE(initial_position(D))
        initial_position=([0.,0.,1.,0.,0.,-1.])

    ELSE
        print*,'wrong choice, im using Harmonic oscillator'
        potential_choice=1
        D=1
        ALLOCATE(initial_position(D))
        initial_position=0.
    END IF

    !file to write the number of alive replicas after every
    !time iteration
    OPEN(unit=2, file='N.txt',ACCESS='append', ACTION='write',status='replace')
    !file to write the refernce energy    
    OPEN(unit=44, file='ref_energy.txt',ACCESS='append', ACTION='write',status='replace')
    
    !file store the ground state wavefunction 
    !first column is x, second psi(x)
    OPEN(unit=1, file='wavefunction.txt',ACCESS='append', ACTION='write',status='replace')

    IF (potential_choice .GT. 3) THEN
        OPEN(unit=88, file='positions.txt',ACCESS='append', ACTION='write',status='replace')
    END IF

    !writing the potential choice
    OPEN(unit=99, file='potential_choice.txt',ACCESS='append', ACTION='write',status='replace')
    WRITE(99,*) potential_choice

    !writing the parameter time_steps on file; useful for trace plots
    OPEN(unit=77, file='time_steps.txt',ACCESS='append', ACTION='write',status='replace')
    WRITE(77,*) tau
    
    !initial reference energy
    E_R=POTENTIAL(initial_position,potential_choice)
    
    !running the simulation
    CALL CPU_TIME(T1)
    CALL DMC(N0,Nmax,tau,delta_tau,D,initial_position,E_R,xmin,xmax,nb,potential_choice)
    CALL CPU_TIME(T2)
    T=T2-T1
    IF (T .GT. 180) THEN
        print*, 'CPU time: ',T/60.,' m'
    ELSE
        print*, 'CPU time: ',T,' s'
    END IF

    DEALLOCATE(initial_position)
    CLOSE(2)
    CLOSE(44)
    CLOSE(1)
    CLOSE(99)
    CLOSE(77)
    IF (potential_choice .gt. 3) THEN
        CLOSE(88)
    END IF

    IF (sound .EQV. .TRUE.) THEN
        call system( "echo simulation completed | festival --tts") 
    END IF
END PROGRAM