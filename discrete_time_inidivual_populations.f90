PROGRAM boombust ! DISCRETE GENERATION COMPETITION MODEL WITH BOOM-BUST POPULATION DYNAMICS
  	IMPLICIT NONE ! general model with arbitrary coefficients
  	DOUBLE PRECISION, ALLOCATABLE :: X(:),X0(:),CORR(:),PAV(:)
  	DOUBLE PRECISION, ALLOCATABLE :: CC(:,:),DC(:)
  	DOUBLE PRECISION, ALLOCATABLE :: B1(:,:),DSEED(:)
  	DOUBLE PRECISION, ALLOCATABLE :: PC(:),PA(:),PA_NEW(:)
  	INTEGER N,K,L,IW,K1,K2,K3,MC,K4,IC,KC,NH,plotcount,nctot,iw0,NM,km
  	INTEGER T, DT, TFIN, KC1,KC2,JE,NE,J1, J2, JL, MC_max, NCORR
  	INTEGER NGAUSS,MCIN,NR,DTMERGE,TMERGE,TPURGE, DTPURGE, TWrite, DTWrite
  	DOUBLE PRECISION ran3,pt,pm, extinct,carr, comp, lambda,b,RCORR,DST
  	DOUBLE PRECISION SA,DTP,TP,POPMIN,CB_mut,CB_mig,POPBIRTH, mut_dist,AA
  	DOUBLE PRECISION DISP,SIG,C1,C2,DIST,dw,lw,d1,mig_density, dy,w, extinct_rad
  	DOUBLE PRECISION GS,XP,GW,mu,TSEED,DTSEED,DCC,DCMIN,DCFIN,DTmigrate, pop_size
  	DOUBLE PRECISION LSEED, CBIN
  	LOGICAL GAUSS 
  	LOGICAL, ALLOCATABLE :: LIVE(:) 


! ecological parameters 
 	lambda=1.2!1.05!1.2 !intrinsic growth rate
  	b = 45.d0! 120.d0! 9.D0  ! Bellows exponent: complexity of discrete map is 1-b(lambda-1)/lambda
  
  	N=2    !SPACE DIMENSIONALITY
  	SA=.5D0   ! WIDTH OF THE GAUSSIAN PART OF COMP. KERN.
  	GAUSS=.true. ! IF TRUE, THE G. PART WITH SA IS IN COMP. KERNEL
  	NGAUSS=N ! NUMBER OF DIMENS. (STARTING FROM 1) WITH GAUSSIAN TERM
        NCORR=100  ! NUMBER OF POINTS IN CORR. HISTOGRAM
        RCORR=2.D0 ! MAX DISTANCE IN CORR. HISTOGRAM
  	POPMIN=1.D-12! THRESHOLD FOR THE POPULATION to go extinct 
  	
! mutation paramneters
	DTSEED=0.005D0 ! probability of seeding a new cluster by mutation per generatihospitalidadon
  	CB_mut=0.01 ! Standard deviation of the Gaussian that is used ot choose the actual distance 
  				  ! BETWEEN THE OFFSPRING AND THE ANSESTOR
  	POPBIRTH=0.01 ! RARE MUTANT POPULATION
  
! merge parameters
 	DTMERGE=1000 ! time interval (in generations) between merging procedures 
	DCMIN=1*CB_mut ! MINIMAL DISTANCE BELOW WHICH THE CLUSTERS ARE MERGED

! purge parameters
  	DTPURGE = 1000 ! time interval (in generations) between purging procedures

! write parameters
  	DTwrite = 1000	! time interval (in generations) between writing to file procedures

! intialization parameters
  	GW=1.D0   ! WIDTH OF A GAUSSIAN DISTRIBUTION FOR COEFFICIENTS
 	CBIN=10*CB_mut ! TYPICAL DISTANCE BETWEEN CLUSTERS ON INITIAL SEEDING
  	MCIN=1 !1   ! INITIAL NUMBER OF CLUSTERS
 
 
  	MC=1000   ! maximal number of clusters
  	MC_max= MCIN ! index of "highest" live cluster  
  	TFIN=20000000 ! length of simulation


  ALLOCATE(X(N*MC),X0(N))
  ALLOCATE(CC(N,MC),LIVE(MC),DSEED(N))
  ALLOCATE(B1(N,N),CORR(NCORR))
  ALLOCATE(PC(MC),PA(MC),pa_new(MC),PAV(MC))
 

  	open(unit=12,file='rand.dat')
  	read(12,*)iw
!	iw=-7434563
  	close(12)
  
  	do k1=1,n
     	do k2=1,n
43     		XP=(RAN3(IW)*2-1.D0)*4*GW ! GAUSSIAN ASSIGNMENT
       		mu=EXP(-XP*XP/(GW*GW*2))
       		IF(mu.LE.RAN3(IW)) GO TO 43
       		B1(K1,K2)=XP/dsqrt(1.d0*N) ! THIS IS INSTEAD OF READING IT FROM coeff*.dat FILE
       		B1(K1,K2)=0
       	end do
  	end do

  	X0=0.D0
  	CC=0.D0
  	LIVE=.FALSE.
  	PA=0.D0

  	DO KC=1,MCIN ! SEEDING CLUSTERS
     	LIVE(KC)=.TRUE.
     	PA(KC)=0.5D0
     	DO K1=1,N
   			CC(K1,KC)=X0(K1)+(2*ran3(iw)-1.d0)*CBIN ! RANDOM CLUSTER PLACING
    	END DO
  	END DO
  	T=1
  	Twrite=DTWrite
  	TMERGE=DTMERGE
  	TPURGE=DTPURGE
  	JL=MCIN ! NUMBER OF LIVE CLUSTERS
  	mc_max=jl ! largest live index
  	OPEN(UNIT=15, FILE='movie.dat') ! writing a movie

         plotcount=0
         PAV=0.D0
  
	DO WHILE (T.LE.TFIN) 
	
                PAV=PAV+PA
		T=T+1
! population dynamics
		pa_new=0.D0
     	DO KC1=1,MC_max  ! COMPUTING THE TABLE OF DEATH RATES 
        	IF(LIVE(KC1)) THEN
        		comp=0.D0
        		carr=0.D0
        		do j1=1,n 
        			carr=carr-(cc(J1,kc1)**4)/4
        		end do
        		carr=EXP(carr)
   				DO KC2=1,MC_max ! affected by the cluster kc2
        			IF(LIVE(KC2)) THEN
           				W=0.D0
        				do j1=1,n ! cluster c1 is affected by c2
      						dy=cc(j1,kc2)-cc(j1,kc1)
     						do j2=1,n
        						W=W-b1(j1,j2)*(cc(j1,kc1)-cc(j1,kc2))*cc(j2,kc2)
        					end do
     					     if(GAUSS.and.j1.le.ngauss) W=W-DY*DY/(2*sa*sa) ! ADDING GAUSSIAN PART TO A PARTICULAR. D
    					end do 
    					comp=comp+EXP(W)*pa(kc2)
    				end if
    			end do
  				pa_new(kc1)=pa(kc1)*lambda/(1+(lambda-1)*(comp/carr)**b)
  			end if
  		end do
  		pa=pa_new

! counting live clusters 
 		jl=0
 		Do kc=1,MC_max
 			IF (live(kc)) then
 				if (pa(kc).lt.popmin) then
 					pa(kc)=0
 					live(kc)=.False.
 				else
 					jl=jl+1
 				end if
 			end if
 		end do
 		
 ! SEEDING A NEW CLUSTER by mutation:
		IF(JL.LT.MC.AND.RAN3(IW).LT.DTSEED) THEN 
			pop_size=0		! choosing the parent
			do kc1=1,MC_max
				if (live(kc1)) then
					pop_size=pop_size+pa(kc1)
				end if
			end do
			mu=ran3(iw)
			kc1=1
			pm=0
			if (live(kc1)) then 
				pm=pm+pa(kc1)/pop_size
			end if
			do while (pm.lt.mu)
				kc1=kc1+1
				if (live(kc1)) then 
					pm=pm+pa(kc1)/pop_size
				end if
			end do	! now kc1 is the parent, chosen proportionally to population size
        	KC2=1
        	DO WHILE (LIVE(KC2)) ! PLACE FOR OFFSPRING
           		KC2=KC2+1
        	END DO
        	LSEED=0.D0
        	DO J1=1,N ! MAKING THE DISTANCE EXACTLY CB_mut
           		DSEED(J1)=2*RAN3(IW)-1
           		LSEED=LSEED+DSEED(J1)*DSEED(J1)
        	END DO
26     		mut_dist=(RAN3(IW)*2-1.D0)*4*CB_mut ! GAUSSIAN ASSIGNMENT of mutation size
       		mu=EXP(-mut_dist*mut_dist/(cb_mut*cb_mut*2))
       		IF(mu.LE.RAN3(IW)) GO TO 26
        	LSEED=DSQRT(LSEED)/mut_dist
        	DO J1=1,N
           		CC(J1,KC2)=CC(J1,KC1)+DSEED(J1)/LSEED
	    	END DO
      		PA(KC2)=POPBIRTH 
      		LIVE(KC2)=.TRUE.
      		Jl=Jl+1;
      		If (KC2.gt.mc_max) mc_max=kc2
 		END IF !end mutation seeding 

   ! MERGING CLOSE CLUSTERS:
		IF(T.GE.TMERGE) THEN ! MERGING CLOSE CLUSTERS
    		TMERGE=T+DTMERGE
     		DO KC=1,MC_max  ! REMOVING CLUSTERS WHICH ARE TOO CLOSE
				IF(live(kc)) THEN ! removing coinciding clusters
              		DO KC2=1,KC-1
               			IF(LIVE(KC2)) THEN
							DCC=0.D0
                			DO J1=1,N
                    			DCC=DCC+(CC(J1,KC)-CC(J1,KC2))**2
							END DO
							DCC=DSQRT(DCC)
                  			IF (DCC.LE.DCMIN) THEN !LIQUIDATING CLUSTER KC
        						DO J1=1,N ! PRESERVING CENTRE OF MASS POSITION
           			 CC(J1,KC2)=(CC(J1,KC)*PA(KC)+CC(J1,KC2)*PA(KC2))/(PA(KC)+PA(KC2)) 
        						END DO
                    			LIVE(KC)=.FALSE.
                    			PA(KC2)=PA(KC2)+PA(KC)
                    			PA(KC)=0.D0
                    			JL=JL-1
                  			END IF
                  		END IF
               		END DO
            	END IF
			END DO
		END IF !END MERGE

! purging of cluster array:
		IF(T.GE.TPURGE) THEN ! "purging" the array (moving all non-zero clusters to the start) 
    		TPURGE=T+DTPURGE
			Do KC=1,MC
				IF(.NOT.LIVE(KC)) THEN
					kc1=kc
					do while (.NOT.LIVE(KC1).and.kc1.lt.mc)
						kc1=kc1+1
					end do
					if (live(kc1)) then 
						do j1=1,n
							cc(j1,kc)=cc(j1,kc1)
						end do
						pa(kc)=pa(kc1)
						pa(kc1)=0.0
						live(kc1)=.False.
						live(kc)=.TRUE.	
					end if
				end if
			END DO
			jl=0
			mc_max=0
			Do kc=1,MC
 				IF (live(kc)) then
 					mc_max=mc_max+1
 					if (pa(kc).lt.popmin) then
 						pa(kc)=0
 						live(kc)=.False.
 					else
 						jl=jl+1
 					end if
 				end if
 			end do
 		end if

		IF(T.GE.TWrite) THEN ! writing to movie file 
    		Twrite=T+DTwrite
                 PAV=PAV/DTwrite
                WRITE(15,*) (plotcount*100,K1=1,4)
     		plotcount=plotcount+1
			DO KC=1,MC
        		IF(LIVE(KC)) THEN
            		WRITE(15,*)(CC(J1,KC),J1=1,2),PA(KC),T
        		END IF
                        END DO
                        DO KC1=1,MC
                           IF(LIVE(KC1)) THEN
                              AA=PA(KC1)-PAV(KC1)
                              DO KC2=1,MC
                                 IF(LIVE(KC2)) THEN
                                    DST=0.D0
                                    DO J1=1,N
                                       DST=DST + (CC(J1,KC1)-CC(J1,KC2))**2
                                    END DO
                                    DST=DSQRT(DST)
                                    J2=DST/RCORR*NCORR+1
                    IF(J2.LE.NCORR) CORR(J2)=CORR(J2)+AA*(PA(KC2)-PAV(KC2))
                                 END IF
                               END DO
                                 END IF
                               END DO
                    PAV=0.d0
                  END IF

	END DO ! end of generation loop

	
		WRITE(*,*) 'FINAL TIME'
  		write(*,*) MC_max
             OPEN(UNIT=21,FILE='corr3.dat')
             DO KC1=1,NCORR
                WRITE(21,*)KC1*RCORR/NCORR,CORR(KC1)
             END DO   
             CLOSE(21)
    
    
  		CLOSE(15)


  open(unit=12, file='rand.dat')
  write(12,*)nint(-ran3(iw)*1.d+6)
  close(12)   


END PROGRAM boombust


!---------------------------------------------------------------------
      function ran3(idum)
      integer idum
      integer mbig,seed,mz
      double precision ran3,fac
      parameter (mbig=1000000000,mseed=161803398,mz=0,fac=1./mbig)
       integer i,iff,ii,inext,inextp,k
      integer mj,mk,ma(55)
      save iff,inext,inextp,ma
      data iff /0/
1      if(idum.lt.0.or.iff.eq.0)then
        iff=1
        mj=mseed-iabs(idum)
        mj=mod(mj,mbig)
        ma(55)=mj
        mk=1
        do 11 i=1,54
          ii=mod(21*i,55)
          ma(ii)=mk
          mk=mj-mk
          if(mk.lt.mz)mk=mk+mbig
          mj=ma(ii)
11      continue
        do 13 k=1,4
          do 12 i=1,55
            ma(i)=ma(i)-ma(1+mod(i+30,55))
            if(ma(i).lt.mz)ma(i)=ma(i)+mbig
12        continue
13      continue
        inext=0
        inextp=31
        idum=1
      endif
      inext=inext+1
      if(inext.eq.56)inext=1
      inextp=inextp+1
      if(inextp.eq.56)inextp=1
      mj=ma(inext)-ma(inextp)
      if(mj.lt.mz)mj=mj+mbig
      ma(inext)=mj
      ran3=mj*fac
      if(ran3.le.0.or.ran3.ge.1) goto 1
      return
      end function ran3

         
