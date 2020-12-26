c#######################################################################
	subroutine Steel_2 (matpar,hstvP,hstv,epsP,sigP,deps,sig,E,ist)
c-----------------------------------------------------------------------
c MENEGOTTO-PINTO STEEL MODEL WITH FILIPPOU ISOTROPIC HARDENING
c             written by MOHD YASSIN (1993)
c           adapted to FEDEAS material library
c     by E. Spacone, G. Monti and F.C. Filippou (1994)
c-----------------------------------------------------------------------
c This steel model uses: 10 material parameters 
c                         8 history variables          
c-----------------------------------------------------------------------
	implicit none
c 
c arguments
	integer ist
	real*8  matpar(10),hstvP(8),hstv(8)
	real*8  epsP,sigP,deps
	real*8  sig,E
c-----------------------------------------------------------------------
c Note :
c   variables with P suffix is from last converged step
c   variables without suffix is from current iterative step
c
c I  matpar : STEEL FIXED PROPERTIES
c    Fy  = matpar(1)  : yield stress
c    E0  = matpar(2)  : initial stiffness
c    b   = matpar(3)  : hardening ratio (Esh/E0)
c    R0  = matpar(4)  : exp transition elastic-plastic
c    cR1 = matpar(5)  : coefficient for changing R0 to R
c    cR2 = matpar(6)  : coefficient for changing R0 to R
c    a1  = matpar(7)  : coefficient for isotropic hardening in compression
c    a2  = matpar(8)  : coefficient for isotropic hardening in compression
c    a3  = matpar(9)  : coefficient for isotropic hardening in tension
c    a4  = matpar(10) : coefficient for isotropic hardening in tension
c I  hstvP : STEEL HISTORY VARIABLES
c    epsminP = hstvP(1) : max eps in compression
c    epsmaxP = hstvP(2) : max eps in tension
c    epsplP  = hstvP(3) : plastic excursion
c    epss0P  = hstvP(4) : eps at asymptotes intersection
c    sigs0P  = hstvP(5) : sig at asymptotes intersection
c    epssrP  = hstvP(6) : eps at last inversion point
c    sigsrP  = hstvP(7) : sig at last inversion point
c    konP    = hstvP(8) : index for loading/unloading
c O  hstv : STEEL HISTORY VARIABLES   
c    same quantities as above without final P
c I  epsP  = strain at previous converged step
c I  sigP  = stress at previous converged step
c I  deps  = strain increment
c O  sig   = stress at current step
c O  E     = stiffness modulus (on return)
c              the type of stiffness depends on variable ist
c
c    ist   = switch that controls read/calculation
c            0: read material data and initialize
c            1: calculate tangent stiffness
c            2: calculate incremental secant stiffness
c            3: calculate total secant stiffness
c-----------------------------------------------------------------------
c material parameters
	real*8  Fy,E0,b,R0,cR1,cR2,a1,a2,a3,a4	  
c history variables
	integer kon
	real*8  epsmin,epsmax,epspl,epss0,sigs0,epsr,sigr  
c local variables
	real*8  eps		! eps = strain at current step
	real*8  shft,xi,R,epsrat,dum1,dum2,Esh,epsy
c-----------------------------------------------------------------------
	if (ist.eq.0) then
c  read material data and identify no of history variables
c
	   hstvP (1) =10	! no of material parameters for material model
	   hstv  (1) = 8	! no of history variables for material model
c
	else
c  material state determination     
c  retrieve concrete fixed material properties
c 
	   Fy   = matpar(1)	! always positive
	   E0   = matpar(2)
	   b    = matpar(3)
	   R0   = matpar(4)
	   cR1  = matpar(5)
	   cR2  = matpar(6)
	   a1   = matpar(7)
	   a2   = matpar(8)   
	   a3   = matpar(9)
	   a4   = matpar(10)  
c 
c	calculate other fixed material properties
c 
	   Esh  = b*E0
	   epsy = Fy/E0		! always positive
c     
c	retrieve steel history variables
c 
	   epsmax = hstvP(1)
	   epsmin = hstvP(2)
	   epss0  = hstvP(3)
	   sigs0  = hstvP(4)
	   epsr   = hstvP(5)
	   sigr   = hstvP(6)   
	   epspl  = hstvP(7)
	   kon    = idint(hstvP(8))
c                       
c	calculate current strain
c 
	   eps    = epsP + deps
c 
c	check for virgin steel: if kon==0 and the strain increment deps is 0,
c	set the tangent modulus, update the history variables and return
c 
	   if (kon.eq.0) then
	      if (deps.eq.0.d0) then
		 E   = E0
		 sig = 0.d0
c
		 hstv(1) = epsmax
		 hstv(2) = epsmin
		 hstv(3) = epss0
		 hstv(4) = sigs0
		 hstv(5) = epsr
		 hstv(6) = sigr  
		 hstv(7) = epspl
		 hstv(8) = kon
		 return
	      else
		 epsmax =  epsy
		 epsmin = -epsy
		 if (deps.lt.0.d0) then   
		    kon   = 2
		    epss0 = epsmin
		    sigs0 = -Fy
		    epspl = epsmin
		 else
		    kon   = 1
		    epss0 = epsmax
		    sigs0 = Fy
		    epspl = epsmax
		 endif
	      endif
	   endif
c
c	in case of load reversal from negative to positive strain increment,
c	update the minimum previous strain, store the last load reversal
c	point and calculate the stress and strain (sigs0 and epss0) at the
c	new intersection between elastic and strain hardening asymptote
c	To include isotropic strain hardening shift the strain hardening
c	asymptote by sigsft before calculating the intersection point
c	Constants a3 and a4 control this stress shift on the tension side 
c 
	   if (kon.eq.2.and.deps.gt.0.d0) then
	      kon    = 1
	      epsr   = epsP
	      sigr   = sigP
	      epsmin = dmin1 (epsP,epsmin)
	      shft   = 1.d0 + a3*(( (epsmax-epsmin)/(2.d0*(a4*epsy)) )**0.8)
	      epss0  = (Fy*shft - Esh*epsy*shft - sigr + E0*epsr)/(E0-Esh)
	      sigs0  =  Fy*shft + Esh*(epss0-epsy*shft)
	      epspl  = epsmax
	   else
c 
c	in case of load reversal from positive to negative strain increment
c	update the maximum previous strain, store the last load reversal
c	point and calculate the stress and strain (sigs0 and epss0) at the
c	new intersection between elastic and strain hardening asymptote
c	To include isotropic strain hardening shift the strain hardening
c	asymptote by sigsft before calculating the intersection point
c	Constants a1 and a2 control this stress shift on compression side 
c 
	      if (kon.eq.1.and.deps.lt.0.d0) then
		 kon    = 2
		 epsr   = epsP
		 sigr   = sigP
		 epsmax = dmax1(epsP,epsmax)
		 shft   = 1.d0 + a1*(( (epsmax-epsmin)/(2.d0*(a2*epsy)) )**0.8)
		 epss0  = (-Fy*shft + Esh*epsy*shft - sigr + E0*epsr)/(E0-Esh)
		 sigs0  =  -Fy*shft + Esh*(epss0+epsy*shft)
		 epspl  = epsmin
	      endif
	   endif   
c 
c	calculate current stress sig and tangent modulus E
c 
	   xi     = dabs((epspl-epss0)/epsy)
	   R      = R0*(1.d0 - (cR1*xi)/(cR2+xi))
	   epsrat = (eps-epsr)/(epss0-epsr)
c
	   dum1  = 1.d0 + dabs(epsrat)**R
	   dum2  = dum1**(1/R)      
	   sig   = b*epsrat +(1.d0-b)*epsrat/dum2
	   sig   = sig*(sigs0-sigr)+sigr
c      
	   E = b + (1.d0-b)/(dum1*dum2)									 
	   E = E*(sigs0-sigr)/(epss0-epsr)
c 
c	calculate appropriate stiffness depending on value of ist
c 
	   if (ist.eq.2.and.(deps.ne.0.d0)) then
	      E = (sig-sigP)/deps
	   else if (ist.eq.3.and.((eps-epsr).ne.0.d0)) then
	      E = (sig-sigr)/(eps-epsr)
	   else
c		add additional cases, if needed
	   endif     
c 
c	transfer all variables to history vector and return
c 
	   hstv(1) = epsmax
	   hstv(2) = epsmin
	   hstv(3) = epss0
	   hstv(4) = sigs0
	   hstv(5) = epsr
	   hstv(6) = sigr  
	   hstv(7) = epspl
	   hstv(8) = kon
c
	endif
	return
c#######################################################################
	end subroutine Steel_2
