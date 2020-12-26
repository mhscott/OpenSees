c#######################################################################
	subroutine Steel_1 (matpar,hstvP,hstv,epsP,sigP,deps,sig,E,ist)
c-----------------------------------------------------------------------
c       BILINEAR STEEL MODEL WITH ISOTROPIC HARDENING
c                by F.C. Filippou (1994)
c-----------------------------------------------------------------------
c This steel model uses: 7 material parameters  
c                        7 history variables          
c-----------------------------------------------------------------------
	implicit none
c 
c arguments
	integer ist
	real*8  matpar(7),hstvP(7),hstv(7)
	real*8  epsP,sigP,deps
	real*8  sig,E
c-----------------------------------------------------------------------
c Note :
c  variables with P suffix is from last converged step
c  variables without suffix is from current iterative step
c
c I  matpar : STEEL FIXED PROPERTIES
c    Fy     = matpar(1)  : yield stress
c    E0     = matpar(2)  : initial stiffness
c    b      = matpar(3)  : hardening ratio (Esh/E0)
c    a1     = matpar(4)  : coefficient for isotropic hardening
c    a2     = matpar(5)  : coefficient for isotropic hardening
c    a3     = matpar(6)  : coefficient for isotropic hardening
c    a4     = matpar(7)  : coefficient for isotropic hardening
c I  hstvP : STEEL HISTORY VARIABLES
c    epsminP = hstvP(1) : max eps in compression
c    epsmaxP = hstvP(2) : max eps in tension
c    kon     = hstvP(3) : index for loading/unloading
c                kon = 1 : positive strain increment
c                kon = 2 : negative strain increment
c    epsr    = hstvP(4) : strain at last reversal
c    sigr    = hstvP(5) : stress at last reversal
c    
c O  hstv : STEEL HISTORY VARIABLES   
c        same quantities as above without final P
c
c I  epsP  = strain at previous converged step
c I  sigP  = stress at previous converged step
c I  deps  = strain increment
c O  sig   = stress at current step
c O  E     = stiffness modulus (on return)
c            the type of stiffness depends on variable ist
c
c    ist   = switch that controls read/calculation
c            0: read material data and initialize
c            1: calculate tangent stiffness
c            2: calculate incremental secant stiffness
c            3: calculate total secant stiffness
c-----------------------------------------------------------------------
c material parameters
	real*8  E0,Fy,b,a1,a2,a3,a4,Esh,epsy  
c history variables
	integer kon
	real*8  epsmin,epsmax,epsr,sigr,shft_p,shft_n  
c local variables
	real*8  eps,c1,c2,c3
c-----------------------------------------------------------------------  
	if (ist.eq.0) then
c
c  read material data and identify no of history variables
c
	   hstvP (1) = 7	! no of material parameters for material model
	   hstv  (1) = 7	! no of history variables for material model
c
	else
c
c  material state determination     
c  retrieve concrete fixed material properties
c 
	   Fy  = matpar(1)
	   E0  = matpar(2)
	   b   = matpar(3)
	   a1  = matpar(4)
	   a2  = matpar(5)
	   a3  = matpar(6)
	   a4  = matpar(7)  
c 
c	calculate other fixed material properties
c 
	   Esh  = b*E0
	   epsy = Fy/E0
c       
c	retrieve steel history variables
c 
	   epsmax = hstvP(1)
	   epsmin = hstvP(2)
	   kon    = idint(hstvP(3))
	   epsr   = hstvP(4)
	   sigr   = hstvP(5)
	   shft_p = hstvP(6)
	   shft_n = hstvP(7)
c       
c	calculate current strain
c 
	   eps   = epsP + deps
c
c       Old -> New
c       ----------
c	Fy -> E0
c	E0 -> E1
c	b -> E2
c	a1 -> eps1
c	a2 -> eps2

c	c1 -> sig1
c	c2 -> sig2
c	c3 -> epstmp

	   c1 = E0*a1
	   c2 = c1 + E0*(a2-a1)
	   c3 = dabs(eps)

c	write (*,*) c3

	   if (c3.lt.a1) then
	      E = Fy
	      sig = E*c3
	   else if (c3.lt.a2) then
	      E = E0
	      sig = c1 + E*(c3-a1)
	   else
	      E = b
	      sig = c2 + E*(c3-a2)
	   endif
	   
	   if (eps.lt.0.0) then
	      sig = -sig
	   endif
c 
c	transfer all variables to history vector and return
c 
	   hstv(1) = epsmax 
	   hstv(2) = epsmin
	   hstv(3) = kon
	   hstv(4) = epsr
	   hstv(5) = sigr
	   hstv(6) = shft_p
	   hstv(7) = shft_n
c
	endif
c       
	return
c#######################################################################
	end subroutine Steel_1
