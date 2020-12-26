!#######################################################################
	subroutine Concrete_2 (matpar,hstvP,hstv,epscP,sigcP,depsc,sigc,Ect,ist)
!-----------------------------------------------------------------------
! concrete model with damage modulus    
!       by MOHD YASSIN (1993)
! adapted to FEDEAS material library
! by D. Sze and Filip C. Filippou in 1994
!
! This concrete model requires: 7 fixed material parameters
!                               2 history variables
!-----------------------------------------------------------------------
! fc    = concrete compression strength           : matpar(1)
! epsc0 = strain at compression strength          : matpar(2)
! fcu   = stress at ultimate (crushing) strain	  : matpar(3)
! epscu = ultimate (crushing) strain              : matpar(4)       
! rat   = ratio between unloading slope at epscu  : matpar(5)
!             and original slope
! ft    = concrete tensile strength               : matpar(6)
! Ets   = tension stiffening slope                : matpar(7)
!                                                               
! Ec0   = initial tangent modulus                 :fixed value
!                    
! ecmin = smallest previous concrete strain
! sigmm = corresponding concrete stress
!
! dept  = shift in tensile strain due to strain softening
!
! epscP      = previous strain
! sigcP,sigc = stress (previous at input, current at output)
! depsc      = strain increment  (input)
!     Ect    = stiffness modulus (on return)
!
!    ist   = switch that controls read/calculation
!            0: read material data and initialize
!            1: calculate tangent stiffness
!            2: calculate incremental secant stiffness
!            3: calculate total secant stiffness
!-----------------------------------------------------------------------  
 !use Var_Precision
	implicit none
!
	integer ist
	real*8  matpar(7),hstvP(2),hstv(2) 
	real*8  epscP,sigcP,depsc,epsc,sigc,Ect
	real*8  fc,epsc0,fcu,epscu,rat,ft,Ets,Ec0
	real*8  ecmin,dept
	real*8  epsR,sigmR,sigmm,Er,ept,epn,sicn,sigmin,sigmax,epstmp
	real*8  dumy
!-----------------------------------------------------------------------
	if (ist==0) then
!
!  read material data and identify no of history variables
!
	   hstvP (1) = 7	! no of material parameters for material model
	   hstv  (1) = 2	! no of history variables for material model
!
	else
!
!  material state determination     
!  retrieve concrete fixed material properties
! 
	   fc    = matpar(1)
	   epsc0 = matpar(2)
	   fcu   = matpar(3)
	   epscu = matpar(4)
	   rat   = matpar(5)               
	   ft    = matpar(6)
	   Ets   = matpar(7)
! 
!	Calculate the other fixed material properties
!       
	   Ec0  = 2.0d0*fc/epsc0
!       
!	retrieve concrete history variables
! 
	   ecmin = hstvP(1)
	   dept  = hstvP(2)
!                    
!	calculate current strain
! 
	   epsc   = epscP+depsc          
! 
!	if the current strain is less than the smallest previous strain
!	call the monotonic envelope in compression and reset minimum strain
!    
	   if (epsc<ecmin) then
	      call Compr_Envlp (epsc,sigc,Ect,matpar)
	      ecmin = epsc
!      
	   else
! 
!		else, if the current strain is between the minimum strain and ept 
!		(which corresponds to zero stress) the material is in the unloading-
!		reloading branch and the stress remains between sigmin and sigmax
! 
!		calculate strain-stress coordinates of point R that determines
!		the reloading slope according to Fig.2.11 in EERC Report
!		(corresponding equations are 2.31 and 2.32
!		the strain of point R is epsR and the stress is sigmR
!       
	      epsR  = (fcu-rat*Ec0*epscu)/(Ec0*(1.d0-rat))
	      sigmR = Ec0*epsR
! 
!		calculate the previous minimum stress sigmm from the minimum
!		previous strain ecmin and the monotonic envelope in compression
!       
	      call Compr_Envlp (ecmin,sigmm,dumy,matpar)
! 
!		calculate current reloading slope Er (Eq. 2.35 in EERC Report)
!		calculate the intersection of the current reloading slope Er
!		with the zero stress axis (variable ept) (Eq. 2.36 in EERC Report)
!               
	      Er  = (sigmm-sigmR)/(ecmin-epsR)
	      ept = ecmin-sigmm/Er
!       
	      if (epsc<=ept) then
		 sigmin = sigmm+Er*(epsc-ecmin)
		 sigmax = 0.5*Er*(epsc-ept)
		 sigc   = sigcP+Ec0*depsc
		 Ect    = Ec0
		 if (sigc<=sigmin) then
		    sigc = sigmin
		    Ect  = Er
		 endif
		 if (sigc>=sigmax) then
		    sigc = sigmax
		    Ect  = 0.5*Er
		 endif
	      else
! 
!			else, if the current strain is between ept and epn
!			(which corresponds to maximum remaining tensile strength) 
!			the response corresponds to the reloading branch in tension
!			Since it is not saved, calculate the maximum remaining tensile 
!			strength sicn (Eq. 2.43 in EERC Report)
! 
!			calculate first the strain at the peak of the tensile stress-strain
!			relation epn (Eq. 2.42 in EERC Report)
!      
		 epn = ept+dept                                            
! 
		 if (epsc<=epn) then   
		    call Tens_Envlp (dept,sicn,Ect,matpar)
		    if (dept/=0.d0) then
		       Ect  = sicn/dept
		    else
		       Ect  = Ec0
		    endif
		    sigc = Ect*(epsc-ept)
		 else
!  
!				else, if the current strain is larger than epn the response
!				corresponds to the tensile envelope curve shifted by ept
! 
		    epstmp = epsc-ept
		    call Tens_Envlp (epstmp,sigc,Ect,matpar)
		    dept   = epsc-ept
		 endif
	      endif      
	   endif
! 
!	calculate appropriate stiffness depending on value of ist
! 
	   if (ist==2.and.depsc/=0.d0) then
	      Ect = (sigc-sigcP)/depsc
	   else if (ist==3.and.epsc/=0.d0)then
	      Ect = sigc/(epsc-ept)
	   else
!		add additional cases, if needed
	   endif         
!      
	   hstv(1) = ecmin
	   hstv(2) = dept
	endif
	return
	end subroutine Concrete_2
!#######################################################################
	subroutine  Tens_Envlp (epsc,sigc,Ect,matpar)
!-----------------------------------------------------------------------
! monotonic envelope of concrete in tension (positive envelope)
!
!   ft    = concrete tensile strength
!   Ec0   = initial tangent modulus of concrete 
!   Ets   = tension softening modulus
!   eps   = strain
!
!   returned variables
!    sig  = stress corresponding to eps
!    Ect  = tangent concrete modulus
!-----------------------------------------------------------------------
	implicit none

	real*8 epsc
	real*8 sigc,Ect
	real*8 matpar(7)
	real*8 fc,epsc0,fcu,epscu,rat,ft,Ets,Ec0
      
! local variables
	real*8 eps0,epsu

	fc    = matpar(1)
	epsc0 = matpar(2)
	fcu   = matpar(3)
	epscu = matpar(4)
	rat   = matpar(5)               
	ft    = matpar(6)
	Ets   = matpar(7)      

	Ec0  = 2.0d0*fc/epsc0
	
	eps0 = ft/Ec0
	epsu = ft*(1.d0/Ets+1.d0/Ec0)
	if (epsc<=eps0) then
	   sigc = epsc*Ec0
	   Ect  = Ec0 
	else
	   if (epsc<=epsu) then 
	      Ect  = -Ets
	      sigc = ft-Ets*(epsc-eps0)
	   else
!       Ect  = 0.d0
	      Ect  = 1.d-10
	      sigc = 0.d0
	   endif
	endif
!
	return
	end subroutine Tens_Envlp
!#######################################################################
	subroutine Compr_Envlp (epsc,sigc,Ect,matpar)
!-----------------------------------------------------------------------
! monotonic envelope of concrete in compression (negative envelope)
!
!   fc    = concrete compressive strength
!   epsc0 = strain at concrete compressive strength
!   fcu   = stress at ultimate (crushing) strain 
!   epscu = ultimate (crushing) strain
!   Ec0   = initial concrete tangent modulus
!   epsc  = strain
!
!   returned variables
!   sigc  = current stress
!   Ect   = tangent concrete modulus
!-----------------------------------------------------------------------
	implicit none

	real*8 epsc,rat
	real*8 sigc,Ect
	real*8 matpar(7)
	real*8 fc,epsc0,fcu,epscu,ratio,ft,Ets,Ec0
! 
! parabolic ascending branch until strain epsc0
!
	fc    = matpar(1)
	epsc0 = matpar(2)
	fcu   = matpar(3)
	epscu = matpar(4)
	ratio = matpar(5) 
	ft    = matpar(6)
	Ets   = matpar(7)      

	Ec0  = 2.0d0*fc/epsc0

	rat = epsc/epsc0    
	if (epsc>=epsc0) then
	   sigc = fc*rat*(2.d0-rat)
	   Ect  = Ec0*(1.d0-rat)
	else
!         
!   linear descending branch between epsc0 and epscu
! 
	   if (epsc>epscu) then   
	      sigc = (fcu-fc)*(epsc-epsc0)/(epscu-epsc0)+fc
	      Ect  = (fcu-fc)/(epscu-epsc0)
	   else   
! 
!   flat friction branch for strains larger than epscu
! 
	      sigc = fcu
	      Ect  = 1.d-10
!       Ect  = 0.d0
	   endif
	endif
!
	return
	end subroutine Compr_Envlp
		    
