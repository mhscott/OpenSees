c#######################################################################
      subroutine Hyster_1(matpar,hstvP,hstv,epsP,sigP,deps,sig,E,ist)
c-----------------------------------------------------------------------
c Hysteretic Moment-Rotation Relation with Pinching and Damage
c implemented in FEDEAS library on November 24, 1994 by F.Filippou
c-----------------------------------------------------------------------
c The relation is based on the following fixed and history variables
c
c Fixed Material Parameters
c   mom1p  = matpar(1) : positive moment at first transition
c   rot1p  = matpar(2) : positive rotation at first transition
c   mom2p  = matpar(3) : ultimate positive moment
c   rot2p  = matpar(4) : ultimate positive rotation
c   mom1n  = matpar(5) : negative moment at first transition
c   rot1n  = matpar(6) : negative rotation at first transition
c   mom2n  = matpar(7) : ultimate negative moment
c   rot2n  = matpar(8) : ultimate negative rotation
c   pinchx = matpar(9) : pinching parameter in x
c   pinchy = matpar(10): pinching parameter in y
c   damfc1 = matpar(11): damage factor for max. deformation ratio
c   damfc2 = matpar(12): damage factor for energy dissipation
c
c History Variables
c   rotmax = hstvP(1): maximum previous rotation (positive)
c   rotmin = hstvP(2): minimum previous rotation (negative)
c   rotpu  = hstvP(3): last rotation for complete positive unloading
c   rotnu  = hstvP(4): last rotation for complete negative unloading
c   enrgyd = hstvP(5): total energy dissipation at end of last step
c   kon    = hstvP(6): unloading indicator
c             (1 if previous rotation increment was positive)
c             (2 if previous rotation increment was negative)
c-----------------------------------------------------------------------
      implicit none
      
      integer kon,ist
      
      real*8 matpar(12),hstvP(6),hstv(6)
      real*8 epsP,sigP,deps,sig,E
c
      real*8 mom1p,rot1p,mom2p,rot2p,mom1n,rot1n,mom2n,rot2n
      real*8 E1p,E2p,E1n,E2n
      real*8 pinchy,pinchx,damfc,damfc1,damfc2
      real*8 rotmax,rotmin,rotpu,rotnu,enrgyd,enrgya
      real*8 rot,maxmom,minmom,rotch
      real*8 rotmp1,rotmp2,tmpmo1,tmpmo2,enrgy,dummy,rotun      
c-----------------------------------------------------------------------
      if (ist.eq.0) then
c
c     read material data and identify no of history variables
c
         hstvP (1) = 12         ! no of material parameters for material model
         hstv  (1) =  6         ! no of history variables for material model
c
      else
c
c     material state determination     
c     retrieve concrete fixed material properties
c 
         mom1p  = matpar (1)
         rot1p  = matpar (2)
         mom2p  = matpar (3)
         rot2p  = matpar (4)
         mom1n  = matpar (5)
         rot1n  = matpar (6)
         mom2n  = matpar (7)
         rot2n  = matpar (8)
         pinchx = matpar (9)
         pinchy = matpar(10)
         damfc1 = matpar(11)
         damfc2 = matpar(12)
c
c     retrieve history variables of hysteretic moment-rotation
c 
         rotmax = hstvP(1)
         rotmin = hstvP(2)
         rotpu  = hstvP(3)
         rotnu  = hstvP(4)
         enrgyd = hstvP(5)
         kon    = idint(hstvP(6))
c
c     calculate different stiffness values
c
         E1p = mom1p/rot1p
         E1n = mom1n/rot1n
         E2p =	(mom2p-mom1p)/(rot2p-rot1p)
         E2n = (mom2n-mom1n)/(rot2n-rot1n)
c     
c     calculate energy absorption at ultimate: enrgya
c 
         enrgya = 0.5*rot1p*mom1p + 0.5*(mom1p+mom2p)*(rot2p-rot1p)
     $        + 0.5*rot1n*mom1n + 0.5*(mom1n+mom2n)*(rot2n-rot1n)
c 
c     set unloading index kon
c 
          if (kon.eq.0) then
             enrgyd = 0.d0
             if (deps.ge.0.d0) then
                kon = 1
             else
                kon = 2
             endif
          endif      
c 
c     update rotation value
c 
          rot = epsP + deps
c 
c     if the current rotation value exceeds the maximum previous positive value,
c     call the positive envelope routine. 
          if (rot.gt.rotmax) then
             call Pos_Envlp (E1p,mom1p,E2p,rot,sig,E)
             rotmax = rot
          else if (rot.lt.rotmin) then
             call Neg_Envlp (E1n,mom1n,E2n,rot,sig,E)
             rotmin = rot
          else
c 
c     for positive rotation increment
             if (deps.ge.0.d0) then
c     
c     If unloading takes place, save rotation and moment value.
c     if unloading takes place from the maximum value, update the extreme
c     rotation on the opposite quadrant with the damage parameter
c 
                if (kon.eq.2) then
                   kon = 1
                   if (sigP.le.0.d0) then
                      rotnu = epsP - sigP/E1n
                      enrgy = enrgyd-0.5*(sigP/E1p)*sigP
                      damfc = 0.d0
                      if (rotmin.lt.rot1n) damfc = damfc2*(enrgy/enrgya)
                      if (rotmin.lt.rot1n.and.(epsP.eq.rotmin))
     $                     damfc  = damfc + damfc1*(rotmax-rot1p)/rot1p   
                      rotmax = rotmax*(1.d0+damfc)
                   endif
                endif
c     Reloading behavior after yielding has taken place in the opposite 
c     loading direction (as demonstrated by the value of variable rotnu)
c     The reloading branch consists of three branches itself
c 
                rotmax = dmax1(rotmax,mom1p/E1p)
                call Pos_Envlp (E1p,mom1p,E2p,rotmax,maxmom,dummy)
                rotmp1 = rotnu + pinchy*(rotmax-rotnu)
                rotmp2 = rotmax-(1.d0-pinchy)*maxmom/E1p
                rotch  = rotmp1 + (rotmp2-rotmp1)*pinchx
c     
c     unloading branch
c 
                if (rot.lt.rotnu) then
c     sig = sigP + E*deps
                   E = E1n
                   sig = dmin1(0.d0,sigP + E*deps) ! changed on 7/18/98
                   if (sig.eq.0.d0) E = 1.d-9
c 
c     first reloading branch; for partial un- and reloading cycle check
c     whether the linear elastic un- and reloading behavior controls 
c 
                else if (rot.ge.rotnu.and.rot.lt.rotch) then
                   E  = maxmom*pinchy/(rotch-rotnu)
                   tmpmo2 = (rot-rotnu)*E
                   tmpmo1 = sigP + E1p*deps
                   sig  = dmin1(tmpmo1,tmpmo2)
                   if (sig.eq.tmpmo1) E = E1p
c     
c     second reloading branch; for partial un- and reloading cycle check
c     whether the linear elastic un- and reloading behavior controls 
c 
                else
                   E = (1.d0-pinchy)*maxmom/(rotmax-rotch)
                   tmpmo2 = pinchy*maxmom+(rot-rotch)*E
                   tmpmo1 = sigP + E1p*deps
                   sig  = dmin1(tmpmo1,tmpmo2)
                   if (sig.eq.tmpmo1) E = E1p
                endif
             else
c     for negative rotation increment
c     
c     if unloading takes place, save rotation and moment value
c     if unloading takes place from the maximum value, update the extreme
c     rotation on the opposite quadrant with the damage parameter
c     
                if (kon.eq.1) then
                   kon = 2
                   if (sigP.gt.0.d0) then
                      rotpu = epsP - sigP/E1p
                      enrgy = enrgyd-0.5*(sigP/E1p)*sigP
                      damfc = 0.d0
                      if (rotmax.gt.rot1p) damfc = damfc2*(enrgy/enrgya)
                      if (rotmax.gt.rot1p.and.(epsP.eq.rotmax))
     $                     damfc = damfc + damfc1*(rotmin-rot1n)/rot1n
                      rotmin = rotmin*(1.d0+damfc)
                   endif
                endif
c 
c     Reloading behavior after yielding has taken place in the opposite 
c     loading direction (as demonstrated by the value of variable rotpu)
c     The reloading branch consists of three branches itself
c 
                rotmin = dmin1(rotmin,mom1n/E1n)
                call Neg_Envlp (E1n,mom1n,E2n,rotmin,minmom,dummy)
                rotmp1 = rotpu + pinchy*(rotmin-rotpu)
                rotmp2 = rotmin-(1.d0-pinchy)*minmom/E1n
                rotch  = rotmp1 + (rotmp2-rotmp1)*pinchx
c 
c     unloading branch
c 
                if (rot.gt.rotpu) then
c     sig = sigP + E*deps
                   E = E1p							      
                   sig = dmax1(0.d0,sigP + E*deps) ! changed on 7/18/98
                   if (sig.eq.0.d0) E = 1.d-9
c     
c     first reloading branch; for partial un- and reloading cycle check
c     whether the linear elastic un- and reloading behavior controls 
c 
                else if (rot.le.rotpu.and.rot.gt.rotch) then
                   E  = minmom*pinchy/(rotch-rotpu)
                   tmpmo2 = (rot-rotpu)*E
                   tmpmo1 = sigP + E1p*deps
                   sig  = dmax1(tmpmo1,tmpmo2)
                   if (sig.eq.tmpmo1) E = E1p
c     
c     second reloading branch; for partial un- and reloading cycle check
c     whether the linear elastic un- and reloading behavior controls 
c     
                else
                   E  = (1.d0-pinchy)*minmom/(rotmin-rotch)
                   tmpmo2 = pinchy*minmom+(rot-rotch)*E
                   tmpmo1 = sigP + E1p*deps
                   sig  = dmax1(tmpmo1,tmpmo2)
                   if (sig.eq.tmpmo1) E = E1p
                endif
             endif
          endif

          enrgyd = enrgyd + 0.5*(sigP+sig)*deps
c 
c	calculate appropriate stiffness depending on value of ist
c 
	if (ist.eq.2) then
           E = (sig-sigP)/deps
	else if (ist.eq.3) then
           if (kon.eq.1) then
              rotun = rotnu
           else
              rotun = rotpu
           endif
           E = sig/(epsP+deps-rotun)
	else
c     add additional cases, if needed
	endif     
c 
c	transfer all variables to history vector and return
c 
	hstv(1) = rotmax 
	hstv(2) = rotmin
	hstv(3) = rotpu
	hstv(4) = rotnu
	hstv(5) = enrgyd 
	hstv(6) = kon
c     
      endif
      
      return
      end subroutine Hyster_1
c#######################################################################
      subroutine Pos_Envlp (E1,My,E2,rot,momnt,stiff)
c-----------------------------------------------------------------------
c     Bilinear elastic-strain hardening envelope (positive)
c     
c     Material parameters:
c     
c     E1 = elastic (pre-yield stiffness)
c     My = yield moment
c     E2 = strain hardening stiffness
c     
c     The subroutine returns the value of moment (momnt) and current
c     stiffness (stiff) for the given rotation value (rot)
c-----------------------------------------------------------------------
      implicit none
c
      real*8  E1,My,E2,rot
      real*8  momnt,stiff
c
c local variables      
      real*8 roty,rotu
c     
      roty = My/E1
c      
      if (rot.lt.0.d0) then
	 momnt = 0.d0
         stiff = 1.d-9
      else if (rot.ge.0.d0.and.rot.le.roty) then
	 stiff = E1
	 momnt = E1*rot
      else
	 if (E2.ge.0.d0) then
            stiff = E2
            momnt = My + E2*(rot-roty)
	 else
            rotu  = roty - My/E2
            if (rot.lt.rotu) then
               stiff = E2
               momnt = My + E2*(rot-roty)
            else
               momnt = 0.d0
               stiff = 1.d-9
            endif
	 endif
      endif
c 
      return
      end subroutine Pos_Envlp
c#######################################################################
      subroutine Neg_Envlp (E1,My,E2,rot,momnt,stiff)
c-----------------------------------------------------------------------
c     Bilinear elastic-strain hardening envelope (negative)
c     
c     Material parameters:
c     
c     E1 = elastic (pre-yield stiffness)
c     My = yield moment
c     E2 = strain hardening stiffness
c     
c     The subroutine returns the value of moment (momnt) and current
c     stiffness (stiff) for the given rotation value (rot)
c-----------------------------------------------------------------------
      implicit none
c
      real*8  E1,My,E2,rot
      real*8  momnt,stiff
c
c local variables 
      real*8  roty,rotu
c     
      roty = My/E1
c      
      if (rot.gt.0.d0) then
	 momnt = 0.d0
         stiff = 1.d-9
      else if (rot.le.0.d0.and.rot.ge.roty) then
         stiff = E1
         momnt = E1*rot
      else
	 if (E2.ge.0.d0) then
            stiff = E2
            momnt = My + E2*(rot-roty)
	 else
            rotu  = roty - My/E2
            if (rot.gt.rotu) then
               stiff = E2
               momnt = My + E2*(rot-roty)
            else
               momnt = 0.d0
               stiff = 1.d-9
            endif
	 endif
      endif
c 
      return
      end subroutine Neg_Envlp
c#######################################################################
