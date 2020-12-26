c#######################################################################
      subroutine Hyster_2(matpar,hstvP,hstv,epsP,sigP,deps,sig,E,ist)
c-----------------------------------------------------------------------
c Hysteretic Stress-Strain Relation with Trilinear Envelope
c implemented in FEDEAS material library on April 21, 1995 by F.Filippou
c-----------------------------------------------------------------------
c The relation is based on the following fixed and history variables
c
c Fixed Material Parameters
c   mom1p   = matpar(1) : positive moment at first transition
c   rot1p   = matpar(2) : positive rotation at first transition
c   mom2p   = matpar(3) : positive moment at second transition
c   rot2p   = matpar(4) : positive rotation at second transition
c   mom3p   = matpar(5) : ultimate positive moment
c   rot3p   = matpar(6) : ultimate positive rotation
c   mom1n   = matpar(7) : negative moment at first transition
c   rot1n   = matpar(8) : negative rotation at first transition
c   mom2n   = matpar(9) : negative moment at second transition
c   rot2n   = matpar(10): negative rotation at second transition
c   mom3n   = matpar(11): ultimate negative moment
c   rot3n   = matpar(12): ultimate negative rotation
c   pinchx  = matpar(13): pinching parameter in x
c   pinchy  = matpar(14): pinching parameter in y
c   damfc1  = matpar(15): damage factor for max. deformation ratio
c   damfc2  = matpar(16): damage factor for energy dissipation
c
c History Variables
c   rotmax = hstvP(1): maximum previous rotation (positive)
c   rotmin = hstvP(2): minimum previous rotation (negative)
c   rotpu  = hstvP(3): last rotation for complete positive unloading
c   rotnu  = hstvP(4): last rotation for complete negative unloading
c   enrgyd = hstvP(5): total energy dissipation at end of last step
c   kon    = hstvP(6): unloading indicator
c                     (1 if previous rotation increment was positive)
c                     (2 if previous rotation increment was negative)
c-----------------------------------------------------------------------
      implicit none
      
      integer kon,ist
      real*8 matpar(16),hstvP(6),hstv(6)
      real*8 epsP,sigP,deps,sig,E
c     
      real*8 mom1p,rot1p,mom2p,rot2p,mom3p,rot3p
      real*8 mom1n,rot1n,mom2n,rot2n,mom3n,rot3n
      real*8 E1p,E2p,E3p,E1n,E2n,E3n
      real*8 pinchy,pinchx,damfc,damfc1,damfc2
      real*8 rotmax,rotmin,rotpu,rotnu,enrgyd,enrgya
      real*8 rot,maxmom,minmom,rotch,rotlim,rotrel
      real*8 rotmp1,rotmp2,tmpmo1,tmpmo2,enrgy,dummy1,dummy2,rotun      
c-----------------------------------------------------------------------
      if (ist.eq.0) then
c
c  read material data and identify no of history variables
c
         hstvP (1) = 16         ! no of material parameters for material model
         hstv  (1) =  6         ! no of history variables for material model
c
      else
c
c  material state determination     
c  retrieve concrete fixed material properties
c 
         mom1p   = matpar (1)
         rot1p   = matpar (2)
         mom2p   = matpar (3)
         rot2p   = matpar (4)
         mom3p   = matpar (5)
         rot3p   = matpar (6)
         mom1n   = matpar (7)
         rot1n   = matpar (8)
         mom2n   = matpar (9)
         rot2n   = matpar(10)
         mom3n   = matpar(11)
         rot3n   = matpar(12)
c
         pinchx  = matpar(13)
         pinchy  = matpar(14)
         damfc1  = matpar(15)
         damfc2  = matpar(16)
c 
c	retrieve history variables of hysteretic moment-rotation
c 
         rotmax = hstvP(1)
         rotmin = hstvP(2)
         rotpu  = hstvP(3)
         rotnu  = hstvP(4)
         enrgyd = hstvP(5)
         kon    = idint(hstvP(6))
c
c  calculate different stiffness values
c
         E1p = mom1p/rot1p
         E1n = mom1n/rot1n
         E2p =	(mom2p-mom1p)/(rot2p-rot1p)
         E2n = (mom2n-mom1n)/(rot2n-rot1n)
         E3p =	(mom3p-mom2p)/(rot3p-rot2p)
         E3n = (mom3n-mom2n)/(rot3n-rot2n)
c 
c	calculate energy absorption at ultimate: enrgya
c 
         enrgya = 0.5*rot1p*mom1p + 0.5*(mom1p+mom2p)*(rot2p-rot1p) +
     $        0.5*(mom2p+mom3p)*(rot3p-rot2p) +
     $        0.5*rot1n*mom1n +	0.5*(mom1n+mom2n)*(rot2n-rot1n) +
     $        0.5*(mom2n+mom3n)*(rot3n-rot2n)
c 
c	set unloading index kon
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
c	update rotation value
c 
         rot = epsP + deps
c 
         if (rot.gt.rotmax) then
c     if the current rotation value exceeds the maximum previous positive value,
c     call the positive envelope routine. 
            call Pos_Envlp2 (E1p,mom1p,E2p,rot2p,E3p,rot,sig,dummy1,E)
            rotmax = rot
         elseif (rot.lt.rotmin) then
c     if the current rotation value exceeds the minimum previous negative value,
c     call the negative envelope routine. 
            call Neg_Envlp2 (E1n,mom1n,E2n,rot2n,E3n,rot,sig,dummy1,E)
            rotmin = rot
         else
c     unload/reloading branch
c     
c     for positive rotation increment
            if (deps.ge.0.d0) then
c     
c     if unloading takes place, save rotation and moment value;
c     if unloading takes place from the maximum value, update the extreme
c     rotation on the opposite quadrant with the damage parameter
               if (kon.eq.2) then
                  kon = 1
                  if (sigP.le.0.d0) then
                     rotnu = epsP - sigP/E1n
                     enrgy = enrgyd-0.5*(sigP/E1p)*sigP
                     damfc = 0.d0
                     if (rotmin.lt.rot1n) damfc = damfc2*(enrgy/enrgya)
                     if (rotmin.lt.rot1n.and.(epsP.eq.rotmin))
     $                    damfc = damfc + damfc1*(rotmax-rot1p)/rot1p   
                     rotmax = rotmax*(1.d0+damfc)
                  endif
               endif
c     
c     reloading behavior consisting of three branches
c     first, recalculate key values
c
               rotmax = dmax1(rotmax,mom1p/E1p)
               call Pos_Envlp2 (E1p,mom1p,E2p,rot2p,E3p,rotmax,maxmom,
     $              dummy1,dummy2)
               call Neg_Envlp2 (E1n,mom1n,E2n,rot2n,E3n,rotmin,dummy1,
     $              rotlim,dummy2)
               rotrel = dmax1(rotlim,rotnu)
               rotmp1 = rotrel + pinchy*(rotmax-rotrel)
               rotmp2 = rotmax-(1.d0-pinchy)*maxmom/E1p
               rotch  = rotmp1 + (rotmp2-rotmp1)*pinchx
c     
c			unloading branch
               if (rot.lt.rotnu) then
c     sig = sigP + E*deps
                  E = E1n
                  sig = dmin1(0.d0,sigP + E*deps) ! changed on 7/18/98
                  if (sig.eq.0.d0) E = 1.d-9
c     
               else if (rot.ge.rotnu.and.rot.lt.rotch) then
c     first reloading branch; for partial un- and reloading cycle check
c     whether the linear elastic un- and reloading behavior controls 
c     
                  if (rot.le.rotrel) then
                     sig = 0.d0
                     E = 0.d0
                  else
                     E  = maxmom*pinchy/(rotch-rotrel)
                     tmpmo2 = (rot-rotrel)*E
                     tmpmo1 = sigP + E1p*deps
                     sig  = dmin1(tmpmo1,tmpmo2)
                     if (sig.eq.tmpmo1) E = E1p
                  endif
                  
               else
c     second reloading branch; for partial un- and reloading cycle check
c     whether the linear elastic un- and reloading behavior controls 
c     
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
               if (kon.eq.1) then
                  kon = 2
                  if (sigP.ge.0.d0) then
                     rotpu = epsP - sigP/E1p
                     enrgy = enrgyd-0.5*(sigP/E1p)*sigP
                     damfc = 0.d0
                     if (rotmax.gt.rot1p) damfc = damfc2*(enrgy/enrgya)
                     if (rotmax.gt.rot1p.and.(epsP.eq.rotmax))  
     $                    damfc = damfc + damfc1*(rotmin-rot1n)/rot1n
                     rotmin = rotmin*(1.d0+damfc)
                  endif
               endif
c     determine key values
c     
               rotmin = dmin1(rotmin,mom1n/E1n)
               call Neg_Envlp2 (E1n,mom1n,E2n,rot2n,E3n,rotmin,minmom,
     $              dummy1,dummy2)
               call Pos_Envlp2 (E1p,mom1p,E2p,rot2p,E3p,rotmax,dummy1,
     $              rotlim,dummy2)
               rotrel = dmin1(rotlim,rotpu)
               rotmp1 = rotrel + pinchy*(rotmin-rotrel)
               rotmp2 = rotmin -(1.d0-pinchy)*minmom/E1n
               rotch  = rotmp1 + (rotmp2-rotmp1)*pinchx
               if (rot.gt.rotpu) then 
c     unloading branch
c     sig = sigP + E*deps
                  E = E1p							      
                  sig = dmax1(0.d0,sigP + E*deps) ! changed on 7/18/98
                  if (sig.eq.0.d0) E = 1.d-9
c     
               else if (rot.le.rotpu.and.rot.gt.rotch) then
c     first reloading branch; for partial un- and reloading cycle check
c     whether the linear elastic un- and reloading behavior controls 
c     
                  if (rot.ge.rotrel) then
                     sig = 0.d0
                     E = 0.d0
                  else
                     E  = minmom*pinchy/(rotch-rotrel)
                     tmpmo2 = (rot-rotrel)*E
                     tmpmo1 = sigP + E1p*deps
                     sig  = dmax1(tmpmo1,tmpmo2)
                     if (sig.eq.tmpmo1) E = E1p
                  endif
c     
c     second reloading branch; for partial un- and reloading cycle check
c     whether the linear elastic un- and reloading behavior controls 
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
c     E = sig/(epsP+deps)
         else
c     add additional cases, if needed
         endif     
c     
c     transfer all variables to history vector and return
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
      end subroutine Hyster_2
c#######################################################################
      subroutine Pos_Envlp2 (E1,M1,E2,rot2,E3,rot,momnt,rotlim,stiff)
c-----------------------------------------------------------------------
c     Bilinear elastic-strain hardening envelope (positive)
c     
c     Material parameters:
c     
c     E1   = initial stiffness
c     M1   = moment at first transition
c     E2   = second stiffness
c     rot2 = rotation at second transition
c     E3   = third stiffness
c     
c     The subroutine returns the value of moment (momnt) and current
c     stiffness (stiff) for the given rotation value (rot)
c-----------------------------------------------------------------------
      implicit none
c     
      real*8 E1,M1,E2,E3,rot2,rot,momnt,rotlim,stiff
      real*8 rot1
c     
      rot1   = M1/E1
      rotlim = 1.d+15
c     
      if (rot.lt.0.d0) then
	 momnt = 0.d0
         stiff = 1.d-9
      else if (rot.le.rot1) then
         momnt = E1*rot
         stiff = E1
      else if (rot.le.rot2) then
	 stiff = E2
         if (E2.ge.0.d0) then
            momnt = M1 + E2*(rot-rot1)
         else
            rotlim = rot1-M1/E2
            if (rot.le.rotlim) then
               momnt = M1 + E2*(rot-rot1)
            else
               momnt = 0.d0
            endif
         endif
      else
         stiff = E3
         if (stiff.ge.0.d0) then
            momnt = M1 + E2*(rot2-rot1)+E3*(rot-rot2)
         else
            rotlim = rot2-(M1 + E2*(rot2-rot1))/E3
            if (rot.le.rotlim) then
               momnt = M1 + E2*(rot2-rot1)+E3*(rot-rot2)
            else
               momnt = 0.d0
               stiff = 0.d0
            endif
         endif
      endif
c     
      return
      end subroutine Pos_Envlp2
c#######################################################################
      subroutine Neg_Envlp2 (E1,M1,E2,rot2,E3,rot,momnt,rotlim,stiff)
c-----------------------------------------------------------------------
c     Bilinear elastic-strain hardening envelope (negative)
c     
c     Material parameters:
c     
c     E1   = initial stiffness
c     M1   = moment at first transition
c     E2   = second stiffness
c     rot2 = rotation at second transition
c     E3   = third stiffness
c     
c     The subroutine returns the value of moment (momnt) and current
c     stiffness (stiff) for the given rotation value (rot)
c-----------------------------------------------------------------------
      implicit none
c     
      real*8 E1,M1,E2,E3,rot2,rot,momnt,rotlim,stiff
      real*8 rot1
c     
      rot1   = M1/E1
      rotlim = -1.d+15
c     
      if (rot.gt.0.d0) then
	 momnt = 0.d0
         stiff = 1.d-9
      else if (rot.ge.rot1) then
         momnt = E1*rot
         stiff = E1
      else if (rot.ge.rot2) then
         stiff = E2
         if (E2.ge.0.d0) then
            momnt = M1 + E2*(rot-rot1)
         else
            rotlim = rot1-M1/E2
            if (rot.ge.rotlim) then
               momnt = M1 + E2*(rot-rot1)
            else
               momnt = 0.d0
            endif
         endif
      else 
         stiff  = E3
         if (stiff.ge.0.d0) then
            momnt = M1 + E2*(rot2-rot1)+E3*(rot-rot2)
         else
            rotlim = rot2-(M1+E2*(rot2-rot1))/E3
            if (rot.ge.rotlim) then
               momnt = M1 + E2*(rot2-rot1)+E3*(rot-rot2)
            else
               momnt = 0.d0
               stiff = 0.d0
            endif
         endif
      endif
c     
      return
      
      end subroutine Neg_Envlp2
c#######################################################################
