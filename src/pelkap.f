	subroutine pelkap(xmom,para)
c***********************************************************************
c*                                                                     *
c*  fortran code written for inclusion in ibm research report rc20525, *
c*  'fortran routines for use with the method of l-moments, version 3' *
c*                                                                     *
c*  j. r. m. hosking                                                   *
c*  ibm research division                                              *
c*  t. j. watson research center                                       *
c*  yorktown heights                                                   *
c*  new york 10598, u.s.a.                                             *
c*                                                                     *
c*  version 3     august 1996                                          *
c*                                                                     *
c***********************************************************************
c
c  parameter estimation via l-moments for the kappa distribution
c
c  parameters of routine:
c  xmom   * input* array of length 4. contains the l-moments lambda-1,
c                  lambda-2, tau-3, tau-4.
c  para   *output* array of length 5. on exit, contains the parameters
c                  in the order xi, alpha, k, h and the ifail flag.
c  ifail           fail flag. on exit, it is set as follows.
c                  0  successful exit
c                  1  l-moments invalid
c                  2  (tau-3, tau-4) lies above the generalized-logistic
c                     line (suggests that l-moments are not consistent
c                     with any kappa distribution with h.gt.-1)
c                  3  iteration failed to converge
c                  4  unable to make progress from current point in
c                     iteration
c                  5  iteration encountered numerical difficulties -
c                     overflow would have been likely to occur
c                  6  iteration for h and k converged, but overflow
c                     would have occurred when calculating xi and alpha
c
c  n.b.  parameters are sometimes not uniquely defined by the first 4
c  l-moments. in such cases the routine returns the solution for which
c  the h parameter is largest.
c
c  other routines used: dlgama,digamd
c
c  the shape parameters k and h are estimated using newton-raphson
c  iteration on the relationship between (tau-3,tau-4) and (k,h).
c  the convergence criterion is that tau-3 and tau-4 calculated from
c  the estimated values of k and h should differ by less than 'eps'
c  from the values supplied in array xmom.
c
      implicit double precision (a-h,o-z)
      double precision xmom(4),para(5)
      integer ifail
      data zero/0d0/,half/0.5d0/,one/1d0/,two/2d0/,three/3d0/,four/4d0/
      data five/5d0/,six/6d0/,twelve/12d0/,twenty/20d0/,thirty/30d0/
      data p725/0.725d0/,p8/0.8d0/
c
c         eps,maxit control the test for convergence of n-r iteration
c         maxsr is the max. no. of steplength reductions per iteration
c         hstart is the starting value for h
c         big is used to initialize the criterion function
c         oflexp is such that dexp(oflexp) just does not cause overflow
c         oflgam is such that dexp(dlgama(oflgam)) just does not cause
c           overflow
c
      data eps/1d-6/,maxit/20/,maxsr/10/,hstart/1.001d0/,big/10d0/
      data oflexp/170d0/,oflgam/53d0/
c
 
      t3=xmom(3)
      t4=xmom(4)
      do 10 i=1,4
   10 para(i)=zero
c
c         test for feasibility
c
      if(xmom(2).le.zero)goto 1000
      if(dabs(t3).ge.one.or.dabs(t4).ge.one)goto 1000
      if(t4.le.(five*t3*t3-one)/four)goto 1000
      if(t4.ge.(five*t3*t3+one)/six )goto 1010
c
c         set starting values for n-r iteration:
c         g is chosen to give the correct value of tau-3 on the
c         assumption that h=1 (i.e. a generalized pareto fit) -
c         but h is actually set to 1.001 to avoid numerical
c         difficulties which can sometimes arise when h=1 exactly
c
      g=(one-three*t3)/(one+t3)
      h=hstart
      z=g+h*p725
      xdist=big
c
c         start of newton-raphson iteration
c
      do 100 it=1,maxit
c
c         reduce steplength until we are nearer to the required
c         values of tau-3 and tau-4 than we were at the previous step
c
      do 40 i=1,maxsr
c
c         - calculate current tau-3 and tau-4
c
c           notation:
c           u.    - ratios of gamma functions which occur in the pwm's
c                   beta-sub-r
c           alam. - l-moments (apart from a location and scale shift)
c           tau.  - l-moment ratios
c
      if(g.gt.oflgam)goto 1020
      if(h.gt.zero)goto 20
      u1=dexp(dlgama(  -one/h-g)-dlgama(  -one/h+one))
      u2=dexp(dlgama(  -two/h-g)-dlgama(  -two/h+one))
      u3=dexp(dlgama(-three/h-g)-dlgama(-three/h+one))
      u4=dexp(dlgama( -four/h-g)-dlgama( -four/h+one))
      goto 30
   20 u1=dexp(dlgama(  one/h)-dlgama(  one/h+one+g))
      u2=dexp(dlgama(  two/h)-dlgama(  two/h+one+g))
      u3=dexp(dlgama(three/h)-dlgama(three/h+one+g))
      u4=dexp(dlgama( four/h)-dlgama( four/h+one+g))
   30 continue
      alam2=u1-two*u2
      alam3=-u1+six*u2-six*u3
      alam4=u1-twelve*u2+thirty*u3-twenty*u4
      if(alam2.eq.zero)goto 1020
      tau3=alam3/alam2
      tau4=alam4/alam2
      e1=tau3-t3
      e2=tau4-t4
c
c         - if nearer than before, exit this loop
c
      dist=dmax1(dabs(e1),dabs(e2))
      if(dist.lt.xdist)goto 50
c
c         - otherwise, halve the steplength and try again
c
      del1=half*del1
      del2=half*del2
      g=xg-del1
      h=xh-del2
   40 continue
c
c         too many steplength reductions
c
      ifail=4
      para(5) = ifail
      return
c
c         test for convergence
c
   50 continue
      if(dist.lt.eps)goto 110
c
c         not converged: calculate next step
c
c         notation:
c         u1g  - derivative of u1 w.r.t. g
c         dl2g - derivative of alam2 w.r.t. g
c         d..  - matrix of derivatives of tau-3 and tau-4 w.r.t. g and h
c         h..  - inverse of derivative matrix
c         del. - steplength
c
      xg=g
      xh=h
      xz=z
      xdist=dist
      rhh=one/(h*h)
      if(h.gt.zero)goto 60
      u1g=-u1*digamd(  -one/h-g)
      u2g=-u2*digamd(  -two/h-g)
      u3g=-u3*digamd(-three/h-g)
      u4g=-u4*digamd( -four/h-g)
      u1h=      rhh*(-u1g-u1*digamd(  -one/h+one))
      u2h=  two*rhh*(-u2g-u2*digamd(  -two/h+one))
      u3h=three*rhh*(-u3g-u3*digamd(-three/h+one))
      u4h= four*rhh*(-u4g-u4*digamd( -four/h+one))
      goto 70
   60 u1g=-u1*digamd(  one/h+one+g)
      u2g=-u2*digamd(  two/h+one+g)
      u3g=-u3*digamd(three/h+one+g)
      u4g=-u4*digamd( four/h+one+g)
      u1h=      rhh*(-u1g-u1*digamd(  one/h))
      u2h=  two*rhh*(-u2g-u2*digamd(  two/h))
      u3h=three*rhh*(-u3g-u3*digamd(three/h))
      u4h= four*rhh*(-u4g-u4*digamd( four/h))
   70 continue
      dl2g=u1g-two*u2g
      dl2h=u1h-two*u2h
      dl3g=-u1g+six*u2g-six*u3g
      dl3h=-u1h+six*u2h-six*u3h
      dl4g=u1g-twelve*u2g+thirty*u3g-twenty*u4g
      dl4h=u1h-twelve*u2h+thirty*u3h-twenty*u4h
      d11=(dl3g-tau3*dl2g)/alam2
      d12=(dl3h-tau3*dl2h)/alam2
      d21=(dl4g-tau4*dl2g)/alam2
      d22=(dl4h-tau4*dl2h)/alam2
      det=d11*d22-d12*d21
      h11= d22/det
      h12=-d12/det
      h21=-d21/det
      h22= d11/det
      del1=e1*h11+e2*h12
      del2=e1*h21+e2*h22
c
c         take next n-r step
c
      g=xg-del1
      h=xh-del2
      z=g+h*p725
c
c         reduce step if g and h are outside the parameter space
c
      factor=one
      if(g.le.-one)factor=p8*(xg+one)/del1
      if(h.le.-one)factor=dmin1(factor,p8*(xh+one)/del2)
      if(z.le.-one)factor=dmin1(factor,p8*(xz+one)/(xz-z))
      if(h.le.zero.and.g*h.le.-one)
     *  factor=dmin1(factor,p8*(xg*xh+one)/(xg*xh-g*h))
      if(factor.eq.one)goto 80
      del1=del1*factor
      del2=del2*factor
      g=xg-del1
      h=xh-del2
      z=g+h*p725
   80 continue
c
c         end of newton-raphson iteration
c
  100 continue
c
c         not converged
c
      ifail=3
      para(5) = ifail
      return
c
c         converged
c
  110 ifail=0
      para(4)=h
      para(3)=g
      temp=dlgama(one+g)
      if(temp.gt.oflexp)goto 1030
      gam=dexp(temp)
      temp=(one+g)*dlog(dabs(h))
      if(temp.gt.oflexp)goto 1030
      hh=dexp(temp)
      para(2)=xmom(2)*g*hh/(alam2*gam)
      para(1)=xmom(1)-para(2)/g*(one-gam*u1/hh)
	para(3)=-g
      para(5) = ifail
      return
c
 1000 ifail=1
      para(5) = ifail
      return
 1010 ifail=2
      para(5) = ifail
      return
 1020 ifail=5
      para(5) = ifail
      return
 1030 ifail=6
      para(5) = ifail
      return
      
c
      end subroutine

c===================================================== digamd.for
      double precision function digamd(x)
c***********************************************************************
c*                                                                     *
c*  fortran code written for inclusion in ibm research report rc20525, *
c*  'fortran routines for use with the method of l-moments, version 3' *
c*                                                                     *
c*  j. r. m. hosking                                                   *
c*  ibm research division                                              *
c*  t. j. watson research center                                       *
c*  yorktown heights                                                   *
c*  new york 10598, u.s.a.                                             *
c*                                                                     *
c*  version 3     august 1996                                          *
c*                                                                     *
c***********************************************************************
c
c  digamma function (euler's psi function) - the first derivative of
c  log(gamma(x))
c
c  based on algorithm as103, appl. statist. (1976) vol.25 no.3
c
      implicit double precision (a-h,o-z)
      data zero/0d0/,half/0.5d0/,one/1d0/
      data small/1d-9/,crit/13d0/
c
c         c1...c7 are the coeffts of the asymptotic expansion of digamd
c         d1 is  -(euler's constant)
c
      data c1,c2,c3,c4,c5,c6,c7,d1/
     *  0.83333 33333 33333 333d-1,  -0.83333 33333 33333 333d-2,
     *  0.39682 53968 25396 825d-2,  -0.41666 66666 66666 666d-2,
     *  0.75757 57575 75757 575d-2,  -0.21092 79609 27960 928d-1,
     *  0.83333 33333 33333 333d-1,  -0.57721 56649 01532 861d 0/
      digamd=zero
      if(x.le.zero)goto 1000
c
c         use small-x approximation if x.le.small
c
      if(x.gt.small)goto 10
      digamd=d1-one/x
      return
c
c         reduce to digamd(x+n) where x+n.ge.crit
c
   10 y=x
   20 if(y.ge.crit)goto 30
      digamd=digamd-one/y
      y=y+one
      goto 20
c
c         use asymptotic expansion if y.ge.crit
c
   30 digamd=digamd+dlog(y)-half/y
      y=one/(y*y)
      sum=((((((c7*y+c6)*y+c5)*y+c4)*y+c3)*y+c2)*y+c1)*y
      digamd=digamd-sum
      return
c
 1000 write(6,7000)x
      return
c
 7000 format(' *** error *** routine digamd :',
     *  ' argument out of range :',d24.16)
      end function
c===================================================== dlgama.for
      double precision function dlgama(x)
c***********************************************************************
c*                                                                     *
c*  fortran code written for inclusion in ibm research report rc20525, *
c*  'fortran routines for use with the method of l-moments, version 3' *
c*                                                                     *
c*  j. r. m. hosking                                                   *
c*  ibm research division                                              *
c*  t. j. watson research center                                       *
c*  yorktown heights                                                   *
c*  new york 10598, u.s.a.                                             *
c*                                                                     *
c*  version 3     august 1996                                          *
c*                                                                     *
c***********************************************************************
c
c  logarithm of gamma function
c
c  based on algorithm acm291, commun. assoc. comput. mach. (1966)
c
      implicit double precision (a-h,o-z)
      data small,crit,big,toobig/1d-7,13d0,1d9,2d36/
c
c         c0 is 0.5*log(2*pi)
c         c1...c7 are the coeffts of the asymptotic expansion of dlgama
c
      data c0,c1,c2,c3,c4,c5,c6,c7/
     *   0.91893 85332 04672 742d 0,  0.83333 33333 33333 333d-1,
     *  -0.27777 77777 77777 778d-2,  0.79365 07936 50793 651d-3,
     *  -0.59523 80952 38095 238d-3,  0.84175 08417 50841 751d-3,
     *  -0.19175 26917 52691 753d-2,  0.64102 56410 25641 026d-2/
c
c         s1 is -(euler's constant), s2 is pi**2/12
c
      data s1/-0.57721 56649 01532 861d 0/
      data s2/ 0.82246 70334 24113 218d 0/
c
      data zero/0d0/,half/0.5d0/,one/1d0/,two/2d0/
      dlgama=zero
      if(x.le.zero)goto 1000
      if(x.gt.toobig)goto 1000
c
c         use small-x approximation if x is near 0, 1 or 2
c
      if(dabs(x-two).gt.small)goto 10
      dlgama=dlog(x-one)
      xx=x-two
      goto 20
   10 if(dabs(x-one).gt.small)goto 30
      xx=x-one
   20 dlgama=dlgama+xx*(s1+xx*s2)
      return
   30 if(x.gt.small)goto 40
      dlgama=-dlog(x)+s1*x
      return
c
c         reduce to dlgama(x+n) where x+n.ge.crit
c
   40 sum1=zero
      y=x
      if(y.ge.crit)goto 60
      z=one
   50 z=z*y
      y=y+one
      if(y.lt.crit)goto 50
      sum1=sum1-dlog(z)
c
c         use asymptotic expansion if y.ge.crit
c
   60 sum1=sum1+(y-half)*dlog(y)-y+c0
      sum2=zero
      if(y.ge.big)goto 70
      z=one/(y*y)
      sum2=((((((c7*z+c6)*z+c5)*z+c4)*z+c3)*z+c2)*z+c1)/y
   70 dlgama=sum1+sum2
      return
c
 1000 write(6,7000)x
      return
c
 7000 format(' *** error *** routine dlgama :',
     *  ' argument out of range :',d24.16)
      end function
