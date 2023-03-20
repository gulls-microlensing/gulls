c
c     global data passed as common blocks
c
      implicit none

c     lens positions and masses
      real*8 z1, z2, m1, m2
      common /lensPass/ z1, z2, m1, m2

c     source position and radius, and limb darkening profile parameter us
      real*8 xs, ys, rs, Gamma, xsCenter, ysCenter
      common /sourcePass/ xs, ys, rs, Gamma, xsCenter, ysCenter
 
c     image properties: positions, magnifications, total magnification
c     number of images
      complex*16 images(5) 
      real*8 mu(5), muTotal
      integer nimage
      real*8 lambda1, lambda2, theta, phixx, phixy, phiyy

      common /imagePass/ images, mu, muTotal, lambda1, lambda2, theta,
     $     phixx, phixy, phiyy, nimage

c     image track and segment properties
      integer bisectMax
      integer nmax              !maximum number of points in an image track
      parameter (bisectMax=18)
      parameter (nmax=2**bisectMax)

      integer nsegmentsMax      !maximum segments in the true images
      parameter (nsegmentsMax=16)

c     some mathematical constants
      real*8 pi
      real*8 dtor
      parameter (pi=3.14159265358979323846264338327950288d0)
      parameter (dtor=pi/180d0)

      logical debug
      common /debugPass/ debug

      real*8 rootEPS
      parameter (rootEPS=1.0d-8)

