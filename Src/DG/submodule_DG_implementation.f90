submodule(DG_structures) DG_implementation
!! The implementation of the interfaces defined in [[DG_structures]] as well
!! as internal routines that should never be called by a user.

  contains

    module procedure init_ref_element

implicit none

      real(wp), dimension(:,:), allocatable :: vinv
      integer(ip) :: n, i

! Set the order of the element.
      init_ref_element%n = order

! Calculate the total number of nodes within the element.
      n = order+1

! Allocate storage for the r coordinates.
      allocate ( init_ref_element%r(n) )

! Calculate the 1D nodal coordinates.
      call JacobiGL ( 0.0_wp, 0.0_wp, order, init_ref_element%r )

! Allocate storage for the integration weigths.
      allocate ( init_ref_element%w(n) )

! Calculate the integration weights.
      call GaussWeigths ( 0.0_wp, 0.0_wp, order, init_ref_element%r, &
                          init_ref_element%w )

! Allocate storage for the 1D Vandermonde matrix.
      allocate ( init_ref_element%v(n,n) )

! Calculate the 1D Vandermonde matrix.
      init_ref_element%v = VanderMonde1D ( order, init_ref_element%r )

! Allocate storage for the derivative matrix.
      allocate ( init_ref_element%dr(n,n) )

! Calculate the derivative matrix.
      init_ref_element%dr = Dmatrix1D ( order, init_ref_element%r, &
                                        init_ref_element%v )

! Allocate storage for the lift matrix.
      allocate (init_ref_element%lift(n,2) )

! Calculate the lift matrix.
      init_ref_element%lift = Lift1D ( order, init_ref_element%v ) 

      if ( present(sorder) .and. present(nc) ) then
        print*,'Setting up filter with order = ', sorder, ' and cutoff = ', nc
! Allocate storage for the filter matrix (and the inverse of v)
        allocate ( init_ref_element%filter(n,n), vinv(n,n) )
! Calculate the inverse of the Vandermonde matrix
        vinv = Inverse ( init_ref_element%v )
! Setup the filter
        init_ref_element%filter = Filter1D ( order, nc, sorder, &
                                             init_ref_element%v, vinv )
        deallocate ( vinv )
        init_ref_element%have_filter = .true.
      else
        init_ref_element%have_filter = .false.
      endif

    end procedure init_ref_element


    module procedure deallocate_ref_element

      if ( allocated(refel%r) ) deallocate ( refel%r )
      if ( allocated(refel%w) ) deallocate ( refel%w )
      if ( allocated(refel%v) ) deallocate ( refel%v )
      if ( allocated(refel%dr) ) deallocate ( refel%dr )
      if ( allocated(refel%lift) ) deallocate ( refel%lift )
      if ( allocated(refel%filter) ) deallocate ( refel%filter )
    end procedure deallocate_ref_element


    module procedure char_flux_real

      implicit none

      integer(ip) :: i, j, k
      real(wp), dimension(2), parameter :: nx = (/ -1.0_wp, 1.0_wp /)
      real(wp), dimension(nvar,nvar) :: lambdaminus, lambdaplus
      real(wp), dimension(2,nvar) :: du, nflux
      
      do j = 1, 2
        lambdaminus = rzero
        lambdaplus = rzero
        do k = 1, nvar
          if ( nx(j) * lambda(j,k) <= 0.0_wp ) then
            lambdaminus(k,k) = nx(j) * lambda(j,k)
          else
            lambdaplus(k,k) = nx(j) * lambda(j,k)
          end if
        end do
        nflux(j,:) = matmul(lambdaplus,matmul(sinv(j,:,:),uint(j,:))) &
                    +matmul(lambdaminus,matmul(sinv(j,:,:),uext(j,:)))
        nflux(j,:) = matmul(s(j,:,:),nflux(j,:))
        du(j,:) =  (nx(j)*flux(j,:) - nflux(j,:))
      end do
      char_flux_real = matmul(this%lift, du)

    end procedure char_flux_real


    module procedure char_flux_complex

      implicit none

      integer(ip) :: i, j, k
      real(wp), dimension(2), parameter :: nx = (/ -1.0_wp, 1.0_wp /)
      real(wp), dimension(nvar,nvar) :: lambdaminus, lambdaplus
      complex(wp), dimension(2,nvar) :: du, nflux
      
      do j = 1, 2
        lambdaminus = czero
        lambdaplus = czero
        do k = 1, nvar
          if ( nx(j) * lambda(j,k) <= 0.0_wp ) then
            lambdaminus(k,k) = nx(j) * lambda(j,k)
          else
            lambdaplus(k,k) = nx(j) * lambda(j,k)
          end if
        end do
        nflux(j,:) = matmul(lambdaplus,matmul(sinv(j,:,:),uint(j,:))) &
                    +matmul(lambdaminus,matmul(sinv(j,:,:),uext(j,:)))
        nflux(j,:) = matmul(s(j,:,:),nflux(j,:))
        du(j,:) =  (nx(j)*flux(j,:) - nflux(j,:))
      end do
      char_flux_complex = matmul(this%lift, du)
      if (debug_output) then
        print*,'lambdaminus = ', lambdaminus
        print*,'lambdaplus = ', lambdaplus
        print*,'nflux = ', nflux
        print*,'du = ', du
        print*,'lift.du = ', matmul(this%lift, du)
      end if

    end procedure char_flux_complex


    subroutine JacobiGQ(alpha, beta, n, x, w)
    !! Compute the n'th order Gauss quadrature points, \(x\) and weights,
    !! \(w\), associated with the Jacobi polynomial of type
    !! \((\alpha,\beta)>-1\).
      real(wp), intent(in) :: alpha
      !! The value of \(\alpha\).
      real(wp), intent(in) :: beta
      !! The value of \(\beta\).
      integer(ip), intent(in) :: n
      !! The order of the Jacobi polynomial.
      real(wp), dimension(n+1), intent(out) :: x
      !! On output the location of the Gauss quadrature points.
      real(wp), dimension(n+1), intent(out) :: w
      !! On output the weights.

      real(wp), dimension(1:n+1) :: h1, jdiag1, ireal
      real(wp), dimension(1:n) :: jdiag2, d1, d2, d3
      real(wp), dimension(max(1,2*n)) :: work
      real(wp), dimension(1:n+1,1:n+1) :: vect
      real(wp) :: lngammaab, lngammaa, lngammab
      integer(ip) :: info, i

      if (n==0) then
        x(1) = ( alpha - beta ) / ( alpha + beta + 2.0_wp )
        w(1) = 2.0_wp
        return
      end if

      forall(i=0:n) ireal(i+1) = real(i,wp)
      h1 = 2.0_wp*ireal+alpha+beta
      where (h1>0.0_wp)
        jdiag1 = -(alpha**2-beta**2)/(h1*(h1+2.0_wp))
      elsewhere
        jdiag1 = 0.0_wp
      end where
      d1 = 2.0_wp/(h1(1:n)+2.0_wp)
      d2 = ireal(2:n+1)*(ireal(2:n+1)+alpha+beta)* &
                        (ireal(2:n+1)+alpha)*(ireal(2:n+1)+beta)
      d3 = 1.0_wp/((h1(1:n)+1)*(h1(1:n)+3))
      jdiag2 = d1*sqrt(d2*d3)
      if (wp == dp ) then
        call dsteqr('I', n+1, jdiag1, jdiag2, vect, n+1, work, info )
      else if (wp == qp ) then
        call qsteqr('I', n+1, jdiag1, jdiag2, vect, n+1, work, info )
      end if
      if (info <0) then
        print*,'Parameter ', i, ' in call to dsteqr has illegal value'
        stop
      end if
      if (info >0) then
        print*, i, ' off-diagonal elements have not converged to zero in call to dsteqr'
        stop
      end if
      x = jdiag1
      lngammaab = log_gamma(alpha+beta+1.0_wp)
      lngammaa  = log_gamma(alpha+1.0_wp)
      lngammab  = log_gamma(beta+1.0_wp)
      w = vect(1,:)**2*2.0_wp**(alpha+beta+1)/(alpha+beta+1)* &
                       exp(lngammaa+lngammab-lngammaab)
      return
    end subroutine JacobiGQ


    subroutine JacobiGL(alpha, beta, n, r)
    !! Compute the n'th order Gauss Lobatto quadrature points, \(x\),
    !! associated with the Jacobi polynomial of type
    !! \((\alpha,\beta)>-1\).
      real(wp), intent(in) :: alpha
      !! The value of \(\alpha\).
      real(wp), intent(in) :: beta
      !! The value of \(\beta\).
      integer(ip), intent(in) :: n
      !! The order of the Jacobi polynomial.
      real(wp), dimension(n+1), intent(out) :: r
      !! On output the location of the Gauss Lobatto quadrature points.
      real(wp), dimension(n-1) :: w

      if ( n<=0 ) then
        print*,'JacobiGL called with n<=0. Aborting'
        stop
      end if
      r = 0.0_wp
      if ( n==1 ) then
        r(1) = -1.0_wp
        r(2) = 1.0_wp
        return
      end if

      call JacobiGQ ( alpha+1.0_wp, beta+1.0_wp, n-2, r(2:n), w )
      r(1) = -1.0_wp
      r(n+1) = 1.0_wp
      return
    end subroutine JacobiGL


    function JacobiP(x, alpha, beta, n)
    !! Function to evaluate Jacobi Polynomial \(P^{(\alpha,\beta)}_n\)
    !! of type \((\alpha,\beta)>-1\) (with \(\alpha+\beta\neq -1\)) at points
    !! \(x\) for order \(n\).
      real(wp), dimension(:), intent(in) :: x
      !! The points \(x\in[-1:1]\) at which to evaluate the Jacobi Polynomial. 
      real(wp), intent(in) :: alpha
      !! The value of \(\alpha\).
      real(wp), intent(in) :: beta
      !! The value of \(\beta\).
      integer(ip), intent(in) :: n
      !! The order of the Jacobi polynomial.
      real(wp), dimension(size(x)) :: JacobiP
      !! Returns \(P^{(\alpha,\beta)}_n(x)\).

      real(wp), dimension(n+1,size(x)) :: pl
      real(wp) :: lngammaab, lngammaa, lngammab, invsqgamma0, gamma0, gamma1
      real(wp) :: fac1, fac2
      real(wp) :: aold, anew, bnew, h1, ireal, irealp1
      integer(ip) :: i

      pl = 0.0_wp

      lngammaab = log_gamma(alpha+beta+1.0_wp)
      lngammaa  = log_gamma(alpha+1.0_wp)
      lngammab  = log_gamma(beta+1.0_wp)

      invsqgamma0 = 2.0_wp**(alpha+beta+1.0_wp)/(alpha+beta+1.0_wp)* &
                             exp(lngammaa+lngammab-lngammaab)
      gamma0 = 1.0_wp/sqrt(invsqgamma0)

      pl(1,:) = gamma0

      if (n==0) then
        JacobiP(:) = pl(1,:)
        return
      end if

      gamma1 = 0.5_wp*sqrt((alpha+beta+3.0_wp) &
                           /((alpha+1.0_wp)*(beta+1.0_wp)))*gamma0
      fac1 = (alpha+beta+2.0_wp)
      fac2 = (alpha-beta)
      pl(2,:) = gamma1 * ( fac1*x(:) + fac2 )

      if (n==1) then
        JacobiP(:) = pl(2,:)
        return
      end if

      aold = 2.0_wp / (2.0_wp+alpha+beta) * &
             sqrt ( (1.0_wp+alpha)*(1.0_wp+beta) / (3.0_wp+alpha+beta) )

      do i = 1, n-1
        ireal = real(i,wp)
        irealp1 = ireal+1.0_wp
        h1 = 2.0_wp*ireal+alpha+beta
        anew = 2.0_wp/(h1+2.0_wp)* &
               sqrt( irealp1*(irealp1+alpha+beta) * &
                     (irealp1+alpha) * (irealp1+beta) / &
                     (h1+1.0_wp)/(h1+3.0_wp) )
        bnew = - (alpha**2-beta**2) / (h1*(h1+2.0_wp) )
        pl(i+2,:) = 1.0_wp / anew * ( -aold*pl(i,:) + (x(:)-bnew)*pl(i+1,:) )
        aold = anew
      end do
      JacobiP(:) = pl(n+1,:)
      return
    end function JacobiP

    subroutine GaussWeigths( alpha, beta, n, r, w )
    !! Routine to calculate the integration weights for the Gauss-Lobatto
    !! Quadrature points.
    !!
    !! Note. I currently don't remember where this equation comes from.
      real(wp), dimension(:), intent(in) :: r
      !! The points for which the weights should be determined.
      real(wp), intent(in) :: alpha
      !! The value of \(\alpha\).
      real(wp), intent(in) :: beta
      !! The value of \(\beta\).
      integer(ip), intent(in) :: n
      !! The order of the Jacobi polynomials.
      real(wp), dimension(:), intent(out) :: w
      !! On return the integration weights.

      w = (2.0_wp*n+1.0_wp)/(n*(n+1.0_wp)) / JacobiP(r, alpha, beta, n)**2
    end subroutine GaussWeigths


    function Vandermonde1D ( n, x )
    !! Initialize the 1D Vandermonde matrix, \(\mathcal{V}_{ij}=P_j(x_i)\) 
    !! for the Legendre-Gauss-Lobatto quadrature points.
      integer(ip), intent(in) :: n
      !! The order of the Jacobi polynomials.
      real(wp), dimension(:), intent(in) :: x
      !! The points \(x\in[-1:1]\) at which to evaluate the Vandermonde matrix. 
      real(wp), dimension(size(x),n+1) :: Vandermonde1D
      !! The return value is the Vandermonde matrix, \(\mathcal{V}_{ij}\)

      integer :: j

      do j=1,n+1
        Vandermonde1D(:,j) = JacobiP(x, 0.0_wp, 0.0_wp, j-1)
      end do
      return
    end function VanderMonde1D


    function GradJacobiP ( x, alpha, beta, n )
    !! Evaluate the derivative of the Jacobi polynomial of type
    !! \((\alpha,\beta)>-1\) at points \(x\) for order \(n\).
      real(wp), dimension(:), intent(in) :: x
      !! The points \(x\in[-1:1]\) at which to evaluate the derivative of the
      !! Jacobi polynomial. 
      real(wp), intent(in) :: alpha
      !! The value of \(\alpha\).
      real(wp), intent(in) :: beta
      !! The value of \(\beta\).
      integer(ip), intent(in) :: n
      !! The order of the Jacobi polynomials.
      real(wp), dimension(size(x)) :: GradJacobiP
      !! The return value is the derivative of the Jacobi polynomial evaluated
      !! at \(x_i\), \(\left .\frac{dP_n^{(\alpha,\beta)}}{dx}
      !! \right |_{x=x_i}\)

      if (n==0) then
        GradJacobiP = 0.0_wp
        return
      end if
      GradJacobiP = sqrt(n*(n+alpha+beta+1))* &
                      JacobiP(x,alpha+1.0_wp,beta+1.0_wp,n-1)
      return
    end function GradJacobiP


    function GradVandermonde1D ( n, x )
    !! Initialize the gradient of the modal basis \(j\) at \(r_i\).
      integer(ip), intent(in) :: n
      !! The order of the Jacobi polynomials.
      real(wp), dimension(:), intent(in) :: x
      !! The points \(x\in[-1:1]\) at which to initialize the gradient of the
      !! modal basis. 
      real(wp), dimension(size(x),n+1) :: GradVandermonde1D
      !! The return value is the derivative of the modal basis,
      !! \(\mathcal{V}_{r,(ij)}\).

      integer :: j

      do j=0,n
        GradVandermonde1D(:,j+1) = GradJacobiP(x, 0.0_wp, 0.0_wp, j)
      end do

      return
    end function GradVandermonde1D


    function Dmatrix1D ( n, x, v )
    !! Initialize the differentiation matrix for a reference element of
    !! order \(n\).
      integer(ip), intent(in) :: n
      !! The order of the Jacobi polynomials.
      real(wp), dimension(n+1), intent(in) :: x
      !! The nodal points \(x\in[-1:1]\) where the differentiation matrix
      !! will provide approximations to the derivative.
      real(wp), dimension(n+1,n+1), intent(in) :: v
      !! The Vandermonde matrix, \(\mathcal{V}_{ij}\).
      real(wp), dimension(n+1,n+1) :: Dmatrix1D
      !! The return value is the differentiation matrix, \(\mathcal{D}_r=
      !! \mathcal{V}_r \mathcal{V}^{-1}\).

      real(wp), dimension(n+1,n+1) :: vr, vrt, vt, drt
      integer(ip), dimension(n+1) :: ipiv
      integer(ip) :: i, j, info

      vt = transpose(v)
      vr = GradVandermonde1D ( n, x )
      vrt = transpose(vr)

      if ( wp == dp ) then
        call dgesv ( n+1, n+1, vt, n+1, ipiv, vrt, n+1, info )
      else if ( wp == qp ) then
        call qgesv ( n+1, n+1, vt, n+1, ipiv, vrt, n+1, info )
      endif

      if (info <0) then
        print*,'Parameter ', i, ' in call to dgesv has illegal value'
        stop
      end if
      if (info >0) then
        print*, 'Matrix in call to dgesv is singular'
        stop
      end if

      Dmatrix1D = transpose(vrt)
    end function Dmatrix1D


    function Lift1D ( n, v )
    !! Initialize the lift matrix, \(L_{ij}\), used to compute surface
    !! integral terms in the Discontinuous Galerkin formulation.
      integer(ip), intent(in) :: n
      !! The order of the Jacobi polynomials.
      real(wp), dimension(n+1,n+1), intent(in) :: v
      !! The Vandermonde matrix, \(\mathcal{V}_{ij}\).
      real(wp), dimension(n+1,2) :: Lift1D
      !! The return value is the lift matrix, \(L= \mathcal{V}
      !! \mathcal{V}^{\mathrm{T}}\mathcal{E}\), where \(\mathcal{E}\) is
      !! a \(n+1\times 2\) array of the form
      !! \[
      !!   \begin{pmatrix}
      !!     1 & 0 & \cdots & 0 & 0 \\
      !!     0 & 0 & \cdots & 0 & 1
      !!   \end{pmatrix}
      !! \]

      real(wp), dimension(n+1,2) :: emat
      integer(ip) :: i, j

      emat = 0.0_wp
      emat(1,1) = 1.0_wp
      emat(n+1,2) = 1.0_wp

      Lift1D = matmul(v,matmul(transpose(v),emat))
      return
    end function Lift1D


    subroutine Jacobian ( n, vx, r, dr, x, rx )
    !! Calculate the Jacobian for transforming derivatives from the
    !! reference element to the physical element.
      integer(ip), intent(in) :: n
      !! The order of the element.
      real(wp), dimension(2), intent(in) :: vx
      !! A 1d real array of size 2 that contains the physical coordinates at
      !! the boundaries of the element.`
      real(wp), dimension(n+1), intent(in) :: r
      !! The node locations within the reference element, \(r_i\).
      real(wp), dimension(:,:), intent(in) :: dr
      !! The derivative matrix for the reference element, \(\mathcal{D}_r\).
      real(wp), dimension(n+1), intent(out) :: x
      !! On output contains the physical coordinates for this physical element,
      !! \(x_i\).
      real(wp), dimension(n+1), intent(out) :: rx
      !! On output contains \(\frac{d r_i}{d x_j}.
      real(wp), dimension(n+1) :: xr
      integer(ip) :: i,j
      
      do i = 1, n+1
        x(i) = vx(1)+0.5_wp*(1.0_wp+r(i))*(vx(2)-vx(1))
      end do
      xr = 0.0_wp
      do i = 1, n+1
        do j = 1, n+1
          xr(i) = xr(i) + dr(i,j) * x(j)
        end do
      end do
      rx = 1.0_wp / xr
    end subroutine Jacobian


    function Inverse ( a )
    !! Helper function that calculates the inverse of a matrix.
    !!
    !! Calls the LAPACK routine GESV.
      real(wp), dimension(:,:), intent(in) :: a
      !! A matrix, \(A\).
      real(wp), dimension(size(a,1),size(a,2)) :: Inverse
      !! On output contains \(A^{-1}\).

      real(wp), dimension(size(a,1),size(a,2)) :: acopy, rhs
      integer(ip), dimension(size(a,1)) :: ipiv
      integer(ip) :: i, n, info

      acopy = a
      n = size(a,1)
      rhs = 0.0_wp
      forall(i=1:n) rhs(i,i) = 1.0_wp
      call dgesv ( n, n, acopy, n, ipiv, rhs, n, info )

      if (info <0) then
        print*,'Parameter ', i, ' in call to dgesv has illegal value'
        stop
      end if
      if (info >0) then
        print*, 'Matrix in call to dgesv is singular'
        stop
      end if

      Inverse = rhs
      return
    end function Inverse


    function Filter1D ( n, nc, s, v, invv )
    !! Create an exponential filter matrix that can be used to filter out
    !! high-frequency noise.
    !!
    !! The filter matrix \(\mathcal{F}\) is defined as \(\mathcal{F}=
    !! \mathcal{V}\Lambda\mathcal{V}^{-1}\) where the diagonal matrix, 
    !! \(\Lambda\) has the entries \(\Lambda_{ii}=\sigma(i-1)\) for
    !! \(i=1,\ldots,n+1\) and the filter function, \(\sigma(i)\) has the form
    !! \[
    !!   \sigma(i) =
    !!     \begin{cases}
    !!         1 & 0\le i\le n_c \\
    !!         e^{-\alpha\left (\frac{i-n_c}{n-n_c}\right )^s} & n_c<i\le n.
    !!     \end{cases}
    !! \]
    !! Here \(\alpha=-\log(\epsilon_M)\), where \(\epsilon_M\) is the machine
    !! precision in working precision, \(n\) is the order of the element,
    !! \(n_c\) is a cutoff, below which the low modes are left untouched and
    !! \(s\) (has to be even) is the order of the filter.
      integer(ip), intent(in) :: n
      !! The order of the element.
      integer(ip), intent(in) :: nc
      !! The cutoff, below which the low modes are left untouched.
      integer(ip), intent(in) :: s
      !! The order of the filter.
      real(wp), dimension(:,:), intent(in) :: v
      !! The Vandermonde matrix, \(\mathcal{V}\).
      real(wp), dimension(:,:), intent(in) :: invv
      !! The inverse of the Vandermonde matric, \(\mathcal{V}^{-1}\).
      real(wp), dimension(n+1,n+1) :: Filter1D
      !! The return value is the filter matrix, \(\mathcal{F}\).
      real(wp), dimension(n+1,n+1) :: tmp

      real(wp) :: alpha
      integer(ip) :: i

      alpha = -log(epsilon(1.0_wp))
      Filter1D = 0.0_wp
      forall(i=1:nc) Filter1D(i,i) = 1.0_wp
      forall(i=nc:n) Filter1D(i+1,i+1) = exp(-alpha*(real(i-nc,wp)/real(n-nc,wp))**s)
!      Filter1D = matmul(Filter1D,invv)
!      Filter1D = matmul(v,Filter1D)

      tmp = matmul(Filter1D,invv)
      Filter1D = matmul(v,tmp)

      return
    end function Filter1D

end submodule DG_implementation
