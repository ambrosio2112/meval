! Modelo de Evaluaciones
! Modelo ordinal para las evaluaciones de marcas
!

module modEvaluaciones
  use randomlib
  implicit none

  integer :: N ! numero de individuos (evaluadores)
  integer :: nmarcas ! numero maximo de marcas evaluadas
  integer, allocatable :: marcaevaluada(:,:) ! marcas evaluadas por cada individuo N x nmarcas
  integer :: E ! Maximo de la escala de evaluacion : 1,..,E
  integer :: nbeta ! numero de atributos evaluados
  integer :: nclass ! numero de clases de individuos (e.g. ciudades)

  real(8), allocatable :: ev(:,:,:) ! evaluaciones por atributo N x nmarcas x nbeta
  integer, allocatable :: y(:,:) ! evaluaciones generales N x nmarcas
  integer, allocatable :: class(:) ! clases de los indivudios dim=N
  integer, allocatable :: iclass(:,:) ! posicion de los individuos en una clase :  nclass x max(bclass)
  integer, allocatable :: bclass(:) ! numero de individuos por clase dim=nclass
  integer, allocatable :: ibj(:) ! numero de evaluaciones en la clase para todas las marcas .. nclass
  integer, allocatable :: missing(:,:) ! valores faltantes N x nmarcas

  real(8), allocatable :: beta(:,:) ! coef. nclass x nbeta
  real(8), allocatable :: tau(:) ! varianza coef evaluacion nclass
  real(8), allocatable :: cutoff(:,:) ! cutoffs para las evaluaciones N x E-1
  real(8), allocatable :: sigma(:) ! varianza de los evaluadores dim = N
  real(8), allocatable :: t(:,:) ! latentes evaluadores por marca N x nmarcas
  real(8), allocatable :: z(:,:) ! latentes marcas N x nmarcas

  real(8) :: alpha, lambda, sigmaMH
  real(8), allocatable :: Sbeta(:), mubeta(:) ! varianzas y medias previas para beta
  real(8), allocatable :: cholxx(:,:,:), diagxx(:,:) ! cholesky X'X para cada clase nclass x nbeta x nbeta y nclass x nbeta
  real(8), allocatable :: cholixx(:,:,:), diagixx(:,:) ! cholesky (X'X)^-1 para cada clase nclass x nbeta x nbeta y nclass x nbeta
  logical adjust
contains

  !--------------------------------------------------------------------------------
  subroutine genBetaTauZ(nc)
    integer :: nc ! class
    real(8) :: b(nbeta), a(nbeta), y(ibj(nc)), s(ibj(nc)), mut(ibj(nc))
    real(8) :: rss
    real(8) :: nu
    real(8) :: ez(ibj(nc)), x(ibj(nc),nbeta)
    real(8) :: c(nbeta)
    real(8) :: g(ibj(nc))
    real(8) :: R
    real(8) :: xb(ibj(nc))
    real(8) :: yy,dd,gg
    integer i,j, iic, m

    iic=0
    do i=1,bclass(nc)
       do m=1,nmarcas
	     if (marcaevaluada(iclass(nc,i), m)==1) then
	        iic=iic+1
            x(iic,:) = ev(iclass(nc,i), m, :)
            y(iic) = z(iclass(nc,i), m)
            s(iic) = sigma(iclass(nc,i))
            mut(iic)= t(iclass(nc,i),m)
	      end if
	   end do
    end do

    yy = sum(y(1:iic)*y(1:iic))
    dd = sum(mubeta*mubeta)

    b = matmul( transpose(x(1:iic,:)) , y(1:iic) ) +  mubeta/Sbeta ! b = X'z + D^-1 d
    call solvechol( cholxx(nc,:,:), diagxx(nc,:), b, a) ! a = (X'X + D^-1)^-1 (X'z + D^-1 d)

    beta(nc,:)=rmvnorm(a, cholixx(nc,:,:), diagixx(nc,:) )*sqrt(tau(nc))

    call solvechol( cholixx(nc,:,:), diagixx(nc,:), beta(nc,:), c ) ! c = (X'X + D^-1) beta

    ! tau
    nu = iic/2. + alpha + nbeta/2.
    r = 1./rgamma( nu )
    rss = yy + dd + sum(beta(nc,:)*( c - 2.*b )) ! beta'(X'X + D^-1) beta - 2 beta'(X'z + D^-1 d)
    tau(nc) = r*(lambda + rss/2.)

    !zeta

    !proposal
    g(1:iic) = rnorm( iic )
    g(1:iic) = g(1:iic)/sqrt( 1. + 1./s(1:iic)) + mut(1:iic)/(1. + s(1:iic))
    gg = sum(g(1:iic)*g(1:iic))

    if(.not.adjust) then
       !acceptance ratio (ln)
       xb(1:iic) = matmul(x(1:iic,:),beta(nc,:)) ! xb = X beta

       R = (sum(  xb(1:iic)*(g(1:iic) - y(1:iic)) ) - gg/2. + yy/2. )/tau(nc) - & ! beta'X'g - beta'X'y
            (alpha + iic/2. ) * log( 2.*lambda + yy + dd - sum(a*b))  ! (X'z + D^-1 d)' {(X'X + D^-1)^-1} (X'z + D^-1 d)

       b = matmul( transpose(x(1:iic,:)) , g(1:iic) ) +  mubeta/Sbeta ! X'g + D^-1 d
       call solvechol( cholxx(nc,:,:), diagxx(nc,:), b, a) ! (X'X + D^-1)^-1(X'g + D^-1 d)
       R = R + (alpha + iic/2. )*log( 2.*lambda + gg + dd - sum(a*b)) ! (X'g + D^-1 d)' {(X'X + D^-1)^-1} (X'g + D^-1 d)

       if (log(runif()) < R ) then
	      iic=0
          do i=1,bclass(nc)
		     do m=1,nmarcas
		        if (marcaevaluada(iclass(nc,i),m)==1 ) then
			       iic=iic+1
                   z(iclass(nc,i),m) = g(iic)
			    end if
			 end do
          end do
       end if
    else
	   iic=0
       do i=1,bclass(nc)
	      do m=1,nmarcas
	         if (marcaevaluada(iclass(nc,i),m)==1 ) then
  		        iic=iic+1
                z(iclass(nc,i),m) = g(iic)
		     end if
		  end do
       end do
    end if

  end subroutine genBetaTauZ

  !--------------------------------------------------------------------------------
  subroutine genSigma()
    real(8) :: r(N)
    integer i,m
    r = 1./rgamma(N, alpha + nmarcas*.5 )
    do i=1,N
       sigma(i) = r(i) * lambda
	   do m=1 , nmarcas
	      if (marcaevaluada(i,m)==1) then
 	         sigma(i) = sigma(i) + r(i) * ( (t(i,m) - z(i,m))**2 )/2.
		  end if
	   end do
    end do
  end subroutine genSigma

  !--------------------------------------------------------------------------------
  subroutine genT()
    integer :: i,j
	real(8) :: r

	logical warn
    do i=1,N
       do j=1,nmarcas
	      if (marcaevaluada(i,j)==1) then
		   warn=.false.
           if ( y(i,j)==1 ) then
              r=rrtruncnorm( z(i,j),sqrt(sigma(i)), cutoff(i, y(i,j)) , warn)
           else if ( y(i,j)==E ) then
		      r=rltruncnorm( z(i,j),sqrt(sigma(i)), cutoff(i, y(i,j)-1) ,warn)
           else
              r=rtruncnorm( z(i,j),sqrt(sigma(i)),cutoff(i, y(i,j)-1), cutoff(i, y(i,j)), warn )
           end if
		   if (.not.warn) t(i,j) = r
		  end if
       end do
    end do
  end subroutine genT

  !--------------------------------------------------------------------------------
  subroutine genCutoff (i)
    integer :: i
    real(8) :: g(E-1)
    integer :: j
    real(8) :: R, sH

    sH=sqrt(sigmaMH)

    ! proposal
    g(1)=0
    do j=2,E-2
       g(j) = rtruncnorm( cutoff(i,j), sH, g(j-1), cutoff(i,j+1))
    end do
    g(E-1) = rltruncnorm( cutoff(i,E-1), sH, g(E-2))

    !acceptance ratio
    R=0.
    do j=1,nmarcas
	   if ( marcaevaluada(i,j)==1 ) then
          if(y(i,j)==1) then
             R = R + log( pnorm( g(y(i,j)), z(i,j) )  - 0. ) - &
                 log( pnorm( cutoff(i,y(i,j)), z(i,j)) - 0. )
          else if(y(i,j)==E) then
             R = R + log( 1. - pnorm( g(y(i,j)-1), z(i,j)) ) - &
                 log( 1. - pnorm( cutoff(i,y(i,j)-1) , z(i,j)) )
          else
             R = R + log( pnorm( g(y(i,j)), z(i,j) )  - pnorm( g(y(i,j)-1), z(i,j)) ) - &
                 log( pnorm( cutoff(i,y(i,j)), z(i,j)) - pnorm( cutoff(i,y(i,j)-1) , z(i,j)) )
          end if
	   end if
    end do

    do j=2,E-2
       R = R + log( pnorm( cutoff(i,j+1),cutoff(i,j),sH ) - pnorm( g(j-1),cutoff(i,j),sH ) ) - &
            log( pnorm( g(j+1),g(j),sH ) - pnorm( cutoff(i,j-1),g(j),sH ) )
    end do

    if (log(runif()) < R ) then
       cutoff(i,:) = g
    end if
  end subroutine genCutoff

  !--------------------------------------------------------------------------------
  subroutine genY(i,m)
    integer :: i,m
    integer :: j

    if (t(i,m) > cutoff(i,E-1) ) then
       y(i,m) = E
    else
       do j=E-1,1, -1
          if (t(i,m) > cutoff(i,j)) then
             y(i,m)=j
             return
          end if
       end do
       y(i,m) = 1
    end if
  end subroutine genY


  !--------------------------------------------------------------------------------
  subroutine predY(y, i,m)
    integer :: i,m
    integer :: j
	integer :: y

    if (t(i,m) > cutoff(i,E-1) ) then
       y = E
    else
       do j=E-1,1, -1
          if (t(i,m) > cutoff(i,j)) then
             y=j
             return
          end if
       end do
       y = 1
    end if
  end subroutine predY

  !--------------------------------------------------------------------------------
  subroutine setCholXX(iwarn)
    integer i,j, id
    logical warn,iwarn
    real(8), allocatable :: x(:,:)
    integer mx,cx

	iwarn=.false.
    ! count and position class of individuals
    mx=0
    do i=1,nclass
       cx=count( class==i )
       if (mx<cx) mx=cx
    end do
    allocate ( iclass(nclass, mx) )
    bclass = 0
    do i=1,N
       bclass(class(i))=bclass(class(i))+1
       iclass(class(i), bclass(class(i)) ) = i
    end do

    do i=1,nclass
       allocate( x(bclass(i)*nmarcas,nbeta) )
	   ibj(i)=0
       do j=1,nmarcas
          do id=1,bclass(i)
		     if (marcaevaluada(iclass(i,id),j)==1 ) then
			    ibj(i)=ibj(i)+1
                x(ibj(i),:)=ev(iclass(i,id),j,:)
			 end if
          end do
       end do
       if (ibj(i)>0) then
          cholxx(i,:,:) = matmul( transpose(x(1:ibj(i),:)), x(1:ibj(i),:) )
          do id=1,nbeta
             cholxx(i,id,id) = cholxx(i,id,id) + 1./Sbeta(id)
          end do
	      warn=.false.
          call chol(cholxx(i,:,:), diagxx(i,:), warn)

          call InvertChol(cholxx(i,:,:), diagxx(i,:), cholixx(i,:,:) )
          warn=.false.
          call chol(cholixx(i,:,:), diagixx(i,:),warn)
          if (warn) iwarn=.true.
	   end if
       deallocate(x)
    end do
  end subroutine setCholXX

  !--------------------------------------------------------------------------------
  subroutine setInitials(tau0)
    real(8) :: tau0
    real(8) :: p
    integer j,i

    tau = tau0
    sigma = tau0

    do j=1,E-1
       p = j/dble(E)
       cutoff(:,j) = qnorm(p)
    end do
    cutoff(:,:) = cutoff(:,:) - qnorm(1./dble(E))

    do i=1,N
       do j=1,nmarcas
	      if (marcaevaluada(i,j)==1) then
            if(missing(i,j)==1) then
               y(i,j) = int(E/2.) ! incializa valores faltantes en E/2
            end if

            if(y(i,j)==1) then
               z(i,j) = cutoff(i,1) - 2
            else if(y(i,j)==E) then
               z(i,j) = cutoff(i,E-1) + 2
            else
               z(i,j) = 0.5*(cutoff(i, y(i,j)) + cutoff(i,y(i,j)-1))
            end if
            t(i,j)=z(i,j)
		  end if
       end do
    end do
  end subroutine setInitials

  !--------------------------------------------------------------------------------
  subroutine saveState(it)
    integer it
    integer i,j
    open (unit=10, file="betas.txt", access='APPEND')
    do i=1,nclass
       write(10,*) beta(i,:)
    end do
    close(10)

    open (unit=10, file="tau.txt", access='APPEND')
    do i=1,nclass
       write(10,*) tau(i)
    end do
    close(10)

    open (unit=10, file="sigma.txt", access='APPEND')
    write(10,*) sigma(:)
    close(10)

    open(unit=10,file="teval.txt", access='APPEND')
    do i=1,N
       write(10,*) t(i,:)
    end do
    close(10)

    open(unit=10,file="status.txt",access='APPEND')
    write(10,*) it
    close(10)

  end subroutine saveState

  !--------------------------------------------------------------------------------
  subroutine iterar(niter,nburn,nstep,nadjust, betaiter, cutoffiter, titer, BA )
    integer :: niter, nburn, nstep, nadjust
    real(8), optional :: betaiter( int( (niter-nburn+1)/nstep ), nclass, nbeta)
    real(8), optional :: cutoffiter( int( (niter-nburn+1)/nstep ), N, nmarcas)
    real(8), optional :: titer( int( (niter-nburn+1)/nstep ), N, nmarcas)
    integer, optional  :: BA(E,E)

    integer it,i ,j, is, yaux

    adjust = .true.
    is = 0
    do it=1,niter
       do i=1,nclass
          call genBetaTauZ(i)
       end do

       do i=1,N
          call genCutoff(i)
          do j=1,nmarcas
		     if ( marcaevaluada(i,j)==1) then
               if(missing(i,j)==1) then
                  call genY(i,j)
               end if
               if (present(BA)) then
                 call predY( yaux, i,j)
                 BA( y(i,j), yaux ) = BA( y(i,j), yaux ) + 1
			   end if
			 end if
          end do
       end do
       call genT()
       call genSigma()

       if (it == nadjust) adjust=.false.

       if( it >= nburn .and. mod(it-nburn,nstep)==0) then
          is = is + 1
		  if (present(betaiter)) then
		    if( is <= int( (niter-nburn+1)/nstep )) then
              betaiter(is,:,:) = beta
			  cutoffiter(is,:,:) = cutoff
			  titer(is, :,:) = t
			end if
		  else
            call saveState(is)
		  end if
       end if

    end do

  end subroutine iterar
end module modEvaluaciones


!=====================================================================================
subroutine evaluaciones(niter,nburn,nstep, nadjust, N_, nmarcas_, E_, nbeta_, nclass_, &
                        marcaevaluada_, &
                        ev_, y_, missing_, class_, alpha_,lambda_,sigmaMH_, tau0, &
                        Sbeta_,mubeta_)

!dec$attributes dllexport :: evaluaciones
  use modEvaluaciones
  implicit none
  integer :: niter,nburn,nstep, nadjust
  integer :: N_, nmarcas_, E_, nbeta_, nclass_
  integer :: marcaevaluada_(N_,nmarcas_)
  real(8) :: ev_(N_,nmarcas_,nbeta_)
  integer :: y_(N_,nmarcas_), class_(N_), missing_(N_,nmarcas_)
  real(8) :: alpha_, lambda_, sigmaMH_
  real(8) :: tau0 ! initial varianza for betas
  real(8) :: Sbeta_(nbeta_), mubeta_(nbeta_)


logical::iwarn
  N=N_
  nmarcas=nmarcas_
  E=E_
  nbeta=nbeta_
  nclass=nclass_

  alpha=alpha_
  lambda=lambda_
  sigmaMH=sigmaMH_

  allocate(ev(N,nmarcas,nbeta), y(N,nmarcas), missing(N,nmarcas), class(N), bclass(nclass), &
           beta(nclass,nbeta), marcaevaluada(N,nmarcas), ibj(nclass), &
           tau(nclass), cutoff(N,E-1), sigma(N), t(N,nmarcas), z(N,nmarcas), &
           cholxx(nclass,nbeta,nbeta), diagxx(nclass, nbeta), &
           cholixx(nclass,nbeta,nbeta), diagixx(nclass, nbeta), &
           Sbeta(nbeta), mubeta(nbeta) &
           ) !iclass allocated in setCholXX
  ev=ev_
  marcaevaluada=marcaevaluada_
  y=y_
  missing=missing_
  class=class_
  Sbeta=Sbeta_
  mubeta=mubeta_

  ! prepara las matrices X'X y (X'X)^-1

  call setCholXX(iwarn)
  if (iwarn) then
  open (unit=10, file="warn.txt", access='APPEND')
  write(10,*) "Problemas en XX o XX-1"
  endif
  ! setup initial values: tau,
  call setInitials(tau0)

  call iterar(niter,nburn,nstep, nadjust)

  deallocate(ev, y, missing, class, iclass, bclass, &
           beta, marcaevaluada, ibj, &
           tau, cutoff, sigma, t, z, &
           cholxx, diagxx, &
           cholixx, diagixx, &
           Sbeta, mubeta &
           ) !iclass allocated in setCholXX

end subroutine evaluaciones

