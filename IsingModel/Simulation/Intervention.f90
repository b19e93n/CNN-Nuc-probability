module common100_XXX
! XXX will be replaced by a number using bash script later, since this Fortran script needs to be ran hundreds of times with different random seed parallelly on computer cluster.
implicit none

public :: initializeMC, computeOffSets, equilibrate, Metropolis, pbc, deltaE, rand

integer, public, parameter :: r8 = selected_real_kind(8)
integer, public, parameter :: i1 = selected_int_kind(1)
integer, public, parameter :: L = 100, N = L*L, R = 10, threshold = 0
! start time and end time is fine tuned to best balance data points with different nucleation probability.
integer, public, parameter :: ntrials = 5, ninter = 50, store_length = 200, start_time = 116, end_time = 168, history = 100
real (kind=r8), public, parameter :: h = 1.17
integer, public, parameter :: up = 1, down = -1
integer (kind = i1), public, dimension(:,:), allocatable :: spin
integer, public, dimension(:), allocatable :: nnx, nny     ! offset vectors
real (kind=r8), public :: E
integer, public :: M, z

contains
	
	subroutine initializeMC(tequil, betaJ, betaH)
		integer, intent(out) :: tequil
		real (kind=r8), intent(out) :: betaJ, betaH
		real (kind=r8) :: T, J, beta, Tc, rnd
		integer, dimension(:), allocatable :: seed
		integer :: i, sizeOfSeed, clock, seed_c
		character(len=1024) :: config, me, pnuc
		character(len=8) :: format
	
		!set random seed to clock time
		call random_seed(size=sizeOfSeed)
		allocate(seed(sizeOfSeed))
		call system_clock(count=clock)
		seed_c = XXX
		seed = seed_c * 100000 + 37*(/ (i - 1, i = 1, sizeOfSeed) /)
		call random_seed(put = seed)
		
	
		if (R == 1) then
			Tc = 2.0/log(1 + sqrt(2.0))
			J = 1.0
		else
			Tc = 4.
			J = 4./(2.0 * z)
		endif
		T = (4./9.) * Tc
		beta = 1./T
		betaJ = beta * J
		betaH = beta * h
		tequil = 50 !duration of run before flip of field
	
		!open record files
		if (L<100) then
			write(config, "(I2,A,I3.3,A)") L, "_config_", seed_c, ".txt"
			write(pnuc, "(I2,A,I3.3,A)") L, "_pnuc_", seed_c, ".txt"
			write(me, "(I2,A,I3.3,A)") L, "_me_", seed_c, ".txt"
		else
			write(config, "(I3,A,I3.3,A)") L, "_config_", seed_c, ".txt"
			write(pnuc, "(I3,A,I3.3,A)") L, "_pnuc_", seed_c, ".txt"
			write(me, "(I3,A,I3.3,A)") L, "_me_", seed_c, ".txt"
		endif
		open(unit = 10, status = 'replace', file = config)
		write (10,*), "L, R, T, h, tequil, ntrials, ninter, seed", L, R, T, h, tequil, ntrials, ninter, seed(0)
		open(unit = 12, status = 'replace', file = me)
		write (12,*), "L, R, T, h, tequil, ntrials, ninter, seed", L, R, T, h, tequil, ntrials, ninter, seed(0)
		open(unit = 14, status = 'replace', file = pnuc)
		write (14,*), "L, R, T, h, tequil, ntrials, ninter, seed", L, R, T, h, tequil, ntrials, ninter, seed(0)

	end subroutine initializeMC
	
	subroutine computeOffSets()
		integer :: x, y, i
	    allocate(spin(L,L))
	    if (R == 1) then
	        z = 2
	    else
	        z = 2*R*(R + 1)
		endif
		allocate(nnx(z))
		allocate(nny(z))
	    nnx = 0
	    nny = 0
		if (R .eq. 1) then
	        nnx(1) = 1
	        nny(1) = 0
	        nnx(2) = 0
	        nny(2) = -1
	    else
	        i = 0
	        do y = -R, R
	            do x = 0, R
	                if (x .ne. 0 .or. y .lt. 0) then
	                    i = i + 1
	                    nnx(i) = x
	                    nny(i) = y
	                end if
	            end do
	        end do
	    end if
	end subroutine computeOffSets

	subroutine equilibrate(betaJ, betaH)
		real (kind=r8), intent(in) :: betaJ, betaH
		real :: rnd
		integer :: x, y, i, r, xn, yn, Heff
		E = 0.
	    spin = 1
	    M = N
	    do y = 1, L
	        do x = 1, L
	            Heff = 0
	            do i = 1, z         ! sum over half of the neighbors of x, y
	                xn = x + nnx(i)
	                yn = y + nny(i)
	                xn = pbc(xn)
	                yn = pbc(yn)
	                Heff = Heff + spin(xn, yn)
	            end do
	            E = E - spin(x,y)*Heff*betaJ - spin(x,y)*betaH
	        end do
	    end do
	end subroutine equilibrate	

	subroutine Metropolis(betaJ, betaH)
		real (kind=r8), intent(in) :: betaJ, betaH
		real (kind=r8) :: prob, de
		integer :: ispin, x, y
	    logical :: flag
	    real :: rnd
	    flag = .false.
	    do ispin = 1, N
			call random_number(rnd)
			x = int(L*rnd) + 1
			call random_number(rnd)
	        y = int(L*rnd) + 1
	        de = deltaE(x,y,betaJ, betaH)
	        if (de > 0.) then
	            call random_number(rnd)
	            prob = exp(-de)      ! de already contains factor of beta
	            if (rnd <= prob) flag = .true.
	        else        ! de <= 0
	            flag = .true.
	        end if
	        if (flag) then
	            spin(x,y) = -spin(x,y)
	            M = M + 2*spin(x,y)
	            E = E + de
	            flag = .false.
	        end if
	    end do
	end subroutine Metropolis
	
	subroutine determineClusterLabels(label, pointer, dir, pb)
	    integer, dimension(:,:), intent(out) :: label
	    integer, dimension(:), intent(out) :: pointer
	    integer, intent(in) :: dir
	    real (kind=r8), intent(in) :: pb
	    integer :: x, y, xn, yn, i, bigLabel, smallLabel, next
	    label = 0
	    pointer = 0
	    next = 1
	    do x = 1, L
	        do y = 1, L
	            if (spin(x,y) .eq. dir) then
	                do i = 1, z
	                    xn = pbc(x+nnx(i))
	                    yn = pbc(y+nny(i))
	                    label(x,y) = locate(pointer,label(x,y))
	                    label(xn,yn) = locate(pointer,label(xn,yn))
	                    if ( (spin(xn,yn) .eq. dir) .and. (rand() <= pb) ) then
	                        if (label(x,y) .eq. 0 .and. label(xn,yn) .eq. 0) then
	                            label(x,y) = next
	                            label(xn,yn) = next
	                            pointer(next) = next
	                            next = next + 1
	                        else if (label(x,y) .eq. 0) then
	                            label(x,y) = label(xn,yn)
	                        else if (label(xn,yn) .eq. 0) then
	                            label(xn,yn) = label(x,y)
	                        else if (label(x,y) .ne. label(xn,yn)) then
	                            bigLabel = max(label(x,y),label(xn,yn))
	                            smallLabel = min(label(x,y),label(xn,yn))
	                            label(x,y) = smallLabel
	                            label(xn,yn) = smallLabel
	                            pointer(bigLabel) = smallLabel
	                        end if
	                    end if
	                end do
	                if (label(x,y) .eq. 0) then
	                    label(x,y) =  next
	                    pointer(next) = next
	                    next = next + 1
	                end if
	            end if
	        end do
	    end do
	end subroutine determineClusterLabels
	
	subroutine determineClusterSizes(label, pointer, dir, nlab)
	    integer, dimension(:,:), intent(out) :: label
	    integer, dimension(:), intent(in) :: pointer
	    integer, intent(in) :: dir
	    integer, dimension(:), intent(out) :: nlab  ! nlab is size of cluster with label k
	    ! determine size of cluster with label k
	    integer :: k, nd     ! cluster label
	    integer :: x, y
	    nlab = 0
	    nd = 0
	     do x = 1, L
	        do y = 1, L
	            if (spin(x,y) .eq. dir) then
	                label(x,y) = locate(pointer,label(x,y))
	                if (label(x,y) .ne. 0) nlab(label(x,y)) = nlab(label(x,y)) + 1
	            end if
	        end do
	    end do
	end subroutine determineClusterSizes

	subroutine determineHistogram(nlab, Hs, sm, km)
	    integer, dimension(:), intent(in) :: nlab
	    integer, dimension(:), intent(out) :: Hs    ! number of clusters of size s
	    integer, intent(inout) :: sm, km
	    integer :: k         ! cluster label
	    integer :: s         ! cluster size
	    Hs = 0
	    sm = 0        ! largest cluster size at this time
	    do k = 1, N/2
	        s = nlab(k)
	        if (s > 0) then
	            Hs(s) = Hs(s) + 1
	        end if
	        if (sm < s) then
	            sm = s
	            km = k
	        end if
	    end do
	end subroutine determineHistogram
	
	subroutine determineCOM(xcom, ycom, sm, km, label) !determine the Center of Mass of the largest droplet
	    integer :: n, x, y
	    real(kind = 8), intent(inout) :: xcom, ycom
	    integer, dimension(:), allocatable :: xclus, yclus
	    integer, intent(in) :: sm, km
	    integer, dimension(:,:), intent(in) :: label
	    real(kind = 8) :: amx, amy, bmx, bmy, ax, ay, bx, by, thetax, thetay, thetamx, thetamy
	    allocate(xclus(sm))
	    allocate(yclus(sm))
	    xcom = 0
	    ycom = 0
	    amx = 0
	    amy = 0
	    bmx = 0
	    bmy = 0
	    n = 1
	    do x = 1, L
	         do y = 1, L
	              if(label(x,y) == km) then
	                    xclus(n) = x
	                    yclus(n) = y
	                    n = n + 1
	              end if
	         end do
	    end do
	    do n = 1, sm
	         thetax = xclus(n) / real(L) * 6.28318
	         thetay = yclus(n) / real(L) * 6.28318
	         ax = sin(thetax)
	         bx = cos(thetax)
	         ay = sin(thetay)
	         by = cos(thetay)
	         amx = amx + ax
	         amy = amy + ay
	         bmx = bmx + bx
	         bmy = bmy + by
	    end do
	    amx = amx / sm
	    amy = amy / sm
	    bmx = bmx / sm
	    bmy = bmy / sm
	    thetamx = atan2(-amx, -bmx) + 3.14159
	    thetamy = atan2(-amy, -bmy) + 3.14159
	    xcom = L * thetamx / 6.28318
	    ycom = L * thetamy / 6.28318    
	end subroutine determineCOM
	
	subroutine determineDistance(xcom, ycom, km, label, maxdist) !determine the distance between the farthest spin in the droplet and COM
	    real(kind = 8), intent(in) :: xcom, ycom
	    real(kind = 8), intent(inout) :: maxdist
	    integer, intent(in) :: km
	    integer, dimension(:,:), intent(in) :: label
	    integer :: x, y
	    real :: dist
	    maxdist = 0
	    dist = 0
	    do x = 1, L
	         do y = 1, L
	              if (label(x,y)==km) then
	                    dist = separation (real(xcom), real(x), real(ycom), real(y))
	                    if (dist >=  maxdist) then 
	                         maxdist=dist
	                    end if
	              end if
	         end do
	     end do
	end subroutine determineDistance
	
	integer function locate(pointer,k)
	    integer, dimension(:), intent(in) :: pointer
	    integer, intent(inout) :: k
	    if (k .ne. 0) then
	        do while (pointer(k) .ne. k)
	            k = pointer(k)
	        end do
	    end if
	    locate = k
	end function locate

	integer function pbc(x)
	    integer, intent(in) :: x
	    if (x .gt. L) then
	        pbc = x - L
	    else if (x .lt. 1) then
	        pbc = x + L
	    else
	        pbc = x
	    end if
	end function pbc
	
	real function separation(x1, x2, y1, y2)
	     real, intent(in):: x1, y1, x2, y2
	     real :: dx, dy
	     if (abs(x1 - x2) .gt. L/2) then
	          dx = L - abs(x1 - x2)
	     else
	          dx = abs(x1 - x2)
	     end if
	     if (abs(y1 - y2) .gt. L/2) then
	          dy = L - abs(y1 - y2)
	     else 
	          dy = abs(y1 - y2)
		end if
	     separation = sqrt(dx ** 2 + dy ** 2)
	end function separation

	real (kind=r8) function deltaE(x,y, betaJ, betaH)
	    integer, intent (in) :: x, y
	    real (kind=r8), intent(in) :: betaJ, betaH
	    integer :: i, r, xn, yn, Heff
	    Heff = 0
	    do i = 1, z
	        do r = 1, -1, -2
	            xn = pbc(x + r*nnx(i))      ! use reflection symmetry
	            yn = pbc(y + r*nny(i))
	            Heff = Heff + spin(xn, yn)
	        end do
	    end do
	    deltaE = 2*spin(x,y)*Heff
	    deltaE = betaJ*deltaE + 2.0*spin(x,y)*betaH
	
	end function deltaE

	real function rand()
	    real :: rnd
	    call random_number(rnd)
	    rand = rnd
	end function rand

end module common100_XXX

program intervention
	use common100_XXX
	implicit none
	real (kind=r8) :: betaJ, betaH, pb, Ndown, rho, xcom, ycom, t_0, t_1
	integer :: dir, tequil, itrial, imcs, tmax, nmcs, success, sm, km, sample_time, sizeofseed
	integer, dimension(:), allocatable :: seed
	integer :: x, y, i, t, i_inter, x_t, y_t
	integer, dimension(store_length) :: M_store
	real (kind=r8), dimension(store_length) :: E_store
	integer(kind = i1), dimension(store_length, L, L) :: spin_store
	integer, dimension(:), allocatable :: pointer, nlab, Hs
	integer, dimension(L,L) :: label
	
	call random_seed(size=sizeOfSeed)
	allocate(seed(sizeOfSeed))
	
	tmax = 1e6
	dir = down
	
	call computeOffSets()
	
    allocate(Hs(N))
    allocate(nlab(N/2))
	allocate(pointer(N))
	
	call initializeMC(tequil, betaJ, betaH)
	
	call cpu_time(t_0)
	do itrial = 1, ntrials
		call equilibrate(betaJ, betaH)
		do imcs = 1, tequil
			call Metropolis(betaJ, betaH)
		end do
		betaH = -betaH
		imcs = 0
		
		do while (M > threshold)
			call Metropolis(betaJ, betaH)
			imcs = imcs + 1
		end do
		betaH = -betaH
		print *, imcs
	end do
	call cpu_time(t_1)
	print *, t_1 - t_0
	
end program intervention