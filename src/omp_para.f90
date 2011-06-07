! This is the general version, for M nodes in a cuboid blob. 

! Extension to general geometry:
! We need to be able to handle a spherical blob, for example.
! Simplifying assumptions:
! (1) The partitioning of the blob is always created by slice boundaries that
!     are perpendicular to the X axis, i.e. parallel to the YZ plane.
! (2) At any x value in every slice the available sites form a convex set in 2D.
!     I.e. for any y value there are either no available sites, or all the
!     sites in the range y1 <= y <= y2 are available.  Similarly for any z.
!     This means that the blob shape is fully described by 2 sets of spans:
!     span(x,y) x=1,NX, y=1,NY => available sites are for span%lo <= z <= span%hi
!     span(z,x) x=1,NX, z=1,NZ => available sites are for span%lo <= y <= span%hi

! Implementing a list of cells:
! In this version the cells in the domain are stored in a list, while the
! occupancy array holds the indices of cells in the list.  When a cell
! leaves the domain (i.e. crosses domain boundary, leaves the paracortex,
! or dies) a gap is created in the list.  The locations of such gaps are stored
! in the gaplist, the total number of such gaps is ngaps.  A cell entering the
! domain is allocated an index from the tail of this list, if ngaps > 0, or
! else it is added to the end of the cell list.

module omp_main_mod
use omp_global
use omp_behaviour
use omp_diffuse
use winsock
!use aviewer 

IMPLICIT NONE

contains

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine rng_initialisation
!integer :: my_seed(8)
integer, allocatable :: zig_seed(:)
integer :: i
!integer :: kpar = 0
integer :: npar, grainsize = 32

!do i = 1,8
!    my_seed(i) = i
!enddo
!call random_seed(size = m)
!write(*,*) 'random_number seed size: ',m
!my_seed(1:2) = seed(1:2)
!call random_seed(put=my_seed(1:m))

npar = Mnodes
allocate(zig_seed(0:npar-1))
do i = 0,npar-1
    zig_seed(i) = seed(1)*seed(2)*(i+1)
enddo
call par_zigset(npar,zig_seed,grainsize)
par_zig_init = .true.
!do i = 1,10
!    write(*,*) i,par_uni(kpar)
!enddo
!stop
end subroutine


!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine omp_initialisation(ok)
logical :: ok
integer :: npr, nth

ok = .true.
if (Mnodes == 1) return
!!DEC$ IF ( DEFINED (_OPENMP) .OR. DEFINED (IBM))
#if defined(OPENMP) || defined(_OPENMP)
write(logmsg,'(a,i2)') 'Requested Mnodes: ',Mnodes
call logger(logmsg)
!close(nflog)
!ok = .false.
!return
npr = omp_get_num_procs()
write(logmsg,'(a,i2)') 'Machine processors: ',npr
call logger(logmsg)

nth = omp_get_max_threads()
write(logmsg,'(a,i2)') 'Max threads available: ',nth
call logger(logmsg)
if (nth < Mnodes) then
    Mnodes = nth
    write(logmsg,'(a,i2)') 'Setting Mnodes = max thread count: ',nth
	call logger(logmsg)
endif

call omp_set_num_threads(Mnodes)
!$omp parallel
nth = omp_get_num_threads()
write(logmsg,*) 'Threads, max: ',nth,omp_get_max_threads()
call logger(logmsg)
!$omp end parallel
#endif
!!DEC$ END IF
!write(*,*) 'Threads: ',nth
!if (nth > npr) then
!    write(*,*) 'Number of threads exceeds CPUs'
!    stop
!endif

!stop

!call set_affinity
!write(*,*) 'set affinity mappings'

!if (track_DCvisits) then
!    write(*,*) 'To track DC visits %DClist(:) must be allocated for cells'
!    stop
!endif
call logger('did omp_initialisation')
end subroutine

!-----------------------------------------------------------------------------------------
! This doesn't seem to help
!-----------------------------------------------------------------------------------------
!subroutine set_affinity
!integer :: tmax, tnum, nproc, ncores
!integer (kind=kmp_affinity_mask_kind) :: mask
!integer :: mask

!tmax = omp_get_max_threads()
!tnum = omp_get_thread_num()
!nproc = omp_get_num_procs()
!ncores = nproc / 2
!if (tmax /= Mnodes) stop

!!!call kmp_create_affinity_mask(mask)

!do i = tnum % ncores; i < tmax; i += ncores) {
!!!do i = 0, tmax-1
!!!    if (kmp_set_affinity_mask_proc(i, mask) /= 0) then
!!!        write(*,*) 'kmp_set_affinity_mask_proc failed for i = ',i
!!!        stop
!!!    endif
!!!enddo

!!!if (kmp_set_affinity(mask) /= 0) then
!!!    write(*,*) 'set_affinity: error: '
!!!    stop
!!!endif
!end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine array_initialisation(ok)
logical :: ok
integer :: x,y,z,k
!integer :: xdc,ydc,zdc,xmin,xmax,ymin,ymax,zmin,zmax
integer :: MAXX, z1, z2
integer :: cog_size !, noncog_size
real :: d, rr(3)
!type(cell_type) :: cell
type(cog_type) :: cog
!type(noncog_type) :: noncog

ok = .false.
call rng_initialisation

cog_size = (sizeof(cog) + 2)/4
!noncog_size = (sizeof(noncog) + 2)/4
!write(*,*) 'noncog_size,cog_size: ',NCSIZE,CSIZE    !noncog_size,cog_size

nsteps_per_min = 1.0/DELTA_T
NY = NX
NZ = NX
if (IN_VITRO) then
	NZ = NZ_IN_VITRO
endif
!max_big_nlist = NX*NY*NZ      ! we need to allow for gaps
ngaps = 0
max_ngaps = 5*NY*NZ
ID_offset = BIG_INT/Mnodes
nlist = 0
MAX_COG = 0.5*NX*NY*NZ
if (use_DC) then
    MAX_DC = 5*(NX*NX*NX)/TC_TO_DC
else
    MAX_DC = 10
endif
allocate(zoffset(0:2*Mnodes))
allocate(zdomain(NZ))
allocate(xoffset(0:2*Mnodes))
allocate(xdomain(NX))
allocate(zrange2D(NX,NY,2))
zrange2D = 0
x0 = (NX + 1.0)/2.                ! global value
y0 = (NY + 1.0)/2.
if (IN_VITRO) then
	z0 = 1
else
	z0 = (NZ + 1.0)/2.
ENDIF
if (IN_VITRO) then
	z0 = 1
endif
if (use_blob) then
    globalvar%Radius = BLOB_RADIUS     ! starting value
    if (2*globalvar%Radius > NX-2) then
        write(logmsg,'(a,i6,f7.1)') 'ERROR: NX too small for RADIUS: ',NX,globalvar%Radius
		call logger(logmsg)
        if (IN_VITRO) then
			logmsg = 'Fix: Reduce IV_WELL_DIAMETER or increase NX'
			call logger(logmsg)
		else
			logmsg = 'Fix: Reduce BLOB_RADIUS or increase NX'
			call logger(logmsg)
		endif
        return
    endif
else
    globalvar%Radius = 0
endif

max_nlist = 1.5*NX*NY*NZ

allocate(occupancy(NX,NY,NZ))
allocate(cellist(max_nlist))
allocate(gaplist(max_ngaps))
allocate(nz_sites(NZ))
allocate(nz_totsites(NZ))
allocate(nz_cells(NZ))
allocate(nz_excess(NZ))
!allocate(globalvar%Inflow(0:Mnodes-1))
!allocate(globalvar%Outflow(0:Mnodes-1))
allocate(cognate_list(MAX_COG))

do k = 1,max_nlist
	nullify(cellist(k)%cptr)
enddo
if (use_DC) then

    allocate(DClist(MAX_DC))
    allocate(DCdeadlist(MAX_DC))
    allocate(DCvisits(0:MAX_DC+1))
    z1 = -DC_RADIUS
    z2 = DC_RADIUS
    if (IN_VITRO) then
		z1 = 0
		z2 = NZ - 1
	endif
    k = 0
    do x = -DC_RADIUS,DC_RADIUS
	    do y = -DC_RADIUS,DC_RADIUS
		    do z = z1, z2
		        rr = (/x,y,z/)
		        d = norm(rr)
			    if (d <= DC_RADIUS) then
				    k = k+1
			    endif
		    enddo
	    enddo
    enddo
    k = k - NDCsites
    nbindmax = min(k,MAX_TC_BIND)
    nbind1 = ABIND1*nbindmax
    nbind2 = ABIND2*nbindmax
!    write(*,*) 'Space for T cells/DC: ',k
!    write(*,*) 'nbind1,nbind2,nbindmax: ',nbind1,nbind2,nbindmax
    if (.not.IN_VITRO) then
	    DCoffset(:,1) = (/ 0,0,0 /)
	    if (NDCsites > 1 .and. NDCsites <= 7) then
	        DCoffset(:,2:NDCsites) = neumann(:,1:NDCsites-1)
	    endif
	else
		! need to create DCoffset() for IN_VITRO case
		if (NDCsites == 9) then
			k = 0
			do x = -1,1
				do y = -1,1
					k = k+1
					DCoffset(:,k) = (/ x, y, 0 /)
				enddo
			enddo
		elseif (NDCsites >= 5 .and. NDCsites <=8) then
			DCoffset(:,1) = (/ 0, 0, 0 /)
			DCoffset(:,2) = (/-1, 0, 0 /)
			DCoffset(:,3) = (/ 1, 0, 0 /)
			DCoffset(:,4) = (/ 0,-1, 0 /)
			DCoffset(:,5) = (/ 0, 1, 0 /)
			if (NDCsites > 5) DCoffset(:,6) = (/ 0, 0, 1 /)
			if (NDCsites > 6) DCoffset(:,7) = (/-1, 0, 1 /)
			if (NDCsites > 7) DCoffset(:,8) = (/ 1, 0, 1 /)
		endif
	endif
    ndeadDC = 0
	DClist(:)%nbound = 0		! moved from initial_binding
	DClist(:)%ncogbound = 0
endif

call make_reldir

nz_excess = 0
Centre = (/x0,y0,z0/)   ! now, actually the global centre (units = grids)
!write(*,*) 'Centre: ',Centre
ncogseed = 0
lastcogID = 0
lastID = 0
k_nonrandom = 0
lastNTcells = 0
nadd_sites = 0
lastbalancetime = 0
localres%dN_EffCogTC = 0
localres%dN_EffCogTCGen = 0
localres%N_EffCogTC = 0
localres%N_EffCogTCGen = 0
totalres%dN_EffCogTC = 0
totalres%dN_EffCogTCGen = 0
totalres%N_EffCogTC = 0
totalres%N_EffCogTCGen = 0
totalres%N_dead = 0
totalres%dN_dead = 0

if (optionA == 1 .and. .not.use_cytokines) then
    write(logmsg,*) 'optionA == 1 => use_cytokines must be true'
    call logger(logmsg)
    stop
endif
if (optionC == 1 .and. .not.use_cytokines) then
    write(logmsg,*) 'optionC == 1 => use_cytokines must be true'
    call logger(logmsg)
    stop
endif
if (track_DCvisits) then
    if (use_DCflux) then
        write(logmsg,*) 'track_DCvisits => use_DCflux must be false'
	    call logger(logmsg)
        stop
    endif
    if (use_cytokines .or. use_diffusion) then
        write(logmsg,*) 'track_DCvisits => use_cytokines, use_diffusion must be false'
	    call logger(logmsg)
        stop
    endif
    nvisits = 0
    nrevisits = 0
    DCvisits = 0
endif

if (evaluate_residence_time) then
    allocate(Tres_dist(int(days*24)))
    Tres_dist = 0
endif

if (use_cytokines) then
    allocate(cytp(NX,NY,NZ,N_CYT))
endif
if (use_diffusion) then
    if (.not.use_cytokines) then
        write(logmsg,*) 'Cannot use_diffusion without use_cytokines'
	    call logger(logmsg)
        stop
    endif
    allocate(xminmax(NY,NZ,2))
    allocate(inblob(NX,NY,NZ))
    MAXX = 1.5*PI*(NX/2)**3/(2*Mnodes)
    allocate(sitelist(MAXX,3,8))
    allocate(neighbours(0:6,MAXX,8))
endif

ok = .true.

end subroutine

!--------------------------------------------------------------------------------
! Generates the arrays wz(), zoffset() and zdomain().
! The domains (slices) are numbered 0,...,2*Mnodes-1
! wz(k) = width of the slice for kth domain
! zoffset(k) = offset of kth domain occupancy array in the occupancy array.
! zdomain(x) = domain that global z lies in.
! The kth domain (slice) extends from z = zoffset(k)+1 to z = zoffset(k+1)
! The idea is to set the domain boundaries such that each domain has roughly the
! same number of available sites.
! This is the initial split, which will continue to be OK if:
! not using a blob, or Mnodes <= 2
! If IN_VITRO, cells lie in the x-y plane (z=1), and the split is based on x value,
! => xoffset(k), xdomain(k)
!--------------------------------------------------------------------------------
subroutine make_split
integer :: k, wsum, kdomain, nsum, Ntot, N, last, x, y, z
integer, allocatable :: scount(:)
integer, allocatable :: wz(:), ztotal(:)
integer :: Mslices
real :: dNT, diff1, diff2
!integer, save :: lastNTcells = 0
!type(cell_type), pointer :: cell
logical :: show = .false.

if (IN_VITRO) then
	xdomain = -1
	if (Mnodes == 1) then
		xdomain = 0
		Mslices = 1
		return
	endif
	Mslices = 2*Mnodes
	! Divide up a disk into Mslices equal pieces
	! First find total number of sites inside
	Ntot = 0
	z = 1
	do x = 1,NX
		do y = 1,NY
			if (occupancy(x,y,z)%indx(1) /= OUTSIDE_TAG) Ntot = Ntot + 1
		enddo
	enddo
	k = 0
	nsum = 0
	do x = 1,NX
		do y = 1,NY
			if (occupancy(x,y,z)%indx(1) /= OUTSIDE_TAG) then
				if (k == 0) then
					xoffset(0) = x - 1
					k = 1
				endif
				nsum = nsum + 1
			endif
		enddo
		if (k > 0 .and. nsum >= (k*Ntot)/Mslices) then
			xoffset(k) = x
			k = k+1
			if (k > Mslices) exit
		endif
	enddo
	do k = 0,Mslices-1
		do x = xoffset(k)+1,xoffset(k+1)
			xdomain(x) = k
		enddo
	enddo
	write(logmsg,'(a,12i4)') 'xoffset: ',xoffset
	call logger(logmsg)
	write(nflog,'(a,12i4)') 'xoffset: ',xoffset
	write(nflog,*) 'xdomain:'
	write(nflog,'(10i5)') xdomain(1:NX)
	return
endif

!write(*,*) 'make_split: Mnodes: ',Mnodes,use_blob
if (Mnodes == 1) then
    Mslices = 1
    zdomain = 0
    return
endif
Mslices = 2*Mnodes
allocate(wz(0:Mslices))
allocate(ztotal(0:Mslices))
allocate(scount(NX))
if (use_blob) then
    dNT = abs(globalvar%NTcells - lastNTcells)/real(lastNTcells+1)
    if (dNT < 0.03) then
!       write(*,*) 'debugging make_split: ',globalvar%NTcells,lastNTcells,dNT
        return
    endif
    lastNTcells = globalvar%NTcells
    if (show) write(*,*) 'make_split: dNT: ',dNT
    nsum = 0
    do z = 1,NZ
        k = 0
        do y = 1,NY
            do x = 1,NX
                if (occupancy(x,y,z)%indx(1) /= OUTSIDE_TAG) then
                    k = k + 1
                endif
            enddo
        enddo
        scount(z) = k
        nsum = nsum + scount(z)
    enddo
    Ntot = nsum
    N = Ntot/Mslices
    nsum = 0
    last = 0
    k = 0
    do z = 1,NZ
        nsum = nsum + scount(z)
        if (nsum >= (k+1)*N) then
            diff1 = nsum - (k+1)*N
            diff2 = diff1 - scount(z)
            if (abs(diff1) < abs(diff2)) then
                wz(k) = z - last
                last = z
            else
                wz(k) = z - last - 1
                last = z - 1
            endif
            k = k+1
            if (k == Mslices-1) exit
        endif
    enddo
    wz(Mslices-1) = NZ - last
    if (show) then
        write(*,*) 'Ntot, N: ',Ntot,N
        write(*,'(10i6)') scount
    endif
else
    wz = NZ/Mslices
    wsum = 0
    do k = 0,Mslices-1
        wsum = wsum + wz(k)
    enddo
    do k = 0,Mslices-1
        if (wsum < NZ) then
            wz(k) = wz(k) + 1
            wsum = wsum + 1
        endif
    enddo
endif
zoffset(0) = 0
do k = 1,Mslices-1
    zoffset(k) = zoffset(k-1) + wz(k-1)
enddo
zoffset(Mslices) = NZ
z = 0
do kdomain = 0,Mslices-1
    do k = 1,wz(kdomain)
        z = z+1
        zdomain(z) = kdomain      ! = kpar with two sweeps
    enddo
enddo
if (show) then
    write(*,*) 'zoffset: ',zoffset
    write(*,*) 'wz:      ',wz
    write(*,*) 'zdomain: '
    write(*,'(10i4)') zdomain
endif
ztotal = 0
do k = 0,2*Mnodes-1
    do z = zoffset(k)+1,zoffset(k+1)
        ztotal(k) = ztotal(k) + scount(z)
    enddo
    if (show) write(*,*) k,ztotal(k)
enddo
deallocate(wz)
deallocate(ztotal)
deallocate(scount)
end subroutine

!--------------------------------------------------------------------------------
! Makes an approximate count of the number of sites of the spherical blob that
! are in the xth slice.  Uses the area of the slice.
! The blob centre is at (x0,y0,z0), and the blob radius is R = globalvar%Radius
!--------------------------------------------------------------------------------
integer function slice_count(x)
integer :: x
real :: r2

r2 = globalvar%Radius**2 - (x-x0)**2
if (r2 < 0) then
    slice_count = 0
else
    slice_count = PI*r2
endif
end function

!--------------------------------------------------------------------------------
!-------------------------------------------------------------------------------- 
subroutine initialise_vascularity
if (VEGF_MODEL == 2) then	! Not used
	VEGF_beta = 0
	VEGF_baserate = 0
    VEGF_decayrate = 0.0012
    vasc_maxrate = 0.001
    globalvar%VEGF = 0
    vasc_beta = 0.00001
    vasc_decayrate = 0.001
else	! VEGF_MODEL = 1
!	VEGF_beta = 4.0e-8
	VEGF_baserate = VEGF_beta*globalvar%NTcells0
!    VEGF_decayrate = 0.002		! delta_G
!    vasc_maxrate = 0.0006		! alpha_V
    globalvar%VEGF = VEGF_baserate/VEGF_decayrate    ! steady-state VEGF level M_G0
    c_vegf_0 = globalvar%VEGF/globalvar%NTcells0	! taking K_V = 1.0
!    vasc_beta = 1.5				! beta_V
    vasc_decayrate = vasc_maxrate*hill(c_vegf_0,vasc_beta*c_vegf_0,vasc_n)	! delta_V
!    write(*,*) 'Vascularity parameters:'
!    write(*,*) 'alpha_G,beta_G,delta_G: ',VEGF_alpha, VEGF_beta, VEGF_decayrate
!    write(*,*) 'alpha_V,beta_V,delta_V: ',vasc_maxrate,vasc_beta,vasc_decayrate
!    write(*,*) 'c_vegf_0,VEGF0,VEGF_baserate: ',c_vegf_0,globalvar%VEGF,VEGF_baserate
endif
globalvar%vascularity = 1.00
!write(*,*) 'VEGF_MODEL, VEGF_baserate: ',VEGF_MODEL, VEGF_baserate
end subroutine

!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
subroutine motility_calibration
integer :: ic,id,k,ds(3),isub
integer :: NTcells, nvar0, nvar, ntime2
integer :: ns = 10000
integer :: npaths = 20
integer :: npos = 301
integer :: nbeta = 30, nrho = 50
integer :: ibeta, irho, kpath
real :: dt, time1 = 5, time2 = 15   ! minutes (used for ICB paper results)
!real :: dt, time1 = 5, time2 = 81   ! minutes  (for writing out cell paths)
real :: tagrad = 10
real :: dbeta, drho
real :: betamin = 0.25, betamax = 0.9		! 0.25 - 0.90 for Model_N, 0.15 - 0.80 for Model_M
real :: rhomin = 0.20, rhomax = 0.85			! 0.20 - 0.85 for Model_N, 0.20 - 0.85 for Model_M
!real :: r2,rsum,r2sum,r2sump,r2max
real :: Cm,speed,ssum,d
!real :: rxy,rxysum,speedxy
!integer, allocatable :: site0(:,:),selectcell(:),selectsite(:,:),path(:,:,:)
integer, allocatable :: tagid(:), tagseq(:), tagsite(:,:,:), pathcell(:)
integer, allocatable :: prevsite(:,:)   ! for mean speed computation
real, allocatable :: Cm_array(:,:), S_array(:,:)
!real, allocatable :: r2mean(:), pathlen(:), totdisp(:)
type(cell_type) :: cell
logical :: ok

    write(*,*) 'motility_calibration'

    NTcells = NX*NY*NZ
    if (motility_save_paths) then
!        allocate(path(3,npaths,0:nvar))
        allocate(pathcell(npaths))
    endif
!    allocate(site0(3,NTcells))
!    allocate(selectcell(ntmax))
!    allocate(selectsite(3,ntmax))
!    allocate(r2mean(ntmax))
!    allocate(pathlen(ntmax))
!    allocate(totdisp(ntmax))

write(*,*) 'rho,beta: ',rho,beta

ntime2 = time2
nvar0 = time1
nvar = time2    ! number of minutes to simulate
dt = DELTA_T*nsteps_per_min
TagRadius = tagrad

if (motility_param_range) then
	drho = (rhomax-rhomin)/(nrho-1)
	dbeta = (betamax-betamin)/(nbeta-1)
elseif (n_multiple_runs > 1) then
	nbeta = n_multiple_runs
	nrho = 1
	betamin = beta
	rhomin = rho
	drho = 0
	dbeta = 0
else
	nbeta = 1
	nrho = 1
	betamin = beta
	rhomin = rho
	drho = 0
	dbeta = 0
endif

allocate(Cm_array(nrho,nbeta))
allocate(S_array(nrho,nbeta))

write(*,*) 'nbeta,nrho: ',nbeta,nrho
do ibeta = 1,nbeta
    do irho = 1,nrho
	    rho = rhomin + (irho-1)*drho
	    beta = betamin + (ibeta-1)*dbeta
	    write(*,'(a,2i3,2f6.2)') ' beta, rho: ',ibeta,irho,BETA,RHO
	    call compute_dirprobs
    	    call place_cells(ok)
    	    if (.not.ok) stop
            call make_split
            if (nlist > 0) then
    	        write(*,*) 'make tag list: NTcells,nlist,ntagged: ',NTcells,nlist,ntagged

                allocate(tagseq(NTcells))
                allocate(tagid(ntagged))
                allocate(tagsite(3,ntagged,0:nvar))
                tagseq = 0
                k = 0
    	        kpath = 0
                do ic = 1,nlist
                    if (cellist(ic)%ctype == TAGGED_CELL) then
                        id = cellist(ic)%ID
                        k = k+1
                        tagid(k) = id
                        tagseq(id) = k
                        tagsite(:,k,0) = cellist(ic)%site
						if (motility_save_paths) then
							if (kpath < npaths) then
								kpath = kpath + 1
								pathcell(kpath) = ic
							endif
						endif
                    endif
                enddo
            endif

	    if (checking > 0) call checker
	    if (checking > 0) call checker
!        do ic = 1,nlist
!            if (cellist(ic)%ID == IDtest) then
!                write(*,*) ic,IDtest,cellist(ic)%site
!            endif
!        enddo

        ns = min(ns,nlist)
        allocate(prevsite(3,ns))
        do ic = 1,ns
            prevsite(:,ic) = cellist(ic)%site
        enddo
        ssum = 0

        !
        ! Now we are ready to run the simulation
        !
        if (motility_save_paths) then
            open(nfpath,file='path.out',status='replace')
            write(nfpath,'(i3,a)') npaths,' paths'
        endif

        do istep = 1,nvar
            do isub = 1,nsteps_per_min
                call mover(ok)
                if (.not.ok) stop
                if (.not.IN_VITRO) then
	                call squeezer(.false.)
	            endif
                do ic = 1,ns
                    ds = cellist(ic)%site - prevsite(:,ic)
                    prevsite(:,ic) = cellist(ic)%site
                    d = sqrt(real(ds(1)*ds(1) + ds(2)*ds(2) + ds(3)*ds(3)))
                    ssum = ssum + d*DELTA_X/DELTA_T
                enddo
                if (motility_save_paths) then
                    k = (istep-1)*nsteps_per_min + isub
                    if (k >= nvar0*nsteps_per_min .and. k < nvar0*nsteps_per_min + npos) then
                        write(nfpath,'(160i4)') (cellist(pathcell(kpath))%site(1:2),kpath=1,npaths)
                    endif
                endif
            enddo
            write(*,*) 'speed: ',ssum/(ns*nsteps_per_min*istep)

            do ic = 1,nlist
                cell = cellist(ic)
                if (cell%ctype == TAGGED_CELL) then
                    id = cell%ID
                    k = tagseq(id)
                    tagsite(:,k,istep) = cell%site
!                   if (k < 20) write(*,*) 'site: ',cell%ID,cell%site
                endif
            enddo
        enddo

        call compute_Cm(tagsite,ntagged,nvar0,nvar,dt,Cm)
        speed = ssum/(ns*nvar*nsteps_per_min)
        write(*,'(a,2f8.2)') 'speed, Cm: ',speed,Cm
        if (allocated(tagid))   deallocate(tagid)
        if (allocated(tagseq))  deallocate(tagseq)
        if (allocated(tagsite)) deallocate(tagsite)
	enddo
enddo
if (allocated(pathcell)) deallocate(pathcell)
deallocate(Cm_array)
deallocate(S_array)
if (motility_save_paths) then
    close(nfpath)
endif

end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine diffusion_calibration
!!!use ifport
integer :: x, y, z, istep
integer, parameter :: Ntimes = 50   !500
real(DP) :: t1, t2, Cnumeric(Ntimes),Canalytic(Ntimes)
logical :: ok

write(*,*) 'diffusion_calibration: '
call init_cytokine

x = NX/4
call analytical_soln(x,Canalytic,Ntimes)
call place_cells(ok)
if (.not.ok) stop
call make_split

!t1 = timef()
call cpu_time(t1)
do x = 1,NX
    do y = 1,NY
        do z = 1,NZ
            occupancy(x,y,z)%indx = 1
        enddo
    enddo
enddo

do istep = 1,Ntimes
    call diffuser
    write(*,*) istep,cyt(NX/4,NY/2,NZ/2,1)
    Cnumeric(istep) = cyt(NX/4,NY/2,NZ/2,1)
enddo
!t2 = timef()
call cpu_time(t2)
write(*,'(a,f10.2)') 'Time: ',t2-t1
write(nfout,*) 'NX, x, NDIFFSTEPS: ',NX,x,NDIFFSTEPS
do istep = 1,Ntimes
    write(nfout,'(i6,f8.2,2f8.4)') istep,istep*DELTA_T,Cnumeric(istep),Canalytic(istep)
enddo

call wrapup
stop

end subroutine

!-----------------------------------------------------------------------------------------
! A(n) = 2*Int[0,L](C0(x).cos(n.pi.x/L).dx)
! with C0(x) = 1  x < L/2
!            = 0  x > L/2
! L = NX*DX
! Site i corresponds to x = (i - 1/2)DX
! With n = 2m-1, the integral for A gives for m = 1,2,3,...:
! B(m) = A(n) = (-1)^m.2/(pi*n), other A(n) = 0
! and B(0) = 0.5
! and solution is the series:
! C(x,t) = B(0) + Sum[m=1,...] B(m).cos(n.pi.x/L).exp(-K.t.(n*pi/L)^2)
!-----------------------------------------------------------------------------------------
subroutine analytical_soln(xsite,C,Ntimes)
integer :: xsite, Ntimes
real(DP) :: C(Ntimes)
integer, parameter :: Nterms = 100
real(DP) :: DX, DT, x, t, csum, bsum, L, xL, Kdiff, dC, tfac, Kfac, B(0:Nterms)
integer :: n, m, k

Kdiff = K_diff(1)
!Kdiff = 2.0e-12         ! m^2/s
write(*,*) 'DELTA_T: ',DELTA_T, DELTA_X, PI, Kdiff
DX = DELTA_X*1.0e-6     ! m
DT = DELTA_T*60         ! s
L = NX*DX
B(0) = 0.5d0
bsum = 0
do m = 1,Nterms
    n = 2*m-1
    B(m) = ((-1)**(m+1))*2/(n*PI)
!    B1 = (2/(n*PI))*sin(n*PI/2)
!    write(*,*) m,n,B(m)-B1
    bsum = bsum + B(m)
enddo

Kfac = Kdiff*PI*PI/(L*L)
write(*,*) 'Bsum: ',B(0),bsum
x = (xsite-0.5)*DX
xL = PI*x/L
do k = 1,Ntimes
    t = k*DT
    csum = B(0)
    do m = 1,Nterms
        n = 2*m-1
        tfac = exp(-Kfac*n*n*t)
        dC = B(m)*cos(n*xL)*tfac
        csum = csum + dC
!        write(*,'(2i4,3e12.4)') m,n,tfac,dC,csum
    enddo
    C(k) = csum
enddo
write(*,'(10f7.4)') C
!write(nfout,*) 'Analytical solution'
!write(nfout,'(10f7.4)') C
end subroutine

!-----------------------------------------------------------------------------------------
! Processes T cell binding and unbinding from DC.
! T cells are allowed to bind to MAX_DC_BIND DC, (MAX_DC_BIND = 1 or 2)
! Note that when a cell is able to bind (less than MAX_DC_BIND bindings) the time of the
! most recent unbinding is stored in %unbindtime(2).
! Note that when there is a single bound DC it is always in slot 1.
! Any DC that has become not "alive" is unbound.
! We do not record a list of T cells bound to each DC, therefore it is not possible to
! handle unbinding by traversing the DC list.
! NOTE: Check that non-cognate cells can bind only one DC
!-----------------------------------------------------------------------------------------
subroutine binder(ok)
logical :: ok
integer :: kcell, nbnd, i, k, site(3), ctype, stage, region, neardc(0:DCDIM-1), idc
integer :: nadd, nsub, nu, nb, nc, nbt, nbc, nt
real :: tnow, bindtime, pMHC, t_travel
logical :: unbound, cognate, bound
type(cell_type), pointer :: cell
type(cog_type), pointer :: cog_ptr
integer :: kpar = 0

nu = 0
nb = 0
nc = 0
nbt = 0
nbc = 0
nt = 0
nadd = 0
nsub = 0
tnow = istep*DELTA_T
do kcell = 1,nlist
    cell => cellist(kcell)
    if (cell%ID == 0) cycle
    nt = nt+1
    if (evaluate_residence_time .or. track_DCvisits) then
        cognate = .false.
    else
        cog_ptr => cell%cptr
        cognate = (associated(cog_ptr))
    endif
    if (cognate) then
		nc = nc+1
		call get_region(cog_ptr,region)
		if (region /= LYMPHNODE) cycle
	endif

    ! Handle unbinding first
    nbnd = 0
    unbound = .false.
    do k = 1,MAX_DC_BIND
        idc = cell%DCbound(k)
        if (idc == 0) cycle
        nbnd = nbnd + 1
        if (tnow >= cell%unbindtime(k) .or. .not.DClist(idc)%capable) then   ! unbind
            DClist(idc)%nbound = DClist(idc)%nbound - 1
            if (DClist(idc)%nbound < 0) then
                write(logmsg,'(a,5i6,3f8.2)') 'Error: binder: DClist(idc)%nbound < 0: ', &
					kcell,nbnd,idc,cell%DCbound,cell%unbindtime,tnow
                call logger(logmsg)
                ok = .false.
                return
            endif
            if (cognate) then
                DClist(idc)%ncogbound = DClist(idc)%ncogbound - 1
                if (DClist(idc)%ncogbound < 0) then
                    write(logmsg,'(a,5i6,3f8.2)') 'Error: binder: DClist(idc)%ncogbound < 0: ', &
						kcell,nbnd,idc,cell%DCbound,cell%unbindtime,tnow
	                call logger(logmsg)
		            ok = .false.
			        return
                endif
            endif
            cell%DCbound(k) = 0
            nbnd = nbnd - 1
            unbound = .true.
!           nsub = nsub + 1
        endif
     enddo
     if (unbound) then  ! ensure that a remaining binding is stored in slot 1
        if (cell%DCbound(1) == 0 .and. cell%DCbound(2) /= 0) then
            cell%DCbound(1) = cell%DCbound(2)
            cell%DCbound(2) = 0
            cell%unbindtime(1) = cell%unbindtime(2)
        endif
        cell%unbindtime(2) = tnow
        if (cognate) nu = nu + 1
    endif

    ! Handle new binding
    bound = .false.
    if (nbnd < MAX_DC_BIND .and. tnow >= cell%unbindtime(2) + DC_BIND_DELAY) then
        if (cognate) then
            if (.not.bindable(cog_ptr)) cycle
!            stage = get_stage(cog_ptr)
			call get_stage(cog_ptr,stage,region)
            if (stage == NAIVE) then    ! on first binding a cognate cell is ready for next stage
                cog_ptr%stagetime = tnow
            endif
        endif
        site = cell%site
        neardc = occupancy(site(1),site(2),site(3))%DC      ! list of DC near this site
        if (neardc(0) /= 0) then
            do k = 1,neardc(0)
                idc = neardc(k)
                if (.not.DClist(idc)%capable) cycle
                if (idc == cell%DCbound(1)) cycle      ! valid for MAX_DC_BIND <= 2
                if (DClist(idc)%ncogbound > MAX_COG_BIND) then      ! ERROR CHECKING
                    write(logmsg,'(a,i6,i16)') 'Error: binder: ncogbound: ',idc,DClist(idc)%ncogbound
			        call logger(logmsg)
		            ok = .false.
	                return
                endif
                if (cognate .and. DClist(idc)%ncogbound == MAX_COG_BIND) cycle
                if (bindDC(idc,kpar)) then
!					if (compute_travel_time .and. cognate .and. nbnd == 1) then
					if (tnow > 2*60 .and. compute_travel_time .and. nbnd == 1) then
						t_travel = tnow - cell%unbindtime(2)
						ntravel = ntravel + 1
						i = t_travel + 0.5
						i = max(i,1)
						i = min(i,N_TRAVEL_DIST)
						travel_dist(k_travel_cog,k_travel_dc,i) = travel_dist(k_travel_cog,k_travel_dc,i) + 1
					endif
                    nbnd = nbnd + 1
                    if (cell%DCbound(nbnd) /= 0) then
                        write(logmsg,'(a,i6)') 'Error: binder: DCbound(nbnd) /= 0: ',cell%DCbound(nbnd)
		                call logger(logmsg)
				        ok = .false.
						return
                    endif
                    bound = .true.
                    cell%DCbound(nbnd) = idc
                    ctype = cell%ctype
                    pMHC = DClist(idc)%density
                    bindtime = get_bindtime(cog_ptr,cognate,ctype,pMHC,kpar)
                    cell%unbindtime(nbnd) = tnow + bindtime
                    nadd = nadd + 1
                    DClist(idc)%nbound = DClist(idc)%nbound + 1
                    if (cognate) then
                        DClist(idc)%ncogbound = DClist(idc)%ncogbound + 1
                    endif
                    if (cell%DCbound(1) == 0 .and. cell%DCbound(2) /= 0) then
                        write(logmsg,'(a,3i6)') 'Error: binder: DCbound order: ',kcell,cell%DCbound
		                call logger(logmsg)
				        ok = .false.
						return
                    endif
					if (track_DCvisits .and. ctype == TAGGED_CELL) then
						if (revisit(kcell,idc)) then
							Nrevisits = Nrevisits + 1
							cell%revisits = cell%revisits + 1
						else
							Nvisits = Nvisits + 1
							cell%visits = cell%visits + 1
							cell%ndclist = cell%ndclist + 1
							cell%dclist(cell%ndclist) = idc
						endif
					endif
                    if (nbnd == MAX_DC_BIND) exit
                endif
            enddo
        endif
    endif
    if (bound .and. cognate) nb = nb+1
    if (cell%dcbound(1) /= 0) then
		nbt = nbt + 1
		if (cognate) nbc = nbc + 1
	endif
enddo
!write(*,*) 'me, nadd, nsub: ',me,nadd,nsub,nadd-nsub
!write(*,'(f6.2,i4,2f6.3)') tnow/60,nc,real(nbc)/nc,real(nbt)/nt
write(nflog,'(f6.2,i4,2f6.3)') tnow/60,nc,real(nbc)/nc,real(nbt)/nt
ok = .true.
end subroutine

!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
logical function revisit(ic,idc)
integer :: ic
integer :: idc
integer :: i,n

revisit = .false.
n = cellist(ic)%ndclist
if (n == 0) return
do i = 1,n
	if (cellist(ic)%dclist(i) == idc) then
		revisit = .true.
		return
	endif
enddo
end function

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine summary
integer :: idc

if (use_DC) then
    do idc = 1,globalvar%NDC
        if (DClist(idc)%alive) then
            write(*,*) idc,DClist(idc)%nbound,real(DClist(idc)%nbound)/nbindmax
        else
            write(*,*) idc,' dead'
        endif
    enddo
endif
end subroutine

!-----------------------------------------------------------------------------------------
! To test add_Tcell()
!-----------------------------------------------------------------------------------------
subroutine add_random_cells(n, ctype, gen, stage, region)
integer :: n, ctype, gen, stage, region
integer :: k, x, y, z, site(3), slots, kcell
integer :: kpar=0
logical :: ok

k = 0
do while (k < n)
    x = random_int(1,NX,kpar)
    y = random_int(1,NY,kpar)
    z = random_int(1,NZ,kpar)
    site = (/x,y,z/)
    slots = getslots(site)
    if (occupancy(x,y,z)%indx(1) >= 0 .and. slots < BOTH) then
        if (dbug) write(*,'(a,7i6)') 'add_random_cells: ',site,occupancy(x,y,z)%indx,slots
        call add_Tcell(site,ctype,gen,stage,region,kcell,ok)
        if (dbug) write(*,'(a,7i6)') 'after add_random_cells: ',site,occupancy(x,y,z)%indx,slots
        call checkslots('add_random_cells: ',site)
        k = k+1
    endif
enddo

end subroutine

!-----------------------------------------------------------------------------------------
! Scans the cell list and builds the counts nz_sites() and nz_excess().
! nz_sites(k)  = number of available sites at the slice z = k
! nz_excess(k) = total excess of cells over sites in the paracortex zone with z >= k
! These values are used in jumper() to adjust the jump probabilities.  The probability
! of jumps in the -z direction is increased by an increment that is proportional to
! nz_excess(k)/nz_sites(k)
! Need to account for DC sites
!-----------------------------------------------------------------------------------------
subroutine scanner
integer :: x, y, z, ns, nc, nst, nct,excess, indx(2), nextra, i, k, idn, imin=0, nz1, nz2
!real, allocatable :: eratio(:), nz_sites0(:)    !, nz_totsites(:)
real :: eratio(NZ), nz_sites0(NZ)
real :: df, dfmin

if (exit_region == EXIT_EVERYWHERE) return

!allocate(eratio(NZ))
!allocate(nz_sites0(NZ))
!allocate(nz_totsites(NZ))
nz_excess = 0
excess = 0
nct = 0
nst = 0
do z = NZ,1,-1
    ns = 0
    nc = 0
    do y = 1,NY
        do x = 1,NX
            indx = occupancy(x,y,z)%indx
            if (indx(1) < 0) cycle       ! OUTSIDE_TAG or DC site (-idc)
            ns = ns+1
            if (indx(1) > 0) nc = nc+1
            if (indx(2) > 0) nc = nc+1
        enddo
    enddo
    nst = nst + ns
    nct = nct + nc
    nz_sites(z) = ns
    nz_cells(z) = nc
    excess = excess + nc - ns
    nz_excess(z) = excess
enddo

!if (mod(istep,100) == 0) then
!    write(*,*) 'nst,nct: ', nst,nct
!    write(*,*) 'nz_cells'
!    write(*,'(10i6)') nz_cells
!    write(*,*) 'nz_sites'
!    write(*,'(10i6)') nz_sites
!endif

! Note: This check fails when DC start to die, because occupancy(:)%indx is not
! updated immediately.  Therefore, we remove the check.
!if (Mnodes == 1) then
!    if (nst - nct /= globalvar%NDCalive*NDCsites - nadd_sites) then
!        write(*,*) 'Site-cell imbalance: ',nst,nct,globalvar%NDCalive*NDCsites-nadd_sites
!        stop
!    endif
!endif
nextra = nst - nct
! This is the imbalance between number of sites and number of T cells
! resulting from (a) DC sites and (b) sites to be added (nadd_sites)
! We need to adjust either nz_sites(:) or nz_cells(:) to bring them into balance,
! so that eratio(:) can be computed correctly to generate a drift.
! The complication arises because we want to spread the adjustment over the slices
! in a way that is proportionate to the number of sites in the slice.
! Choose to adjust nz_sites(:) (arbitrarily).

if (nextra /= 0) then
    if (nextra > 0) then
        idn = -1
    else
        idn = 1
    endif
    nz_sites0 = nz_sites    ! This conveys approx the shape of the blob - we want to maintain this
    do k = 1,abs(nextra)    ! we need to remove/add this many sites from/to nz_sites(:)
        dfmin = 1.0e10
        do i = 1,NZ
            if (nz_sites(i) == 0) cycle
            df = abs(real(nz_sites(i) + idn)/nz_sites0(i) - 1)
            if (df < dfmin) then
                dfmin = df
                imin = i
            endif
        enddo
        nz_sites(imin) = nz_sites(imin) + idn
        nst = nst + idn
    enddo
    nz_excess = 0
    excess = 0
    do z = NZ,1,-1
        excess = excess + nz_cells(z) - nz_sites(z)
        nz_excess(z) = excess
    enddo
endif
nz1 = NZ
nz2 = 1
do z = NZ,1,-1
    if (z == NZ) then
        nz_totsites(z) = nz_sites(z)
    else
        nz_totsites(z) = nz_sites(z) + nz_totsites(z+1)
    endif
    if (nz_sites(z) > 0) then
!        eratio(z) = nz_excess(z)/real(nz_sites(z))
        eratio(z) = 100*nz_excess(z)/real(nz_totsites(z))
        nz1 = min(z,nz1)
        nz2 = max(z,nz2)
    else
        eratio(z) = 0
    endif
enddo
!write(*,*) 'scanner, istep: ',istep
!nz1 = nz1 + 3
!nz2 = nz2 - 3
!if (mod(istep,100) == 0) then
!    er = eratio((nz1+nz2)/2)     ! actually value at the midline
!    ave_er = (nave*ave_er + er)/(nave+1)
!    nave = nave + 1
!    write(*,'(a,2f6.2)') 'scanner: eratio, average (%): ',er,ave_er
!endif
!deallocate(eratio)
!deallocate(nz_sites0)
!deallocate(nz_totsites)
if (constant_efactor) then
    excess_factor = efactor
else
!    excess_factor = (efactor - excess_factor0 - excess_factor1*globalvar%NTcells)/excess_factor2
    excess_factor = efactor + (1-efactor)*exp(-etheta*(globalvar%NTcells - 50000))
endif
end subroutine

!-----------------------------------------------------------------------------------------
! EGRESS_SUPPRESSION_TIME1 is the start of the ramp down to 0
! EGRESS_SUPPRESSION_TIME2 is the start of the ramp up to 1
!-----------------------------------------------------------------------------------------
real function egressFraction(tnow)
real :: tnow

if (tnow/60 < EGRESS_SUPPRESSION_TIME1) then
	egressFraction = 1
elseif (tnow/60 < EGRESS_SUPPRESSION_TIME1 + EGRESS_SUPPRESSION_RAMP) then
	egressFraction = 1 - (tnow/60 - EGRESS_SUPPRESSION_TIME1)/EGRESS_SUPPRESSION_RAMP
elseif (tnow/60 < EGRESS_SUPPRESSION_TIME2) then
	egressFraction = 0
elseif (tnow/60 < EGRESS_SUPPRESSION_TIME2 + EGRESS_SUPPRESSION_RAMP) then
	egressFraction = (tnow/60 - EGRESS_SUPPRESSION_TIME2)/EGRESS_SUPPRESSION_RAMP
else
	egressFraction = 1
endif
end function


!-----------------------------------------------------------------------------------------
! The number of cells that leave exactly matches the number that enter.
! ZIN_FRACTION = fraction of radius over which inflow traffic occurs
!-----------------------------------------------------------------------------------------
subroutine traffic(ok)
logical :: ok
integer :: x, y, z, k, kcell, indx(2), ctype, gen, stage, region, site(3), n, slot
integer :: zin_min, zin_max, node_inflow, node_outflow, add(3), net_inflow, ihr
real(DP) :: R, df, prob
real :: tnow, exfract
logical :: left
integer :: kpar=0

if (.not.use_blob) then
    write(*,*) 'ERROR: traffic: only for blobs'
    stop
endif
ok = .true.
!write(*,*) 'traffic'
tnow = istep*DELTA_T
region = LYMPHNODE
node_inflow = globalvar%InflowTotal
node_outflow = globalvar%OutflowTotal
df = globalvar%InflowTotal - node_inflow
R = par_uni(kpar)
if (dbug) write(nfres,'(a,2f10.6)') 'traffic: in: ',df,R
if (R < df) then
    node_inflow = node_inflow + 1
endif

exfract = 1
if (suppress_egress) then
	exfract = egressFraction(tnow)
endif
node_outflow = exfract*node_outflow

if (steadystate) then
    node_outflow = node_inflow
else
    df = globalvar%OutflowTotal - node_outflow
    R = par_uni(kpar)
    if (dbug) write(nfres,'(a,2f10.6)') 'traffic: out: ',df,R
    if (R < df) then
        node_outflow = node_outflow + 1
    endif
endif
net_inflow = node_inflow - node_outflow
if (dbug) write(nfres,*) 'traffic: in,out: ',node_inflow,node_outflow
gen = 1
add = 0
zin_max = min(z0 + globalvar%Radius,real(NZ))
zin_min = max(zin_max - ZIN_FRACTION*globalvar%Radius,1.0)

! Inflow
k = 0
do while (k < node_inflow)
    R = par_uni(kpar)
    x = 1 + R*NX
    if (dbug) write(nfres,'(a,i4,f15.9)') 'in x R: ',x,R
    R = par_uni(kpar)
    y = 1 + R*NY
    if (dbug) write(nfres,'(a,i4,f15.9)') 'in y R: ',y,R
    if (exit_region == EXIT_EVERYWHERE) then    ! currently => entry everywhere
        R = par_uni(kpar)
        z = 1 + R*NZ
    if (dbug) write(nfres,'(a,i4,f15.9)') 'in z R: ',z,R
    else
        z = random_int(zin_min,zin_max,kpar)
    endif
    indx = occupancy(x,y,z)%indx
    if (dbug) write(nfres,*) 'site,indx: ',x,y,z,indx
    if (indx(1) < 0) cycle      ! OUTSIDE_TAG or DC
    if (indx(1) == 0 .or. indx(2) == 0) then
        site = (/x,y,z/)
        if (evaluate_residence_time) then
            ! This is for measuring residence time
            if (istep > istep_res1 .and. istep <= istep_res2) then
                ctype = RES_TAGGED_CELL   ! Use ctype = 2 = COG_TYPE_TAG for residence time tagging, to have cptr%entrytime
                ninflow_tag = ninflow_tag + 1
            else
                ctype = 1
            endif
        elseif (track_DCvisits) then
            if (ntagged < ntaglimit) then
                ctype = TAGGED_CELL
            else
                ctype = NONCOG_TYPE_TAG
            endif
        else
            ctype = select_cell_type(kpar)
        endif
        if (ctype /= NONCOG_TYPE_TAG) then
            ncogseed = ncogseed + 1
        endif
        call add_Tcell(site,ctype,gen,NAIVE,region,kcell,ok)
        if (dbug) then
            write(nfres,'(a,5i4,i6)') 'added cell: ',k,site,ctype,kcell
         endif
        if (.not.ok) return
        k = k+1
        cycle
    endif
enddo

! Outflow
k = 0
do while (k < node_outflow)
    R = par_uni(kpar)
    if (dbug) write(nfres,'(a,f10.6)') 'out x R: ',R
    x = 1 + R*NX
    R = par_uni(kpar)
    if (dbug) write(nfres,'(a,f10.6)') 'out y R: ',R
    y = 1 + R*NY
    R = par_uni(kpar)
    if (dbug) write(nfres,'(a,f10.6)') 'out z R: ',R
    z = 1 + R*NZ        ! any z is OK to exit
    if (exit_region == EXIT_EVERYWHERE) then
        ! accept it
    elseif (exit_region == EXIT_LOWERHALF) then
        if (z > z0) cycle
    endif
    indx = occupancy(x,y,z)%indx
    if (indx(1) < 0) cycle      ! OUTSIDE_TAG or DC 
    if (indx(2) > 0) then
        slot = 2
    elseif (indx(1) > 0) then
        slot = 1
    else
        cycle
    endif
    kcell = indx(slot)
    if (SIMULATE_PERIPHERY) then
		prob = 1/PERI_PROBFACTOR	! prob of allowing this cell to exit 
		if (associated(cellist(kcell)%cptr)) then
			gen = get_generation(cellist(kcell)%cptr)
			call get_stage(cellist(kcell)%cptr,stage,region)
			if (gen >= PERI_GENERATION .and. globalvar%NDCcapable == 0) then
				prob = 1
			elseif (gen == 1) then	! suppress egress for undivided cognate cells.  
				prob = 0 
			endif
		endif
		R = par_uni(kpar)
		if (R > prob) cycle
	endif
    site = (/x,y,z/)
    call cell_exit(kcell,slot,site,left)
    if (.not.left) cycle
    k = k+1

    if (evaluate_residence_time) then
        if (cellist(kcell)%ctype == RES_TAGGED_CELL) then
            noutflow_tag = noutflow_tag + 1
            restime_tot = restime_tot + tnow - cellist(kcell)%entrytime
            ihr = (tnow - cellist(kcell)%entrytime)/60. + 1
            Tres_dist(ihr) = Tres_dist(ihr) + 1
        endif
    endif
    if (track_DCvisits .and. cellist(kcell)%ctype == TAGGED_CELL) then
        n = cellist(kcell)%visits
        DCvisits(n) = DCvisits(n) + 1
        DCvisits(globalvar%NDC+1) = DCvisits(globalvar%NDC+1) + cellist(kcell)%revisits
    endif
enddo
nadd_sites = nadd_sites + net_inflow
globalvar%NTcells = globalvar%NTcells + net_inflow
globalvar%NTcellsPer = globalvar%NTcellsPer + node_outflow
end subroutine

!-----------------------------------------------------------------------------------------
! Determine whether a cognate T cell is licensed to exit the DCU.
! Possible exit rules:
! (1) gen >= NGEN_EXIT
! (2) gen > 1 and act < EXIT_THRESHOLD
! (3) Allow exit to any cognate cell that gets there (was gen > 1 and CD69 < CD69_threshold)
!-----------------------------------------------------------------------------------------
logical function exitOK(p,ctype)
type(cog_type), pointer :: p
integer :: ctype
integer :: gen, cd4_8
real :: act
!integer :: kpar = 0

gen = get_generation(p)

exitOK = .false.
if (exit_rule == 1) then
    if (gen >= NGEN_EXIT) then
        exitOK = .true.
    endif
elseif (exit_rule == 2) then
    if (gen > 1) then
        act = get_activation(p)
        cd4_8 = ctype - 1
        if (act < EXIT_THRESHOLD(cd4_8)) then
            exitOK = .true.
        endif
    endif
elseif (exit_rule == 3) then
!    if (p%CD69 == 0) then
!        exitOK = .true.
!    else
!    if (p%CD69 < CD69_threshold) then
!        R = par_uni(kpar)
!        if (R < 1 - p%CD69/CD69_threshold) then
            exitOK = .true.
!        endif
!    endif
endif
end function

!-----------------------------------------------------------------------------------------
! This version aims to determine outflow on the basis of the number of exit sites (which
! depends on NTcells) and on how many cells arrive at the exits, which also involves the
! calibration parameter chemo_K_exit and possibly another.
! The number of exits is greater than 1/2 the desired outflow, and the calibration
! parameter exit_prob is used to adjust the mean outflow.
! PROBLEM:
! It seems that CHEMO_K_EXIT has no effect on the exit rate
! For some reason, CHEMO_K_EXIT = 0 is NOT the same as USE_EXIT_CHEMOTAXIS = false,
! but it SHOULD be.
!-----------------------------------------------------------------------------------------
subroutine portal_traffic(ok)
logical :: ok
integer :: iexit, esite(3), site(3), k, slot, indx(2), kcell, ne, ipermex, ihr, nv
integer :: x, y, z, ctype, gen, region, node_inflow, node_outflow, net_inflow, npossible, iloop, nloops
real(DP) :: R, df
logical :: left
integer, allocatable :: permex(:)
integer :: kpar = 0
real :: tnow, ssfract, exfract
real :: exit_prob   ! 12 hr => 0.2, 24 hr => 0.1

ok = .true.
region = LYMPHNODE
!exit_prob = ep_factor/residence_time
!exit_prob = 0.1	! for Tres = 24
!exit_prob = 0.2	! for Tres = 12

!exit_prob = 0.1*24/RESIDENCE_TIME	! factor of 0.1 is OK for CHEMO_K_EXIT = 1.0
exit_prob = 1.0

node_inflow = globalvar%InflowTotal
node_outflow = globalvar%OutflowTotal
df = globalvar%InflowTotal - node_inflow
R = par_uni(kpar)
if (R < df) then
    node_inflow = node_inflow + 1
endif

if (steadystate) then
    ssfract = real(globalvar%NTcells)/globalvar%NTcells0
    if (ssfract > 1.01) then
        node_inflow = node_inflow - 1
    elseif (ssfract < 0.99) then
        node_inflow = node_inflow + 1
    endif
    exit_prob = ssfract*exit_prob
else
    df = globalvar%OutflowTotal - node_outflow
    R = par_uni(kpar)
    if (R < df) then
        node_outflow = node_outflow + 1
    endif
endif
if (computed_outflow) then
	nloops = 5
else
	nloops = 2		! 1
endif

tnow = istep*DELTA_T
exfract = 1
if (suppress_egress) then
	exfract = egressFraction(tnow)
endif

! Inflow
gen = 1
k = 0
do while (k < node_inflow)
    R = par_uni(kpar)
    x = 1 + R*NX
    R = par_uni(kpar)
    y = 1 + R*NY
    if (exit_region /= EXIT_LOWERHALF) then    ! => entry everywhere
        R = par_uni(kpar)
        z = 1 + R*NZ
    else
!        z = random_int(zin_min,zin_max,kpar)
        write(*,*) 'EXIT_LOWERHALF not simulated with chemotaxis'
        stop
    endif
    indx = occupancy(x,y,z)%indx
    if (indx(1) < 0) cycle      ! OUTSIDE_TAG or DC
    if (indx(1) == 0 .or. indx(2) == 0) then
        site = (/x,y,z/)
        if (evaluate_residence_time) then     ! try cancelling this, use a different method
            ! This is for measuring residence time
            if (istep > istep_res1 .and. istep <= istep_res2) then
                ctype = RES_TAGGED_CELL   ! Use ctype = 2 = COG_TYPE_TAG for residence time tagging, to have cptr%entrytime
                                          ! This is no longer necessary, %entrytime has been shifted.
                ninflow_tag = ninflow_tag + 1
            else
                ctype = NONCOG_TYPE_TAG
            endif
        elseif (track_DCvisits) then
            if (ntagged < ntaglimit) then
                ctype = TAGGED_CELL
            else
                ctype = NONCOG_TYPE_TAG
            endif
        else
            ctype = select_cell_type(kpar)
        endif
        if (ctype /= NONCOG_TYPE_TAG) then
            ncogseed = ncogseed + 1
        endif

        call add_Tcell(site,ctype,gen,NAIVE,region,kcell,ok)
        if (.not.ok) return
        k = k+1
        cycle
    endif
enddo

! Outflow
npossible = 0
ne = 0
if (globalvar%lastexit > max_exits) then
	write(logmsg,'(a,2i4)') 'Error: portal_traffic: globalvar%lastexit > max_exits: ',globalvar%lastexit, max_exits
	call logger(logmsg)
	ok = .false.
	return
endif
allocate(permex(globalvar%lastexit))
do k = 1,globalvar%lastexit
	permex(k) = k
enddo
call permute(permex,globalvar%lastexit,kpar)
do iloop = 1,nloops
do ipermex = 1,globalvar%lastexit
	iexit = permex(ipermex)
	if (exitlist(iexit)%ID == 0) cycle
	if (exfract /= 1) then
	    R = par_uni(kpar)
	    if (R > exfract) cycle
	endif
	esite = exitlist(iexit)%site
	indx = occupancy(esite(1),esite(2),esite(3))%indx
	do slot = 2,1,-1
		kcell = indx(slot)
		if (kcell > 0) npossible = npossible + 1
		if (kcell > 0 .and. par_uni(kpar) < exit_prob) then
			call cell_exit(kcell,slot,esite,left)
			if (left) then
				ne = ne + 1
				if (evaluate_residence_time) then
					if (cellist(kcell)%ctype == RES_TAGGED_CELL) then
						noutflow_tag = noutflow_tag + 1
						restime_tot = restime_tot + tnow - cellist(kcell)%entrytime
						ihr = (tnow - cellist(kcell)%entrytime)/60. + 1
						Tres_dist(ihr) = Tres_dist(ihr) + 1
!                        write(*,*) 'entry: ',cellist(kcell)%entrytime
					endif
				endif
				if (track_DCvisits .and. cellist(kcell)%ctype == TAGGED_CELL) then
					nv = cellist(kcell)%visits
					DCvisits(nv) = DCvisits(nv) + 1
					DCvisits(globalvar%NDC+1) = DCvisits(globalvar%NDC+1) + cellist(kcell)%revisits
				endif
				if (ne == node_outflow .and. computed_outflow) exit
			endif
		endif
	enddo
	if (ne == node_outflow .and. computed_outflow) exit
enddo
if (ne == node_outflow .and. computed_outflow) exit
enddo
deallocate(permex)

!if (npossible /= 2*globalvar%lastexit) then
!	write(*,*) 'npossible: ',npossible,2*globalvar%lastexit
!endif

!if (ne > node_outflow) then
!	write(*,*) 'Excess egress: ',ne - node_outflow
!else
!	write(*,*) 'Deficit outflow: ',node_outflow - ne
!endif
net_inflow = node_inflow - ne
!if (net_inflow < 0) then
!	write(*,*) 'Net outflow: ',-net_inflow,node_inflow,ne,globalvar%Vascularity
!endif
nadd_sites = nadd_sites + net_inflow
total_in = total_in + node_inflow
total_out = total_out + ne
globalvar%NTcells = globalvar%NTcells + net_inflow
globalvar%NTcellsPer = globalvar%NTcellsPer + node_outflow
globalvar%Radius = (globalvar%NTcells*3/(4*PI))**0.33333
!write(*,'(a,4i8)') 'portal_traffic: ',node_inflow,ne,net_inflow,globalvar%NTcells 
end subroutine

!-----------------------------------------------------------------------------------------
! The cell kcell in slot at esite(:) is a candidate to exit.  If it meets the criteria,
! left = .true.
! Currently, when a T cell exits the LN its ID is set to 0, and a gap is recorded in cellist(:).
! To keep track of a T cell in the periphery, we need to keep it in cellist(:), but somehow
! flag it to indicate that it is not in the LN.  Use region.
!-----------------------------------------------------------------------------------------
subroutine cell_exit(kcell,slot,esite,left)
integer :: kcell, slot, esite(3)
logical :: left
integer :: x, y, z, ctype, gen, stage, region
real :: tnow
logical :: cognate, activated
type(cog_type), pointer :: p

tnow = istep*DELTA_T
left = .false.
if (cellist(kcell)%DCbound(1) /= 0) return     ! MUST NOT BE BOUND TO A DC!!!!!!!!!!!!!!!
if (evaluate_residence_time .or. track_DCvisits) then
    cognate = .false.
elseif (associated(cellist(kcell)%cptr)) then
    cognate = .true.
    p => cellist(kcell)%cptr
!    stage = get_stage(p)
	call get_stage(p,stage,region)
    gen = get_generation(p)
!        tin = tnow - p%entrytime
!        if (tin < min_transit_time) cycle     ! prevent immediate exit of newly arrived cell
!        if (stage /= NAIVE .and. gen < NGEN_EXIT) cycle
!        if (stage /= NAIVE) then
!    if (use_DC .and. gen <= 1) return            ! prevents immediate exit of newly arrived cell, in fact of undivided cell
!    if (stage == ACTIVATED) then   ! this allows cells to exit with gen < NGEN_EXIT (unlike DCU paper runs)
        ctype = cellist(kcell)%ctype
        if (.not.exitOK(p,ctype)) then
            write(*,*) 'cell_exit: cognate exit suppressed: ',kcell
            return
        endif
!    endif
else
    cognate = .false.
endif
! For initial testing, remove cells that leave the LN
cellist(kcell)%DCbound = 0
x = esite(1)
y = esite(2)
z = esite(3)
if (slot == 2) then
    occupancy(x,y,z)%indx(2) = 0
else
    if (occupancy(x,y,z)%indx(2) == 0) then
        occupancy(x,y,z)%indx(1) = 0
    else    ! shift cell on indx(2) to indx(1)
        occupancy(x,y,z)%indx(1) = occupancy(x,y,z)%indx(2)
        occupancy(x,y,z)%indx(2) = 0
    endif
endif
if (cognate) then
	if (stage > CLUSTERS) then
		activated = .true.
	else
		activated = .false.
	endif
    if (.not.evaluate_residence_time .and. activated) then
		call efferent(p,ctype)
	endif
!	if (activated) then
!		write(logmsg,'(a,i4,2f8.1)') "activated cognate cell left: stage: ",stage,cellist(kcell)%entrytime,tnow
!		call logger(logmsg)
!	else 
!		call logger("non-activated cognate cell left") 
!    endif
	if (SIMULATE_PERIPHERY .and. activated) then 
		region = PERIPHERY
		call set_stage_region(p,stage,region)
	else
		ngaps = ngaps + 1
		gaplist(ngaps) = kcell
		cellist(kcell)%ID = 0
	    cognate_list(p%cogID) = 0
	endif
else
	ngaps = ngaps + 1
	gaplist(ngaps) = kcell
	cellist(kcell)%ID = 0
endif
left = .true.
end subroutine

!-----------------------------------------------------------------------------------------
! Vascularity responds to the VEGF level.  The rate of production of VEGF is proportional to
! either:
! (VEGF_MODEL_1) a specified inflammation signal
! or
! (VEGF_MODEL_2) the total DC activity level (i.e. antigen load).
!
! In Model 1
!	VEGFsignal depends on the inflammation level
! In Model 2
!   VEGFsignal depends on total DC activity = total antigen density (normalized)
!
!   VEGF_baserate = constitutive rate of VEGF production (per LN volume, i.e. per T cell)
!   dVEGFdt = current rate of VEGF secretion (artificial units)
!   VEGF = current total mass of VEGF (artificial units)
!   c_vegf = current VEGF concentration (artificial units)
!   vasc_beta, vasc_n = Hill function parameters for dependence of vascularity growth
!                    on VEGF concentration
!   vasc_maxrate = maximum rate constant for growth of relative vascularity
!   vasc_decayrate = rate constant for decline of relative vascularity
!   Vascularity = current relative vascularity
!   Note: if inflammation = 0, i.e. VEGFsignal = 0, there should be no change in vascularity,
!   i.e. we should have dVdt = 0.
!-----------------------------------------------------------------------------------------
subroutine vascular
real :: VEGFsignal=0, dVEGFdt, c_vegf, dVdt

!write(*,*) 'vascular'
if (.not.vary_vascularity) then
    globalvar%Vascularity = 1.0
    return
endif
if (VEGF_MODEL == 1) then
    VEGFsignal = get_inflammation() ! Rate of secretion of VEGF is proportional to inflammation 
elseif (VEGF_MODEL == 2) then
    VEGFsignal = get_DCactivity()   ! Rate of secretion of VEGF is proportional to total DC antigen activity
endif
dVEGFdt = VEGFsignal*VEGF_alpha + VEGF_baserate - VEGF_decayrate*globalvar%VEGF
! Mass of VEGF is augmented by rate, and is subject to decay
globalvar%VEGF = globalvar%VEGF + dVEGFdt*DELTA_T
c_vegf = globalvar%VEGF/globalvar%NTcells   ! concentration (proportional to, anyway)
if (VEGF_MODEL == 2) then
    dVdt = vasc_maxrate*hill(c_vegf,vasc_beta,vasc_n)*globalvar%Vascularity - vasc_decayrate*(globalvar%Vascularity - 1)
else	! VEGF_MODEL = 1
    dVdt = vasc_maxrate*hill(c_vegf,vasc_beta*c_vegf_0,vasc_n)*globalvar%Vascularity - vasc_decayrate*globalvar%Vascularity
endif
globalvar%Vascularity = globalvar%Vascularity + dVdt*DELTA_T
!if (mod(istep,10) == 0) then
!	write(logmsg,'(a,i6,2e10.3,3f8.3)') 'vasc: ',istep/240,VEGFsignal,c_vegf,globalvar%VEGF, &
!		globalvar%Vascularity,real(globalvar%NTcells)/globalvar%NTcells0
!	call logger(logmsg)
!endif
!write(*,*) 'dVEGFdt, c_vegf, dVdt: ',dVEGFdt, c_vegf, dVdt, globalvar%Vascularity
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
real function hill(x,b,n)
real :: x, b
integer :: n
hill = x**n/(x**n + b**n)
end function

!-----------------------------------------------------------------------------------------
! The total number of T cells added to the blob since the last balancing is nadd_sites.
! (Note that this can be negative, indicating a net loss of T cells). 
! A balancing is triggered either when this count exceeds a limit nadd_limit, which is a
! specified fraction of the original T cell population NTcells0, or when the time since
! the last balancing exceeds BALANCER_INTERVAL, and an adjustment to site count is needed.
! Sites needed or made available by a change to the DC population are accounted for.
! Either new sites are added (made available) or existing sites are removed (made OUTSIDE).
!-----------------------------------------------------------------------------------------
subroutine balancer(ok)
logical :: ok
integer :: nadd_total, nadd_limit, n, nadded
integer :: k, idc, naddDC, naddex, nremex, Nexits0, dexit
real :: tnow
!real, save :: lasttime = 0
integer :: kpar = 0
logical :: blob_changed

ok = .true.
!write(*,*) 'balancer: ',istep
if (IN_VITRO) then
	if (ndeadDC > 0) then
		do k = 1,ndeadDC
			idc = DCdeadlist(k)
			call clearDC(idc)
		enddo
		ndeadDC = 0
	endif
	return
endif
blob_changed = .false.
tnow = istep*DELTA_T
nadd_limit = 0.01*globalvar%NTcells0
nadd_total = nadd_sites

if ((abs(nadd_total) > nadd_limit) .or. (tnow > lastbalancetime + BALANCER_INTERVAL)) then
    if (dbug) write(nflog,*) 'balancer: nadd_total: ',nadd_total,nadd_limit,lastbalancetime,BALANCER_INTERVAL
    if (dbug) write(nflog,*) 'call squeezer'
	call squeezer(.false.)
	if (dbug) write(nflog,*) 'did squeezer'
    ! adding/removal of sites on occupancy
    if (ndeadDC > 0) then
		if (dbug) write(nflog,*) 'remove DCs'
        do k = 1,ndeadDC
            idc = DCdeadlist(k)
            call clearDC(idc)
        enddo
        nadd_total = nadd_total - ndeadDC*NDCsites
        ndeadDC = 0
        if (dbug) write(nflog,*) 'balancer: removed DCs: nadd_total: ',nadd_total
    endif
    naddDC = 0
    if (use_DC .and. use_DCflux) then    ! need to account for incoming DC
        naddDC = DCinflux(lastbalancetime,tnow,kpar)
        if (naddDC > 0) then
!			write(*,*) 'add DCs: ',naddDC
            call place_DCs(naddDC,nadded)
            nadd_total = nadd_total + nadded*NDCsites
            if (dbug) write(nflog,*) 'balancer: added DCs: NDCtotal: ',NDCtotal
        endif
    endif
    if (dbug) write(nflog,*) 'nadd_total: ',nadd_total
    if (nadd_total > 0) then
        n = nadd_total
	    if (dbug) write(nflog,*) 'call add_sites'
        call addSites(n,ok)
        if (.not.ok) return
	    if (dbug) write(nflog,*) 'did add_sites'
        blob_changed = .true.
    elseif (nadd_total < 0) then
        n = -nadd_total
	    if (dbug) write(nflog,*) 'call removeSites'
        call removeSites(n,ok)
        if (.not.ok) return
	    if (dbug) write(nflog,*) 'did removeSites'
        blob_changed = .true.
    else
        n = 0
    endif
    if (n /= 0) then
        write(logmsg,'(a,4i6)') 'Error: balancer: remaining sites: ',n, &
			globalvar%NTcells,globalvar%NDCalive,globalvar%Nsites
		call logger(logmsg)
        ok = .false.
        return
    endif
!    if (DC_motion) then
!        call test_moveDC
!    endif
   if (naddDC > 0) then
		if (dbug) write(nflog,*) 'call reassign_DC'
        call reassign_DC(kpar,ok)
        if (.not.ok) return
		if (dbug) write(nflog,*) 'did reassign_DC'
    endif
    call growDC
! The cognate list is maintained at the time that a cell arrives or leaves
    if (globalvar%NTcells+globalvar%NDCalive*NDCsites /= globalvar%Nsites) then
	    write(logmsg,'(a,3i6)') 'Error: balancer: cells /= sites: ', &
			globalvar%NTcells,globalvar%NDCalive,globalvar%Nsites
	    call logger(logmsg)
	    ok = .false.
	    return
	endif
    lastbalancetime = tnow
    nadd_sites = 0
    if (blob_changed) then
		if (dbug) write(nflog,*) 'call make_split'
        call make_split
        if (use_diffusion) then
            call setup_minmax
        endif
        if (dbug) write(nflog,'(a,2i6)') 'balancer: nadd_total, radius: ',nadd_total,int(globalvar%radius)
    endif
    if (USE_PORTAL_EGRESS) then
		call adjustExits
	endif
else
    call set_globalvar
endif
if (use_portal_egress .and. use_traffic) then
	if (globalvar%NTcells < globalvar%NTcells0 .and. .not.SURFACE_PORTALS) then
		! Set Nexits = Nexits0 for steady-state maintenance.  To wrap up the end of the response.
!        Nexits0 = exit_fraction*globalvar%NTcells0
		Nexits0 = requiredExitPortals(globalvar%NTcells0)
		if (globalvar%Nexits < Nexits0) then
			do k = 1,Nexits0 - globalvar%Nexits
				call addExitPortal()
			enddo
!			write(*,*) '-------------------------------------'
!			write(*,*) 'Set Nexits to base steady-state level'
!			write(*,*) '--------------------------'
!			write(*,*) 'Nexits: ',globalvar%Nexits
!			write(*,*) '--------------------------'
		elseif (globalvar%Nexits > Nexits0) then
			call removeExits(globalvar%Nexits - Nexits0)
!			write(*,*) '-------------------------------------'
!			write(*,*) 'Set Nexits to base steady-state level'
!			write(*,*) '--------------------------'
!			write(*,*) 'Nexits: ',globalvar%Nexits
!			write(*,*) '--------------------------'
		endif
		return
	endif
!    dexit = exit_fraction*globalvar%NTcells - globalvar%Nexits 
    dexit = requiredExitPortals(globalvar%NTcells) - globalvar%Nexits
    if (dexit > 0) then
        naddex = dexit
!       write(*,*) '--------------------------'
!       write(*,*) 'Nexits: ',globalvar%Nexits
!       write(*,*) '--------------------------'
!       write(nflog,*) 'added: Nexits: ',naddex, globalvar%Nexits
        write(logmsg,*) 'add: Nexits: ',naddex, globalvar%NTcells, globalvar%Nexits, requiredExitPortals(globalvar%NTcells)
        call logger(logmsg) 
        do k = 1,naddex
            call addExitPortal()
        enddo
    elseif (dexit < -1) then
        nremex = -dexit-1
!        write(nflog,*) 'balancer: dexit < 0: ',requiredExitPortals(globalvar%NTcells) 
!            write(*,*) '--------------------------'
!            write(*,*) 'Nexits: ',globalvar%Nexits
!            write(*,*) '--------------------------'
!	        write(nflog,*) 'removed: Nexits: ',nremex,globalvar%Nexits
!        write(logmsg,*) 'remove: Nexits: ',nremex,requiredExitPortals(globalvar%NTcells),globalvar%Nexits,globalvar%NTcells
!        call logger(logmsg)
        call removeExits(nremex)
    endif
endif
end subroutine

!-----------------------------------------------------------------------------------------
! Sites are added to occupancy().  The general idea is to preserve the shape of the 
! T cell zone, i.e. a spherical blob remains roughly spherical.
!-----------------------------------------------------------------------------------------
subroutine addSites(n,ok)
integer :: n
logical :: ok
integer :: maxblist,x,y,z,i,k,nb,nadd,idc,kdc,site0(3),site(3)
real :: x2, y2, z2, r2, r2min, r2max, r(3)
logical :: bdryflag
integer, allocatable :: t(:), bdrylist(:,:)
real, allocatable :: r2list(:)

ok = .true.
!write(logmsg,'(a,i5,i7,f8.2)') 'add_sites: ',n,globalvar%Nsites,globalvar%Radius
!call logger(logmsg)
r2 = globalvar%Radius*globalvar%Radius
maxblist = 4*PI*r2*0.1*globalvar%Radius
allocate(t(maxblist))
allocate(bdrylist(3,maxblist))
allocate(r2list(maxblist))

r2max = 1.3*r2
r2min = 0.8*r2

!write(*,*) 'add_sites: r2min,r2max: ',r2min,r2max
k = 0
do z = 2,NZ-1
    z2 = (z-z0)*(z-z0)
    if (z2 > r2max) cycle
    do y = 2,NY-1
        y2 = (y-y0)*(y-y0)
        if (y2+z2 > r2max) cycle
        do x = 2,NX-1
            x2 = (x-x0)*(x-x0)
            r2 = x2 + y2 + z2
            if (r2 > r2min .and. r2 < r2max) then
                if (occupancy(x,y,z)%indx(1) == OUTSIDE_TAG) then
                    bdryflag = .false.
                    site0 = (/x,y,z/)
                    do i = 1,6
                        site = site0 + neumann(:,i)
                        if (occupancy(site(1),site(2),site(3))%indx(1) >= 0) then
                            bdryflag = .true.
                            exit
                        endif
                    enddo
                    if (bdryflag .and. k < maxblist) then
                        k = k+1
                        bdrylist(:,k) = site0
                        r2list(k) = r2
                    endif
                endif
            endif
        enddo
    enddo
enddo
nb = k
!write(*,*) 'bdry sites: ',nb,maxblist
if (nb < n) then
    write(logmsg,'(a,2i4)') 'Error: add_sites: insufficient candidate sites: ',nb,n
    call logger(logmsg)
    ok = .false.
    return
endif
do i = 1,nb
    t(i) = i
enddo
call qsort(r2list,nb,t)     ! sort in increasing order
nadd = min(n,nb)
do i = 1,nadd
    site = bdrylist(:,t(i))
    occupancy(site(1),site(2),site(3))%indx = 0
    occupancy(site(1),site(2),site(3))%DC = 0
    if (use_cytokines .and. use_diffusion) then
        cyt(site(1),site(2),site(3),:) = 0
    endif
!    call smear_cyt(site)
    ! Now need to determine %DC(0:DCDIM-1)     ! DC(0) = number of near DCs, DC(k) = ID of kth near DC
    if (use_DC) then
        do idc = 1,globalvar%NDC
            r = DClist(idc)%site - site
            if (norm(r) <= DC_RADIUS) then
	            kdc = occupancy(site(1),site(2),site(3))%DC(0)
		        if (kdc < 0) cycle   ! don't touch a DC site
	            if (kdc < DCDIM-1) then     ! can add to the list
	                kdc = kdc + 1
	                occupancy(site(1),site(2),site(3))%DC(0) = kdc
	                occupancy(site(1),site(2),site(3))%DC(kdc) = idc
		        endif
            elseif (norm(r) <= chemo_radius) then
	            kdc = occupancy(site(1),site(2),site(3))%cDC(0)
		        if (kdc < 0) cycle   ! don't touch a DC site
	            if (kdc < cDCDIM-1) then     ! can add to the list
	                kdc = kdc + 1
	                occupancy(site(1),site(2),site(3))%cDC(0) = kdc
	                occupancy(site(1),site(2),site(3))%cDC(kdc) = idc
		        endif
	        endif
        enddo
    endif
    ! Check for adjacent exit site (portal), and if necessary move it.
!    call adjustExit(site)
enddo
n = n - nadd
globalvar%Nsites = globalvar%Nsites + nadd
deallocate(t)
deallocate(bdrylist)
deallocate(r2list)
end subroutine

!-----------------------------------------------------------------------------------------
!   removeExit(site)
!   placeExitPortal(iexit,site) 
!-----------------------------------------------------------------------------------------
subroutine adjustExits
integer :: x, y, z, dx, dy, dz, iexit, dir, site1(3), site2(3)
real :: u(3)
logical :: ok

!write(*,*) 'adjustExits'
do iexit = 1,globalvar%lastexit
    if (exitlist(iexit)%ID == 0) cycle
	site1 = exitlist(iexit)%site
	ok = .false.
	do dir = 1,27
		if (dir == 14) cycle
		x = site1(1) + jumpvec(1,dir)
		if (x < 1 .or. x > NX) cycle
		y = site1(2) + jumpvec(2,dir)
		if (y < 1 .or. y > NY) cycle
		z = site1(3) + jumpvec(3,dir)
		if (z < 1 .or. z > NZ) cycle
		if (occupancy(x,y,z)%indx(1) < 0) then
			ok = .true.		! site1 is on the blob boundary 
			exit
		endif
	enddo
	if (ok) cycle
	! The exit portal has been enclosed, need to move it to site2 on the boundary if possible 
	u = site1 - Centre
	u = u/sqrt(dot_product(u,u))
	call getBoundarySite(u,site2,ok)
	if (.not.ok) cycle
	if (use_DC) then
		if (toonearDC(site2,globalvar%NDC,exit_DCprox)) cycle    ! exit_DCprox is min distance in sites
	endif
	call removeExitPortal(site1)
!	write(logmsg,*) 'adjustExits: did removeExitPortal: ',iexit,site1
!	call logger(logmsg)
!	call checkExits
	
	! Note: since we are reusing an existing exit index, no need to increment %lastexit
	globalvar%Nexits = globalvar%Nexits + 1
	if (globalvar%Nexits > max_exits) then
		write(*,*) 'Error: adjustExit: too many exits: need to increase max_exits: ',max_exits
		stop
	endif
	call placeExitPortal(iexit,site2)
!	write(logmsg,*) 'adjustExits: did placeExitPortal: ',iexit,site2
!	call logger(logmsg)
!	call checkExits
enddo
end subroutine

!-----------------------------------------------------------------------------------------
! Try to remove n sites from the boundary of the blob (occupancy).
! A boundary site is one with at least 3 Neumann neighbours outside the blob.
!-----------------------------------------------------------------------------------------
subroutine removeSites(n,ok)
integer :: n
logical :: ok
integer :: maxblist,x,y,z,i,nout,k,nb,nr,count,site0(3),site(3),indx(2),kcell
integer :: freeslot
real :: x2, y2, z2, r2, r2min, r2max
logical :: moved
integer, allocatable :: t(:), bdrylist(:,:)
real, allocatable :: r2list(:)

!write(*,'(a,i6)') 'removeSites: ',n
ok = .true.
r2 = globalvar%Radius*globalvar%Radius
maxblist = 4*PI*r2*0.1*globalvar%Radius
allocate(t(maxblist))
allocate(bdrylist(3,maxblist))
allocate(r2list(maxblist))
r2max = 1.3*r2
r2min = 0.8*r2

! Make a list of boundary sites.
k = 0
do z = 1,NZ
    z2 = (z-z0)*(z-z0)
    if (z2 > r2max) cycle
    do y = 1,NY
        y2 = (y-y0)*(y-y0)
        if (y2+z2 > r2max) cycle
        do x = 1,NX
            x2 = (x-x0)*(x-x0)
            r2 = x2 + y2 + z2
            if (r2 > r2min .and. r2 < r2max) then
                indx = occupancy(x,y,z)%indx
                if (indx(1) >= 0) then
                    if (indx(1) > 0 .and. indx(2) > 0) cycle    ! consider only sites with 0 or 1 cell
                    site0 = (/x,y,z/)
                    nout = 0
                    do i = 1,6
                        site = site0 + neumann(:,i)
                        if (site(1) < 1 .or. site(1) > NX) cycle
                        if (site(2) < 1 .or. site(2) > NY) cycle
                        if (site(3) < 1 .or. site(3) > NZ) cycle
                        if (occupancy(site(1),site(2),site(3))%indx(1) == OUTSIDE_TAG) then
                            nout = nout + 1
                        endif
                    enddo
                    if (nout >= 1 .and. k < maxblist) then
                        k = k+1
                        bdrylist(:,k) = site0
                        if (r2 > r2max) then
                            write(logmsg,'(a,4i6,f8.2)') 'Error: removeSites: bad r2: ',k,x,y,z,r2
                            call logger(logmsg)
                            ok = .false.
                            return
                        endif
                        r2list(k) = r2
                    endif
                endif
            endif
        enddo
    enddo
enddo
nb = k
!write(*,*) 'bdry sites: ',nb,maxblist,r2max
do i = 1,nb
    t(i) = i
enddo
call qsort(r2list,nb,t)     ! sort in increasing order 
nr = min(nb,n)
count = 0
do i = nb,1,-1
    if (count == nr) exit
    site0 = bdrylist(:,t(i))
    indx = occupancy(site0(1),site0(2),site0(3))%indx
    if (indx(1) > 0 .and. indx(2) > 0) cycle    ! consider only sites with 0 or 1 cell
    if (indx(1) < 0) then
        write(logmsg,'(a,i6)') 'Error: removeSites: site has a DC: ',indx(1)
        call logger(logmsg)
        ok = .false.
        return
    endif
    moved = .false.
    do k = 1,2
        kcell = indx(k)
        if (kcell > 0) then
            ! This cell needs to be moved to a nearby site
            call get_free_slot(occupancy,NX,site0,site,freeslot)
            if (freeslot == 0) cycle
            occupancy(site(1),site(2),site(3))%indx(freeslot) = kcell
            cellist(kcell)%site = site
            occupancy(site0(1),site0(2),site0(3))%indx = OUTSIDE_TAG
            occupancy(site0(1),site0(2),site0(3))%DC = 0
            if (occupancy(site0(1),site0(2),site0(3))%exitnum < 0) then
!				write(*,*) 'removeSites:  need to move exit'
                call removeExitPortal(site0)
                ! Let the balancer add a site back again (if needed) 
            endif
            globalvar%Nsites = globalvar%Nsites - 1
            count = count + 1
            moved = .true.
        elseif (kcell < 0) then
            write(logmsg,'(a,5i6)') 'Error: removeSites: negative indx: ',site0,indx
            call logger(logmsg)
            ok = .false.
            return
        endif
        if (moved) exit
    enddo
enddo
n = n - count
deallocate(t)
deallocate(bdrylist)
deallocate(r2list)
end subroutine

!--------------------------------------------------------------------------------
! When a new empty site is added, the cytokines in the neighbouring sites are
! mixed in.
! NOT USED
!--------------------------------------------------------------------------------
subroutine smear_cyt(site0)
integer :: site0(3)
real :: csum(MAX_CYT)
integer :: k, n, site(3), jlist(MAXRELDIR)

n = 0
csum(1:Ncytokines) = 0
do k = 1,27
    if (k == 14) cycle
    site = site0 + jumpvec(:,k)
    if (occupancy(site(1),site(2),site(3))%indx(1) >= 0) then
        n = n+1
        jlist(n) = k
        csum(1:Ncytokines) = csum(1:Ncytokines) + cyt(site(1),site(2),site(3),1:Ncytokines)
    endif
enddo
csum = csum/(n+1)
do k = 1,n
    site = site0 + jumpvec(:,jlist(k))
    cyt(site(1),site(2),site(3),1:Ncytokines) = csum(1:Ncytokines)
enddo
cyt(site0(1),site0(2),site0(3),:) = csum(1:Ncytokines)

end subroutine

!--------------------------------------------------------------------------------
! Determines the average radius of the blob, which is needed for exitter().
! Note: this code only applies to the spherical blob case.
! Get range by averaging over a small number of grid sites near the centre.
! NOT USED
!--------------------------------------------------------------------------------
subroutine sizer(nb,r2list)
integer :: nb
real :: r2list(:)
integer :: k
real :: rmean

rmean = 0
do k = 1,nb
    rmean = rmean + sqrt(r2list(k))
enddo
rmean = rmean/nb
write(*,*) 'sizer: nb,rmean: ',nb,rmean
end subroutine

!--------------------------------------------------------------------------------
! Counts efferent cognate cells, and records their generation distribution.
! Only activated cells (stage >= CLUSTERS) are counted
!--------------------------------------------------------------------------------
subroutine efferent(p,ctype)
type (cog_type), pointer :: p
integer :: ctype, gen, i
real :: avid

gen = get_generation(p)
if (ctype > NCTYPES) then
    write(*,*) 'efferent: bad cell type:', ctype
    stop
endif
!localres%dN_EffCogTC(ctype)  = localres%dN_EffCogTC(ctype) + 1
!localres%dN_EffCogTCGen(gen) = localres%dN_EffCogTCGen(gen) + 1
!localres%N_EffCogTC(ctype)   = localres%N_EffCogTC(ctype) + 1
!localres%N_EffCogTCGen(gen)  = localres%N_EffCogTCGen(gen) + 1
totalres%dN_EffCogTC(ctype)  = totalres%dN_EffCogTC(ctype) + 1
totalres%dN_EffCogTCGen(gen) = totalres%dN_EffCogTCGen(gen) + 1
totalres%N_EffCogTC(ctype)   = totalres%N_EffCogTC(ctype) + 1
totalres%N_EffCogTCGen(gen)  = totalres%N_EffCogTCGen(gen) + 1

if (log_results) then
    ! Record avidity statistics for exiting cells
    avid = p%avidity
    if (avid_count%logscale) then
        avid = log10(avid)
    endif
    if (avid_count%nbins == 1) then
        i = 1
    else
        !i = (avid-avidity_min)*1.01/avidity_step + 1
        i = (avid-avid_count%binmin)/avid_count%binstep + 1.5
        i = max(i,1)
        !i = min(i,avidity_nlevels)
        i = min(i,avid_count%nbins)
    endif
    avid_count%ndist(i) = avid_count%ndist(i) + 1
endif
end subroutine

!-----------------------------------------------------------------------------------------
! Using the complete list of cells, cellist(), extract info about the current state of the
! paracortex.  This info must be supplemented by counts of cells that have died and cells that
! have returned to the circulation.
!-----------------------------------------------------------------------------------------
subroutine show_snapshot(ok)
logical :: ok
integer :: kcell, ctype, stype, ncog, noncog, ntot, stage, region, i, iseq, error
integer :: gen, ngens, neffgens, teffgen, dNdead, Ndead
real :: stim(FINISHED), IL2sig(FINISHED), tgen, tnow, fac, act, cyt_conc, mols_pM
type (cog_type), pointer :: p
integer :: nst(FINISHED)
integer, allocatable :: gendist(:)
integer, allocatable :: div_gendist(:)  ! cells that are capable of dividing
character*(6) :: numstr
character*(256) :: msg

ok = .true.
allocate(gendist(TC_MAX_GEN))
allocate(div_gendist(TC_MAX_GEN))
tnow = istep*DELTA_T
noncog = 0
ncog = 0
ntot = 0
nst = 0
stim = 0
IL2sig = 0
gendist = 0
div_gendist = 0
do kcell = 1,nlist
    if (cellist(kcell)%ID == 0) cycle
    p => cellist(kcell)%cptr
    ntot = ntot + 1
    ctype = cellist(kcell)%ctype
    stype = struct_type(ctype)
    if (stype == COG_TYPE_TAG) then
        ncog = ncog + 1
!        stage = get_stage(p)
		call get_stage(p,stage,region)
        nst(stage) = nst(stage) + 1
        stim(stage) = stim(stage) + p%stimulation
        IL2sig(stage) = IL2sig(stage) + get_IL2store(p)
        gen = get_generation(p)
        if (gen < 0 .or. gen > TC_MAX_GEN) then
            write(logmsg,'(a,2i6)') 'show_snapshot: bad gen: ',kcell,gen
            call logger(logmsg)
            ok = .false.
            return
        endif
        gendist(gen) = gendist(gen) + 1
        if ((gen == 1 .and. p%stimulation > FIRST_DIVISION_THRESHOLD(1)) .or. &
			(gen > 1 .and. p%stimulation > DIVISION_THRESHOLD(1))) then
			div_gendist(gen) = div_gendist(gen) + 1
        endif
        max_TCR = max(p%stimulation,max_TCR)
    elseif (stype == NONCOG_TYPE_TAG) then
        noncog = noncog + 1
    else
        write(*,*) 'ERROR: show_snapshot: bad stype: ',ctype,stype
        stop
    endif
enddo
do i = 1,FINISHED
    if (nst(i) > 0) then
        stim(i) = stim(i)/nst(i)
        IL2sig(i) = IL2sig(i)/nst(i)
    else
        stim(i) = 0
        IL2sig(i) = 0
    endif
enddo
tgen = sum(gendist)
do i = TC_MAX_GEN,1,-1
    if (gendist(i) /= 0) exit
enddo
ngens = i

teffgen = sum(totalres%N_EffCogTCGen(1:TC_MAX_GEN))
do i = TC_MAX_GEN,1,-1
    if (totalres%N_EffCogTCGen(i) /= 0) exit
enddo
neffgens = i

dNdead = totalres%dN_Dead
Ndead = totalres%N_Dead
act = get_DCactivity()
mols_pM = L_um3*M_pM/(globalvar%NTcells*Vc*Navo)

if (teffgen > 0) then
    fac = 1/real(teffgen)
else
    fac = 0
endif
if (.not.use_TCP .and. use_cognate) then
write(*,'(a)') '----------------------------------------------------------------------'
write(*,*) 'use_cognate: ',use_cognate
write(*,'(a,i6,3i8,a,2i8)') 'snapshot: ',istep,ntot,ncogseed,ncog,'     dead: ',dNdead,Ndead
write(*,'(a,7i7)')   '# in stage:  ',nst
write(*,'(a,7f7.0)') 'stimulation: ',stim
write(*,'(a,7f7.0)') 'IL2 signal:  ',IL2sig
write(*,'(a,2i8,4x,i8)') 'Recent efferent: ',totalres%dN_EffCogTC(2:3),sum(totalres%dN_EffCogTC)
write(*,'(a,2i8,4x,i8)') 'Total efferent:  ',totalres%N_EffCogTC(2:3),teffgen
write(*,'(a,10i6)')   'gen dist: ',(i,i=1,10)
write(*,'(a)')        'In node:  '
write(*,'(10x,10f6.3)') gendist(1:ngens)/tgen
write(*,'(a)')        'Efferent: '
write(*,'(10x,10f6.3)') fac*totalres%N_EffCogTCGen(1:neffgens)
write(*,'(a,i6,a,f6.0)') 'Live DC: ',globalvar%NDCalive,'  max_TCR: ',max_TCR
if (use_cytokines) then
    do iseq = 1,Ncytokines
        if (use_diffusion) then
            cyt_conc = cyt_mean(iseq)
        else
            cyt_conc = cyt_mols(iseq)*mols_pM   ! -> conc in pM
        endif
        write(*,'(3a,f8.4)') 'Mean cytokine conc: ',cyt_name(cyt_tag(iseq)),'  ',cyt_conc
    enddo
endif
write(*,'(a)') '----------------------------------------------------------------------'

!call check_cognate_list
!kcog = 1
!kcell = cognate_list(kcog)
!if (kcell > 0) then
!    write(*,'(2i6,f8.4,f8.1)') kcog,kcell,cellist(kcell)%cptr%stimrate,cellist(kcell)%cptr%CD69
!endif
write(*,*) '========= Average time to IL2 threshold: ',nIL2thresh,tIL2thresh/max(1,nIL2thresh)
write(*,'(a)') '----------------------------------------------------------------------'
endif

!call get_cognate_dist(ncog1,ncog2)

write(nfout,'(2(i8,f8.0),6i8,25f7.4)') istep,tnow/60,globalvar%NDCalive,act,ntot,ncogseed,ncog,dNdead,Ndead,teffgen, &
    fac*totalres%N_EffCogTCGen(1:TC_MAX_GEN)
if (use_tcp) then
!    if (.not.awp_1%is_open) then
!        call logger("in show_snapshot: awp_1 is not open")
!    endif
!    write(msg,'(2(i6,f8.0),5i8)') istep,tnow,globalvar%NDCalive,act,ntot,ncogseed,ncog,Ndead,teffgen
!    call winsock_send(awp_1,msg,len_trim(msg),error)
    msg = ''
    do i = 1,ngens
		write(numstr,'(i6)') gendist(i)
		msg = trim(msg)//trim(adjustl(numstr))
		msg = trim(msg)//'('
		write(numstr,'(i6)') div_gendist(i)
		msg = trim(msg)//trim(adjustl(numstr))
		msg = trim(msg)//')-'
	enddo
    call logger(msg)
endif

! To plot outflow variation with time
!write(nfout,'(2f8.2)') tnow/60,globalvar%OutflowTotal

!write(nfout,'(a,7f7.0)') 'stimulation: ',stim
!write(nfout,'(a,7f7.0)') 'IL2 signal:  ',IL2sig
deallocate(gendist)
totalres%dN_EffCogTC = 0
totalres%dN_EffCogTCGen = 0
totalres%dN_Dead = 0

if (save_DCbinding) then
    ncog = sum(dcbind(0:MAX_COG_BIND))
    write(*,'(21i6)') dcbind(0:MAX_COG_BIND)
    write(*,*) 'ncog: ',ncog
    write(nfdcbind,'(f8.2,i6,30f7.4)') tnow,ncog,dcbind(0:MAX_COG_BIND)/real(ncog)
    dcbind = 0
endif
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine compute_stim_dist
integer :: kcell, ctype, stype, ncog, noncog
real :: s
type (cog_type), pointer :: p

ncog = 0
noncog = 0
do kcell = 1,nlist
    if (cellist(kcell)%ID == 0) cycle
    p => cellist(kcell)%cptr
    ctype = cellist(kcell)%ctype
    stype = struct_type(ctype)
    if (stype == COG_TYPE_TAG) then
        ncog = ncog + 1
        s = p%stimulation
    elseif (stype == NONCOG_TYPE_TAG) then
        noncog = noncog + 1
    else
        write(*,*) 'ERROR: compute_stim_dist: bad stype: ',ctype,stype
        stop
    endif
enddo
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine check_cognate_list
integer :: k, kcell

do k = 1,lastcogID
    kcell = cognate_list(k)
    if (kcell > nlist) then
        write(*,*) 'check_cognate_list: bad cognate_list entry > nlist: ',lastcogID,k,kcell,nlist
        stop
    endif
enddo
end subroutine

!-----------------------------------------------------------------------------------------
! The aim is to display the distribution of cognate cells in the z direction, in order
! to investigate failure of cognate cell proliferation to scale with NTcells0.
! For now, just look at cognate fractions in upper and lower hemispheres.
!-----------------------------------------------------------------------------------------
subroutine get_cognate_dist(ncog1,ncog2)
integer :: z,kcell,ntot1,ntot2,ncog1,ncog2
logical :: cognate

ntot1 = 0
ntot2 = 0
ncog1 = 0
ncog2 = 0
do kcell = 1,nlist
    if (cellist(kcell)%ID == 0) cycle
    z = cellist(kcell)%site(3)
    cognate = (associated(cellist(kcell)%cptr))
    if (z < z0) then
        ntot1 = ntot1 + 1
        if (cognate) ncog1 = ncog1 + 1
    else
        ntot2 = ntot2 + 1
        if (cognate) ncog2 = ncog2 + 1
    endif
enddo
write(*,*) 'Lower: ',ntot1,ncog1,real(ncog1)/ntot1
write(*,*) 'Upper: ',ntot2,ncog2,real(ncog2)/ntot2
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine testqsort(n)
integer :: n
integer :: i
integer :: kpar = 0
real, allocatable :: a(:)
integer, allocatable :: t(:)

call rng_initialisation

allocate(a(n))
allocate(t(n))
do i = 1,n
    a(i) = par_uni(kpar)
    t(i) = i
enddo
call qsort(a,n,t)
end subroutine

!-----------------------------------------------------------------------------------------
! To compare with the Matlab code: Lauff_a.m
!-----------------------------------------------------------------------------------------
subroutine testupdater
real :: state(IL2_NP), statep(IL2_NP)
integer :: ctype = COG_CD4_TAG, ID = 1
real :: t, dt, th, S_TCR, C_IL2, mrate, crate, dm, dc
logical :: producing = .true., dbgflag = .false.
logical :: first, ok
integer :: i, nsub = 1, nhours = 24
real ::  C0 = 100

call read_cell_params(ok)
first = .true.
!state = 0  ! for Lauff_a.m comparison
state = (/1200,1300,300,1400,4000,90000/)
t = 0
crate = 0.01
S_TCR = 5000
C_IL2 = C0
dt = DELTA_T/nsub
!write(*,*)
!write(*,'(6x,9a8)') 'Min ','CD25_S','Comp_S','CD25_I','Comp_I','IL2_I','S_Rec','IL-2','rate'
do i = 1,nhours*60*4*nsub
    t = t + dt
    th = t/60

    C_IL2 = C0

!    C_IL2 = 20.0*(1 - 0.5*th/nhours)
!    if (th > 24) producing = .false.
!    call IL2_update(ctype,t,S_TCR,state,producing,C_IL2,dt,mrate,dbgflag) ! old IL2_update, no xcase
    call IL2_update(ID,ctype,t,S_TCR,state,statep,first,producing,C_IL2,Vc,dt,mrate,dbgflag)
    dm = (mrate/Vc)*DELTA_T*M_pM*L_um3/Navo
    dc = crate*DELTA_T*M_pM*L_um3/Navo
!    if (mod(i,60*4*nsub) == 0) then
!        write(*,'(i6,7f8.0,2f8.1)') i,t,state,C_IL2,mrate
!    endif
    first = .false.
enddo

end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine vascular_test
integer :: nsteps = 10*24*60/DELTA_T
real :: inflow0, act, tnow, exfract

globalvar%VEGF = 0
globalvar%vascularity = 1.00
globalvar%NTcells0 = 100000
use_exit_chemotaxis = .false.

globalvar%NDC = DC_FACTOR*globalvar%NTcells0/TC_TO_DC
globalvar%NTcells = globalvar%NTcells0
inflow0 = globalvar%NTcells0*DELTA_T/(residence_time*60)
write(*,*) TC_TO_DC,globalvar%NDC,nsteps,globalvar%NTcells

do istep = 1,nsteps
	tnow = istep*DELTA_T
    call vascular
    call generate_traffic(inflow0)
    if (suppress_egress) then
		exfract = egressFraction(tnow)
		globalvar%OutflowTotal = exfract*globalvar%OutflowTotal
	endif
    globalvar%NTcells = globalvar%NTcells + globalvar%InflowTotal - globalvar%OutflowTotal
    if (mod(istep,60) == 0) then
        act = get_DCactivity()
        write(nfout,'(2f8.3,e12.3,f8.3,i8,2f6.1)') istep*DELTA_T/60,act,globalvar%VEGF/globalvar%NTcells, &
            globalvar%vascularity,globalvar%NTcells,globalvar%InflowTotal,globalvar%OutflowTotal
    endif
enddo

end subroutine

subroutine testrnor
integer :: n = 100000
integer :: k, j
integer :: kpar = 0
real :: r, rmin=1.0e10, rmax = -1.0e10

do k = 1,n
    do j = 1,n
        r = par_rnor(kpar)
        rmin = min(r,rmin)
        rmax = max(r,rmax)
    enddo
    write(*,'(i12,2e12.4)') k,rmin,rmax
enddo
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine testrandom
integer :: n = 100000
integer :: k, j, ncog
integer :: kpar = 0

ncog = 0
do k = 1,n
    j = select_cell_type(kpar)
    if (j == COG_CD4_TAG) ncog = ncog + 1
    if (mod(k,1000) == 0) then
        write(*,*) ncog,real(ncog)/k
    endif
enddo
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine old_process_command_line
integer :: i, cnt, len, status
character :: c*(64), b*(256)
character*(64) :: progname

call get_command (b, len, status)
if (status .ne. 0) then
    write (*,*) 'get_command failed with status = ', status
    stop
end if
!write (*,*) 'command line = ', b (1:len)
call get_command_argument (0, c, len, status)
if (status .ne. 0) then
    write (*,*) 'Getting command name failed with status = ', status
    stop
end if
!write (*,*) 'command name = ', c (1:len)
progname = c(1:len)
cnt = command_argument_count ()
!write (*,*) 'number of command arguments = ', cnt
if (cnt < 1) then
    write(*,*) 'Use: ',trim(progname),' num_cpu'
    stop
!    Mnodes = 4      ! for profiling
!    write(*,*) 'Ruuning with Mnodes = ',Mnodes
endif
do i = 1, cnt
    call get_command_argument (i, c, len, status)
    if (status .ne. 0) then
        write (*,*) 'get_command_argument failed: status = ', status, ' arg = ', i
        stop
    end if
!    write (*,*) 'command arg ', i, ' = ', c (1:len)
    if (i == 1) then
!!!        read(c(1:len),'(i)') Mnodes
        read(c(1:len),*) Mnodes
        write(*,*) 'Requested threads: ',Mnodes
    elseif (i == 2) then
        inputfile = c(1:len)
        write(*,*) 'Input file: ',inputfile
    elseif (i == 3) then
        outputfile = c(1:len)
        write(*,*) 'Output file: ',outputfile
!    elseif (i == 4) then
!        resultfile = c(1:len)
!        write(*,*) 'Result file: ',resultfile
    endif
end do
write (*,*) 'command line processed'
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine write_Tres_dist
integer :: k, kmax

do k = int(days*24),1,-1
    if (Tres_dist(k) /= 0) then
        kmax = k
        exit
    endif
enddo
write(nfout,'(a,i6)') 'Tres_dist ',kmax
do k = 1,kmax
    write(nfout,'(f6.1,e12.4)') k-0.5,Tres_dist(k)/noutflow_tag
enddo
end subroutine

!-----------------------------------------------------------------------------------------
! The distributions gendist, tcrdist and aviddist are all for the current DCU population.
! avid_count is for the recent efferent population.
!-----------------------------------------------------------------------------------------
subroutine write_results
integer :: kcell, ctype, stype, ncog, ntot, i
integer :: gen
real :: tcr, avid, dtcr, hour
type (cog_type), pointer :: p
integer :: gendist(TC_MAX_GEN),aviddist(MAX_AVID_LEVELS),tcrdist(tcr_nlevels)
character*(60) :: fmtstr = '(f6.2,2i8,4x,15f7.4,4x,10f7.4,4x,10f7.4,4x,10i7)'

write(fmtstr(14:15),'(i2)') TC_MAX_GEN
write(fmtstr(24:25),'(i2)') tcr_nlevels
write(fmtstr(34:35),'(i2)') avidity_nlevels
write(fmtstr(44:45),'(i2)') avidity_nlevels
hour = istep*DELTA_T/60
dtcr = TCR_limit/TCR_nlevels
gendist = 0
aviddist = 0
tcrdist = 0
ntot = 0
ncog = 0
do kcell = 1,nlist
    if (cellist(kcell)%ID == 0) cycle
    p => cellist(kcell)%cptr
    ntot = ntot + 1
    ctype = cellist(kcell)%ctype
    stype = struct_type(ctype)
    if (stype == COG_TYPE_TAG) then
        ncog = ncog + 1
        ! TCR stimulation distribution
        tcr = p%stimulation
        i = tcr/dtcr + 1
        i = min(i,TCR_nlevels)
        tcrdist(i) = tcrdist(i) + 1
        ! T cell generation distribution
        gen = get_generation(p)
        gendist(gen) = gendist(gen) + 1
        ! T cell avidity distribution
        avid = p%avidity
        if (avidity_logscale) then
            avid = log10(avid)
        endif
        if (avidity_nlevels == 1) then
            i = 1
        else
            i = (avid-avidity_min)/avidity_step + 1.5
!           write(nfout,'(a,2f8.4,3i7)') 'Count: ',p%avidity,avid,i,kcell,p%cogID
            i = max(i,1)
            i = min(i,avidity_nlevels)
        endif
        aviddist(i) = aviddist(i) + 1
    endif
enddo
if (fix_avidity) then
    write(nfres,fmtstr) hour,ntot,ncog,real(gendist)/ncog,real(tcrdist)/ncog,real(aviddist)/ncog, &
        avid_count%ndist
else
    write(nfres,fmtstr) hour,ntot,ncog,real(gendist)/ncog,real(tcrdist)/ncog
endif
avid_count_total%ndist = avid_count_total%ndist + avid_count%ndist
avid_count%ndist = 0
!write(nfout,'(8i6)') aviddist
end subroutine

!-----------------------------------------------------------------------------------------
! A different approach to tagging cells for residence time computation
!-----------------------------------------------------------------------------------------
subroutine tag_cells
integer :: kcell
type(cell_type),pointer :: cell

do kcell = 1,nlist
    cell => cellist(kcell)
    if (cell%ID == 0) cycle
    cell%ctype = RES_TAGGED_CELL
enddo
end subroutine

!-----------------------------------------------------------------------------------------
! Various logging counters are initialized here.
!-----------------------------------------------------------------------------------------
subroutine init_counters

ninflow_tag = 0
noutflow_tag = 0
restime_tot = 0
if (log_results) then
    if (.not.use_cognate) then
        write(*,*) 'No use logging results with no cognate cells'
        stop
    endif
    avid_count%nbins = avidity_nlevels
    allocate(avid_count%ndist(avid_count%nbins))
    avid_count%period = ntres
    avid_count%logscale = avidity_logscale
    avid_count%binmin = avidity_min
    avid_count%binstep = avidity_step
    avid_count%ndist = 0
    avid_count%total = 0

    allocate(avid_count_total%ndist(avid_count%nbins))
    avid_count_total = avid_count
!    avid_count_total%nbins = avid_count%nbins
!    avid_count_total%period = avid_count%period
!    avid_count_total%logscale = avid_count%logscale
!    avid_count_total%binmin = avidity_min
!    avid_count_total%binstep = avidity_step
!    avid_count_total%ndist = 0
!    avid_count_total%total = 0
endif
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine check_chemoactivity(cave)
real :: cave
type (cell_type), pointer :: cell
integer :: kcell

cave = 0
do kcell = 1,nlist
    if (cellist(kcell)%ID == 0) cycle
    cell => cellist(kcell)
    cave = cave + chemo_active_exit(cell)
enddo
cave = cave/globalvar%NTcells
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine test_cum_prob
real :: m = 30, s = 2.0
real :: p1, p2, a
integer :: i

p1 = log(m)
p2 = log(s)
do i = 1,30
    a = 2*i
    write(*,'(i4,2f8.3)') i,a,1-cum_prob_lognormal(a,p1,p2)
enddo
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine check_pause
!!!use ifport
logical :: paused

inquire(file=pausefile,exist=paused)
if (paused) then
	call logger('Pause order received')
	do while (paused)
!!!		call sleepqq(100)
        call sleeper(1)   ! Too coarse!
		inquire(file=pausefile,exist=paused)
	enddo
	call logger('Resuming ...')
endif
end subroutine

!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
subroutine get_dimensions(NX_dim,NY_dim,NZ_dim) BIND(C)
!DEC$ ATTRIBUTES DLLEXPORT :: get_dimensions
use, intrinsic :: iso_c_binding
integer(c_int) :: NX_dim,NY_dim,NZ_dim

NX_dim = NX
NY_dim = NY
NZ_dim = NZ
end subroutine

!-----------------------------------------------------------------------------------------
! Using the complete list of cells, cellist(), extract info about the current state of the
! paracortex.  This info must be supplemented by counts of cells that have died and cells that
! have returned to the circulation.
! We now store stim() and IL2sig() for cells in the periphery.
!-----------------------------------------------------------------------------------------
subroutine get_summary(summaryData) BIND(C)
!DEC$ ATTRIBUTES DLLEXPORT :: get_summary
use, intrinsic :: iso_c_binding
integer(c_int) :: summaryData(*)
logical :: ok
integer :: kcell, ctype, stype, ncog(2), noncog, ntot, nbnd, stage, region, i, iseq, error
integer :: gen, ngens, neffgens, teffgen, dNdead, Ndead, nact
real :: stim(2*STAGELIMIT), IL2sig(2*STAGELIMIT), tgen, tnow, fac, act, cyt_conc, mols_pM
type (cog_type), pointer :: p
integer :: nst(FINISHED)
integer, allocatable :: gendist(:)
integer, allocatable :: div_gendist(:)  ! cells that are capable of dividing
character*(6) :: numstr
character*(256) :: msg

ok = .true.
allocate(gendist(TC_MAX_GEN))
allocate(div_gendist(TC_MAX_GEN))
tnow = istep*DELTA_T
noncog = 0
ncog = 0
ntot = 0
nbnd = 0
nst = 0
stim = 0
IL2sig = 0
gendist = 0
div_gendist = 0
do kcell = 1,nlist
    if (cellist(kcell)%ID == 0) cycle
    p => cellist(kcell)%cptr
    if (associated(p)) then
		call get_stage(p,stage,region)
	else
		stage = 0
		region = LYMPHNODE
	endif
	if (region == LYMPHNODE) then
	    ntot = ntot + 1
	endif
    ctype = cellist(kcell)%ctype
    stype = struct_type(ctype)
    if (stype == COG_TYPE_TAG) then
        ncog(region) = ncog(region) + 1
		if (cellist(kcell)%DCbound(1) > 0 .or. cellist(kcell)%DCbound(2) > 0) then
			nbnd = nbnd + 1
		endif
        nst(stage) = nst(stage) + 1
        stim(stage) = stim(stage) + p%stimulation
        IL2sig(stage) = IL2sig(stage) + get_IL2store(p)
        gen = get_generation(p)
        if (gen < 0 .or. gen > TC_MAX_GEN) then
            write(logmsg,'(a,2i6)') 'get_summary: bad gen: ',kcell,gen
            call logger(logmsg)
            ok = .false.
            return
        endif
        gendist(gen) = gendist(gen) + 1
        if ((gen == 1 .and. p%stimulation > FIRST_DIVISION_THRESHOLD(1)) .or. &
			(gen > 1 .and. p%stimulation > DIVISION_THRESHOLD(1))) then
			div_gendist(gen) = div_gendist(gen) + 1
        endif
        max_TCR = max(p%stimulation,max_TCR)
    elseif (stype == NONCOG_TYPE_TAG) then
        noncog = noncog + 1
    else
        write(*,*) 'ERROR: show_snapshot: bad stype: ',ctype,stype
        stop
    endif
enddo
do i = 1,FINISHED
    if (nst(i) > 0) then
        stim(i) = stim(i)/nst(i)
        IL2sig(i) = IL2sig(i)/nst(i)
    else
        stim(i) = 0
        IL2sig(i) = 0
    endif
enddo
tgen = sum(gendist)
do i = TC_MAX_GEN,1,-1
    if (gendist(i) /= 0) exit
enddo
ngens = i

teffgen = sum(totalres%N_EffCogTCGen(1:TC_MAX_GEN))
do i = TC_MAX_GEN,1,-1
    if (totalres%N_EffCogTCGen(i) /= 0) exit
enddo
neffgens = i

dNdead = totalres%dN_Dead
Ndead = totalres%N_Dead
act = get_DCactivity()
mols_pM = L_um3*M_pM/(globalvar%NTcells*Vc*Navo)

if (teffgen > 0) then
    fac = 1/real(teffgen)
else
    fac = 0
endif
if (.not.use_TCP .and. use_cognate) then
write(*,'(a)') '----------------------------------------------------------------------'
write(*,'(a,i6,4i8,a,2i8)') 'snapshot: ',istep,ntot,ncogseed,ncog,'     dead: ',dNdead,Ndead
write(*,'(a,7i7)')   '# in stage:  ',nst
write(*,'(a,7f7.0)') 'stimulation: ',stim
write(*,'(a,7f7.0)') 'IL2 signal:  ',IL2sig
write(*,'(a,2i8,4x,i8)') 'Recent efferent: ',totalres%dN_EffCogTC(2:3),sum(totalres%dN_EffCogTC)
write(*,'(a,2i8,4x,i8)') 'Total efferent:  ',totalres%N_EffCogTC(2:3),teffgen
write(*,'(a,10i6)')   'gen dist: ',(i,i=1,10)
write(*,'(a)')        'In node:  '
write(*,'(10x,10f6.3)') gendist(1:ngens)/tgen
write(*,'(a)')        'Efferent: '
write(*,'(10x,10f6.3)') fac*totalres%N_EffCogTCGen(1:neffgens)
write(*,'(a,i6,a,f6.0)') 'Live DC: ',globalvar%NDCalive,'  max_TCR: ',max_TCR
if (use_cytokines) then
    do iseq = 1,Ncytokines
        if (use_diffusion) then
            cyt_conc = cyt_mean(iseq)
        else
            cyt_conc = cyt_mols(iseq)*mols_pM   ! -> conc in pM
        endif
        write(*,'(3a,f8.4)') 'Mean cytokine conc: ',cyt_name(cyt_tag(iseq)),'  ',cyt_conc
    enddo
endif
write(*,'(a)') '----------------------------------------------------------------------'

write(*,*) '========= Average time to IL2 threshold: ',nIL2thresh,tIL2thresh/max(1,nIL2thresh)
write(*,'(a)') '----------------------------------------------------------------------'
endif

write(nfout,'(2(i8,f8.0),7i8,25f7.4)') istep,tnow/60,globalvar%NDCalive,act,ntot,ncogseed,ncog,dNdead,Ndead,teffgen, &
    fac*totalres%N_EffCogTCGen(1:TC_MAX_GEN)

if (use_tcp) then
    msg = ''
    do i = 1,ngens
		write(numstr,'(i6)') gendist(i)
		msg = trim(msg)//trim(adjustl(numstr))
		msg = trim(msg)//'('
		write(numstr,'(i6)') div_gendist(i)
		msg = trim(msg)//trim(adjustl(numstr))
		msg = trim(msg)//')-'
	enddo
    call logger(msg)
endif
deallocate(gendist)

totalres%dN_EffCogTC = 0
totalres%dN_EffCogTCGen = 0
totalres%dN_Dead = 0

if (save_DCbinding) then
    ncog = sum(dcbind(0:MAX_COG_BIND))
    write(*,'(21i6)') dcbind(0:MAX_COG_BIND)
    write(*,*) 'ncog: ',ncog
    write(nfdcbind,'(f8.2,i6,30f7.4)') tnow,ncog,dcbind(0:MAX_COG_BIND)/real(ncog)
    dcbind = 0
endif

nact = 100*act
summaryData(1:12) = (/istep,globalvar%NDCalive,nact,ntot,ncogseed,ncog,Ndead,teffgen,nbnd,int(globalvar%InflowTotal),globalvar%Nexits/)
write(nflog,*) 'ndivisions = ',ndivisions
end subroutine
!-------------------------------------------------------------------------------- 
!--------------------------------------------------------------------------------
subroutine get_scene(nTC_list,TC_list,nDC_list,DC_list,nbond_list,bond_list) BIND(C)
!DEC$ ATTRIBUTES DLLEXPORT :: get_scene
use, intrinsic :: iso_c_binding
!integer(c_int) :: nmono_list, mono_list(*), ncap_list, npit_list, nclast_list
!real(c_float) :: cap_list(*), pit_list(*), clast_list(*)
integer(c_int) :: nDC_list, nTC_list, nbond_list, DC_list(*), TC_list(*), bond_list(*)
integer :: k, kc, kcell, site(3), j, jb, idc, dcsite(3)
real :: dcstate
integer :: idcstate, itcstate, stype, ctype, stage, region
real :: Tcell_diam = 0.9
real :: DC_diam = 1.8
integer :: gen, bnd(2)

! DC section
if (globalvar%NDC > 0) then
	k = 0
    do kcell = 1,globalvar%NDC
        if (DClist(kcell)%alive) then
			k = k+1
			j = 5*(k-1)
            site = DClist(kcell)%site
            dcstate = min(1.0,DClist(kcell)%density/DC_ANTIGEN_MEDIAN)
            idcstate = 100*dcstate
            ! Need dcstate to convey antigen density level (normalized to 0-1)
!            write(nfpos,'(a2,i4,3i4,f4.1,f5.2)') 'D ',k-1, site, DC_diam, dcstate
			DC_list(j+1) = kcell-1
			DC_list(j+2:j+4) = site
			DC_list(j+5) = idcstate
        endif
    enddo
    nDC_list = k
endif

if (.not.IV_SHOW_NONCOGNATE) then
	! T cell section
	k = 0
	do kc = 1,lastcogID
		kcell = cognate_list(kc)
		if (kcell > 0) then
			call get_stage(cellist(kcell)%cptr,stage,region)
			if (region /= LYMPHNODE) cycle
			k = k+1
			j = 5*(k-1)
			site = cellist(kcell)%site
	!        tcstate = mod(kcell,2) + 1
			gen = get_generation(cellist(kcell)%cptr)
			bnd = cellist(kcell)%DCbound
!			if (get_stage(cellist(kcell)%cptr) == NAIVE) then
			if (stage == NAIVE) then
				itcstate = 0
			else
				if (bnd(1) == 0 .and. bnd(2) == 0) then
	!				tcstate = (gen-1.0)/(TC_MAX_GEN-1.0)*spectrum_max*spectrum_freefraction
	!				tcstate = (gen/TC_MAX_GEN)*spectrum_max*spectrum_freefraction
					itcstate = gen
				else
	!				tcstate = spectrum_max
					itcstate = 99
				endif
			endif
			! Need tcstate to convey non-activated status, i.e. 0 = non-activated
!			write(nfpos,'(a2,i6,3i4,f4.1,i3)') 'T ',kc-1, site, Tcell_diam, itcstate
			TC_list(j+1) = kc-1
			TC_list(j+2:j+4) = site
			TC_list(j+5) = itcstate
		endif
	enddo
    nTC_list = k
	! Bond section
	k = 0
	do kc = 1,lastcogID
		kcell = cognate_list(kc)
		if (kcell > 0) then
			site = cellist(kcell)%site
			do jb = 1,2
				idc = cellist(kcell)%DCbound(jb)
				if (idc /= 0) then
					if (DClist(idc)%capable) then
						k = k+1
						j = 2*(k-1)
						dcsite = DClist(idc)%site
	!                    dcstate = mod(idc,2) + 1
						dcstate = 1
	!                    write(nfpos,'(a,i8,3i4,4i4)') 'Node: ',nd, site, dcsite-site, dcstate
!						write(nfpos,'(a2,2i5)') 'B ',kc-1,idc-1
						bond_list(j+1) = kc-1
						bond_list(j+2) = idc-1
					endif
				endif
			enddo
		endif
	enddo
	nbond_list = k
else
	! T cell section
	k = 0
	do kcell = 1,nlist
		if (cellist(kcell)%ID == 0) cycle  ! gap
		k = k+1
		j = 5*(k-1)
		site = cellist(kcell)%site
		bnd = cellist(kcell)%DCbound
		ctype = cellist(kcell)%ctype
		stype = struct_type(ctype)
		if (stype == NONCOG_TYPE_TAG) then
			itcstate = -1
		else
			gen = get_generation(cellist(kcell)%cptr)
!			if (get_stage(cellist(kcell)%cptr) == NAIVE) then
			call get_stage(cellist(kcell)%cptr,stage,region)
			if (stage == NAIVE) then
				itcstate = 0
			else
				if (bnd(1) == 0 .and. bnd(2) == 0) then
					itcstate = gen
				else
					itcstate = 99
				endif
			endif
		endif
		! Need tcstate to convey non-activated status, i.e. 0 = non-activated
!		write(nfpos,'(a2,i6,3i4,f4.1,i3)') 'T ',kcell-1, site, Tcell_diam, itcstate
		TC_list(j+1) = kc-1
		TC_list(j+2:j+4) = site
		TC_list(j+5) = itcstate
	enddo
	nTC_list = k
	! Bond section
	k = 0
	do kcell = 1,nlist
		if (cellist(kcell)%ID == 0) cycle  ! gap
		site = cellist(kcell)%site
		do jb = 1,2
			idc = cellist(kcell)%DCbound(jb)
			if (idc /= 0) then
				if (DClist(idc)%capable) then
					k = k+1
					j = 2*(k-1)
					dcsite = DClist(idc)%site
!                    dcstate = mod(idc,2) + 1
					dcstate = 1
!					write(nfpos,'(a2,2i5)') 'B ',kcell-1,idc-1
					bond_list(j+1) = kcell-1
					bond_list(j+2) = idc-1
				endif
			endif
		enddo
	enddo
	nbond_list = k
endif

end subroutine

!-----------------------------------------------------------------------------------------
! Now this is used only to set use_TCP = .false.
! The lines
!    call get_command (b, len, status)
!    call get_command_argument (0, c, len, status)
! were failing with gfortran (don't know why), but in any case there was no need to
! get the command line arguments in this way.
!-----------------------------------------------------------------------------------------
subroutine process_command_line(ncpu,infile,outfile)
!DEC$ ATTRIBUTES DLLEXPORT :: process_command_line
!DEC$ ATTRIBUTES STDCALL, REFERENCE, MIXED_STR_LEN_ARG, ALIAS:"PROCESS_COMMAND_LINE" :: process_command_line
integer :: i, cnt, len, status
integer :: ncpu
character :: c*(64), b*(256)
character*(64) :: infile,outfile
character*(64) :: progname

!write(*,*) 'process_command_line'
use_TCP = .false.   ! because this is called from para_main()							! --> use_TCP

return

ncpu = 3
infile = 'omp_para.inp'
outfile = 'omp_para.out'
!resfile = 'result.out'
!runfile = ' '

call get_command (b, len, status)
if (status .ne. 0) then
    write (logmsg,'(a,i4)') 'get_command failed with status = ', status
    call logger(logmsg)
    stop
end if
call logger('command: ')
call logger(b)
c = ''
call get_command_argument (0, c, len, status)
if (status .ne. 0) then
    write (*,*) 'Getting command name failed with status = ', status
    write(*,*) c
    stop
end if
progname = c(1:len)
cnt = command_argument_count ()
if (cnt < 1) then
    write(*,*) 'Use: ',trim(progname),' num_cpu'
    stop
endif

do i = 1, cnt
    call get_command_argument (i, c, len, status)
    if (status .ne. 0) then
        write (*,*) 'get_command_argument failed: status = ', status, ' arg = ', i
        stop
    end if
    if (i == 1) then
!        read(c(1:len),'(i)') ncpu
        read(c(1:len),*) ncpu															! --> ncpu
        write(*,*) 'Requested threads: ',ncpu
    elseif (i == 2) then
        infile = c(1:len)																! --> infile
        write(*,*) 'Input file: ',infile
    elseif (i == 3) then
        outfile = c(1:len)																! --> outfile
        write(*,*) 'Output file: ',outfile
!    elseif (i == 4) then
!        resfile = c(1:len)																! --> resfile
!        write(*,*) 'Result file: ',resfile
    endif
end do

end subroutine



!-----------------------------------------------------------------------------------------
! Until I find a time function in gfortran
!----------------------------------------------------------------------------------------- 
!real(8) function timef()
!timef = 0
!end function

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine simulate_step(res) BIND(C)
!DEC$ ATTRIBUTES DLLEXPORT :: simulate_step 
use, intrinsic :: iso_c_binding
integer(c_int) :: res
logical :: ok

res = 0
dbug = .false.
ok = .true.
istep = istep + 1
if (dbug) then
	write(logmsg,*) 'simulate_step: ',istep
	call logger(logmsg)
endif

if (mod(istep,240) == 0) then
!	call checkExits
    globalvar%Radius = (globalvar%NTcells*3/(4*PI))**0.33333
    if (log_traffic) then
        write(nftraffic,'(5i8,3f8.3)') istep, globalvar%NTcells, globalvar%Nexits, total_in, total_out, &
                globalvar%InflowTotal, globalvar%vascularity
    endif
    total_in = 0
    total_out = 0
endif
if (use_cytokines) then
    if (use_diffusion) then
        call diffuser
    else
        call molsynch
    endif
endif

if (use_traffic .and. mod(istep,SCANNER_INTERVAL) == 0) then
	if (dbug) write(nflog,*) 'call scanner'
	call scanner
	if (dbug) write(nflog,*) 'did scanner'
endif
if (dbug) write(nflog,*) 'call mover'
call mover(ok)
if (dbug) write(nflog,*) 'did mover'
if (.not.ok) then
	call logger("mover returned error")
	res = 1
	return
endif

if (use_DC .and. globalvar%NDCalive > 0) then
	if (dbug) write(nflog,*) 'call binder'
    call binder(ok)
    if (.not.ok) then
		call logger('binder returned error')
		res = 1
		return
	endif
	if (dbug) write(nflog,*) 'did binder'
    if (.not.track_DCvisits .and. .not.evaluate_residence_time) then
		if (dbug) write(nflog,*) 'call update_DCstate'
        call update_DCstate(ok)
        if (.not.ok) then
			call logger('update_DCstate returned error')
			res = 1
			return
		endif
		if (dbug) write(nflog,*) 'did update_DCstate'
    endif
endif
if (.not.track_DCvisits .and. .not.evaluate_residence_time) then
	if (dbug) write(nflog,*) 'call updater'
    call updater(ok)
	if (dbug) write(nflog,*) 'did updater'
    if (.not.ok) then
		call logger('updater returned error')
		res = 1
		return
	endif
endif

if (use_traffic) then
    if (vary_vascularity) then	! There is a problem with this system
        call vascular
    endif
    if (use_portal_egress) then
		if (dbug) write(nflog,*) 'call portal_traffic'
        call portal_traffic(ok)
	    if (.not.ok) then
			call logger('portal_traffic returned error')
			res = 1
			return
		endif
		if (dbug) write(nflog,*) 'did portal_traffic'
    else
		if (dbug) write(nflog,*) 'call traffic'
        call traffic(ok)
		if (dbug) write(nflog,*) 'did traffic'
	    if (.not.ok) then
			call logger('traffic returned error')
			res = 1
			return
		endif
    endif
endif
if (dbug) call check_xyz(3)

if (dbug) write(nflog,*) 'call balancer'
call balancer(ok)
if (dbug) write(nflog,*) 'did balancer'
if (.not.ok) then
	call logger('balancer returned error')
	res = 1
	return
endif
!if (dbug) write(nflog,*) 'call set_globalvar'
call set_globalvar
!if (dbug) write(nflog,*) 'did set_globalvar'

end subroutine

!-----------------------------------------------------------------------------------------
! Runs to compute the travel time distributions are based on varying two parameters:
!	TC_COGNATE_FRACTION	! fraction of T cells that are cognate initially
!	TC_TO_DC                ! number of T cells for every DC
! The first takes N_TRAVEL_COG values, the second N_TRAVEL_DC values, while the distribution
! is evaluated at N_TRAVEL_TIME values (intervals of a minute).
! The case to be simulated corresponds to (k_travel_cog, k_travel_dc)
!-----------------------------------------------------------------------------------------
subroutine set_travel_params

!TC_COGNATE_FRACTION = travel_cog(k_travel_cog)
TC_COGNATE_FRACTION = 0
TC_TO_DC = travel_dc(k_travel_dc)
ntravel = 0
travel_dist(k_travel_cog,k_travel_dc,:) = 0
use_traffic = .false.
use_DCflux = .false.
use_cognate = .false.
write(*,*) 'TC_COGNATE_FRACTION: ',TC_COGNATE_FRACTION
write(*,*) 'TC_TO_DC: ',TC_TO_DC
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine connection(awp,port,error)
TYPE(winsockport) :: awp
integer :: port, error
integer :: address = 0
!!!character*(64) :: ip_address = "127.0.0.1"C      ! need a portable way to make a null-terminated C string
character*(64) :: host_name = "localhost"

if (.not.winsock_init(1)) then
    call logger("winsock_init failed")
    stop
endif
!write(nftemp,*) 'did winsock_init'

awp%handle = 0
awp%host_name = host_name
awp%ip_port = port
awp%protocol = IPPROTO_TCP
call Set_Winsock_Port (awp,error)

if (.not.awp%is_open) then
    write(nflog,*) 'Error: connection: awp not open: ',port
endif
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine connecter(ok)
logical :: ok
integer :: error

! Main connection
ok = .true.
call connection(awp_0,TCP_PORT_0,error)
if (awp_0%handle < 0 .or. error /= 0) then
    write(logmsg,'(a)') 'TCP connection to TCP_PORT_0 failed'
    call logger(logmsg)
    ok = .false.
    return
endif
if (.not.awp_0%is_open) then
	write(logmsg,'(a)') 'No connection to TCP_PORT_0'
    call logger(logmsg)
    ok = .false.
    return
endif
write(logmsg,'(a)') 'Connected to TCP_PORT_0  '
call logger(logmsg)

if (use_CPORT1) then
	call connection(awp_1,TCP_PORT_1,error)
	if (awp_1%handle < 0 .or. error /= 0) then
		write(logmsg,'(a)') 'TCP connection to TCP_PORT_1 failed'
		call logger(logmsg)
		ok = .false.
		return
	endif
	if (.not.awp_1%is_open) then
		write(logmsg,'(a)') 'No connection to TCP_PORT_1'
		call logger(logmsg)
		ok = .false.
		return
	endif
	write(logmsg,'(a)') 'Connected to TCP_PORT_1  '
	call logger(logmsg)
endif
! Allow time for completion of the connection
call sleeper(2)
end subroutine


!-----------------------------------------------------------------------------------------
! This subroutine is called to initialize a simulation run.
! ncpu = the number of processors to use
! infile = file with the input data
! outfile = file to hold the output
! runfile = file to pass info to the master program (e.g. Python) as the program executes.
!-----------------------------------------------------------------------------------------
subroutine setup(ncpu,infile,outfile,ok)
integer :: ncpu
character*(*) :: infile, outfile
logical :: ok
character*(64) :: msg
integer :: error

ok = .true.
initialized = .false.
par_zig_init = .false.
Mnodes = ncpu
inputfile = infile
outputfile = outfile
!resultfile = resfile
write(logmsg,*) 'ncpu: ',ncpu 
call logger(logmsg)

#if defined(OPENMP) || defined(_OPENMP)
    call logger("OPENMP defined")
    call omp_initialisation(ok)
    if (.not.ok) return
#else
    call logger("OPENMP NOT defined")
    if (Mnodes > 1) then
        write(logmsg,'(a)') 'No OpenMP, using one thread only'
        call logger(logmsg)
        Mnodes = 1
    endif
#endif

call logger("read_cell_params")
call read_cell_params(ok)
if (.not.ok) return

if (compute_travel_time .and. .not.IN_VITRO) then
	call set_travel_params
endif
Fcognate = TC_COGNATE_FRACTION

ndivisions = 0
call CD69_setparameters(K1_S1P1,K2_S1P1,K1_CD69,K2_CD69)
call array_initialisation(ok)
if (.not.ok) return

if (calibrate_motility) then
	use_DC = .false.
	call motility_calibration
	stop
endif

if (log_traffic) then
    open(nftraffic,file='traffic.out',status='replace')
endif

call place_cells(ok)

if (.not.ok) return

if (use_portal_egress) then
    call placeExits
else
    globalvar%Nexits = 0
    globalvar%lastexit = 0
endif

if (use_exit_chemotaxis) then
	call chemo_setup
endif

if (vary_vascularity) then
	call initialise_vascularity
endif
if (inflammation_level == 0) then
	suppress_egress = .false.
elseif (EGRESS_SUPPRESSION_TIME2 - EGRESS_SUPPRESSION_TIME1 <  EGRESS_SUPPRESSION_RAMP) then
	write(logmsg,*) 'ERROR: EGRESS_SUPPRESSION_RAMP too big'
	call logger(logmsg)
	stop
endif
!write(*,*) 'NTcells, NDCalive, NTsites: ',globalvar%NTcells, globalvar%NDCalive, globalvar%Nsites
!write(*,*) 'nlist: ',nlist

if (use_DC) then
	call update_DCstate(ok)
	if (.not.ok) return
endif
if (use_cytokines) then
    call init_cytokine
!    write(*,*) 'did init_cytokine'
!    write(*,*) 'Ncytokines: ',Ncytokines
endif
call set_globalvar
call make_split
if (.not.IN_VITRO) then
	call scanner
endif
call init_counters

!    call arrayview(NX)
!    read(*,*)
!    stop
if (save_input) then
    call save_inputfile(inputfile)
    call save_parameters
	call save_inputfile(fixedfile)
endif

initialized = .true.

write(logmsg,'(a,i6)') 'Startup procedures have been executed: initial T cell count: ',globalvar%NTcells0
call logger(logmsg)

end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine wrapup
integer :: ierr
logical :: isopen

call logger('doing wrapup ...')
ierr = 0
if (allocated(xoffset)) deallocate(xoffset)
if (allocated(zoffset)) deallocate(zoffset)
if (allocated(xdomain)) deallocate(xdomain)
if (allocated(zdomain)) deallocate(zdomain)
if (allocated(zrange2D)) deallocate(zrange2D)
if (allocated(occupancy)) deallocate(occupancy)
if (allocated(cellist)) deallocate(cellist,stat=ierr)
if (ierr /= 0) then
    write(*,*) 'cellist deallocate error: ',ierr
    stop
endif
ierr = 0
if (allocated(gaplist)) deallocate(gaplist,stat=ierr)
if (allocated(DClist)) deallocate(DClist)
if (allocated(DCdeadlist)) deallocate(DCdeadlist)
if (allocated(DCvisits)) deallocate(DCvisits)
if (allocated(nz_sites)) deallocate(nz_sites)
if (allocated(nz_totsites)) deallocate(nz_totsites)
if (allocated(nz_cells)) deallocate(nz_cells)
if (allocated(nz_excess)) deallocate(nz_excess)
if (allocated(cognate_list)) deallocate(cognate_list)
if (allocated(life_dist)) deallocate(life_dist)
if (allocated(divide_dist)) deallocate(divide_dist)
if (allocated(exitlist)) deallocate(exitlist)
if (allocated(chemo_r)) deallocate(chemo_r)
if (allocated(chemo_p)) deallocate(chemo_p)
if (allocated(cytp)) deallocate(cytp)
if (allocated(xminmax)) deallocate(xminmax)
if (allocated(inblob)) deallocate(inblob)
if (allocated(sitelist)) deallocate(sitelist)
if (allocated(neighbours)) deallocate(neighbours)

! Close all open files
inquire(unit=nfout,OPENED=isopen)
if (isopen) then
	close(nfout)
	call logger('closed nfout')
endif
inquire(nfres,OPENED=isopen)
if (isopen) close(nfres)
inquire(nftraffic,OPENED=isopen)
if (isopen) close(nftraffic)
inquire(nfdcbind,OPENED=isopen)
if (isopen) close(nfdcbind)

if (par_zig_init) then
	call par_zigfree
endif

end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine terminate_run(res) BIND(C)
!DEC$ ATTRIBUTES DLLEXPORT :: terminate_run
use, intrinsic :: iso_c_binding
integer(c_int) :: res
character*(8), parameter :: quit = '__EXIT__'
integer :: error, i

call wrapup

if (res == 0) then
	call logger(' Execution successful!')
else
	call logger('  === Execution failed ===')
	call sleeper(1)
endif
close(nflog)

if (use_TCP) then
	if (stopped) then
	    call winsock_close(awp_0)
	    if (use_CPORT1) call winsock_close(awp_1)
	else
	    call winsock_send(awp_0,quit,8,error)
	    call winsock_close(awp_0)
!	    call logger("closed PORT_0")
		if (use_CPORT1) then
			call winsock_send(awp_1,quit,8,error)
			call winsock_close(awp_1)
!			call logger("closed PORT_1")
		endif
	endif
endif

end subroutine

!-----------------------------------------------------------------------------------------
! This is the DLL procedure that can be called from an external non-Fortran program to
! make a simulation run.
! Called from Python with:
!     mydll.EXECUTE(byref(ncpu),infile,n1,outfile,n2,resfile,n3,runfile,n4)
! Note that the arguments n1,n2,n3,n4, the lengths of the filename strings, are
! hidden arguments (they are not explicitly in the Fortran subroutine argument list).
! Every Python string needs to be followed by a hidden length parameter, and unless
! the declared length of a Fortran string matches the actual length of that passed from
! Python, the form character*(*) must be used.
!-----------------------------------------------------------------------------------------
subroutine execute(ncpu,infile_array,inbuflen,outfile_array,outbuflen) BIND(C)
!!DEC$ ATTRIBUTES DLLEXPORT :: EXECUTE
!!DEC$ ATTRIBUTES C, REFERENCE, MIXED_STR_LEN_ARG, ALIAS:"EXECUTE" :: execute
!DEC$ ATTRIBUTES DLLEXPORT :: execute
use, intrinsic :: iso_c_binding
character(c_char) :: infile_array(128), outfile_array(128)
integer(c_int) :: ncpu, inbuflen, outbuflen
character*(128) :: infile, outfile
logical :: ok, success
integer :: i, res

use_CPORT1 = .false.	! DIRECT CALLING FROM C++
infile = ''
do i = 1,inbuflen
	infile(i:i) = infile_array(i)
enddo
outfile = ''
do i = 1,outbuflen
	outfile(i:i) = outfile_array(i)
enddo

open(nflog,file='para.log',status='replace')
awp_0%is_open = .false.
awp_1%is_open = .false.

#ifdef GFORTRAN
    write(logmsg,'(a)') 'Built with GFORTRAN'
	call logger(logmsg)
#endif

logmsg = 'OS??'
#ifdef LINUX
    write(logmsg,'(a)') 'OS is Linux'
#endif
#ifdef OSX
    write(logmsg,'(a)') 'OS is OS-X'
#endif
#ifdef _WIN32
    write(logmsg,'(a)') 'OS is Windows'
#endif
#ifdef WINDOWS
    write(logmsg,'(a)') 'OS is Windows'
#endif
call logger(logmsg)

!#ifdef OPENMP
#if defined(OPENMP) || defined(_OPENMP)
    write(logmsg,'(a)') 'Executing with OpenMP'
	call logger(logmsg)
#endif

write(logmsg,*) 'inputfile:  ', infile
call logger(logmsg)
write(logmsg,*) 'outputfile: ', outfile 
call logger(logmsg)
!write(logmsg,*) 'resultfile: ', resfile 
!call logger(logmsg)
if (use_tcp) then
	call connecter(ok)
	if (.not.ok) then
		call logger('Failed to make TCP connections')
		return
	endif
endif
call setup(ncpu,infile,outfile,ok)
if (ok) then
	clear_to_send = .true.
	simulation_start = .true.
	istep = 0
else
	call logger('=== Setup failed ===')
endif
if (ok) then
	res = 0
else
	res = 1
endif
if (test_vascular) then
	call vascular_test
	stop
endif
return
!call terminate_run(res)

end subroutine

end module

