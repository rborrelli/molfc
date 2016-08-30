module parameters

    implicit none

    public

    integer, parameter :: i4b = selected_int_kind(4)
    !integer, parameter :: dp = kind(1.0D0)
    integer, parameter :: DFTI_DPKP = SELECTED_REAL_KIND(15,307)
    integer, parameter :: DP = SELECTED_REAL_KIND(15,307)

    integer, parameter :: sp = kind(1.0)
    integer, parameter :: spc = kind((1.0,1.0))
    !integer, parameter :: dpc = kind((1.0D0,1.0D0))
    integer, parameter :: DPC = SELECTED_REAL_KIND(15,307)

    integer, parameter :: STRLEN = 100

    integer, parameter :: NDIM  = 3000,            &
    NLA   = 150,             &
    MAXAT = 200,             &
    NUMATM= 200,             &
    MAXIT = 2000,            &
    MAXEXACT = 3000,         &
    IDLEN = 20,              &
    NISOT = 1844,            &
    NELE  = 29,              &
    ATLEN = 8,               &
    SYMLEN = 4,              &
    MAXFC = 200000000,       &
    MAXNINT = 600,           &
    MAXNI = 600,             &
    MAXNC = 200,             &
    IFLEN = 30,              &		
    FINAL_STATE = 1,         &
    INITIAL_STATE = 2

    real(kind=dp), parameter ::         &
    PI=3.141592653589793238462643383279502884197_dp,     &
    PIO2=1.57079632679489661923132169163975144209858_dp, &
    TWOPI=6.283185307179586476925286766559005768394_dp,  &
    hplanck  = 6.62602957d-27,     & ! erg * s
    htplanck = 6.62602957d-27/(TWOPI),     &
    clgt     = 2.99792458d+10,  &
    clgt0    = 2.99792458_dp,   &
    htacm    = 5.3087224d-12,   &
    hpacm    = 3.3355684d-11,   &
    hacm     = 219474.6_dp,     &
    ev       = 8065.5_dp,       &
    cal      = 349.75_dp,       &
    ph       = 6.62602957d-34,  &  ! J * s
    uma      = 1.660538921d-24, &
    cfac     = (TWOPI*clgt*uma/htplanck)*1.0d-16, & ! = 0.0296596899
    one      = 1.0_dp,            &
    zero     = 0.0_dp,            &
    bohr     = 0.529177248999_dp, &
    nav      = 0.60221_dp,        &
    bk       = 1.3806503e-16_dp,  & ! erg/K 
    bkcm     = 0.6950515002_dp,   &
    cmtohr   = TWOPI*clgt0

    real(kind=dp), parameter  :: epsla = 1.0d-15, &
    epslin = 1.0d-8, &
    tm4 = 1.0d-4,    &
    tolpsi = zero

    integer, parameter ::  fout = 6, &
    finp = 10, &
    flog = 99, &
    fdos = 9, &
    fdusch = 11, &
    ftransf = 12, &
    fnst = 30, &
    fivb = 31, &
    fiup = 14, &
    fizg = 15, &
    ffcpt = 32, &
    ffcwd = 33, &
    iatc = 25, &
    inpf = 50, &
    iuspec = 19, &
    foct = 8, &
    fjm = 22, &
    fkm = 23, &
    tmpseq = 98

    character(len=1), parameter :: BLANK=" "
    character(len=1), parameter :: axslbl(1:3)=(/'x', 'y', 'z'/)

end module parameters
