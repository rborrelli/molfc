#add_definitions(-DMAXMEM=200000000)

IF(Enable_LAPACK)
  set(ITRANSF_SRC itransf-lapack.f90)
  set(NATTRANSF_SRC nattransf-lapack.f90)
ELSE(Enable_LAPACK)
  set(ITRANSF_SRC itransf-mkl.f90)
  set(NATTRANSF_SRC nattransf-mkl.f90)
ENDIF(Enable_LAPACK)

IF(Enable_FFTW3)
  set(KUBO_SRC kubo-fftw.f90)
ELSE(Enable_FFTW3)
  set(KUBO_SRC kubo-mkl.f90)
ENDIF(Enable_FFTW3)

add_executable(molfc3
active.f90
agric.f90
dos_t.f90
errors.f90
fc_t.f90
hypfun.f90
input.f90
intc.f90
iofiles.f90
iogeom.f90
iomatrix.f90
iomodes.f90
isotopes.f90
job_t.f90
matfun.f90
messages.f90
mkl_dfti.f90
molfc.F90
output.F90
parameters.f90
proc.f90
proc_t.f90
sysop.f90
sys_t.f90
tmderivative.f90
transf_t.f90
util.f90
vecm.f90
vibut.f90
xmvar.f90

fc/axsw.f90
fc/ctransf.f90
fc/fc.f90
fc/fc_h.f90
fc/fcint.f90
fc/fcint_h.f90
fc/fcio.f90
fc/${ITRANSF_SRC}
fc/mtransf.f90
fc/${NATTRANSF_SRC}
fc/spectrum.f90
fc/${KUBO_SRC}
fc/reorg.f90
fc/fcpt.f90

dos/dos.f90
)

link_directories(${CMAKE_BUILD_DIR}/lib)

IF(Enable_MKL) 
  message("Linking with MKL")
  # For 64bit machines
  target_link_libraries(molfc3 slatec fox_dom mkl_lapack95_lp64 mkl_intel_lp64 mkl_intel_thread mkl_core iomp5 pthread)
ELSE(Enable_MKL)
  message("Linking without MKL")
  target_link_libraries(molfc3 slatec fftw3 fox_dom lapack blas pthread)
ENDIF(Enable_MKL) 

install(TARGETS molfc3 DESTINATION /Users/lello/Applications/bin)

