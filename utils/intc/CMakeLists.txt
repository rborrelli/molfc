add_subdirectory(colib)

add_executable(tm2fc
tm2fc.f90
)

add_executable(int2fc
int2fc.f90
)

add_executable(intc
intc1.f 
intc2.f 
intc3.f 
intc4.f 
intc5.f 
intc6.f
)

target_link_libraries(intc colib -lmkl_core -lmkl_intel -lmkl_sequential)

#install(TARGETS intc DESTINATION /Users/lello/Applications/bin)

