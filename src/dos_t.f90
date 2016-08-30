MODULE dos_type

	use parameters
	
	implicit none
	
	private
	
	type, public :: dos_t
		character(len=STRLEN) :: method
		character(len=STRLEN) :: state
		character(len=STRLEN) :: molecule
		character(len=STRLEN) :: file
		real(kind=dp) :: emin, emax, egrain
	end type dos_t

	type, public :: fcwd_t
		character(len=STRLEN) :: bra, ket
		character(len=STRLEN) :: molecule
		character(len=STRLEN) :: file
		character(len=STRLEN) :: method
		real(kind=dp) :: emin, emax, egrain
	end type fcwd_t

END MODULE dos_type
