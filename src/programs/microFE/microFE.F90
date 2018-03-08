PROGRAM microFE
!------------------------------------------------------------------------------
!   program xx15:  finite strain elasto-plastic analysis with Newton-Raphson
!------------------------------------------------------------------------------

  USE choose_solvers
  USE parafem_petsc
  USE PRECISION
  USE GLOBAL_VARIABLES
  USE mp_interface
  USE INPUT
  USE OUTPUT
  USE GATHER_SCATTER
  USE PARTITION
  USE MATHS
  USE TIMING
  USE LARGE_STRAIN

  USE shared_data
  USE setup
  USE core
  
  IMPLICIT NONE

  !----------------------------------------------------------------------------
  ! 1. Input and initialisation
  !----------------------------------------------------------------------------

  timest = 0.0_iwp

  timest(1) = elap_time()
  timest(2) = elap_time()

  CALL FIND_PE_PROCS(numpe,npes) ! mp_interface

  peak_memory_use = p_memory_peak() ! parafem_petsc
  IF (numpe == 1) THEN
    WRITE(*,'(A,F7.2,A)') "peak memory use at start:           ", &
      peak_memory_use, " GB"
  ENDIF

  CALL parse_user_input ! setup.F90

  CALL minimal_init ! setup.F90

  ! Setup parafem specific data structures
  CALL parafem_init ! setup.F90 

  ! Setup PETSc as needed
  IF (solvers == petsc_solvers) THEN
    CALL petsc_init ! setup.F90
  ENDIF

  CALL setup_bcs ! setup.F90

  ! Allocate and initialise equation arrays now we know how many equations
  ! this rank is dealing with.
  CALL allocate_equation_storage ! setup.F90

  peak_memory_use = p_memory_peak()
  IF (numpe == 1) WRITE(*,'(A,F7.2,A)') &
    "peak memory use after setup:        ", peak_memory_use, " GB"


  ! Limit for Newton-Raphson iterations
  limit_2=10
  
  ! Limit for the PCG iteration
  limit=12500

  ! Establish a maximum number of increments
  max_inc=100

  ! Establish the minimum time increment
  min_inc=0.001_iwp

  ! Establish the number of times the output is printed
  number_output=1
  counter_output=1
  
  ! Set the initial guess (normalized)
  initial_guess  = 0.05_iwp
  lambda         = initial_guess
  lambda_total   = 0.0_iwp
  lambda_prev    = 0.0_iwp
  next_output = 1.0_iwp/number_output
  
  peak_memory_use = p_memory_peak()

  ! Run the solver
  CALL solver ! core.F90

  IF (numpe==1) THEN
    IF (solvers == parafem_solvers) THEN
      WRITE(log_unit,'(A)') "ParaFEM revision "//REVISION
    ELSE IF (solvers == petsc_solvers) THEN
      WRITE(log_unit,'(A)') "ParaFEM revision "//REVISION//"; "//p_version()
    END IF
    WRITE(log_unit,'(A,I5,A)') "This job ran on ",npes," processors"
    WRITE(log_unit,'(A,3(I8,A))') "There are ",n_nodes," nodes",n_elem, &
      " elements and ",neq," equations"
    WRITE(log_unit,'(A,F10.4)') "Time to read input:       ", timest(30)
    WRITE(log_unit,'(A,F10.4)') "Time for setup:           ", timest(31)
    WRITE(log_unit,'(A,F10.4)') "Time for matrix assemble: ", timest(33)
    WRITE(log_unit,'(A,F10.4)') "Time for linear solve:    ", timest(34)
    WRITE(log_unit,'(A,F10.4)') "Other time in load loop:  ", timest(32)
    WRITE(log_unit,'(A,F10.4)') "Time to write results:    ", timest(35)
    WRITE(log_unit,'(A)')       "                          ----------"
    WRITE(log_unit,'(A,F10.4)') "Total:                    ", SUM(timest(30:35))
    WRITE(log_unit,'(A)')       "                          ----------"
    WRITE(log_unit,'(A,F10.4)') "This analysis took:       ", elap_time()-timest(1)
    WRITE(log_unit,*)
    WRITE(log_unit,'(A,F10.2,A)') "Peak memory use: ",peak_memory_use," GB"
    WRITE(log_unit,*)
    ! REVISION is substituted with a string like "2108" (including the double
    ! quotes) by the preprocessor, see the makefile.
    WRITE(log_unit,'(2A,2(I0,A),4(F0.2,A),F0.2)',advance='no') &
      REVISION,tab, &
      npes,tab, &
      neq,tab, &
      timest(31),tab, &
      timest(33),tab, &
      timest(34),tab, &
      timest(31)+timest(33)+timest(34),tab, &
      peak_memory_use
    IF (solvers == petsc_solvers) THEN
      WRITE(log_unit,'(2A)') tab,p_version()
    END IF
!   CALL FLUSH(log_unit)
    CLOSE(log_unit)
  END IF

  IF (numpe==1) THEN
    WRITE(*,*) 'The simulation is finished'
  END IF

!---------------------------------- shutdown ----------------------------------

  IF (solvers == petsc_solvers) THEN
    CALL p_shutdown
  END IF

  CALL SHUTDOWN()

CONTAINS
  
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

END PROGRAM microFE
