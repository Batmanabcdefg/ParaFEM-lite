MODULE setup

  ! from ParaFEM/PETSc
  USE parafem_petsc
  USE input
  USE mp_interface 
  USE gather_scatter
  USE large_strain
  USE timing

  ! From microFE
  USE shared_data

  PRIVATE

  PUBLIC :: allocate_equation_storage, minimal_init, parse_user_input, &
    parafem_init, petsc_init, setup_bcs

CONTAINS

  SUBROUTINE allocate_equation_storage
    ! Allocate and initialise equation arrays now we know 
    ! how many equations this rank is dealing with.

  ALLOCATE(r_pp(0:neq_pp), xnew_pp(0:neq_pp), diag_precon_pp(0:neq_pp))
  ALLOCATE(res_pp(0:neq_pp), deltax_pp(0:neq_pp), fint_pp(0:neq_pp)) 
  ALLOCATE(deltax_pp_temp(0:neq_pp), xnew_pp_previous(0:neq_pp))
  ALLOCATE(xnewprev_pp(0:neq_pp))
  ALLOCATE(diag_precon_tmp(ntot,nels_pp))

  r_pp             = 0.0_iwp
  xnew_pp          = 0.0_iwp
  diag_precon_pp   = 0.0_iwp
  res_pp           = 0.0_iwp
  deltax_pp        = 0.0_iwp
  fint_pp          = 0.0_iwp
  xnew_pp_previous = 0.0_iwp
  xnewprev_pp      = 0.0_iwp

  END SUBROUTINE allocate_equation_storage


  SUBROUTINE minimal_init

    IF (n_elem < npes) THEN
      IF (numpe==1) THEN
        WRITE(*,*)"Error: fewer elements than processors"
      END IF
      CALL SHUTDOWN
      STOP
    END IF

    solvers = get_solvers()
    IF (.NOT. solvers_valid(solvers)) THEN
      CALL SHUTDOWN
    END IF


    IF(numpe==1) THEN
      fname = fname_base(1:INDEX(fname_base," ")-1) // ".res"
      OPEN (newunit=log_unit, file=fname, status='replace', action='write')
    END IF

    IF (n_nodpel==8) THEN
      element='hexahedron'
      dimH=8
    ELSE IF (n_nodpel==4) THEN
      element='tetrahedron'
      dimH=4
    END IF
    
    ntot = n_nodpel * n_dof
  
  END SUBROUTINE minimal_init


  SUBROUTINE parse_user_input
  
    INTEGER :: argc

  ! Assert filename has been given
  argc = COMMAND_ARGUMENT_COUNT()
  IF (argc < 1) THEN
    IF (numpe == 1) THEN
      WRITE(*,*) "Need name of filename_base!!"
    END IF
    CALL SHUTDOWN
    STOP
  END IF

  ! Get the passed filename
  CALL GET_COMMAND_ARGUMENT(1,fname_base)

  fname = fname_base(1:INDEX(fname_base," ")-1) // ".dat"
  CALL READ_DATA_XX7(fname,numpe,n_elem,n_nodes,n_restrained,loaded_nodes,fixed_nodes, &
                     n_intp,limit,tol,e,v,n_nodpel,num_load_steps,jump,tol2,partitioner)

  END SUBROUTINE parse_user_input

  
  SUBROUTINE parafem_init

    !----------------------------------------------------------------------------
    ! 1. Calculate subdivision over available ranks
    !----------------------------------------------------------------------------

    CALL CALC_NODES_PP(n_nodes,npes,numpe,node_end,node_start,nodes_pp)

    !------------------------------------------------------------------------------
    ! 2. Get integration Gauss points and weights in the element
    !------------------------------------------------------------------------------
      
    ALLOCATE(points(ndim,n_intp), weights(n_intp))
    CALL GET_GAUSS_POINTS(element,points,weights)

    !------------------------------------------------------------------------------
    ! 3. Import and distribute mesh
    !------------------------------------------------------------------------------

    CALL CALC_NELS_PP(fname_base,n_elem,npes,numpe,partitioner,nels_pp)
    
    ALLOCATE(g_num_pp(n_nodpel, nels_pp)) 

    fname = fname_base(1:INDEX(fname_base," ")-1) // ".d"
    CALL READ_G_NUM_PP(fname_base,iel_start,n_nodes,npes,numpe,g_num_pp)
    ! read_elements_2() does not work for METIS partitions
    !CALL READ_ELEMENTS_2(fname,npes,n_nodes,numpe,g_num_pp)

    CALL CALC_NN_PP(g_num_pp,nn_pp,nn_start)

    ALLOCATE(g_coord_pp(ndim, nn_pp))
    CALL READ_NODES(fname,n_nodes,nn_start,numpe,g_coord_pp)

    ALLOCATE(rest(n_restrained,n_dof+1))
    fname = fname_base(1:INDEX(fname_base," ")-1) // ".bnd"
    CALL READ_RESTRAINTS(fname,numpe,rest)

    timest(30) = timest(30) + elap_time()-timest(2) ! 30 = read
    timest(2) = elap_time()

  !------------------------------------------------------------------------------
  ! 4. Allocate dynamic arrays used in the main program  
  !------------------------------------------------------------------------------

    ALLOCATE(coord(n_nodpel,ndim))
    ALLOCATE(upd_coord(n_nodpel,ndim))
    ALLOCATE(bee(nst,ntot))
    ALLOCATE(num(n_nodpel))
    ALLOCATE(g_g_pp(ntot,nels_pp))
    ALLOCATE(load_value(ndim,loaded_nodes))
    ALLOCATE(load_node(loaded_nodes))
    ALLOCATE(nf_pp(n_dof,nn_pp))
    IF (solvers == parafem_solvers) THEN
      ALLOCATE(storekm_pp(ntot,ntot,nels_pp))
    END IF
    ALLOCATE(kmat_elem(ntot,ntot))
    ALLOCATE(xnewel_pp(ntot,nels_pp))
    ALLOCATE(xnewel_pp_previous(ntot,nels_pp))
    ALLOCATE(comp(n_nodpel,ndim))
    ALLOCATE(jacF(ndim,ndim))
    ALLOCATE(jacFinc(ndim,ndim))
    ALLOCATE(lnstrainelas_mat_con(nels_pp,n_intp,nst))
    ALLOCATE(lnstrainelas(nst))
    ALLOCATE(statev(statevar_num))
    ALLOCATE(statevar_con(nels_pp,n_intp,statevar_num))
    ALLOCATE(sigma(ndim,ndim))
    ALLOCATE(sigma1C_mat_con(nels_pp,n_intp,nst))
    ALLOCATE(auxm(n_nodpel,ndim))
    ALLOCATE(jacFinv(ndim,ndim))
    ALLOCATE(derivFtran(n_nodpel,ndim))
    ALLOCATE(derivF(ndim,n_nodpel))
    ALLOCATE(beeF(nst,ntot))
    ALLOCATE(geeF((ndim*ndim),ntot))
    ALLOCATE(geeFT(ntot,(ndim*ndim)))
    ALLOCATE(sigma1C(nst))
    ALLOCATE(storefint_pp(ntot,nels_pp))
    ALLOCATE(deeF((ndim*ndim),(ndim*ndim)))
    ALLOCATE(principal(ndim))
    ALLOCATE(value_shape(n_nodpel))
    ALLOCATE(xnewelinc_pp(ntot,nels_pp))
    ALLOCATE(auxm_inc(n_nodpel,ndim))
    ALLOCATE(auxm_previous(n_nodpel,ndim))
    ALLOCATE(stiffness_mat_con(nels_pp,n_intp,(ndim*ndim),(ndim*ndim)))
    ALLOCATE(unload_pp(nels_pp,n_intp))
    ALLOCATE(rm_vec(npes))
    ALLOCATE(km(ntot,ntot))

    ! Vector comp to compute F (gradient of deformation)
    DO i = 1,n_nodpel
      DO j = 1,ndim
        comp(i,j) = (i-1)*ndim + j
      END DO
    END DO

    ! Initialise the converged stresses, elastic strains and the state variables
    sigma1C_mat_con      = 0.0_iwp
    lnstrainelas_mat_con = 0.0_iwp
    statevar_con         = 0.0_iwp
    stiffness_mat_con    = 0.0_iwp

    !------------------------------------------------------------------------------
    ! 5. Loop the elements to find the steering array and the number of 
    !    equations to solve. 
    !------------------------------------------------------------------------------

    CALL REST_NF(nn_start,rest,nf_pp)

    ! Assign global equation numbers to each dof of each element
    DO iel = 1,nels_pp
      CALL NUM_TO_G2(g_num_pp(:,iel), nf_pp, g_g_pp(:,iel), nn_start)
    END DO

    CALL CALC_NEQ(n_nodes,rest,neq)

    !------------------------------------------------------------------------------
    ! 6. Create interprocessor communications tables
    !------------------------------------------------------------------------------  
    
    CALL CALC_NPES_PP(npes,npes_pp)

    CALL CALC_NEQ_PP

    CALL MAKE_GGL(npes_pp,npes,g_g_pp)

    timest(31) = timest(31) + elap_time()-timest(2) ! 31 = setup
    timest(2) = elap_time()

  END SUBROUTINE parafem_init


  SUBROUTINE petsc_init
    ! Initialize PETSc and create the needed equation system

    CALL p_initialize(fname_base,error)
    IF (error) THEN
      CALL shutdown
    END IF
    CALL p_create(ntot,g_g_pp,error)
    IF (error) THEN
      CALL p_finalize
      CALL shutdown
    END IF

  END SUBROUTINE petsc_init


  SUBROUTINE setup_bcs

    CALL setup_fixed_bcs
    CALL setup_load_bcs

  END SUBROUTINE setup_bcs


  SUBROUTINE setup_fixed_bcs

    !------------------------------------------------------------------------------
    ! 1. Read and distribute essential boundary conditions
    !------------------------------------------------------------------------------

    numfix_pp = 0

    IF (fixed_nodes>0) THEN

      ALLOCATE(fixed_node(fixed_nodes))
      ALLOCATE(fixed_dof(fixed_nodes))
      ALLOCATE(fixed_value(fixed_nodes))

      CALL READ_FIXED(fname_base,numpe,fixed_node,fixed_dof,fixed_value)

      IF (element=='hexahedron') THEN
        fixdim=4
      ELSE IF (element=='tetrahedron') THEN
        fixdim=20
      END IF

      ALLOCATE(fixelem_pp(fixdim*fixed_nodes))
      ALLOCATE(fixdof_pp(fixdim*fixed_nodes))
      ALLOCATE(fixval_pp(fixdim*fixed_nodes))
      ALLOCATE(fixvalpiece_pp(fixdim*fixed_nodes))
      ALLOCATE(fixvalprev_pp(fixdim*fixed_nodes))
      ALLOCATE(fixvaltot_pp(fixdim*fixed_nodes))
      
      fixelem_pp      = 0
      fixdof_pp       = 0
      fixval_pp       = zero
      fixvalpiece_pp  = zero
      fixvalprev_pp   = zero
      fixvaltot_pp    = zero
      
      DO i = 1,fixed_nodes
        DO iel = 1,nels_pp
          DO j = 1,n_nodpel
            IF (fixed_node(i)==g_num_pp(j,iel)) THEN
              numfix_pp = numfix_pp + 1
              fixelem_pp(numfix_pp) = iel
              fixdof_pp(numfix_pp) = (j-1)*ndim + fixed_dof(i)
              fixval_pp(numfix_pp) = fixed_value(i)
            END IF
          END DO
        END DO
      END DO
      
      DEALLOCATE(fixed_node, fixed_dof, fixed_value)

      IF (solvers == petsc_solvers) THEN
        ! storekm_pp is only used by the ParaFEM solvers but would also be needed
        ! to convert the fixed displacements to a force.  Use a (hopefully small)
        ! copy of the elements needed for the fixed displacements.  It would
        ! better to use PETSc matrix-vector multiplies but that would require all
        ! the handling of fixed dofs in ParaFEM to be changed (right down to
        ! changes in gather and scatter I think) to pass them through as zeroes in
        ! the matrix.
        ALLOCATE(fixkm_pp(ntot,ntot,numfix_pp))
      END IF

    END IF

  END SUBROUTINE setup_fixed_bcs


  SUBROUTINE setup_load_bcs

    !----------------------------------------------------------------------------
    ! 2. Read and distribute natural boundary conditions
    !----------------------------------------------------------------------------

    ALLOCATE(fextpiece_pp(0:neq_pp))
    ALLOCATE(fext_pp(0:neq_pp))
    ALLOCATE(fextprev_pp(0:neq_pp))
    
    fextpiece_pp = zero
    fext_pp      = zero
    fextprev_pp  = zero

    IF (loaded_nodes>0) THEN

      CALL READ_LOADS(fname_base,numpe,load_node,load_value)

      CALL LOAD_2(nn_start,g_num_pp,load_node,load_value,nf_pp,fext_pp(1:))
     
      DEALLOCATE(load_node,load_value)

    END IF

    timest(30) = timest(30) + elap_time()-timest(2) ! 30 = read
    timest(2) = elap_time()

  END SUBROUTINE setup_load_bcs

END MODULE setup
