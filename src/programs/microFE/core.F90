MODULE core

  ! ParaFEM
  USE parafem_petsc
  USE timing
  USE mp_interface
  USE gather_scatter
  USE large_strain
  USE output

  ! microFE
  USE shared_data
  USE linalg
  USE plastic

  IMPLICIT NONE
  
  PUBLIC:: solver

  PRIVATE

CONTAINS

  SUBROUTINE solver
  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  !----------------------- Start Load LOOP --------------------------------------
  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------

    timest(31) = timest(31) + elap_time()-timest(2) ! 31 = setup
    timest(2) = elap_time()

    DO iload = 1,max_inc
    
      converged = .FALSE.

      timest(2) = elap_time()

      ! Increase the stage control by 50% if the number of iterations of the two 
      ! previously converged full increments (not affected by the output  
      ! requests) are below 6
      IF ((iter<6).AND.(prev_iter<6).AND.(iload>2)) THEN
        IF(numpe==1) THEN
          WRITE(*,*) 'The load increment is increased by 50%'
        END IF
        lambda=1.5*lambda
      END IF

      timest(32) = timest(32) + elap_time()-timest(2) ! 32 = other work in load loop
      timest(2) = elap_time()

      100 CONTINUE

      timest(2) = elap_time()

      ! Special case: If lambda is larger than 0.25, the simulation is finished
      IF (lambda>=(1.00_iwp-tol_increment)) THEN
  !!$    ! If lambda is larger than unity, the simulation is finished
  !!$    IF (lambda_total>=(1._iwp-tol_increment)) THEN

        timest(32) = timest(32) + elap_time()-timest(2) ! 32 = other work in load loop
        timest(2) = elap_time()
        
        exit_iload = .TRUE.
        GOTO 400
      END IF

      ! Update the load increments
      lambda_prev=lambda_total
      lambda_total=lambda_total+lambda
      
      ! Special case: If the total lambda is larger than 0.25, cut it off
      IF (lambda_total>=(1.00_iwp-tol_increment)) THEN
        lambda_total=1.00_iwp
  !!$    ! If the total lambda is larger than 1.0_iwp, cut it off
  !!$    IF (lambda_total>=(1._iwp-tol_increment)) THEN
  !!$      lambda_total=1._iwp
        lambda=lambda_total-lambda_prev
      END IF

      ! Exit if the time increment is less that the specified minimum
      IF (lambda<min_inc) THEN
        IF(numpe==1) THEN
          WRITE(*,*) 'The load increment is too small'
        END IF

        timest(32) = timest(32) + elap_time()-timest(2) ! 32 = other work in load loop
        timest(2) = elap_time()
        
        exit_iload = .TRUE.
        GOTO 400
      END IF

      ! Display the incremental, total and previous load increments
      IF(numpe==1) THEN
        WRITE(*,*) 'For the increment ',iload,' the data is as following:'
        WRITE(*,*) 'The current load increment is ',lambda
        WRITE(*,*) 'The previous total load multiplier is ',lambda_prev
        WRITE(*,*) 'The total load multiplier is ',lambda_total
      END IF

      ! Calculate the displacement increment, total displacement, and previously 
      ! converged displacement
      fixvalpiece_pp(1:)=fixval_pp(1:)*lambda
      fixvaltot_pp(1:)=fixval_pp(1:)*lambda_total
      fixvalprev_pp(1:)=fixval_pp(1:)*lambda_prev

  !------------------------------------------------------------------------------
  !----------------------- Start Newton-Raphson iterations ----------------------
  !------------------------------------------------------------------------------

      deltax_pp_temp = 0.0_iwp
      inewton = 0
      
      timest(32) = timest(32) + elap_time()-timest(2) ! 32 = other work in load loop
      timest(2) = elap_time()

      iterations: DO

        timest(2) = elap_time()

        inewton = inewton + 1

        storefint_pp = 0.0_iwp

        CALL GATHER(xnew_pp(1:),xnewel_pp)
        CALL GATHER(deltax_pp_temp(1:),xnewelinc_pp)
        CALL GATHER(xnew_pp_previous(1:),xnewel_pp_previous)

        DO i = 1,numfix_pp
          xnewel_pp_previous(fixdof_pp(i),fixelem_pp(i))=fixvalprev_pp(i)
        END DO

        IF (inewton>1 .AND. numfix_pp>0) THEN
          DO i = 1,numfix_pp
            xnewel_pp(fixdof_pp(i),fixelem_pp(i)) = fixvaltot_pp(i)
          END DO
          DO i = 1,numfix_pp
            xnewelinc_pp(fixdof_pp(i),fixelem_pp(i)) = fixvalpiece_pp(i)
          END DO
        END IF

        IF (iload>1 .AND. inewton==1 .AND. numfix_pp>0) THEN
          DO i = 1,numfix_pp
            xnewel_pp(fixdof_pp(i),fixelem_pp(i))=fixvalprev_pp(i)
          END DO
        END IF      

  !------------------------------------------------------------------------------
  ! 11. Element stiffness integration and storage
  !------------------------------------------------------------------------------
        
        timest(32) = timest(32) + elap_time()-timest(2) ! 32 = other work in load loop
        timest(2) = elap_time()
        
        IF (solvers == parafem_solvers) THEN
          storekm_pp = 0.0_iwp
        ELSE IF (solvers == petsc_solvers) THEN
          CALL p_zero_matrix
        END IF
        
        DO iel = 1,nels_pp

          km = 0.0_iwp

          DO i = 1,n_nodpel
            num(i) = g_num_pp(i,iel) - nn_start + 1
          END DO

          DO i = 1,ndim
            auxm(:,i) = xnewel_pp(comp(:,i),iel)
          END DO

          DO i = 1,ndim
            auxm_inc(:,i) = xnewelinc_pp(comp(:,i),iel)
          END DO

          DO i = 1,ndim
            auxm_previous(:,i) = xnewel_pp_previous(comp(:,i),iel)
          END DO 
          
          coord = TRANSPOSE(g_coord_pp(:,num))
          upd_coord=coord+auxm_previous
          
          DO igauss = 1,n_intp

            ! Initialise the state variables to the same converged value of     
            ! the last time increment
            statev(:)=statevar_con(iel,igauss,:)
            lnstrainelas(:)=lnstrainelas_mat_con(iel,igauss,:)

            ! Calculates the total and incremental deformation gradients
            CALL DEFGRA(igauss,auxm,coord,points,det,detF,beeF,geeF,jacF,ndim, &
             n_nodpel)
            CALL DEFGRAINC(igauss,auxm_inc,upd_coord,points,jacFinc,ndim,n_nodpel)
            
            CALL PLASTICITY(deeF,jacF,jacFinc,sigma1C,statev,lnstrainelas,sigma, &
                            detF,statevar_num,iel,igauss,noncon_flag, &
                            umat_quadric_linear_hardening)
             
            ! Exit if the local Newton-Raphson has not converged
            IF (noncon_flag) THEN
              GOTO 200
            END IF
             
            ! During the first Newton-Raphson iteration, retrieve the previous  
            ! converged stress and stiffness tensor
            IF ((iload>1).AND.(inewton==1)) THEN
              deeF(:,:)=stiffness_mat_con(iel,igauss,:,:)
              sigma1C(:)=sigma1C_mat_con(iel,igauss,:)
            END IF
            
            dw = det*weights(igauss)
            
            ! Calculate the internal force vector of the element
            storefint_pp(:,iel)=storefint_pp(:,iel) + MATMUL(TRANSPOSE(beeF), &
             sigma1C)*dw
             
            ! Calculate the stiffness tensor of the element
            geeFT = TRANSPOSE(geeF)
            km = km + (MATMUL(MATMUL(geeFT,deeF),geeF)*dw)

          END DO

          IF (solvers == parafem_solvers) THEN
            storekm_pp(:,:,iel) = km
          ELSE IF (solvers == petsc_solvers) THEN
            CALL p_add_element(g_g_pp(:,iel),km)
            FORALL (i = 1:numfix_pp, fixelem_pp(i) == iel) fixkm_pp(:,:,i) = km
          END IF
        END DO

        200 CONTINUE

        peak_memory_use = p_memory_peak()
        IF (numpe == 1) WRITE(*,'(A,F7.2,A)') &
          "peak memory use after add elements: ", peak_memory_use, " GB"

        ! When using PETSc, even if local convergence fails, the matrix has to be
        ! assembled to be in the correct state for p_0.0_iwp_matrix().
        IF (solvers == petsc_solvers) THEN
          CALL p_assemble
        END IF

        peak_memory_use = p_memory_peak()
        IF (numpe == 1) WRITE(*,'(A,F7.2,A)') &
          "peak memory use after assemble:     ", peak_memory_use, " GB"

        timest(33) = timest(33) + elap_time()-timest(2) ! 33 = matrix assemble
        timest(2) = elap_time()
        
        ! This code deals with the fail of local convergence
        ! ---------------------------------------------------------------------------

        ! Call a barrier to ensure all processes are synchronised
        CALL MPI_BARRIER(MPI_COMM_WORLD,ier)
        
        IF (noncon_flag) THEN
          noncon_flag = .FALSE.
          rm_vec(numpe) = 1
        ELSE
          rm_vec(numpe) = 0
        END IF
              
        rm_sum=0
        CALL MPI_ALLREDUCE(rm_vec(numpe),rm_sum,1,MPI_INTEGER,MPI_SUM, &
         MPI_COMM_WORLD,ier)
        rm_vec(numpe)=0

        IF (rm_sum.GE.1) THEN
          
          IF (numpe==1) WRITE(*,*) 'Fail of local convergence - The load ', &
           'increment is cut in 0.5_iwp'
               
          IF (print_output) THEN
            lambda_total=lambda_total-temp_inc
          ELSE
            lambda_total=lambda_total-lambda
          END IF
              
          lambda=0.5_iwp*lambda
          xnew_pp=xnewprev_pp
          flag=.FALSE.
              
          IF (print_output) THEN
            counter_output=counter_output-1
            next_output=FLOAT(counter_output)/number_output
            print_output=.FALSE.
          END IF
          
          ! Set the number of iterations to a number bigger than six to reset the  
          ! iteration counter
          prev_iter=7
          iter=7

          timest(33) = timest(33) + elap_time()-timest(2) ! 33 = matrix assemble
          timest(2) = elap_time()

          GOTO 100
        END IF
        ! ---------------------------------------------------------------------------

        timest(33) = timest(33) + elap_time()-timest(2) ! 33 = matrix assemble
        timest(2) = elap_time()
        
  !------------------------------------------------------------------------------
  ! 12. Build and invert the preconditioner (ParaFEM only)
  !------------------------------------------------------------------------------

        IF (solvers == parafem_solvers) THEN
          diag_precon_tmp = 0.0_iwp
          
          DO iel = 1,nels_pp
            DO k = 1,ntot 
              diag_precon_tmp(k,iel)=diag_precon_tmp(k,iel) + storekm_pp(k,k,iel)
            END DO
          END DO
          
          diag_precon_pp(:) = 0.0_iwp
          CALL SCATTER(diag_precon_pp(1:),diag_precon_tmp)
          
          diag_precon_pp(1:) = 1.0_iwp/diag_precon_pp(1:)
          diag_precon_pp(0)  = 0.0_iwp
        END IF

        timest(34) = timest(34) + elap_time()-timest(2) ! 34 = solve
        timest(2) = elap_time()

  !------------------------------------------------------------------------------
  ! 13. Initialize PCG
  !------------------------------------------------------------------------------
        
        ! During the first Newton-Raphson iteration, the incremental 
        ! displacements are applied through a linear mapping of these 
        ! displacements to the internal force vector
        IF (inewton==1 .AND. numfix_pp>0) THEN
          DO i = 1,numfix_pp
            DO j = 1,ntot
              IF (g_g_pp(j,fixelem_pp(i))>0) THEN
                IF (solvers == parafem_solvers) THEN
                  storefint_pp(j,fixelem_pp(i)) = storefint_pp(j,fixelem_pp(i)) &
                    + fixvalpiece_pp(i)*storekm_pp(j,fixdof_pp(i),fixelem_pp(i))
                ELSE IF (solvers == petsc_solvers) THEN
                  storefint_pp(j,fixelem_pp(i)) = storefint_pp(j,fixelem_pp(i)) &
                    + fixvalpiece_pp(i)*fixkm_pp(j,fixdof_pp(i),i)
                END IF
              END IF
            END DO
          END DO
        END IF
        
        fint_pp = 0.0_iwp
        CALL SCATTER(fint_pp(1:),storefint_pp)

        r_pp(1:) = fext_pp(1:) - fint_pp(1:)       !Residual
        r_pp(0) = 0.0_iwp

        ! Compute maxdiff of residual 
        maxdiff =  MAXABSVAL_P(r_pp(1:))

        ! Normalise residual vector and stiffness matrix for pcg
        IF (maxdiff == 0.0_iwp) THEN
          IF(numpe==1) THEN
            WRITE(*,*) "maxdiff = 0.0_iwp and now exiting loop"
          END IF

          timest(32) = timest(32) + elap_time()-timest(2) ! 32 = other work in load loop
          timest(2) = elap_time()
        
          EXIT
        END IF

  !------------------------------------------------------------------------------
  !----------------- Solve using preconditioned Krylov solver -------------------
  !------------------------------------------------------------------------------
        
        deltax_pp = 0.0_iwp
        res_pp    = r_pp

        timest(32) = timest(32) + elap_time()-timest(2) ! 32 = other work in load loop
        timest(2) = elap_time()

        IF (solvers == parafem_solvers) THEN
          CALL PCG_VER1(inewton,limit,tol,storekm_pp,r_pp(1:), &
                        diag_precon_pp(1:),rn0,deltax_pp(1:),iters)
        ELSE IF (solvers == petsc_solvers) THEN
          CALL p_use_solver(1,error)
          IF (error) THEN
            CALL p_shutdown
            CALL shutdown
          END IF
          CALL p_solve(r_pp(1:),deltax_pp(1:))
        END IF

        timest(34) = timest(34) + elap_time()-timest(2) ! 34 = solve
        timest(2) = elap_time()

        IF (solvers == petsc_solvers) THEN
          CALL p_print_info(1,log_unit)
        END IF

        timest(32) = timest(32) + elap_time()-timest(2) ! 32 = other work in load loop
        timest(2) = elap_time()

        IF (numpe==1) THEN
          WRITE(91,*)iload,inewton,iters
          CALL FLUSH(91)
        END IF

        timest(35) = timest(35) + elap_time()-timest(2) ! 35 = write
        timest(2) = elap_time()

        ! Total displacements
        xnew_pp(1:) = xnew_pp(1:) + deltax_pp(1:)
        xnew_pp(0) = 0.0_iwp
        
        ! Incremental displacements corresponding to the time increment
        deltax_pp_temp(1:) = deltax_pp_temp(1:) + deltax_pp(1:)
        deltax_pp_temp(0) = 0.0_iwp

  !---------------------------------------------------------------------------
  !---------------------------------------------------------------------------
  !---------------------------------------------------------------------------

  ! Check convergence for Newton-Raphson iterations 
        energy_prev_prev=energy_prev
        energy_prev=energy
        energy = ABS(DOT_PRODUCT_P(res_pp(1:),deltax_pp(1:)))
        IF (inewton==1) THEN
          energy1 = energy
        END IF
        IF (numpe==1) THEN
          WRITE(*,energy_fmt) ' Energy  ',iload,inewton,energy,energy/energy1
        END IF
        IF (inewton>1) THEN
          IF ((energy/energy1)<=tol2) THEN
            converged = .TRUE.
          END IF 
        END IF

        ! The time increment is cut in 0.5_iwp if 1.0_iwp of the following situations
        ! happen:
        ! - The maximum number of iterations has been met
        ! - The current iteration diverged (the current energy is larger than
        ! the two previous energies)
        ! - There is some numerical instability which causes some NaN values
        IF ((inewton==limit_2).OR.(((energy_prev_prev/energy1)<(energy_prev/ &
         energy1)).AND.((energy_prev/energy1)<(energy/energy1)).AND. &
          (inewton>2)).OR.(ISNAN(energy))) THEN
          
          IF (numpe==1) THEN
            WRITE(*,*) 'The load increment is cut in 0.5_iwp'
          END IF
          IF (print_output) THEN
            lambda_total=lambda_total-temp_inc
          ELSE
            lambda_total=lambda_total-lambda
          END IF
          lambda=0.5_iwp*lambda
          xnew_pp=xnewprev_pp
          flag=.FALSE.
          IF (print_output) THEN
            counter_output=counter_output-1
            next_output=FLOAT(counter_output)/number_output
            print_output=.FALSE.
          END IF
          
          ! Set the number of iterations to a number bigger than six to reset the  
          ! iteration counter
          prev_iter=7
          iter=7

          timest(32) = timest(32) + elap_time()-timest(2) ! 32 = other work in load loop
          timest(2) = elap_time()

          GOTO 100
        END IF
   
        ! After convergence, a last "iteration" is needed to calculate the fully
        ! converged values (some FE algorithms omit this step, but we perform
        ! it because we are interested in a fully precise solution)
        IF (converged) THEN 

          ! If the current increment is not an output request increment, save the  
          ! current and previous increment number of iterations to allow for   
          ! for lambda increment to be increased
          IF (.NOT.(lambda_inc)) THEN
            IF (iload>1) THEN
              prev_iter=iter
            END IF
            iter=inewton
            lambda_inc=.TRUE.
          END IF

          xnewprev_pp=xnew_pp

          CALL GATHER(xnew_pp(1:),xnewel_pp)
          CALL GATHER(deltax_pp_temp(1:),xnewelinc_pp)
          CALL GATHER(xnew_pp_previous(1:),xnewel_pp_previous)

          DO i = 1,numfix_pp
            xnewel_pp_previous(fixdof_pp(i),fixelem_pp(i))=fixvalprev_pp(i)
          END DO
          
          DO i = 1,numfix_pp
            xnewel_pp(fixdof_pp(i),fixelem_pp(i)) = fixvaltot_pp(i)
          END DO

          DO i = 1,numfix_pp
            xnewelinc_pp(fixdof_pp(i),fixelem_pp(i)) = fixvalpiece_pp(i)
          END DO
          
          DO iel = 1,nels_pp

            DO i = 1,n_nodpel
              num(i) = g_num_pp(i,iel) - nn_start + 1
            END DO

            DO i = 1,ndim
              auxm(:,i) = xnewel_pp(comp(:,i),iel)
            END DO

            DO i = 1,ndim
              auxm_inc(:,i) = xnewelinc_pp(comp(:,i),iel)
            END DO

            DO i = 1,ndim
              auxm_previous(:,i) = xnewel_pp_previous(comp(:,i),iel)
            END DO 
            
            coord = TRANSPOSE(g_coord_pp(:,num))
            upd_coord=coord+auxm_previous

            DO igauss = 1,n_intp

              statev(:)=statevar_con(iel,igauss,:)
              lnstrainelas(:)=lnstrainelas_mat_con(iel,igauss,:)

              CALL DEFGRA(igauss,auxm,coord,points,det,detF,beeF,geeF,jacF,ndim, &
               n_nodpel)

              CALL DEFGRAINC(igauss,auxm_inc,upd_coord,points,jacFinc,ndim,n_nodpel)

              CALL PLASTICITY(deeF,jacF,jacFinc,sigma1C,statev,lnstrainelas, &
                              sigma,detF,statevar_num,iel,igauss,noncon_flag, &
                              umat_quadric_linear_hardening)
               
              ! Exit if the local Newton-Raphson has not converged
              IF (noncon_flag) THEN
                GOTO 300
              END IF
               
              ! Save the variables
              statevar_con(iel,igauss,:)=statev(:)
              lnstrainelas_mat_con(iel,igauss,:)=lnstrainelas(:)
              sigma1C_mat_con(iel,igauss,:)=sigma1C(:)
              stiffness_mat_con(iel,igauss,:,:)=deeF(:,:)
            END DO
          END DO
          
          300 CONTINUE

          ! This code deals with the fail of local convergence
          ! ---------------------------------------------------------------------------
        
          ! Call a barrier to ensure all processes are synchronised
          CALL MPI_BARRIER(MPI_COMM_WORLD,ier)
        
          IF (noncon_flag) THEN
            noncon_flag = .FALSE.
            rm_vec(numpe) = 1
          ELSE
            rm_vec(numpe) = 0
          END IF
              
          rm_sum=0
          CALL MPI_ALLREDUCE(rm_vec(numpe),rm_sum,1,MPI_INTEGER,MPI_SUM, &
           MPI_COMM_WORLD,ier)
          rm_vec(numpe)=0

          IF (rm_sum.GE.1) THEN
          
            IF (numpe==1) WRITE(*,*) 'Fail of local convergence - The load ', &
             'increment is cut in 0.5_iwp'
               
            IF (print_output) THEN
              lambda_total=lambda_total-temp_inc
            ELSE
              lambda_total=lambda_total-lambda
            END IF
              
            lambda=0.5_iwp*lambda
            xnew_pp=xnewprev_pp
            flag=.FALSE.
              
            IF (print_output) THEN
              counter_output=counter_output-1
              next_output=FLOAT(counter_output)/number_output
              print_output=.FALSE.
            END IF
          
            ! Set the number of iterations to a number bigger than six to reset the  
            ! iteration counter
            prev_iter=7
            iter=7

            timest(32) = timest(32) + elap_time()-timest(2) ! 32 = other work in load loop
            timest(2) = elap_time()
                
            GOTO 100
          END IF
          ! ---------------------------------------------------------------------------
          
          ! Save the converged displacement
          xnew_pp_previous=xnew_pp

          timest(32) = timest(32) + elap_time()-timest(2) ! 32 = other work in load loop
          timest(2) = elap_time()

          EXIT
        END IF

        timest(32) = timest(32) + elap_time()-timest(2) ! 32 = other work in load loop
        timest(2) = elap_time()

      END DO iterations

      timest(2) = elap_time()

  !------------------------------------------------------------------------------
  !------------------------- End Newton-Raphson iterations ----------------------
  !------------------------------------------------------------------------------

      IF (numpe==1) THEN
        WRITE(log_unit,'(a,i3,a,f12.4,a,i4,a)') "Time after load step ",iload,": ", &
        ELAP_TIME() - timest(1),"     (",inewton," iterations )"
        FLUSH(log_unit)
      END IF

      !IF (iload==1) THEN
      !  ALLOCATE(disp_pp(ndim,nn_pp))
      !END IF

      !IF (iload==num_load_steps) THEN
      !  DEALLOCATE(diag_precon_tmp)
      !END IF

  !------------------------------------------------------------------------------
  !-----------------------------print out results -------------------------------
  !------------------------------------------------------------------------------

      IF (numpe==1) THEN
        IF (iload==1) THEN
          ! Displacement and total strain
          fname = fname_base(1:INDEX(fname_base, " ")-1) // "_dis.res"
          OPEN(24, file=fname, status='replace', action='write')
          fname = fname_base(1:INDEX(fname_base, " ")-1) // "_est.res"
          OPEN(27, file=fname, status='replace', action='write')
          
          WRITE(ch,'(I6.6)') numpe
          OPEN(112,file=fname_base(1:INDEX(fname_base, " ")-1)//".ensi.DISPL-" &
                        //ch,status='replace',action='write')
          WRITE(112,'(A)') "Alya Ensight Gold --- Vector per-node variable file"
          WRITE(112,'(A/A/A)') "part", "     1","coordinates"

          !fname = fname_base(1:INDEX(fname_base, " ")-1) // "_str.res"
          !OPEN(25, file=fname, status='replace', action='write')
          !fname = fname_base(1:INDEX(fname_base, " ")-1) // "_rea.res"
          !OPEN(26, file=fname, status='replace', action='write')
          !fname = fname_base(1:INDEX(fname_base, " ")-1) // "_fint.res"
          !OPEN(28, file=fname, status='replace', action='write')
                  
          ! Homogenized stress and strain
          !fname = fname_base(1:INDEX(fname_base, " ")-1) // "_hom_stress.res"
          !OPEN(29, file=fname, status='replace', action='write')
          !fname = fname_base(1:INDEX(fname_base, " ")-1) // "_hom_strain.res"
          !OPEN(30, file=fname, status='replace', action='write')
          
          ! Load and displacement
          !fname = fname_base(1:INDEX(fname_base, " ")-1) // "_disp_load.res"
          !OPEN(31, file=fname, status='replace', action='write')
        END IF
      END IF
      
      ! Calculate the load
      !IF (iload==1) THEN
      !  ALLOCATE(reacnodes_pp(nodes_pp*n_dof))
      !END IF
      
      !CALL SCATTER_NODES(npes,n_nodes,nels_pp,g_num_pp,n_nodpel,n_dof,nodes_pp,node_start, &
      ! node_end,storefint_pp,reacnodes_pp,0)
       
      !max_disp_inc = max_disp*lambda_total
       
      !CALL WRITE_LOAD_DISP(31,nodes_pp,npes,numpe,reacnodes_pp,max_disp_inc, &
      ! load_node,iload)
       
      IF (numpe==1) THEN
        !CALL FLUSH(29)
        !CALL FLUSH(30)
        !CALL FLUSH(31)
      END IF

      IF (iload==max_inc) THEN
        exit_iload = .TRUE.
        GOTO 400
      END IF
      
  !-----print out displacements, stress, principal stress and reactions -------
      GOTO 500 ! remove this to print at every load step; but this prints an
               ! extra output at the end
  400 CONTINUE 
      !IF (print_output) THEN
        writetimes = writetimes + 1
        
        ALLOCATE(xnewnodes_pp(nodes_pp*n_dof))
        ALLOCATE(shape_integral_pp(n_nodpel,nels_pp))
        !ALLOCATE(stress_integral_pp(n_nodpel*nst,nels_pp))
        ALLOCATE(strain_integral_pp(n_nodpel*nst,nels_pp))
        !ALLOCATE(stressnodes_pp(nodes_pp*nst))
        ALLOCATE(strainnodes_pp(nodes_pp*nst))
        !ALLOCATE(reacnodes_pp(nodes_pp*n_dof))
        
        CALL GATHER(xnew_pp(1:),xnewel_pp)
        IF (numfix_pp > 0) THEN
          DO i = 1,numfix_pp
            xnewel_pp(fixdof_pp(i),fixelem_pp(i)) = fixvaltot_pp(i)
          END DO
        END IF

        shape_integral_pp     = 0.0_iwp
        !stress_integral_pp    = 0.0_iwp
        strain_integral_pp    = 0.0_iwp

        DO iel = 1, nels_pp
          DO i = 1, n_nodpel
            num(i) = g_num_pp(i,iel) - nn_start + 1
          END DO

          !coord = TRANSPOSE(g_coord_pp(:,num))
          !DO i = 1,ndim
          !  auxm(:,i) = xnewel_pp(comp(:,i),iel)
          !END DO

          DO igauss = 1,n_intp
            CALL SHAPE_FUNCTIONS(igauss,points,value_shape)
            dw = det * weights(igauss)
  !         DO i = 1,n_nodpel
              shape_integral_pp(:,iel) = shape_integral_pp(:,iel) + &
                                         value_shape*dw
  !         END DO
          END DO

          !DO igauss = 1,n_intp
          !  CALL SHAPE_FUNCTIONS(igauss,points,value_shape)
          !  dw = det * weights(igauss)
          !  DO i = 1,n_nodpel
          !    idx1 = (i-1)*nst 
          !    DO j = 1,nst
          !      stress_integral_pp(idx1+j,iel) = stress_integral_pp(idx1+j,iel) + &
          !                value_shape(i)*sigma1C_mat_con(j,igauss,iel)*dw
          !    END DO
          !  END DO
          !END DO

          DO igauss = 1,n_intp
            CALL SHAPE_FUNCTIONS(igauss,points,value_shape)
            dw = det * weights(igauss)
            DO i = 1,n_nodpel
              idx1 = (i-1)*nst 
              DO j=1,nst
                strain_integral_pp(idx1+j,iel) = strain_integral_pp(idx1+j,iel) + &
                          value_shape(i)*lnstrainelas_mat_con(iel,igauss,j)*dw
              END DO
            END DO
          END DO
        END DO

        text = "*DISPLACEMENT"
        CALL SCATTER_NODES(npes,n_nodes,nels_pp,g_num_pp,n_nodpel,n_dof,nodes_pp, &
                node_start,node_end,xnewel_pp,xnewnodes_pp,1)
        CALL WRITE_NODAL_VARIABLE(text,24,iload,nodes_pp,npes,numpe,n_dof, &
                                  xnewnodes_pp)

        ALLOCATE(temp(n_dof,nodes_pp))
        temp = RESHAPE(xnewnodes_pp,(/n_dof,nodes_pp/))
        DO i = 1, n_dof
          CALL dismsh_ensi_p(112,iload,nodes_pp,npes,numpe,1,temp(i,:))
        END DO
        DEALLOCATE(temp)

        DEALLOCATE(xnewnodes_pp)

  !      text = "*STRESS"
        !CALL NODAL_PROJECTION(npes,n_nodes,nels_pp,g_num_pp,n_nodpel,nst,nodes_pp, &
        ! node_start,node_end,shape_integral_pp,stress_integral_pp,stressnodes_pp)
        !CALL WRITE_NODAL_VARIABLE(text,25,iload,nodes_pp,npes,numpe,nst, &
        !                          stressnodes_pp)
        !DEALLOCATE(stress_integral_pp,stressnodes_pp)

  !      text = "*NODAL REACTIONS"
        !CALL SCATTER_NODES(npes,n_nodes,nels_pp,g_num_pp,n_nodpel,n_dof,nodes_pp, &
        !        node_start,node_end,storefint_pp,reacnodes_pp,0)
        !CALL WRITE_NODAL_VARIABLE(text,26,iload,nodes_pp,npes,numpe,n_dof, &
        !                          reacnodes_pp)
        !DEALLOCATE(reacnodes_pp)

        text = "*ELASTIC STRAIN"
        CALL NODAL_PROJECTION(npes,n_nodes,nels_pp,g_num_pp,n_nodpel,nst,nodes_pp, &
         node_start,node_end,shape_integral_pp,strain_integral_pp,strainnodes_pp)
        CALL WRITE_NODAL_VARIABLE(text,27,iload,nodes_pp,npes,numpe,nst, &
                                  strainnodes_pp)
                                  
        DEALLOCATE(strain_integral_pp,strainnodes_pp)
        DEALLOCATE(shape_integral_pp)

        ! Outputting unloading
        !CALL MPI_ALLREDUCE(yield_ip_pp,yield_tot,1,MPI_INT,MPI_SUM, &
        ! MPI_COMM_WORLD,ier) 

        !CALL MPI_ALLREDUCE(unload_ip_pp,unload_tot,1,MPI_INT,MPI_SUM, &
        ! MPI_COMM_WORLD,ier) 
        
        !IF (numpe==1) THEN
        !  WRITE(*,*) 'The total number of yielded integration points is ', &
        !   yield_tot
        !  WRITE(*,*) 'The total number of unloaded integration points is ', &
        !   unload_tot
        !  WRITE(*,*) 'The percentage of yielded over total is ', &
        !   FLOAT(yield_tot)/FLOAT(n_elem*8)
        !  WRITE(*,*) 'The percentage of unloaded over total is ', &
        !   FLOAT(unload_tot)/FLOAT(n_elem*8)
        !  WRITE(*,*) 'The percentage of unloaded over yielded is ', &
        !   FLOAT(unload_tot)/FLOAT(yield_tot*8)
        !END IF
        
        print_output=.false.

      !END IF  !printing

  500 CONTINUE
      
      timest(35) = timest(35) + elap_time()-timest(2) ! 35 = write
      timest(2) = elap_time()

      IF (exit_iload) THEN
        EXIT
      END IF

    END DO !iload
    
    timest(2) = elap_time()

    IF(numpe==1) THEN
      CLOSE(24)
      !CLOSE(25)
      !CLOSE(26)
      CLOSE(27)
      !CLOSE(29)
      !CLOSE(30)
      !CLOSE(31)
      CLOSE(112)
    END IF

  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  !------------------------- End Load LOOP --------------------------------------
  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------

  !!$            ALLOCATE(etype(n_elem),nf(n_dof,n_nodes),oldlds(n_nodes*ndim)) 
  !!$
  !!$            etype  = 1   ! Only 1.0_iwp material type in this mesh
  !!$            nf     = 0
  !!$            oldlds = 0.0_iwp
  !!$
  !!$            nstep=1; npri=1; dtim=1.0; solid=.true. 
  !!$
  !!$            k=0 
  !!$            DO j=1,loaded_freedoms
  !!$              k=k+1
  !!$              found=.false.
  !!$              DO i=1,n_nodes
  !!$                IF(i==no(k)) THEN
  !!$                  l=i*3
  !!$                  oldlds(l)=val(k)
  !!$                  found=.true.
  !!$                END IF
  !!$                IF(found)CYCLE
  !!$              END DO
  !!$            END DO
  !!$
  !!$            CALL rest_to_nf(rest,nf)
  !!$
  !!$            CALL mesh_ensi(TRIM(fname_base),LEN(TRIM(fname_base),g_coord,g_num,element,etype,nf, &
  !!$                           oldlds(1:),nstep,npri,dtim,solid)
  !!$

  END SUBROUTINE solver

  SUBROUTINE umat_quadric_linear_hardening(stress,statev,ddsdde,stran,dstran, &
                                           ntens,statevar_num,iel,igauss, &
                                           noncon_flag)
    ! This subroutine returns the updated stress, strain and the tangent 
    ! operator for the Eccentric-Ellipsoid model with (linear) isotropic 
    ! hardening and associative plastic flow rule
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: ntens
    INTEGER, INTENT(IN) :: statevar_num, iel, igauss
    REAL(iwp), INTENT(OUT) :: stress(:), ddsdde(:,:)
    REAL(iwp), INTENT(INOUT) :: stran(:), statev(:), dstran(:)
    LOGICAL, INTENT(INOUT)   :: noncon_flag

    REAL(iwp) :: scalar_term, eqplas, e, nu, yield_t, yield_c, f_lower_0, &
     f_upper_0, syield, hard, sfs, zeta, stress_eq, plastic_mul, norm_flow, &
     yield, ran_scalar, norm_solution, nen, eqplas_trial, res_piv, alpha, &
     mprod, mprev, mder, alpha1, alpha2
    REAL(iwp) :: unit_tensor(6), ixi(6,6), ixi_sym(6,6), f_4(6,6), f_2(6), &
     flow_dir(7), fs(6), res_strain(6), dnds(6,6), unit_7(7,7), &
     strain_trial(6), residuals(8), results(8), jacobian(8,8), fd(7), &
     inv_jacobian(8,8), e_comp(6,6), nxn(6,6), en(6), g_mat(7,7), &
     flow_dir_stress(6), grad_flow(7,7), dir(8), results0(8)
    REAL(iwp), PARAMETER :: tol=0.000001_iwp, tol_nr=0.00000001_iwp, beta=0.0001_iwp, &
     ls_const=0.1_iwp, four=4._iwp
    INTEGER :: i, j, iter, ls_iter
    INTEGER, PARAMETER :: max_nr_iter=50, max_ls_iter=25
     
    ! Assign material properties (user defined)
    ! Yield properties obtained from Wolfram et al. 2012
    e=12700._iwp
    nu=0.3_iwp
    yield_t=52._iwp !(0.41%)
    yield_c=105._iwp !(0.83%)
    !zeta=0.2_iwp
    zeta=0.49_iwp
    !hard=0.001 ! Corresponds to 0% of the elastic slope
    !hard=0.296_iwp ! Corresponds to 5% of the elastic slope in tension
    hard=0.038_iwp ! Corresponds to 5% of the elastic slope in compression
     
    ! Assign derived material properties
    f_upper_0=(yield_t+yield_c)/(2.0_iwp*yield_t*yield_c)
    f_lower_0=(1.0_iwp/2.0_iwp)*((1.0_iwp/yield_t)-(1.0_iwp/yield_c))
       
    ! Recover the previous equivalent plastic strain
    eqplas_trial=statev(1)

    ! Initializing variables
    ddsdde=0.0_iwp
    stress=0.0_iwp

    ! Compute the linear elastic isotropic stiffness matrix
    scalar_term=e/((1.0_iwp+nu)*(1.0_iwp-2.0_iwp*nu))

    ddsdde(1,1)=1.0_iwp-nu
    ddsdde(2,2)=1.0_iwp-nu
    ddsdde(3,3)=1.0_iwp-nu
    ! previously (1.0_iwp-2.0_iwp*nu)/2.0_iwp
    ddsdde(4,4)=0.5_iwp-nu
    ddsdde(5,5)=0.5_iwp-nu
    ddsdde(6,6)=0.5_iwp-nu
    ddsdde(1,2)=nu
    ddsdde(1,3)=nu
    ddsdde(2,1)=nu
    ddsdde(2,3)=nu
    ddsdde(3,1)=nu
    ddsdde(3,2)=nu

    ddsdde=scalar_term*ddsdde
    
    ! Calculate the predictor strain and stress
    DO i=1,ntens
      DO j=1,ntens
        stress(i)=stress(i)+ddsdde(i,j)*stran(j)
      END DO
    END DO
    
    ! Define the unit tensor, IxI and IxI_sym
    ixi=0.0_iwp
    ixi_sym=0.0_iwp
    
    DO i=1,3
      unit_tensor(i)=1.0_iwp
      unit_tensor(i+3)=0.0_iwp
      DO j=1,3
        ixi(i,j)=1.0_iwp
      END DO
      
      ixi_sym(i,i)=1.0_iwp
      ixi_sym(i+3,i+3)=0.5_iwp
    END DO

    ! Define the fourth order tensor F_4 and the second order tensor F_2
    f_4=-zeta*(f_upper_0**2)*ixi+(zeta+1.0_iwp)*(f_upper_0**2)*ixi_sym
    f_2=f_lower_0*unit_tensor
        
    ! Calculate the equivalent stress
    fs=MATMUL(f_4,stress)
    fs(4:6)=four*fs(4:6)
    sfs=DOT_PROD(stress,fs,6)
    sfs=SQRT(sfs)
    stress_eq=sfs+DOT_PROD(f_2,stress,6)
        
    ! Calculate the equivalent yield stress
    syield=1.0_iwp+hard*eqplas_trial
    
    ! Determine if there is yielding
    ! =========================================================================
    IF ((stress_eq-syield)>=tol) THEN

      ! This material point is yielding, proceed with the return-mapping
      ! =======================================================================

      ! Initialise some matrices
      g_mat=0.0_iwp
      g_mat(1:6,1:6)=ddsdde
      g_mat(7,7)=hard
      
      unit_7=0.0_iwp
      DO i=1,7
        unit_7(i,i)=1.0_iwp
      END DO
      
      ! Initialize variables before the local Newton-Raphson loop
      plastic_mul=0.0_iwp
      strain_trial=stran
      
      DO iter=1,max_nr_iter
      
        ! Warn if the maximum number of iterations has been reached
        IF (iter==max_nr_iter) THEN
          WRITE(*,*) 'Maximum local Newton-Raphson iterations have been reached'
          noncon_flag=.TRUE.
          EXIT
        END IF
        
        ! Calculate the flow direction
        IF (iter.EQ.1) THEN
          flow_dir_stress=(fs/sfs)+f_2
      
          ! Compute the residuals 
          yield=sfs+DOT_PROD(f_2,stress,6)-(1.0_iwp+hard*(plastic_mul+eqplas_trial))
         
          ! Assemble the residual and the results vector
          residuals=0.0_iwp
          results(1:6)=stran(1:6)  
          
          residuals(7)=0.0_iwp
          results(7)=-eqplas_trial
          
          residuals(8)=yield
          results(8)=0.0_iwp
        END IF
        
        ! Assemble the flow direction
        flow_dir(1:6)=flow_dir_stress(1:6)
        flow_dir(7)=1.0_iwp
         
        ! Calculate the Jacobian
        ! =====================================================================
         
        ! Calculate the derivative of the flow vector with respect to stress
        fs=MATMUL(f_4,stress)
        fs(4:6)=four*fs(4:6)
        
        dnds=(f_4/sfs)
        dnds(4:6,4:6)=four*dnds(4:6,4:6)
        
        dnds=dnds-TENSOR_PRODUCT_22(fs,(fs/(sfs**3)))
         
        grad_flow=0.0_iwp
        grad_flow(1:6,1:6)=dnds
                 
        ! Assemble the jacobian
        jacobian(1:7,1:7)=unit_7+plastic_mul*MATMUL(grad_flow,g_mat)
        fd=MATMUL(flow_dir,g_mat)
        
        jacobian(8,1:7)=fd(1:7)
        jacobian(1:7,8)=flow_dir(1:7)        
        jacobian(8,8)=0.0_iwp
        
        ! Invert the Jacobian
        CALL INVERSE(jacobian,inv_jacobian,8)
        ! =====================================================================
         
        ! Compute direction of advance
        dir=-MATMUL(inv_jacobian,residuals)
         
        ! Line search scheme
        ! =====================================================================
        
        ! Set up some initial results
        alpha=1.0_iwp
        mprod=0.5_iwp*DOT_PROD(residuals,residuals,8)
        mprev=mprod
        mder=-2.0_iwp*mprod
        results0=results
        
        DO ls_iter=1,max_ls_iter
        
          ! Update new results
          results=results0+alpha*dir
          plastic_mul=results(8)
          eqplas=-results(7)
          stran(1:6)=results(1:6)
          stress=MATMUL(ddsdde,stran)
          
          fs=MATMUL(f_4,stress)
          fs(4:6)=four*fs(4:6)
          sfs=DOT_PROD(stress,fs,6)
          sfs=SQRT(sfs)
          flow_dir_stress=(fs/sfs)+f_2
          
          res_strain=stran-strain_trial+plastic_mul*flow_dir_stress
          res_piv=-eqplas+eqplas_trial+plastic_mul
          yield=sfs+DOT_PROD(f_2,stress,6)-(1.0_iwp+hard*(plastic_mul+eqplas_trial))
         
          residuals(1:6)=res_strain(1:6)
          residuals(7)=res_piv
          residuals(8)=yield
          
          mprod=0.5_iwp*DOT_PROD(residuals,residuals,8)
          
          IF (mprod<=((1.0_iwp-2.0_iwp*ls_const*beta)*mprev)) THEN
            EXIT
          ELSE
            alpha1=ls_const*alpha
            alpha2=(-(alpha**2)*mder)/(2.0_iwp*(mprod-mprev-alpha*mder))
            
            IF (alpha1>=alpha2) THEN
              alpha=alpha1
            ELSE
              alpha=alpha2
            END IF
          END IF
        END DO
        ! =====================================================================
         
        ! Exit if convergence
        ! =====================================================================
        norm_solution=0.0_iwp
        DO i=1,8
          norm_solution=norm_solution+(residuals(i)**2)
        END DO
        norm_solution=SQRT(norm_solution)
         
        IF (norm_solution<=tol_nr) THEN
          EXIT
        END IF
        ! =====================================================================
        
      END DO
      
      ! Update equivalent plastic strain     
      eqplas=eqplas_trial+plastic_mul

      ! Assemble tangent operator
      CALL INVERSE(ddsdde,e_comp,6)
      dnds=e_comp+plastic_mul*dnds
      CALL INVERSE(dnds,ddsdde,6)
      en=MATMUL(ddsdde,flow_dir_stress)
      nxn=tensor_product_22(en,en)
      
      DO i=4,6
        flow_dir(i)=0.5_iwp*flow_dir(i)
      END DO
      
      nen=double_contraction_22(flow_dir,en)+hard
      ddsdde=ddsdde-(1.0_iwp/nen)*nxn
      
      ! Update state variables
      statev(1)=eqplas

    END IF
    
    ! End of yielding
    ! =========================================================================
    
    RETURN
  END SUBROUTINE umat_quadric_linear_hardening

END MODULE core
