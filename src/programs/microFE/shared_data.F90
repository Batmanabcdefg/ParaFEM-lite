MODULE shared_data

  USE precision
  USE choose_solvers

!   Formats
  CHARACTER(len=*), PARAMETER :: energy_fmt = '(a,i3,1p,i3,1p,e25.15,1p,e25.5)'

  INTEGER :: log_unit

  INTEGER :: n_elem, n_nodes, n_restrained, n_intp, n_dof=3, n_nodpel, nst=6, loaded_nodes, nn_pp, &
   nf_start, fmt=1, i, j, k, l, ndim=3, iters, limit, iel, nn_start, &
   num_load_steps, iload, igauss, dimH, inewton, jump, npes_pp, partitioner, &
   limit_2, fixed_nodes, numfix_pp, fixdim, writetimes=0, &
   nodes_pp, node_start, node_end, idx1, idx2, rm_sum

  REAL(iwp) :: det, tol, maxdiff, tol2, detF, detFinc, energy, energy1, rn0, &
   timest(40), detF_mean, initial_guess, lambda, lambda_total, lambda_prev, &
   energy_prev, energy_prev_prev, min_inc, next_output, temp_inc, dw, e, v

  REAL(iwp), PARAMETER :: tol_increment=0.000001_iwp, tol_val=0.00000001_iwp
  REAL(iwp), PARAMETER :: zero=0.0_iwp, one=1._iwp, half=0.5_iwp
  INTEGER, PARAMETER :: statevar_num=10
  
  REAL(iwp) :: E_mat(3,3), S_mat(3,3), C_mat(3,3), jacFtransinv(3,3), &
   jacFtrans(3,3), sum_strain(3,3), sum_stress(3,3), sum_strain_pp(3,3), &
   sum_stress_pp(3,3), sum_stress_prev(3,3), sum_strain_prev(3,3)
  
  INTEGER :: iter, prev_iter, max_inc, number_output, counter_output

  CHARACTER(len=15) :: element
  CHARACTER(len=50) :: text, fname_base, fname
  CHARACTER(LEN=6)  :: ch 
 
  LOGICAL :: converged, timewrite=.TRUE., flag=.FALSE., print_output=.FALSE., &
   tol_inc=.FALSE., lambda_inc=.TRUE., noncon_flag=.FALSE.,exit_iload=.false.

  CHARACTER(len=choose_solvers_string_length) :: solvers
  LOGICAL              :: error
  CHARACTER, PARAMETER :: tab = ACHAR(9)
  REAL                 :: peak_memory_use
  
  !-------------------------- dynamic arrays-----------------------------------
  REAL(iwp), ALLOCATABLE:: points(:,:), coord(:,:), weights(:), xnew_pp(:), &
   diag_precon_pp(:), r_pp(:), bee(:,:), load_value(:,:), g_coord_pp(:,:), &
   diag_precon_tmp(:,:), storekm_pp(:,:,:), disp(:,:), g_coord(:,:), &
   val_pp(:), disp_pp(:,:), res_pp(:), upd_coord(:,:), fext_pp(:), auxm(:,:), &
   fextpiece_pp(:), deltax_pp(:), fint_pp(:), kmat_elem(:,:), jacF(:,:), &
   xnewel_pp(:,:),   derivFtran(:,:), derivF(:,:), storefint_pp(:,:), &
   deeF(:,:), jacFinv(:,:), beeF(:,:), sigma1C(:), fixed_value(:), &
   fixval_pp(:), fixvalpiece_pp(:), jacFinc(:,:), statev(:), sigma(:,:), &
   lnstrainelas_mat(:,:,:), lnstrainelas(:), statevar(:,:,:), xnewnodes_pp(:),&
   sigma1C_mat(:,:,:), xnewelinc_pp(:,:), auxm_inc(:,:), sigma1C_temp(:,:,:), &
   deltax_pp_temp(:), statevar_con(:,:,:), sigma1C_mat_con(:,:,:), &
   lnstrainelas_mat_con(:,:,:), detFincmod_vec(:), geeF(:,:), geeFT(:,:), &
   xnew_pp_previous(:), xnewel_pp_previous(:,:), auxm_previous(:,:), &
   fextprev_pp(:), fixvalprev_pp(:), fixvaltot_pp(:), xnewprev_pp(:), &
   value_shape(:), shape_integral_pp(:,:), stress_integral_pp(:,:), &
   stressnodes_pp(:), strain_integral_pp(:,:), strainnodes_pp(:), &
   principal_integral_pp(:,:), princinodes_pp(:), principal(:), &
   reacnodes_pp(:), stiffness_mat_con(:,:,:,:), km(:,:), fixkm_pp(:,:,:), &
   temp(:,:)

  INTEGER, ALLOCATABLE  :: num(:), g_num(:,:), g_num_pp(:,:), g_g_pp(:,:), &
   load_node(:), rest(:,:), nf_pp(:,:), no_pp(:), comp(:,:), fixed_node(:), &
   fixed_dof(:), fixelem_pp(:), fixdof_pp(:), unload_pp(:,:), rm_vec(:)

END MODULE shared_data
