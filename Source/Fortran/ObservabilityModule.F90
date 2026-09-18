!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> Opt-in MatrixMultiply observability and dense-branch safety guard.
MODULE ObservabilityModule
  USE DataTypesModule, ONLY : NTREAL, NTLONG, MPINTINTEGER, MPINTLONG, MPINTREAL
  USE NTMPIModule
#ifdef _OPENMP
  USE OMP_LIB, ONLY : OMP_GET_MAX_THREADS
#endif
  IMPLICIT NONE
  PRIVATE

  LOGICAL, SAVE :: configured = .FALSE.
  LOGICAL, SAVE :: enabled = .FALSE.
  LOGICAL, SAVE :: abort_before_dense = .FALSE.
  INTEGER(NTLONG), SAVE :: dense_memory_limit_bytes = 0_NTLONG
  INTEGER, SAVE :: active_comm = MPI_COMM_WORLD
  INTEGER, SAVE :: active_rank = -1
  INTEGER, SAVE :: call_id = 0
  INTEGER, SAVE :: local_gemm_count = 0
  INTEGER, SAVE :: local_dense_count = 0
  INTEGER, SAVE :: local_sparse_count = 0
  INTEGER(NTLONG), SAVE :: input_a_nnz = 0_NTLONG
  INTEGER(NTLONG), SAVE :: input_b_nnz = 0_NTLONG
  INTEGER, SAVE :: matrix_dimension = 0
  REAL(NTREAL), SAVE :: multiply_threshold = 0.0_NTREAL
  REAL(NTREAL), SAVE :: max_operand_occupancy = 0.0_NTREAL
  INTEGER(NTLONG), SAVE :: max_dense_bytes = 0_NTLONG
  REAL(NTREAL), SAVE :: start_time = 0.0_NTREAL

  PUBLIC :: BeginObservedMatrixMultiply
  PUBLIC :: EndObservedMatrixMultiply
  PUBLIC :: ObservabilityEnabled
  PUBLIC :: RecordLocalGemm

CONTAINS

  SUBROUTINE ConfigureObservability()
    CHARACTER(LEN=128) :: value
    INTEGER :: length, status

    IF (configured) RETURN
    configured = .TRUE.

    CALL GET_ENVIRONMENT_VARIABLE("NTPOLY_OBSERVABILITY", value, length, status)
    IF (status .EQ. 0 .AND. length .GT. 0) THEN
       enabled = (value(1:1) .EQ. "1" .OR. value(1:1) .EQ. "T" .OR. &
            & value(1:1) .EQ. "t" .OR. value(1:1) .EQ. "Y" .OR. &
            & value(1:1) .EQ. "y")
    END IF

    CALL GET_ENVIRONMENT_VARIABLE("NTPOLY_ABORT_BEFORE_DENSE", value, &
         & length, status)
    IF (status .EQ. 0 .AND. length .GT. 0) THEN
       abort_before_dense = (value(1:1) .EQ. "1" .OR. &
            & value(1:1) .EQ. "T" .OR. value(1:1) .EQ. "t" .OR. &
            & value(1:1) .EQ. "Y" .OR. value(1:1) .EQ. "y")
    END IF

    CALL GET_ENVIRONMENT_VARIABLE("NTPOLY_DENSE_MEMORY_LIMIT_BYTES", value, &
         & length, status)
    IF (status .EQ. 0 .AND. length .GT. 0) THEN
       READ(value(1:length), *, IOSTAT=status) dense_memory_limit_bytes
       IF (status .NE. 0) dense_memory_limit_bytes = 0_NTLONG
    END IF
  END SUBROUTINE ConfigureObservability

  FUNCTION ObservabilityEnabled() RESULT(is_enabled)
    LOGICAL :: is_enabled
    CALL ConfigureObservability()
    is_enabled = enabled
  END FUNCTION ObservabilityEnabled

  SUBROUTINE BeginObservedMatrixMultiply(comm, rank, dimension, a_nnz, b_nnz, &
       & threshold)
    INTEGER, INTENT(IN) :: comm, rank, dimension
    INTEGER(NTLONG), INTENT(IN) :: a_nnz, b_nnz
    REAL(NTREAL), INTENT(IN) :: threshold

    CALL ConfigureObservability()
    IF (.NOT. enabled) RETURN
    active_comm = comm
    active_rank = rank
    call_id = call_id + 1
    matrix_dimension = dimension
    input_a_nnz = a_nnz
    input_b_nnz = b_nnz
    multiply_threshold = threshold
    local_gemm_count = 0
    local_dense_count = 0
    local_sparse_count = 0
    max_operand_occupancy = 0.0_NTREAL
    max_dense_bytes = 0_NTLONG
    start_time = MPI_WTIME()
  END SUBROUTINE BeginObservedMatrixMultiply

  SUBROUTINE RecordLocalGemm(rows_a, columns_a, nnz_a, rows_b, columns_b, &
       & nnz_b, rows_c, columns_c, is_dense, element_bytes)
    INTEGER, INTENT(IN) :: rows_a, columns_a, nnz_a
    INTEGER, INTENT(IN) :: rows_b, columns_b, nnz_b
    INTEGER, INTENT(IN) :: rows_c, columns_c, element_bytes
    LOGICAL, INTENT(IN) :: is_dense
    REAL(NTREAL) :: occupancy_a, occupancy_b
    INTEGER(NTLONG) :: dense_bytes, concurrent_dense_bytes
    INTEGER :: dense_concurrency
    INTEGER :: ierr

    CALL ConfigureObservability()
    IF (.NOT. enabled .AND. .NOT. abort_before_dense) RETURN
    occupancy_a = REAL(nnz_a, NTREAL) / &
         & REAL(MAX(1_NTLONG, INT(rows_a, NTLONG)*INT(columns_a, NTLONG)), &
         & NTREAL)
    occupancy_b = REAL(nnz_b, NTREAL) / &
         & REAL(MAX(1_NTLONG, INT(rows_b, NTLONG)*INT(columns_b, NTLONG)), &
         & NTREAL)
    dense_bytes = INT(element_bytes, NTLONG) * &
         & (INT(rows_a, NTLONG)*INT(columns_a, NTLONG) + &
         &  INT(rows_b, NTLONG)*INT(columns_b, NTLONG) + &
         &  INT(rows_c, NTLONG)*INT(columns_c, NTLONG))
    dense_concurrency = 1
#ifdef _OPENMP
    dense_concurrency = MAX(1, OMP_GET_MAX_THREADS())
#endif
    concurrent_dense_bytes = dense_bytes * INT(dense_concurrency, NTLONG)

    IF (enabled) THEN
!$OMP CRITICAL(NTPOLY_OBSERVABILITY_RECORD)
       local_gemm_count = local_gemm_count + 1
       IF (is_dense) THEN
          local_dense_count = local_dense_count + 1
       ELSE
          local_sparse_count = local_sparse_count + 1
       END IF
       max_operand_occupancy = MAX(max_operand_occupancy, occupancy_a, &
            & occupancy_b)
       IF (is_dense) max_dense_bytes = MAX(max_dense_bytes, &
            & concurrent_dense_bytes)
!$OMP END CRITICAL(NTPOLY_OBSERVABILITY_RECORD)
    END IF

    IF (is_dense .AND. abort_before_dense .AND. &
         & dense_memory_limit_bytes .GT. 0_NTLONG .AND. &
         & concurrent_dense_bytes .GT. dense_memory_limit_bytes) THEN
       WRITE(*,'(A,I0,A,I0,A,I0)') &
            & "NTPOLY_DENSE_GUARD rank=", active_rank, &
            & " predicted_bytes=", concurrent_dense_bytes, &
            & " limit_bytes=", dense_memory_limit_bytes
       CALL MPI_ABORT(active_comm, 86, ierr)
    END IF
  END SUBROUTINE RecordLocalGemm

  SUBROUTINE EndObservedMatrixMultiply(output_nnz)
    INTEGER(NTLONG), INTENT(IN) :: output_nnz
    INTEGER :: global_gemm_count, global_dense_count, global_sparse_count
    INTEGER(NTLONG) :: global_max_dense_bytes
    REAL(NTREAL) :: global_max_occupancy, elapsed, global_elapsed
    REAL(NTREAL) :: denom
    INTEGER :: ierr

    IF (.NOT. enabled) RETURN
    elapsed = MPI_WTIME() - start_time
    CALL MPI_Allreduce(local_gemm_count, global_gemm_count, 1, MPINTINTEGER, &
         & MPI_SUM, active_comm, ierr)
    CALL MPI_Allreduce(local_dense_count, global_dense_count, 1, MPINTINTEGER, &
         & MPI_SUM, active_comm, ierr)
    CALL MPI_Allreduce(local_sparse_count, global_sparse_count, 1, &
         & MPINTINTEGER, MPI_SUM, active_comm, ierr)
    CALL MPI_Allreduce(max_operand_occupancy, global_max_occupancy, 1, &
         & MPINTREAL, MPI_MAX, active_comm, ierr)
    CALL MPI_Allreduce(max_dense_bytes, global_max_dense_bytes, 1, MPINTLONG, &
         & MPI_MAX, active_comm, ierr)
    CALL MPI_Allreduce(elapsed, global_elapsed, 1, MPINTREAL, MPI_MAX, &
         & active_comm, ierr)

    IF (active_rank .EQ. 0) THEN
       denom = REAL(MAX(1, matrix_dimension), NTREAL)**2
       WRITE(*,'(A,I0,A,I0,A,ES14.6,A,I0,A,ES14.6,A,I0,A,ES14.6,' // &
            & 'A,I0,A,ES14.6,A,I0,A,I0,A,I0,A,I0,A,ES14.6,A,ES14.6)') &
            & "NTPOLY_OBS call=", call_id, " dim=", matrix_dimension, &
            & " threshold=", multiply_threshold, " a_nnz=", input_a_nnz, &
            & " a_occ=", REAL(input_a_nnz, NTREAL)/denom, &
            & " b_nnz=", input_b_nnz, &
            & " b_occ=", REAL(input_b_nnz, NTREAL)/denom, &
            & " c_nnz=", output_nnz, &
            & " c_occ=", REAL(output_nnz, NTREAL)/denom, &
            & " local_gemm_sum=", global_gemm_count, &
            & " dense_sum=", global_dense_count, &
            & " sparse_sum=", global_sparse_count, &
            & " max_dense_bytes=", global_max_dense_bytes, &
            & " max_operand_occ=", global_max_occupancy, &
            & " walltime_max_s=", global_elapsed
    END IF
  END SUBROUTINE EndObservedMatrixMultiply

END MODULE ObservabilityModule
