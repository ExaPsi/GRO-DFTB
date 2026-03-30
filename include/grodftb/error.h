/*
 * SDD:specs.md:§18.1 — Error codes for the GRO-DFTB driver library.
 *
 * All grodftb_* functions returning int use these codes.
 * 0 = success; nonzero = error.
 */

#ifndef GRODFTB_ERROR_H
#define GRODFTB_ERROR_H

#ifdef __cplusplus
extern "C" {
#endif

/* SDD:specs.md:§18.1 — Error codes (numeric values match spec exactly) */
#define GRODFTB_SUCCESS                     0
#define GRODFTB_ERR_NULL_HANDLE             1
#define GRODFTB_ERR_NOT_INITIALIZED         2
#define GRODFTB_ERR_ALREADY_INIT            3
#define GRODFTB_ERR_INVALID_ARGUMENT        4
#define GRODFTB_ERR_SIZE_MISMATCH           5
#define GRODFTB_ERR_FILE_NOT_FOUND          6
#define GRODFTB_ERR_HSD_PARSE               7
#define GRODFTB_ERR_DFTB_INIT               8
#define GRODFTB_ERR_SCC_NOT_CONVERGED       9
#define GRODFTB_ERR_NO_RESULTS              10
#define GRODFTB_ERR_EXCITED_NOT_CONFIGURED  11
#define GRODFTB_ERR_EXCITED_FAILED          12
#define GRODFTB_ERR_NAC_UNAVAILABLE         13
#define GRODFTB_ERR_CP_SOLVE_FAILED         14
#define GRODFTB_ERR_ALLOC_FAILED            15
#define GRODFTB_ERR_DFTB_INTERNAL           16

/* Forward declaration — full type in driver.h */
typedef struct grodftb_context *grodftb_handle_t;

/**
 * Get human-readable error message for the last error.
 * Thread-local; valid until the next grodftb_* call on the same handle.
 *
 * @param handle  The handle that produced the error (may be NULL)
 * @return  Static or handle-local error string
 */
const char *grodftb_last_error(grodftb_handle_t handle);

/**
 * Get string name for an error code.
 *
 * @param errcode  One of the GRODFTB_ERR_* codes
 * @return  Static string (never NULL)
 */
const char *grodftb_error_string(int errcode);

#ifdef __cplusplus
}
#endif

#endif /* GRODFTB_ERROR_H */
