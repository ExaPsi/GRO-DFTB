/*
 * SDD:specs.md:§5.2 — Smoke test for DFTB+ C API availability.
 * Verifies that the installed libdftbplus exposes dftbp_api() and
 * reports version 0.4.0.
 */
#include <stdio.h>
#include <stdlib.h>
#include "dftbplus.h"

int main(void) {
    int major, minor, patch;
    dftbp_api(&major, &minor, &patch);
    printf("DFTB+ API version: %d.%d.%d\n", major, minor, patch);
    if (major != 0 || minor != 4 || patch != 0) {
        fprintf(stderr, "ERROR: Expected API version 0.4.0, got %d.%d.%d\n",
                major, minor, patch);
        return 1;
    }
    printf("API version check: PASS\n");
    return 0;
}
