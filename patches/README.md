# GROMACS Patches for GRO-DFTB

GRO-DFTB requires modifications to the GROMACS source tree to add the
`DFTBForceProvider` backend. These changes are distributed as a patch
file that applies to GROMACS 2027.0-dev (tag `v2026.0`, commit
`4ac0b1a023`).

## Files Modified/Added

| File | Change |
|---|---|
| `src/gromacs/applied_forces/qmmm/CMakeLists.txt` | Add `dftbcheckpointstate.cpp` to build |
| `src/gromacs/applied_forces/qmmm/dftbforceprovider.cpp` | `DFTBForceProvider::calculateForces()` implementation |
| `src/gromacs/applied_forces/qmmm/dftbforceprovider.h` | Class declaration with PME, embedding, checkpoint members |
| `src/gromacs/applied_forces/qmmm/dftbforceprovider_stub.cpp` | Stub for builds without `GMX_DFTB` |
| `src/gromacs/applied_forces/qmmm/qmmm.cpp` | DFTB+ registration in `QMMMModule` |
| `src/gromacs/applied_forces/qmmm/tests/CMakeLists.txt` | Add checkpoint state test |
| `src/gromacs/applied_forces/qmmm/dftbcheckpointstate.cpp` | **New**: KVT checkpoint serialization |
| `src/gromacs/applied_forces/qmmm/dftbcheckpointstate.h` | **New**: Checkpoint state class |
| `src/gromacs/applied_forces/qmmm/tests/dftbcheckpointstate.cpp` | **New**: Checkpoint GTest |

## Applying the Patch

```bash
# From the GROMACS source directory:
cd external/gromacs
git checkout v2026.0    # or commit 4ac0b1a023
git apply ../../patches/gromacs-dftb.patch

# Then build GROMACS with DFTB+ support:
cmake -B _build -S . \
  -DCMAKE_INSTALL_PREFIX=_install \
  -DGMX_DFTB=ON \
  -Dlibgrodftb_DIR=/path/to/grodftb/build
cmake --build _build --target install
```

## Verifying the Patch

After building, the GROMACS QM/MM tests should include the DFTB+ backend tests:

```bash
cd _build
ctest -R qmmm
```

## Compatibility

This patch is tested against GROMACS tag `v2026.0` (commit `4ac0b1a023`).
Applying to other GROMACS versions may require manual conflict resolution.
