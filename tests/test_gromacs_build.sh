#!/usr/bin/env bash
# US-025 + US-026 + US-027 + US-028 + US-029 + US-030 + US-031 + US-032: Build verification tests for GMX_DFTB CMake flag, DFTB+ enum, MDP params, DFTBForceProvider skeleton, calculateForces, species mapping, preprocessing reuse, and energy bookkeeping
# SDD:specs.md:S3.4 -- Verify GROMACS builds with GMX_DFTB=OFF (stub) and GMX_DFTB=ON (placeholder)
# US-026: Verify DFTBPLUS enum, string mapping, dispatch stub, and CP2K guard in GROMACS source
# US-027: Verify DFTB+-specific MDP parameters (enums, fields, registration, validation)
#
# Usage: bash tests/test_gromacs_build.sh
# Requires: cmake, make/ninja, C++17 compiler

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
GROMACS_SRC="$REPO_ROOT/external/gromacs"

PASS=0
FAIL=0
TOTAL=0

report() {
    local name="$1" result="$2"
    TOTAL=$((TOTAL + 1))
    if [ "$result" -eq 0 ]; then
        PASS=$((PASS + 1))
        echo "PASS: $name"
    else
        FAIL=$((FAIL + 1))
        echo "FAIL: $name"
    fi
}

# Detect number of cores for parallel build
NPROC=$(nproc 2>/dev/null || echo 4)

# --------------------------------------------------------------------------
# Test 1: GMX_DFTB=OFF (stub compiled)
# --------------------------------------------------------------------------
test_gromacs_build_dftb_off() {
    local build_dir
    build_dir=$(mktemp -d)
    trap "rm -rf $build_dir" RETURN

    cmake -S "$GROMACS_SRC" -B "$build_dir" \
        -DGMX_DFTB=OFF \
        -DGMX_CP2K=OFF \
        -DGMX_GPU=off \
        -DGMX_MPI=OFF \
        -DGMX_OPENMP=OFF \
        -DGMX_SIMD=None \
        -DBUILD_TESTING=OFF \
        -DGMXAPI=OFF \
        -DGMX_INSTALL_NBLIB_API=OFF \
        -DCMAKE_BUILD_TYPE=Release \
        > "$build_dir/cmake_out.log" 2>&1

    cmake --build "$build_dir" --target libgromacs -j "$NPROC" \
        > "$build_dir/build_out.log" 2>&1
    return $?
}

# --------------------------------------------------------------------------
# Test 2: GMX_DFTB=ON (placeholder compiled)
# --------------------------------------------------------------------------
test_gromacs_build_dftb_on() {
    local build_dir
    build_dir=$(mktemp -d)
    trap "rm -rf $build_dir" RETURN

    cmake -S "$GROMACS_SRC" -B "$build_dir" \
        -DGMX_DFTB=ON \
        -DGRODFTB_INCLUDE_DIR="$REPO_ROOT/include" \
        -DGMX_CP2K=OFF \
        -DGMX_GPU=off \
        -DGMX_MPI=OFF \
        -DGMX_OPENMP=OFF \
        -DGMX_SIMD=None \
        -DBUILD_TESTING=OFF \
        -DGMXAPI=OFF \
        -DGMX_INSTALL_NBLIB_API=OFF \
        -DCMAKE_BUILD_TYPE=Release \
        > "$build_dir/cmake_out.log" 2>&1

    cmake --build "$build_dir" --target libgromacs -j "$NPROC" \
        > "$build_dir/build_out.log" 2>&1
    return $?
}

# --------------------------------------------------------------------------
# Test 3: Verify stub contains the expected error string
# --------------------------------------------------------------------------
test_stub_throws_error() {
    local stub_file="$GROMACS_SRC/src/gromacs/applied_forces/qmmm/dftbforceprovider_stub.cpp"
    grep -q "DFTB+ support is not available. Reconfigure GROMACS with -DGMX_DFTB=ON." "$stub_file"
    return $?
}

# --------------------------------------------------------------------------
# Test 4 (US-026): DFTBPLUS enum exists in qmmmtypes.h
# --------------------------------------------------------------------------
test_dftb_enum_exists() {
    local types_file="$GROMACS_SRC/src/gromacs/applied_forces/qmmm/qmmmtypes.h"
    grep -q 'DFTBPLUS' "$types_file"
    return $?
}

# --------------------------------------------------------------------------
# Test 5 (US-026): DFTBPLUS string mapping exists in qmmmtypes.h
# --------------------------------------------------------------------------
test_dftb_enum_string_mapping() {
    local types_file="$GROMACS_SRC/src/gromacs/applied_forces/qmmm/qmmmtypes.h"
    grep -q '"DFTBPLUS"' "$types_file"
    return $?
}

# --------------------------------------------------------------------------
# Test 6 (US-026): DFTBForceProvider dispatch stub in qmmm.cpp
# --------------------------------------------------------------------------
test_dftb_dispatch_stub() {
    local qmmm_file="$GROMACS_SRC/src/gromacs/applied_forces/qmmm/qmmm.cpp"
    grep -q 'DFTBForceProvider' "$qmmm_file"
    return $?
}

# --------------------------------------------------------------------------
# Test 7 (US-026): DFTBPLUS guard in qmmmoptions.cpp (skip CP2K input)
# --------------------------------------------------------------------------
test_dftb_skips_cp2k_input() {
    local options_file="$GROMACS_SRC/src/gromacs/applied_forces/qmmm/qmmmoptions.cpp"
    grep -q 'DFTBPLUS' "$options_file"
    return $?
}

# --------------------------------------------------------------------------
# Test 8 (US-027): QMMMDFTBEmbedding enum exists with Cutoff and Pme
# --------------------------------------------------------------------------
test_dftbplus_embedding_enum_exists() {
    local types_file="$GROMACS_SRC/src/gromacs/applied_forces/qmmm/qmmmtypes.h"
    grep -q 'enum class QMMMDFTBEmbedding' "$types_file" \
        && grep -q 'Cutoff' "$types_file" \
        && grep -q 'Pme' "$types_file"
    return $?
}

# --------------------------------------------------------------------------
# Test 9 (US-027): Embedding enum string mapping contains CUTOFF and PME
# --------------------------------------------------------------------------
test_dftbplus_embedding_string_mapping() {
    local types_file="$GROMACS_SRC/src/gromacs/applied_forces/qmmm/qmmmtypes.h"
    grep -q '"CUTOFF"' "$types_file" && grep -q '"PME"' "$types_file"
    return $?
}

# --------------------------------------------------------------------------
# Test 10 (US-027): QMMMDFTBChargeRedistribution enum exists with None, Shift, Zero
# --------------------------------------------------------------------------
test_dftbplus_charge_redist_enum_exists() {
    local types_file="$GROMACS_SRC/src/gromacs/applied_forces/qmmm/qmmmtypes.h"
    grep -q 'enum class QMMMDFTBChargeRedistribution' "$types_file" \
        && grep -q 'None' "$types_file" \
        && grep -q 'Shift' "$types_file" \
        && grep -q 'Zero' "$types_file"
    return $?
}

# --------------------------------------------------------------------------
# Test 11 (US-027): Charge redistribution string mapping contains NONE, SHIFT, ZERO
# --------------------------------------------------------------------------
test_dftbplus_charge_redist_string_mapping() {
    local types_file="$GROMACS_SRC/src/gromacs/applied_forces/qmmm/qmmmtypes.h"
    grep -q '"NONE"' "$types_file" \
        && grep -q '"SHIFT"' "$types_file" \
        && grep -q '"ZERO"' "$types_file"
    return $?
}

# --------------------------------------------------------------------------
# Test 12 (US-027): dftbplus-hsd tag registered in qmmmoptions.cpp
# --------------------------------------------------------------------------
test_dftbplus_hsd_param_registered() {
    local options_file="$GROMACS_SRC/src/gromacs/applied_forces/qmmm/qmmmoptions.cpp"
    grep -q 'dftbplus-hsd' "$options_file"
    return $?
}

# --------------------------------------------------------------------------
# Test 13 (US-027): dftbplus-cutoff tag registered in qmmmoptions.cpp
# --------------------------------------------------------------------------
test_dftbplus_cutoff_param_registered() {
    local options_file="$GROMACS_SRC/src/gromacs/applied_forces/qmmm/qmmmoptions.cpp"
    grep -q 'dftbplus-cutoff' "$options_file"
    return $?
}

# --------------------------------------------------------------------------
# Test 14 (US-027): Four DFTB+ fields exist in QMMMParameters struct
# --------------------------------------------------------------------------
test_dftbplus_params_in_struct() {
    local types_file="$GROMACS_SRC/src/gromacs/applied_forces/qmmm/qmmmtypes.h"
    grep -q 'dftbplusHsd_' "$types_file" \
        && grep -q 'dftbplusEmbedding_' "$types_file" \
        && grep -q 'dftbplusCutoff_' "$types_file" \
        && grep -q 'dftbplusChargeRedistribution_' "$types_file"
    return $?
}

# --------------------------------------------------------------------------
# Test 15 (US-027): Embedding default is Cutoff
# --------------------------------------------------------------------------
test_dftbplus_embedding_default_cutoff() {
    local types_file="$GROMACS_SRC/src/gromacs/applied_forces/qmmm/qmmmtypes.h"
    grep -q 'dftbplusEmbedding_.*=.*QMMMDFTBEmbedding::Cutoff' "$types_file"
    return $?
}

# --------------------------------------------------------------------------
# Test 16 (US-027): Cutoff default is 1.2
# --------------------------------------------------------------------------
test_dftbplus_cutoff_default_1_2() {
    local types_file="$GROMACS_SRC/src/gromacs/applied_forces/qmmm/qmmmtypes.h"
    grep -q 'dftbplusCutoff_.*=.*1\.2' "$types_file"
    return $?
}

# --------------------------------------------------------------------------
# Test 17 (US-027): Charge redistribution default is Shift
# --------------------------------------------------------------------------
test_dftbplus_charge_redist_default_shift() {
    local types_file="$GROMACS_SRC/src/gromacs/applied_forces/qmmm/qmmmtypes.h"
    grep -q 'dftbplusChargeRedistribution_.*=.*QMMMDFTBChargeRedistribution::Shift' "$types_file"
    return $?
}

# --------------------------------------------------------------------------
# Test 18 (US-027): DFTB+ params conditional MDP output guarded by DFTBPLUS
# --------------------------------------------------------------------------
test_dftbplus_params_conditional_mdp_output() {
    local options_file="$GROMACS_SRC/src/gromacs/applied_forces/qmmm/qmmmoptions.cpp"
    # buildMdpOutput must contain a DFTBPLUS guard near dftbplus-hsd output
    grep -A15 'QMMMQMMethod::DFTBPLUS' "$options_file" | grep -q 'c_dftbplusHsdTag_'
    return $?
}

# --------------------------------------------------------------------------
# Test 19 (US-027): HSD validation guard (GMX_THROW for empty HSD)
# --------------------------------------------------------------------------
test_dftbplus_hsd_validation_guard() {
    local options_file="$GROMACS_SRC/src/gromacs/applied_forces/qmmm/qmmmoptions.cpp"
    grep -q 'dftbplusHsd_.empty()' "$options_file"
    return $?
}

# --------------------------------------------------------------------------
# Test 20 (US-027): Cutoff validation guard (positivity check)
# --------------------------------------------------------------------------
test_dftbplus_cutoff_validation_guard() {
    local options_file="$GROMACS_SRC/src/gromacs/applied_forces/qmmm/qmmmoptions.cpp"
    grep -q 'dftbplusCutoff_.*<=.*0' "$options_file"
    return $?
}

# --------------------------------------------------------------------------
# Test 21 (US-028): DFTBForceProvider inherits IForceProvider
# --------------------------------------------------------------------------
test_028_dftb_inherits_iforceprovider() {
    local header="$GROMACS_SRC/src/gromacs/applied_forces/qmmm/dftbforceprovider.h"
    grep -q 'class DFTBForceProvider.*public IForceProvider' "$header"
    return $?
}

# --------------------------------------------------------------------------
# Test 22 (US-028): calculateForces() override present in header
# --------------------------------------------------------------------------
test_028_calculate_forces_override() {
    local header="$GROMACS_SRC/src/gromacs/applied_forces/qmmm/dftbforceprovider.h"
    grep -q 'calculateForces.*ForceProviderInput.*ForceProviderOutput.*override' "$header"
    return $?
}

# --------------------------------------------------------------------------
# Test 23 (US-028): writeCheckpointData() present in header
# --------------------------------------------------------------------------
test_028_write_checkpoint_data() {
    local header="$GROMACS_SRC/src/gromacs/applied_forces/qmmm/dftbforceprovider.h"
    grep -q 'writeCheckpointData.*MDModulesWriteCheckpointData' "$header"
    return $?
}

# --------------------------------------------------------------------------
# Test 24 (US-028): Constructor takes 5 correct parameters
# --------------------------------------------------------------------------
test_028_constructor_signature() {
    local header="$GROMACS_SRC/src/gromacs/applied_forces/qmmm/dftbforceprovider.h"
    grep -q 'DFTBForceProvider(const QMMMParameters&' "$header" \
        && grep -q 'LocalAtomSet&.*localQMAtomSet' "$header" \
        && grep -q 'LocalAtomSet&.*localMMAtomSet' "$header" \
        && grep -q 'PbcType' "$header" \
        && grep -q 'MDLogger&' "$header"
    return $?
}

# --------------------------------------------------------------------------
# Test 25 (US-028): Stub calculateForces throws InternalError
# --------------------------------------------------------------------------
test_028_stub_calculate_forces_throws() {
    local stub="$GROMACS_SRC/src/gromacs/applied_forces/qmmm/dftbforceprovider_stub.cpp"
    grep -A5 'DFTBForceProvider::calculateForces' "$stub" | grep -q 'GMX_THROW.*InternalError'
    return $?
}

# --------------------------------------------------------------------------
# Test 26 (US-028): Stub writeCheckpointData throws InternalError
# --------------------------------------------------------------------------
test_028_stub_checkpoint_throws() {
    local stub="$GROMACS_SRC/src/gromacs/applied_forces/qmmm/dftbforceprovider_stub.cpp"
    grep -A5 'DFTBForceProvider::writeCheckpointData' "$stub" | grep -q 'GMX_THROW.*InternalError'
    return $?
}

# --------------------------------------------------------------------------
# Test 27 (US-028): dftbDescription() free function declared in header
# --------------------------------------------------------------------------
test_028_dftb_description_exists() {
    local header="$GROMACS_SRC/src/gromacs/applied_forces/qmmm/dftbforceprovider.h"
    grep -q 'std::string dftbDescription()' "$header"
    return $?
}

# --------------------------------------------------------------------------
# Test 28 (US-028): qmmm.cpp dispatches to DFTBForceProvider for DFTBPLUS
# --------------------------------------------------------------------------
test_028_dispatch_creates_dftb_provider() {
    local qmmm="$GROMACS_SRC/src/gromacs/applied_forces/qmmm/qmmm.cpp"
    grep -A10 'QMMMQMMethod::DFTBPLUS' "$qmmm" | grep -q 'make_unique<DFTBForceProvider>'
    return $?
}

# --------------------------------------------------------------------------
# Test 29 (US-028): addForceProvider called on DFTB path
# --------------------------------------------------------------------------
test_028_add_force_provider_dftb() {
    local qmmm="$GROMACS_SRC/src/gromacs/applied_forces/qmmm/qmmm.cpp"
    grep -q 'addForceProvider(dftbForceProvider_' "$qmmm"
    return $?
}

# --------------------------------------------------------------------------
# Test 30 (US-028→US-029): calculateForces calls grodftb_compute (real implementation)
# --------------------------------------------------------------------------
test_028_calculate_forces_calls_grodftb() {
    local impl="$GROMACS_SRC/src/gromacs/applied_forces/qmmm/dftbforceprovider.cpp"
    grep -q 'grodftb_compute' "$impl"
}

# --------------------------------------------------------------------------
# Test 31 (US-029): ThDD annotations present in calculateForces implementation
# --------------------------------------------------------------------------
test_029_thdd_annotations() {
    local impl="$GROMACS_SRC/src/gromacs/applied_forces/qmmm/dftbforceprovider.cpp"
    local count=0
    for tag in 'ThDD:T-US-029-2.3' 'ThDD:T-US-029-3.1' 'ThDD:T-US-029-4.1' \
               'ThDD:T-US-029-5.1' 'ThDD:T-US-029-6.1' 'ThDD:06_theory:Eq3.2'; do
        if grep -q "$tag" "$impl"; then
            count=$((count + 1))
        fi
    done
    [ "$count" -ge 6 ]
}

# --------------------------------------------------------------------------
# Test 32 (US-029): Sign convention — grodftb_get_forces returns already-negated forces
# --------------------------------------------------------------------------
test_029_sign_convention() {
    local impl="$GROMACS_SRC/src/gromacs/applied_forces/qmmm/dftbforceprovider.cpp"
    local body
    body="$(grep -B5 -A5 'grodftb_get_forces' "$impl")"
    echo "$body" | grep -qi 'already negated'
}

# --------------------------------------------------------------------------
# Test 33 (US-029): Force accumulation uses += (not bare =)
# --------------------------------------------------------------------------
test_029_force_accumulation() {
    local impl="$GROMACS_SRC/src/gromacs/applied_forces/qmmm/dftbforceprovider.cpp"
    grep -q 'force_\[.*\]\[.*\] +=' "$impl"
}

# --------------------------------------------------------------------------
# Test 34 (US-029): Energy added only on main rank
# --------------------------------------------------------------------------
test_029_main_rank_guard() {
    local impl="$GROMACS_SRC/src/gromacs/applied_forces/qmmm/dftbforceprovider.cpp"
    local body
    body="$(grep -A5 'isMainRank' "$impl")"
    echo "$body" | grep -q 'QuantumMechanicalRegionEnergy'
}

# --------------------------------------------------------------------------
# Test 35 (US-029): Uses GROMACS unit constants (c_bohr2Nm, c_hartree2Kj, c_hartreeBohr2Md)
# --------------------------------------------------------------------------
test_029_gromacs_constants() {
    local impl="$GROMACS_SRC/src/gromacs/applied_forces/qmmm/dftbforceprovider.cpp"
    grep -q 'c_bohr2Nm' "$impl" \
        && grep -q 'c_hartree2Kj' "$impl" \
        && grep -q 'c_hartreeBohr2Md' "$impl"
}

# --------------------------------------------------------------------------
# Test 36 (US-029): No GRO-DFTB internal constants in GROMACS code
# --------------------------------------------------------------------------
test_029_no_grodftb_constants() {
    local impl="$GROMACS_SRC/src/gromacs/applied_forces/qmmm/dftbforceprovider.cpp"
    if grep -q 'GRODFTB_BOHR_TO_NM' "$impl"; then
        return 1
    fi
    return 0
}

# --------------------------------------------------------------------------
# Run tests
# --------------------------------------------------------------------------
echo "=== US-025 + US-026 + US-027 + US-028 + US-029 + US-030 + US-031: GMX_DFTB Build Verification Tests ==="
echo ""

echo "Running: test_stub_throws_error"
test_stub_throws_error
report "test_stub_throws_error" $?

echo "Running: test_gromacs_build_dftb_off (this may take several minutes)"
if test_gromacs_build_dftb_off; then
    report "test_gromacs_build_dftb_off" 0
else
    report "test_gromacs_build_dftb_off" 1
fi

echo "Running: test_gromacs_build_dftb_on (this may take several minutes)"
if test_gromacs_build_dftb_on; then
    report "test_gromacs_build_dftb_on" 0
else
    report "test_gromacs_build_dftb_on" 1
fi

echo "Running: test_dftb_enum_exists"
if test_dftb_enum_exists; then
    report "test_dftb_enum_exists" 0
else
    report "test_dftb_enum_exists" 1
fi

echo "Running: test_dftb_enum_string_mapping"
if test_dftb_enum_string_mapping; then
    report "test_dftb_enum_string_mapping" 0
else
    report "test_dftb_enum_string_mapping" 1
fi

echo "Running: test_dftb_dispatch_stub"
if test_dftb_dispatch_stub; then
    report "test_dftb_dispatch_stub" 0
else
    report "test_dftb_dispatch_stub" 1
fi

echo "Running: test_dftb_skips_cp2k_input"
if test_dftb_skips_cp2k_input; then
    report "test_dftb_skips_cp2k_input" 0
else
    report "test_dftb_skips_cp2k_input" 1
fi

# US-027 tests
for t in \
    test_dftbplus_embedding_enum_exists \
    test_dftbplus_embedding_string_mapping \
    test_dftbplus_charge_redist_enum_exists \
    test_dftbplus_charge_redist_string_mapping \
    test_dftbplus_hsd_param_registered \
    test_dftbplus_cutoff_param_registered \
    test_dftbplus_params_in_struct \
    test_dftbplus_embedding_default_cutoff \
    test_dftbplus_cutoff_default_1_2 \
    test_dftbplus_charge_redist_default_shift \
    test_dftbplus_params_conditional_mdp_output \
    test_dftbplus_hsd_validation_guard \
    test_dftbplus_cutoff_validation_guard \
; do
    echo "Running: $t"
    if "$t"; then
        report "$t" 0
    else
        report "$t" 1
    fi
done

# US-028 tests
for t in \
    test_028_dftb_inherits_iforceprovider \
    test_028_calculate_forces_override \
    test_028_write_checkpoint_data \
    test_028_constructor_signature \
    test_028_stub_calculate_forces_throws \
    test_028_stub_checkpoint_throws \
    test_028_dftb_description_exists \
    test_028_dispatch_creates_dftb_provider \
    test_028_add_force_provider_dftb \
    test_028_calculate_forces_calls_grodftb \
; do
    echo "Running: $t"
    if "$t"; then
        report "$t" 0
    else
        report "$t" 1
    fi
done

# US-029 tests
for t in \
    test_029_thdd_annotations \
    test_029_sign_convention \
    test_029_force_accumulation \
    test_029_main_rank_guard \
    test_029_gromacs_constants \
    test_029_no_grodftb_constants \
; do
    echo "Running: $t"
    if "$t"; then
        report "$t" 0
    else
        report "$t" 1
    fi
done

# US-030 tests
# --------------------------------------------------------------------------
# Test 37 (US-030): Species built from atomNumbers_, not hardcoded
# --------------------------------------------------------------------------
test_030_species_from_atomic_numbers() {
    local impl="$GROMACS_SRC/src/gromacs/applied_forces/qmmm/dftbforceprovider.cpp"
    grep -q 'atomNumbers_' "$impl"
}

# --------------------------------------------------------------------------
# Test 38 (US-030): ThDD annotations for species mapping equations
# --------------------------------------------------------------------------
test_030_thdd_annotations() {
    local impl="$GROMACS_SRC/src/gromacs/applied_forces/qmmm/dftbforceprovider.cpp"
    local count=0
    for tag in 'ThDD:T-US-030-1.2' 'ThDD:T-US-030-1.3' 'ThDD:T-US-030-1.4' 'ThDD:T-US-030-1.5'; do
        if grep -q "$tag" "$impl"; then
            count=$((count + 1))
        fi
    done
    [ "$count" -ge 4 ]
}

# --------------------------------------------------------------------------
# Test 39 (US-030): No hardcoded species placeholder remains
# --------------------------------------------------------------------------
test_030_no_hardcoded_species() {
    local impl="$GROMACS_SRC/src/gromacs/applied_forces/qmmm/dftbforceprovider.cpp"
    if grep -q 'species_\.resize(nQMAtoms_, 0)' "$impl"; then
        return 1
    fi
    return 0
}

# --------------------------------------------------------------------------
# Test 40 (US-030): Species mapping logic in constructor, not calculateForces
# --------------------------------------------------------------------------
test_030_species_build_in_constructor() {
    local impl="$GROMACS_SRC/src/gromacs/applied_forces/qmmm/dftbforceprovider.cpp"
    # atomNumbers_ should appear before calculateForces definition
    local first_atomnum_line
    first_atomnum_line=$(grep -n 'atomNumbers_' "$impl" | head -1 | cut -d: -f1)
    local calc_forces_line
    calc_forces_line=$(grep -n 'DFTBForceProvider::calculateForces' "$impl" | head -1 | cut -d: -f1)
    [ "$first_atomnum_line" -lt "$calc_forces_line" ]
}

for t in \
    test_030_species_from_atomic_numbers \
    test_030_thdd_annotations \
    test_030_no_hardcoded_species \
    test_030_species_build_in_constructor \
; do
    echo "Running: $t"
    if "$t"; then
        report "$t" 0
    else
        report "$t" 1
    fi
done

# US-031 tests
# --------------------------------------------------------------------------
# Test 41 (US-031): modifyQMMMTopology() has no DFTBPLUS guard
# VC-1: Preprocessing must run for all QM methods including DFTBPLUS
# --------------------------------------------------------------------------
test_031_no_dftbplus_guard_in_modify_topology() {
    local options_file="$GROMACS_SRC/src/gromacs/applied_forces/qmmm/qmmmoptions.cpp"
    # Extract the modifyQMMMTopology function body (from definition to next function)
    local body
    body="$(sed -n '/^void QMMMOptions::modifyQMMMTopology/,/^void QMMMOptions::/p' "$options_file" | head -n -1)"
    # Verify the function exists
    [ -n "$body" ] || return 1
    # Verify no DFTBPLUS guard inside the function body
    if echo "$body" | grep -q 'DFTBPLUS'; then
        return 1
    fi
    return 0
}

# --------------------------------------------------------------------------
# Test 42 (US-031): DFTBForceProvider does not call any preprocessing functions
# VC-2: No runtime preprocessing duplication
# --------------------------------------------------------------------------
test_031_no_runtime_preprocessing() {
    local impl="$GROMACS_SRC/src/gromacs/applied_forces/qmmm/dftbforceprovider.cpp"
    for func in splitEmbeddedBlocks removeEmbeddedClassicalCharges \
                addEmbeddedNBExclusions modifyEmbeddedTwoCenterInteractions \
                modifyEmbeddedThreeCenterInteractions modifyEmbeddedFourCenterInteractions \
                buildEmbeddedAtomNumbers buildLinkFrontier checkConstrainedBonds; do
        if grep -q "$func" "$impl"; then
            return 1
        fi
    done
    return 0
}

# --------------------------------------------------------------------------
# Test 43 (US-031): DFTBForceProvider accesses preprocessed data
# VC-3: Constructor uses atomNumbers_ and qmIndices_ from QMMMParameters
# --------------------------------------------------------------------------
test_031_preprocessed_data_accessible() {
    local impl="$GROMACS_SRC/src/gromacs/applied_forces/qmmm/dftbforceprovider.cpp"
    local header="$GROMACS_SRC/src/gromacs/applied_forces/qmmm/dftbforceprovider.h"
    # Constructor receives QMMMParameters by const ref
    grep -q 'const QMMMParameters&' "$header" || return 1
    # Accesses atomNumbers_ and qmIndices_
    grep -q 'parameters_.atomNumbers_' "$impl" || return 1
    grep -q 'parameters_.qmIndices_' "$impl" || return 1
    return 0
}

# --------------------------------------------------------------------------
# Test 44 (US-031): processCoordinates DFTBPLUS guard is separate from modifyQMMMTopology
# VC-4: The early return in processCoordinates does NOT skip preprocessing
# --------------------------------------------------------------------------
test_031_process_coords_separate_from_topology() {
    local options_file="$GROMACS_SRC/src/gromacs/applied_forces/qmmm/qmmmoptions.cpp"
    # processCoordinates and modifyQMMMTopology must be separate functions
    local pc_line mt_line
    pc_line=$(grep -n 'QMMMOptions::processCoordinates' "$options_file" | head -1 | cut -d: -f1)
    mt_line=$(grep -n 'QMMMOptions::modifyQMMMTopology' "$options_file" | head -1 | cut -d: -f1)
    # Both must exist
    [ -n "$pc_line" ] && [ -n "$mt_line" ] || return 1
    # They must be different functions (different line numbers)
    [ "$pc_line" -ne "$mt_line" ] || return 1
    # processCoordinates has a DFTBPLUS reference (early return for CP2K input)
    local pc_body
    pc_body="$(sed -n "${pc_line},${mt_line}p" "$options_file")"
    echo "$pc_body" | grep -q 'DFTBPLUS' || return 1
    return 0
}

# --------------------------------------------------------------------------
# Test 45 (US-031): ThDD annotation for preprocessing reuse
# --------------------------------------------------------------------------
test_031_thdd_annotations() {
    local impl="$GROMACS_SRC/src/gromacs/applied_forces/qmmm/dftbforceprovider.cpp"
    grep -q 'ThDD:T-US-031-1.3' "$impl"
}

for t in \
    test_031_no_dftbplus_guard_in_modify_topology \
    test_031_no_runtime_preprocessing \
    test_031_preprocessed_data_accessible \
    test_031_process_coords_separate_from_topology \
    test_031_thdd_annotations \
; do
    echo "Running: $t"
    if "$t"; then
        report "$t" 0
    else
        report "$t" 1
    fi
done

# ==========================================================================
# US-032: QM Energy Bookkeeping on Main Rank (Tests 46-51)
# ==========================================================================

# --------------------------------------------------------------------------
# Test 46 (US-032): Energy written to QuantumMechanicalRegionEnergy term
# --------------------------------------------------------------------------
test_032_energy_correct_term() {
    local impl="$GROMACS_SRC/src/gromacs/applied_forces/qmmm/dftbforceprovider.cpp"
    grep -q 'QuantumMechanicalRegionEnergy.*+=.*energyKjMol' "$impl"
}

# --------------------------------------------------------------------------
# Test 47 (US-032): Energy write guarded by isMainRank()
# --------------------------------------------------------------------------
test_032_main_rank_guard() {
    local impl="$GROMACS_SRC/src/gromacs/applied_forces/qmmm/dftbforceprovider.cpp"
    # Extract the block around QuantumMechanicalRegionEnergy and verify isMainRank guard
    local energy_line
    energy_line=$(grep -n 'QuantumMechanicalRegionEnergy' "$impl" | head -1 | cut -d: -f1)
    [ -n "$energy_line" ] || return 1
    # Check that isMainRank() appears within 5 lines before the energy write
    local start=$((energy_line - 5))
    [ "$start" -lt 1 ] && start=1
    sed -n "${start},${energy_line}p" "$impl" | grep -q 'isMainRank()'
}

# --------------------------------------------------------------------------
# Test 48 (US-032): Energy conversion uses GROMACS constants
# --------------------------------------------------------------------------
test_032_conversion_constants() {
    local impl="$GROMACS_SRC/src/gromacs/applied_forces/qmmm/dftbforceprovider.cpp"
    # Verify c_hartree2Kj and c_avogadro used in energy conversion
    grep -q 'c_hartree2Kj.*c_avogadro' "$impl" || return 1
    return 0
}

# --------------------------------------------------------------------------
# Test 49 (US-032): Energy accumulated with += (not =)
# --------------------------------------------------------------------------
test_032_energy_accumulation() {
    local impl="$GROMACS_SRC/src/gromacs/applied_forces/qmmm/dftbforceprovider.cpp"
    # Must use += for QuantumMechanicalRegionEnergy
    grep 'QuantumMechanicalRegionEnergy' "$impl" | grep -q '+='
}

# --------------------------------------------------------------------------
# Test 50 (US-032): Pattern matches CP2K provider
# --------------------------------------------------------------------------
test_032_matches_cp2k_pattern() {
    local dftb="$GROMACS_SRC/src/gromacs/applied_forces/qmmm/dftbforceprovider.cpp"
    local cp2k="$GROMACS_SRC/src/gromacs/applied_forces/qmmm/qmmmforceprovider.cpp"
    # Both must use isMainRank() guard
    grep -q 'isMainRank()' "$dftb" || return 1
    grep -q 'isMainRank()' "$cp2k" || return 1
    # Both must write to QuantumMechanicalRegionEnergy
    grep -q 'QuantumMechanicalRegionEnergy' "$dftb" || return 1
    grep -q 'QuantumMechanicalRegionEnergy' "$cp2k" || return 1
    # Both must use += accumulation on that term
    grep 'QuantumMechanicalRegionEnergy' "$dftb" | grep -q '+=' || return 1
    grep 'QuantumMechanicalRegionEnergy' "$cp2k" | grep -q '+=' || return 1
    # Both must use c_hartree2Kj * c_avogadro
    grep -q 'c_hartree2Kj.*c_avogadro' "$dftb" || return 1
    grep -q 'c_hartree2Kj.*c_avogadro' "$cp2k" || return 1
    return 0
}

# --------------------------------------------------------------------------
# Test 51 (US-032): ThDD annotations for energy bookkeeping
# --------------------------------------------------------------------------
test_032_thdd_annotations() {
    local impl="$GROMACS_SRC/src/gromacs/applied_forces/qmmm/dftbforceprovider.cpp"
    grep -q 'ThDD:T-US-032-1.3' "$impl" || return 1
    grep -q 'ThDD:T-US-032-1.4' "$impl" || return 1
    return 0
}

for t in \
    test_032_energy_correct_term \
    test_032_main_rank_guard \
    test_032_conversion_constants \
    test_032_energy_accumulation \
    test_032_matches_cp2k_pattern \
    test_032_thdd_annotations \
; do
    echo "Running: $t"
    if "$t"; then
        report "$t" 0
    else
        report "$t" 1
    fi
done

echo ""
echo "=== Results: $PASS/$TOTAL passed, $FAIL failed ==="
[ "$FAIL" -eq 0 ]
