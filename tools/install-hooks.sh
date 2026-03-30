#!/usr/bin/env bash
# Install git hooks for GRO-DFTB local CI
# Usage: tools/install-hooks.sh

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
HOOKS_DIR="${REPO_ROOT}/.git/hooks"

if [ ! -d "$HOOKS_DIR" ]; then
    echo "Error: .git/hooks/ not found. Are you in a git repository?" >&2
    exit 1
fi

# Install pre-push hook
cat > "${HOOKS_DIR}/pre-push" << 'HOOK'
#!/usr/bin/env bash
# GRO-DFTB pre-push hook: runs local CI before allowing push.
# Bypass with: git push --no-verify

REPO_ROOT="$(git rev-parse --show-toplevel)"
CI_SCRIPT="${REPO_ROOT}/tools/ci.sh"

if [ ! -x "$CI_SCRIPT" ]; then
    echo "Warning: tools/ci.sh not found or not executable — skipping CI check" >&2
    exit 0
fi

echo "Running local CI before push..."
if "$CI_SCRIPT"; then
    echo "CI passed — push allowed."
    exit 0
else
    echo "CI failed — push blocked. Fix issues and retry." >&2
    echo "Bypass with: git push --no-verify" >&2
    exit 1
fi
HOOK

chmod +x "${HOOKS_DIR}/pre-push"
echo "Installed pre-push hook at ${HOOKS_DIR}/pre-push"
