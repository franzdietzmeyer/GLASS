#!/usr/bin/env bash
# First-time setup for GLASS: Python env + JDK 25 + Nextflow
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
VENV_DIR="$REPO_ROOT/.venv"
JDK_BASE="$REPO_ROOT/jdk"
NF_BIN="$VENV_DIR/bin/nextflow"

echo "[GLASS] Step 1/3: Creating Python environment and installing dependencies..."
uv venv
uv sync

echo "[GLASS] Step 2/3: Installing JDK 25..."
if ls "$JDK_BASE"/jdk-* > /dev/null 2>&1; then
    echo "[GLASS] JDK already present at $JDK_BASE, skipping."
else
    "$VENV_DIR/bin/python" -c "import jdk; jdk.install('25', path='$JDK_BASE')"
fi

JAVA_HOME="$("$VENV_DIR/bin/python" -c "
import jdk, os, glob
d = sorted(glob.glob(os.path.join('$JDK_BASE', 'jdk-*')))[0]
system=jdk.OS
print(os.path.join(d, 'Contents', 'Home') if 'mac' in system else d)
")"

export JAVE_HOME=$JAVA_HOME
echo $JAVA_HOME
echo "[GLASS] Step 3/3: Installing Nextflow..."
if [ -x "$NF_BIN" ]; then
    echo "[GLASS] Nextflow already present at $NF_BIN, skipping."
else
    JAVA_HOME="$JAVA_HOME" \
        bash -c "cd \"$VENV_DIR/bin\" && curl -fsSL https://get.nextflow.io | bash"
fi

echo "[GLASS] Setup complete. Run: ./run_glass_nextflow.sh local"
