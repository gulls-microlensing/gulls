#!/bin/bash
# Basic simulation example runner

set -e  # Exit on any error

echo "Running basic Gulls simulation example..."

# Check if Gulls is built
if [ ! -f "../../build/bin/gulls_std.x" ]; then
    echo "Error: Gulls not built. Run 'cmake --build build' from the main directory."
    exit 1
fi

# Check if ESPL table exists
if [ ! -f "../../src/ESPL.tbl" ]; then
    echo "Error: ESPL.tbl not found. Copy it to src/ directory."
    exit 1
fi

# Set environment variables
export GULLS_BASE_DIR="$(pwd)/../../"
export GULLS_STARS_DIR="$(pwd)/../../"

# Run the simulation
echo "Running simulation..."
../../build/bin/gulls_std.x basic_simulation.prm

echo "Simulation complete!"
echo "Check the output files:"
ls -la *.out *.lc *.png 2>/dev/null || echo "No output files found"
