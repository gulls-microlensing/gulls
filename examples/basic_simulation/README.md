# Basic Simulation Example

This example demonstrates how to run a simple microlensing simulation with Gulls.

## Files

- `basic_simulation.prm` - Parameter file
- `run.sh` - Script to run the simulation
- `expected_output/` - Expected output files

## Quick Start

1. Build Gulls (see main README)
2. Run the simulation:
   ```bash
   ./run.sh
   ```

## What This Example Does

This example simulates a single microlensing event with:
- A single source star
- A single lens object
- Basic observing conditions
- Standard output (lightcurve, detection statistics)

## Expected Output

After running, you should see:
- `basic_simulation_0_0.out` - Event summary
- `basic_simulation_0_0_0.all.lc` - Lightcurve data
- `basic_simulation_0_0_0.all_plot.png` - Lightcurve plot

## Customization

Edit `basic_simulation.prm` to:
- Change the number of events (`SUBRUNSIZE`)
- Modify observing conditions
- Adjust physical parameters
