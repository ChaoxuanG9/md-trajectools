#!/bin/bash

# Run residence time analysis using a config file
# Make sure you're in the top-level md-trajectools/ directory

echo "Running residence time analysis..."
python residence_time/residence_time.py --config examples/residence_time/config.yaml

echo "Done. Results saved to 'stickytime/' (default output directory)."
