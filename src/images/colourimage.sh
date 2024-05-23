#!/bin/bash

if [ $# -ne 3 ]; then
  echo "Usage: $0 <detectorlist> <field> <output root>"
  exit
fi

~/Documents/Roman/gulls/src/images/colourImage $1 0 ~/Documents/Roman/input/faint/out-$2 0.000020 ~/Documents/Roman/input/moderate/out-$2 0.000020 ~/Documents/Roman/input/bright/out-$2 0.001000 $3
