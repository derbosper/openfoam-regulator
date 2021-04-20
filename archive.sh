#!/bin/bash

if [ $# -eq 0 ]
  then
    echo "Usage: ./archive.sh arc_name.tar.xz target"
fi

tar -v -c -I "xz -8 -T0" -f "$@"