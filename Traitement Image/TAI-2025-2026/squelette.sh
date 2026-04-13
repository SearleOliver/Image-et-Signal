#!/bin/bash

if [ $# -ne 2 ]; then
    echo "Usage: $0 input.pbm output.pbm"
    exit 1
fi

INPUT=$1
OUTPUT=$2
CURRENT="__current_$$.pbm"
PREV="__prev_$$.pbm"

cp "$INPUT" "$CURRENT"

while true; do
    cp "$CURRENT" "$PREV"
    
    for i in 1 2 3 4 5 6 7 8; do
        ./thinning "$CURRENT" "b${i}.mx" "__tmp_$$.pbm"
        mv "__tmp_$$.pbm" "$CURRENT"
    done
    
    if ./comparison "$CURRENT" "$PREV"; then
        break
    fi
done

cp "$CURRENT" "$OUTPUT"
rm -f "$CURRENT" "$PREV"

echo "Skeleton saved to $OUTPUT"