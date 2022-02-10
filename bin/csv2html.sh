#!/bin/bash

input=$1

echo "<table>"
sed "s/,/<\/td><td>/g" $input | sed "s/^/<tr><td>/g" | sed "s/Capture/Capture<\/td><\/tr>/g" | sed "s/feature_type/feature_type<\/td><\/tr>/g" | sed '1 s/td/th/g'
echo "</table>"


