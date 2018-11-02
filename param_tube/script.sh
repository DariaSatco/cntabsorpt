#!/bin/bash
filename="tube.chiralities"
while read -r n m
do
    mv tube.param.* tube.param.$n$m
    echo "$n$m" > fort.input
    echo "300.00" > tube.param.$n$m
    echo "$n $m" >> tube.param.$n$m
    echo "100" >> tube.param.$n$m
    echo "0.0" >> tube.param.$n$m
    echo "0.0" >> tube.param.$n$m
    echo "1.3" >> tube.param.$n$m
    echo "1.0" >> tube.param.$n$m
    echo "501" >> tube.param.$n$m
    echo "0.0 4.0" >> tube.param.$n$m
    echo "0.05" >> tube.param.$n$m
    echo "501 -10.0 15.0" >> tube.param.$n$m
    echo "90.0" >> tube.param.$n$m
    ../cntabsorpt.out < fort.input > fort.output
done < "$filename"
