#!/bin/sh

g++ raytrace.cpp -o raytrace

dir=testFiles/
ls testFiles/ > B

cat B | while read line
do
    ./raytrace $dir$line
done

#g++ raytrace.cpp -o raytracer




