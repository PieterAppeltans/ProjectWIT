# Generate then display mesh (minimum area 7 square units, minimum angle 30 degrees)
exec ./../triangle/triangle -p -a0.000005 -q30 ../triangle/circle
exec ./../triangle/showme ../triangle/pear.1
