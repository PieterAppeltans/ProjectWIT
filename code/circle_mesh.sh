# Generate then display mesh (minimum area 7 square units, minimum angle 30 degrees)
exec ./../triangle/triangle -p -a0.5 -q30 ../triangle/pear
exec ./../triangle/showme ../triangle/pear.1
