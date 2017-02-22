# Generate then display mesh (minimum area 7 square units, minimum angle 30 degrees)
exec ./triangle -p -a7 -q30 pear
exec ./showme pear.1.ele
