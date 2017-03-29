# Generate then display mesh (minimum area -a in square units, minimum angle -q in degrees)
exec ./triangle -p -a3 -q30 pear.poly
exec ./showme pear.1.ele
