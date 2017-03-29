# Generate then display mesh (minimum area -a in square units, minimum angle -q in degrees)
exec ./triangle -p -a0.00001 -q30 circle.poly
exec ./showme circle.1.ele
