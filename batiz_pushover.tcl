model basic -ndm 2 -ndf 3
file mkdir _results
set width    450
set height   200
node  1       0.0     0.0
node  2    $width     0.0
node  3       0.0 $height
node  4    $width $height
fix   1     1    1    1
fix   2     1    1    1

# Core concrete (confined)  tag -- fc -- ec0 --- fcu ---- ecu
uniaxialMaterial Concrete01  1  -300.0  -0.002   -300.0     -0.004
# Cover concrete (unconfined)
uniaxialMaterial Concrete01  2  -300.0   -0.002   -310.0     -0.004

#                        tag  fy Young    hardening coefficient
uniaxialMaterial Steel01  3  5000. 2.04e6 0.01

set colWidth 30.0
set colDepth 40.0
set cover  2.0
set area_steel    2.85
# some variables derived from the parameters
set y1 [expr $colDepth/2.0]
set z1 [expr $colWidth/2.0]

section Fiber 1 {

    # Create the concrete core fibers
    patch rect 1 10 1 [expr $cover-$y1] [expr $cover-$z1] [expr $y1-$cover] [expr $z1-$cover]

    # Create the concrete cover fibers (top, bottom, left, right)
    patch rect 2 10 1  [expr -$y1] [expr $z1-$cover] $y1 $z1
    patch rect 2 10 1  [expr -$y1] [expr -$z1] $y1 [expr $cover-$z1]
    patch rect 2  2 1  [expr -$y1] [expr $cover-$z1] [expr $cover-$y1] [expr $z1-$cover]
    patch rect 2  2 1  [expr $x1-$cover] [expr $cover-$z1] $y1 [expr $z1-$cover]

    # Create the reinforcing fibers (left, middle, right)
    layer straight 3 3 $area_steel [expr $y1-$cover] [expr $z1-$cover] [expr $y1-$cover] [expr $cover-$z1]
    layer straight 3 2 $area_steel 0.0 [expr $z1-$cover] 0.0 [expr $cover-$z1]
    layer straight 3 3 $area_steel [expr $cover-$y1] [expr $z1-$cover] [expr $cover-$y1] [expr $cover-$z1]

}

geomTransf PDelta 1

# Number of integration points along length of element
set np 5

# Create the coulumns using Beam-column elements
#               e            tag ndI ndJ nsecs secID transfTag
set eleType forceBeamColumn
element $eleType  1   1   3   $np    1       1
element $eleType  2   2   4   $np    1       1

# Define beam elment
# -----------------------------
geomTransf Linear 2
# Create the beam element
#                          tag ndI ndJ     A       E    Iz   transfTag
element elasticBeamColumn   3   3   4    360e6    2.04e6  133333.  2

set P 0.0

pattern Plain 1 "Linear" {

	#    nd    FX          FY  MZ
	load  3   0.0  [expr -$P] 0.0
	load  4   0.0  [expr -$P] 0.0
}
system BandGeneral
constraints Transformation
numberer RCM
test NormDispIncr 1.0e-12  10 3
algorithm Newton
integrator LoadControl 0.1
analysis Static
analyze 10

print node 3 4
print ele 1

puts "Gravity Analysis Completed"

loadConst -time 0.0
set H 10.0;		# Reference lateral load
pattern Plain 2 "Linear" {
     # Create nodal loads at nodes 3 & 4
     #    nd    FX  FY  MZ
     load 3 $H 0.0 0.0
     load 4 $H 0.0 0.0
}

set disp_incr 0.1;	        # Displacement increment

#                             node dof init Jd min max
integrator DisplacementControl  3   1   $disp_incr  1 $disp_incr $disp_incr

recorder Node -file _results/disps.out -node 3 -dof 1 disp
recorder Element -file _results/col_forces.out -ele 1 -dof 1 force
set max_disp 10.0
set currentDisp 0.0;
set ok 0

while {$ok == 0 && $currentDisp < $max_disp} {

	set ok [analyze 1]

	if {$ok != 0} {
	    puts "modified newton"
	    test NormDispIncr 1.0e-12  1000
	    algorithm ModifiedNewton -initial
	    set ok [analyze 1]
	    if {$ok == 0} {puts "regular newton"}
	    test NormDispIncr 1.0e-12  10
	    algorithm Newton
	}

	set currentDisp [nodeDisp 3 1]
}


if {$ok == 0} {
  puts "Pushover SUCCESSFUL";
} else {
  puts "Pushover FAILED";
}

print node 3
