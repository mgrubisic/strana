# http://opensees.berkeley.edu/wiki/index.php/Moment_Curvature_Example
# 3 same bars per side, total 8 bars
model BasicBuilder -ndm 2 -ndf 3
# Core concrete (confined) tag fc ec fcu ecu
uniaxialMaterial Concrete01  1  -22.0e3  -0.003   -22e3    -0.0038
uniaxialMaterial Concrete01  2  -22.0e3  -0.003   -22e3    -0.0038
# Cover concrete (unconfined)
set axial_load -0e3
set peralte 50.0
set ancho 50.0
set recubrimiento 12.0
set area_steel 1.27
set fy 420.0
set Young 2.039e5
uniaxialMaterial Steel01  3  $fy $Young 0.01
set y2 [expr $peralte/2.0]
set y1 [expr $ancho/2.0]

section Fiber 1 {
    # Create the concrete core fibers
    patch rect 1 10 1 [expr $recubrimiento-$y2] [expr $recubrimiento-$y1] [expr $y2-$recubrimiento] [expr $y1-$recubrimiento]
    # Create the concrete cover fibers (top, bottom, left, right)
    patch rect 2 10 1  [expr -$y2] [expr $y1-$recubrimiento] $y2 $y1
    patch rect 2 10 1  [expr -$y2] [expr -$y1] $y2 [expr $recubrimiento-$y1]
    patch rect 2  2 1  [expr -$y2] [expr $recubrimiento-$y1] [expr $recubrimiento-$y2] [expr $y1-$recubrimiento]
    patch rect 2  2 1  [expr $y2-$recubrimiento] [expr $recubrimiento-$y1] $y2 [expr $y1-$recubrimiento]
    # Create the reinforcing fibers (left, middle, right)
    layer straight 3 3 $area_steel [expr $y2-$recubrimiento] [expr $y1-$recubrimiento] [expr $y2-$recubrimiento] [expr $recubrimiento-$y1]
    layer straight 3 2 $area_steel 0.0 [expr $y1-$recubrimiento] 0.0 [expr $recubrimiento-$y1]
    layer straight 3 3 $area_steel [expr $recubrimiento-$y2] [expr $y1-$recubrimiento] [expr $recubrimiento-$y2] [expr $recubrimiento-$y1]

}

set yield_curvature [expr $fy/(0.7*($peralte - $recubrimiento))/$Young]
set fardis_curvature [expr 1.7*$fy/$peralte/$Young]
puts "Estimated yield curvature: $yield_curvature"
puts "Fardis yield curvature: $fardis_curvature"
set target_ductility 8

proc Concrete_rect {tag b h top_cover side_cover top_steel_ratio bot_steel_ratio confined_tag unconfined_tag steel_tag} \
{
    
}

proc MomentCurvature {section_tag load max_curvature {increment 100}} {
    node 1 0. 0. 0.
    node 2 0. 0. 0.
    fix 1 1 1 1
    fix 2 0 1 0
    element zeroLengthSection 1 1 2 $section_tag
    recorder Node -file mk.results.out -time -node 2 -dof 3 disp
    integrator LoadControl 0.0
    system BandGeneral
    test NormUnbalance 1.0e-9 10
    numberer Plain
    algorithm Newton
    analysis Static
    constraints Plain
    pattern Plain 1 "Constant" <-fact cFactor> {
        load 2 $load 0. 0.
    }
    analyze 1
    pattern Plain 2 "Linear" {
        load 2 0. 0. 1.0
    }
    set curvature [expr $max_curvature/$increment]
    integrator DisplacementControl 2 3 $curvature
    analyze $increment
}

MomentCurvature 1 $axial_load [expr $yield_curvature*$target_ductility]
recorder Plot moment_curv_results.txt "moment-curvature" 10 10 400 400 -columns 1 2
