# CDMX 2016 Carlo Ruiz iiUNAM
# procedures for building regular simple rc frames in OpenSees
# SI units

proc Nodebuilder {storeys bays} {
    # storeys & bays format ; [list 0 len1 len2 ... lenn]
    # builds tags for nodes (inelastic 2d rigid-diaph frames) in dict: nodes{st0 {j0 {c:tag1 r:tag2.. u l d} j1 {...} st1 ...}}
    # num convention: down-up, left-right , 0-n..
    # master node: center (center)
    # slaves : r,u,l,d
    #       up
    # left  (center) right
    #       down
    # returns: dict
    # nodesdict(st0 {j0(c u l..) j1(c u l...)} st1 {j..})
    # dict with storeys containing a list of dicts with tags for center, u d l nodes
    # no need to remember tags as plain numbers, just access the jointdict and use a key to get the node's tag
    set tag 0
    for {set st 0} {$st < [llength $storeys]} {incr st} {
        for {set bay 0} {$bay < [llength $bays]} {incr bay} {
            # set x [lindex $storeys $st]
            # set y [lindex $bays $bay]
            set x [expr ([join [lrange $storeys 0 $st] +])]
            set y [expr ([join [lrange $bays 0 $bay] +])]
            puts "$x  :  $y"
            set j$bay [dict create c [incr tag]]
            node $tag $x $y

            if {$st == 0} {
                set directions [list u]
            } elseif {$st == [expr [llength $storeys] - 1]} {
                if {$bay == 0} {
                    set directions [list r d]
                } elseif {$bay == [expr [llength $bays] - 1]} {
                    set directions [list l d]
                } else {
                    set directions [list r l d]
                }
            } else {
            if {$bay == 0} {
                set directions [list r u d]
            } elseif {$bay == [expr [llength $bays] - 1]} {
                set directions [list u l d]
            } else {
                set directions [list r u l d]
            }
            }

            foreach dir $directions {
                dict append j$bay $dir [incr tag]
                node $tag $x $y
                equalDOF [dict get [set j$bay] c] [dict get [set j$bay] $dir] 1 2
            }
            if {$bay>0} {
                equalDOF [dict get $j0 c] [dict get [set j$bay] c] 1
            }
            dict lappend nodes st$st [set j$bay]
        }
    }
    # return [array get nodes]
    return $nodes

}

# now we need a proc for building a generic frame out of ibarra
# modified IMK model (uniaxialMaterial ModIMKPinching), consisting in two options:
# computing alpha, theta, lambda values for every type of column using the equations derived from lignos or
# read database as list for every element and pass the arguments as list elements and assign to $alpha, etc for every rot spring

# ideas for proc simple
# return
# beams, cols and spring tags should be retained in three separate dicts or single dict
# with same coords as Nodebuilder cols-dict{st0[c0,c1,c2]}
# beamdict{st1{b0,b1,b2..}} springdict{st0{b0{left right}}}
# coldict{st0[c0{tag spring-down, spring-up}]}
# beamdict{st0[b1{tag, spring-left, spring-right}]}

# input should could be coldict[st0{c0 c1 c2...}] with properties on the elements
# if [coldict lindex counter] is empty , copy col-default (c0) into other PROPERTIEs
# if that property is empty then use default-col or ask for col-dimentiosn into prompt
# ideally we would want /smart-builder-proc/ called for every column(analysis results)

# we wish a dict[cols] with all the info for the model, same for beams_dict

# -uniaxial material
# mat tag
# elastic mod
# elastic mod material
# column - beam
# -area
# -youngs modulus
# -shear mod
# -torsional inertia
# -Iz
# -Iy
# shear area y axis
# shear area z axisp
# ptransformation tag
# mass density


# zero length
# uniaxial material
# ij nodes


# span
# shear_span
# aspect_ratio
# stirr_yield_str
# conc_compr_str
# steel_elast_mod
# tension_rebar_ratio
# comp_rebar_ratio
# conc_cover_4dir
# axial_load
# stirr_slippage [0,1]


# #computed

# beam_width
# rebar_tension_area
# rebar_comp_area
# conc_gross_area
# conc_elast_mod
# steel_area (goal)
# steel_area_real
# tension_rebar_ratio_real
# eff_depth
# stirr_area
# stirr_ratio_trans
# stirr_ratio_long
# stirr_ratio_eff
# axial_load_ratio
# tension_rebar_diameter
# comp_rebar_diameter
# stirr_diameter
# bars_tension
# bar_size_tension
# bars_comp
# bar_size_tension
# bars_stirrs
# bar_size_stirrs
# eff_width
# rebar_centroid_tension
# rebar_centroid_comp
# rebar_spacing
# delta_prime (fardis)
# unbalanced_reinf_factor

# #output
# rebar_ratio_error
# moment_error
# moment_yield
# moment_cap
# rot_yield
# rot_plastic_precap
# rot_postcap
# rot_ultimate
# elast_stiff_?? real
# brute_conc_stiff
# hardening_stiff
# softening_stiff
# coeff_hardening_stiff
# coeff_softening_stiff
# confinement_eff_factor
# energy_dissipation_capacity
# shear_capacity_conc

# proc conc_beam_builder {rect_or_squared desired_inertia} {
# }
# proc conc_simple_beam_build {rect_or_squared desired_cracked_inertia} {
# }

proc Steelbuilder {tensionArea tensionStr tensionElastMod compArea compStr compElastMod webArea webStr webElastMod} {
    # dict with rebar properties

    dict append steeldict tensionArea $tensionArea
    dict append steeldict tensionStr $tensionStr
    dict append steeldict tensionElastMod $tensionElastMod

    dict append steeldict compArea $compArea
    dict append steeldict compStr $compStr
    dict append steeldict compElastMod $compElastMod

    dict append steeldict webArea $webArea
    dict append steeldict webStr $webStr
    dict append steeldict webElastMod $webElastMod

}

proc FardisHaseltonprop {coldict materialdict webSpacing axial_load slippageYESNO} {

}

proc Colbuilder {height width tensionCentroid compCentroid} {
    # returns dict with geometrical properties
    # 2 is around an axis parallel to height
    # 3 is around an axis parallel to width
    set BruteArea [expr $height * $width]
    set BruteInertia_3 [expr { $width * pow($height, 3) / 12 }]
    set depth [expr {$height - $tensionCentroid}]
    set effArea [expr {$depth*$width}]
    set effInertia_3 [expr {$width * pow($depth,3)/12}]
    set BruteInertia_2 [expr {$height * pow($width, 3)/12}]
    set effInertia_2 [expr {$depth * pow($width, 3)/12}]

    dict append coldict BruteArea $BruteArea
    dict append coldict effArea $effArea

    dict append coldict BruteInertia_3 $BruteInertia_3
    dict append coldict BruteInertia_2 $BruteInertia_2
    dict append coldict effInertia_3 $effInertia_3
    dict append coldict effInertia_2 $effInertia_2

    dict append coldict height $height
    dict append coldict width $width
    dict append coldict depth $depth
    dict append coldict tensionCentroid $tensionCentroid
    dict append coldict compCentroid $compCentroid

    return $coldict
}

proc Framebuilder {nodebuilder_dict coldict beamdict mass_list} {
    # builds a frame with uniform storey col and beam sections
    # includes p-Delta effects
    # nodebuilder_dict is from Nodebuilder proc
    # \coldims_list : {{width-col1,height-col1},{wc2,hc2},...{wc-n,hc-n}} frame having same storey config
    # on all bays, beams follow the same format
    # masses follow same format [list m1 m2 m3 m_ndof]
    # conc_comp_str is nominal f'c for concrete

    if {[expr [llength $coldims_list] + 1] != [dict size $nodebuilder_dict]} {
        puts "procERROR: num of storeys in col_list != dof"
    }

    if {[expr [llength $beamdims_list] + 1] != [dict size $nodebuilder_dict]} {
        puts "procERROR: num of storeys in beam_list != dof"
    }
    if {[expr [llength $mass_list] + 1] != [dict size $nodebuilder_dict]} {
        puts "procERROR: number of masses != dof"
    }
    set pdeltatag 0
    set tag 0
    set mattag 0
    set conc_elast_mod [expr 14000 * sqrt($conc_comp_str)]
    uniaxialMaterial Elastic $mattag $conc_elast_mod
    geomTransf PDelta $pdeltatag
    for {set st 0} {$st < [dict size $nodebuilder_dict]} {incr st} {

        set conc_col [Conc_col $conc_comp_str [lindex [lindex $coldims_list $st] 1] [lindex [lindex $coldims_list $st] 0]]
        set conc_beam [Conc_col $conc_comp_str [lindex [lindex $beamdims_list $st] 1] [lindex [lindex $beamdims_list $st] 0]]
        puts "elements calculated"
        mass [dict get [lindex [dict get $nodebuilder_dict st[expr $st + 1]] 0] c] [lindex $mass_list $st] 0 0
        puts "mass created and assigned"

        for {set col 0} {$col < [llength $coldims_list]} {incr col} {
            set cnode_i [dict get [lindex [dict get $nodebuilder_dict st$st] $col] c]
            set cnode_j [dict get [lindex [dict get $nodebuilder_dict st[expr $st + 1]] $col] c]
            set dnode [dict get [lindex [dict get $nodebuilder_dict st[expr $st + 1]] $col] d]
            set unode [dict get [lindex [dict get $nodebuilder_dict st$st] $col] u]
            puts "joint nodal tags obtained successfully."
            element zeroLength [incr tag] $cnode_i $unode -mat $mattag -dir 6 ;# [dict get [lindex [dict get $nodebuilder_dict st$st] j$col] c]
            element zeroLength [incr tag] $cnode_j $dnode -mat $mattag -dir 6
            element elasticBeamColumn [incr tag] [expr $tag - 2] [expr $tag -1] [dict get $conc_col area] [dict get $conc_col elast_mod] [dict get $conc_col inertia] $pdeltatag
            # possible conceptual error in dicts here. reread how dict should be arranged
            # set current_col [lindex [dict get $nodebuilder_dict st$st] $col]
            # dict append $current_col col $tag
            # dict append $current_col zerolen_i [expr $tag - 2]
            # dict append $current_col zerolen_j [expr $tag - 1]
            # dict lappend cols_dict b$0 $tag
            dict set [dict lappend colsdict st$st col$col] col $tag
            dict set [dict lappend colsdict st$st col$col] zerolen_i [expr $tag - 2]
            dict set [dict lappend colsdict st$st col$col] zerolen_j [expr $tag - 1]
            puts "col : $col created with zerolen elements  zl_i : [expr $tag - 2] zl_j: [expr $tag - 1]"
            puts $colsdict
        }

            for {set beam 0} {$beam < [llength [dict get $nodebuilder_dict [lindex [dict keys $nodebuilder_dict 0]]]]} {incr beam} {
            set cnode_i [dict get [lindex [dict get $nodebuilder_dict st[expr $st + 1]] $beam] c]
            set cnode_j [dict get [lindex [dict get $nodebuilder_dict st[expr $st + 1]] [expr $beam + 1]] c]
            set lnode [dict get [lindex [dict get $nodebuilder_dict st[expr $st + 1]] [expr $beam + 1]] l]
            set rnode [dict get [lindex [dict get $nodebuilder_dict st[expr $st + 1]] $beam] r]

                element zeroLength [incr tag] $cnode_i $rnode -mat $mattag -dir 6
                element zeroLength [incr tag] $cnode_j $lnode -mat $mattag -dir 6
            element elasticBeamColumn [incr tag] [expr $tag - 2] [expr $tag - 1] [dict get $conc_beam area] [dict get $conc_beam conc_elast_mod] [dict get $conc_beam inertia] $pdeltatag

                dict set [dict lappend beamsdict st$st beam$beam] beam $tag
                dict set [dict lappend beamsdict st$st beam$beam] zerolen_i [expr $tag - 2]
                dict set [dict lappend beamsdict st$st beam$beam] zerolen_j [expr $tag - 1]
            }
    }
    return [list $colsdict $beamsdict]
}


proc ConcRect {concStr width height tensionCentroid compCentroid} {
    # concStr in kg/cm2
    # builds a dictionary with all parameters as kwords
    # assumes ultimate linear strain of conc = 0.003
    # units: MPa, mm.
    set concElastmod [expr 14e3 * sqrt($concStr) * 9.81/100]; # 14000 f'c -> MPa
    set concStr [expr $concStr * 9.81/100]
    set depth [expr $height - $tensionCentroid]
    set depthTension [expr $height - $compCentroid]
    set bruteArea [expr $height * $width]
    set effArea [expr $depth * $width]
    set bruteInertia [expr $width * pow($height, 3)/12]
    set effInertia [expr $width * pow($depth, 3)/12]
    dict append concDict concElastmod $concElastmod
    dict append concDict concStr $concStr
    dict append concDict bruteArea $bruteArea
    dict append concDict bruteInertia [expr $bruteInertia]
    dict append concDict effArea $effArea
    dict append concDict effInertia $effInertia
    dict append concDict height $height
    dict append concDict width $width
    dict append concDict depth $depth
    dict append concDict tensionCentroid $tensionCentroid
    dict append concDict compCentroid $compCentroid
    return $concDict
}



proc ConcFardis {concDict steelStr tensionArea compArea stirrupArea} {
    # builds upon given concDict with Fardis RC params
    # todo: axial force in expression
    # steel parameters for given section, assumes uniform yield str for every steel stirrup (including shear)
    # usually fy=412 MPa for steel is used
    set steelElastmod [expr 2e6 * 9.81/100]
    set compCentroid [dict get $concDict compCentroid]
    set depth [dict get $concDict depth]
    set concElastmod [dict get $concDict concElastmod]
    set width [dict get $concDict width]
    set effArea [dict get $concDict effArea]
    set r [expr $tensionArea/$effArea]
    set rp [expr $compArea/$effArea]
    set rv [expr $stirrupArea/$effArea]
    set n [expr $steelElastmod/$concElastmod]
    # puts "find centroid $compCentroid"
    # puts "find depth $depth"
    set dp [expr $compCentroid / $depth]
    # puts "find delta prime $dp"
    set A [expr $r + $rp + $rv]
    set B [expr $r + $rp * $dp + 0.5 * $rv * (1 + $dp)]
    set ky [expr sqrt(pow($n, 2) * pow($A, 2) + 2 * $n * $B) - $n * $A]
    set phiy [lindex [lsort [list [expr $steelStr/$steelElastmod/$depth/(1-$ky)] [expr 1.8*[dict get $concDict concStr]/$depth/$concElastmod/$ky]]] 0]
    set My [expr $width*pow($depth,3)*$phiy*($concElastmod*pow($ky,2)/2*(0.5*(1+$dp) - $ky/3) + (1-$dp)*$steelElastmod/2*((1-$ky)*$r + $rp*($ky-$dp) + $rv/6*(1-$dp)))]

    dict append concDict steelStr $steelStr
    dict append concDict steelElastmod $steelElastmod
    dict append concDict tensionArea $tensionArea
    dict append concDict compArea $compArea
    dict append concDict stirrupArea $stirrupArea
    dict append concDict tensionRatio $r
    dict append concDict shearRatio $rv
    dict append concDict compRatio $rp
    dict append concDict yieldMoment $My
    dict append concDict yieldCurv $phiy

    # puts "tensionratio $r"
    # puts "compratio $rp"
    # puts "Fardis delta prime $dp"
    # puts "shearratio $rv"
    # puts "Fardis n is $n"
    # puts "Fardis A is $A "
    # puts "Fardis B is $B "
    # puts "Fardis ky is $ky"
    # puts "my yieldcurvis $phiy"
    # set yieldMoment [expr $yieldMoment/1]
    # puts "my yieldmoment is [expr $My/9810/1000] t-m"
    # puts "####################################################################### fs"
    # set tensionRatio [expr $tensionArea/[dict get $concDict effArea]]
    # set compRatio [expr $compArea/[dict get $concDict effArea]]
    # set shearRatio [expr $stirrupArea/[dict get $concDict effArea]]
    # set yieldCurv [expr 2.12 * $steelStr / $steelElastmod / [dict get $concDict height]]
    # set Fardis_n [expr $steelElastmod / [dict get $concDict concElastmod]]
    # set FardisDeltaPrime [expr [dict get $concDict cover] / [dict get $concDict depth]]
    # set FardisA [expr $tensionRatio + $compRatio + $shearRatio]
    # set FardisB [expr $tensionRatio + $compRatio * $FardisDeltaPrime + 0.5 * $shearRatio*(1 + $FardisDeltaPrime)]
    # not normalized
    # set Fardis_ky [expr ( sqrt(pow($Fardis_n,2) * pow($FardisA,2) + 2 * $Fardis_n * $FardisB) - $Fardis_n * $FardisA)]
    #minimum of next two yieldcurves
    # set yieldCurv [1.7 $steelStr / $steelElastmod/[dict get $concDict height]]
    # set yieldCurv [expr 1.8 * [dict get $concDict concStr]/([dict get $concDict concElastmod] * $Fardis_ky * [dict get $concDict depth])]
    # set yieldCurv [lindex [lsort [list [expr 1.8 * [dict get $concDict concStr]/([dict get $concDict concElastmod] * $Fardis_ky * [dict get $concDict depth])] [expr $steelStr / $steelElastmod / (1 - $Fardis_ky)/[dict get $concDict depth]] ]] 0]

    # set yieldMoment [expr pow([dict get $concDict depth], 3) * [dict get $concDict width] * $yieldCurv * (($steelElastmod * pow($Fardis_ky,2)/2*(0.5*(1 + $FardisDeltaPrime)) - $Fardis_ky/3) + $steelElastmod/2*(1 - $FardisDeltaPrime)*((1 - $Fardis_ky) * $tensionRatio + ($Fardis_ky - $FardisDeltaPrime) * $compRatio + $shearRatio/6*(1 - $FardisDeltaPrime)))]
    #
return $concDict
}

proc bar2area {number barsize} {
    # returns areas in mm2 given 1/8s in.
    # bar2area ::INt [a] -> a
    set areamm2 [expr {$number * pow($barsize * 25.4/8, 2) * 3.14159/4}]
    # puts "$number bars number $barsize are $areamm2"
    return $areamm2
}

proc addbars {listoflist} {
    # adds bar areas if more than 1 type of bar
    set total 0
    if {[llength [lindex $listoflist 0]] == 1} {
        set total [bar2area [lindex $listoflist 0] [lindex $listoflist 1]]
    } else {
    foreach elem $listoflist {
        set total [expr {$total + [bar2area [lindex $elem 0] [lindex $elem 1]]}]
    }
    }
    return $total
}
#------------------------------------------------------------------------------
# periods  = [1.447 1.391 1.070 0.541 ...]
set cols {{1000 1000} {1000 1000} {900 900} {900 900} {800 800} {800 800} {800 800} {700 700} {700 700}}
puts $outputdir "columns sizes are $cols (width height) mm."
# axial loads
#--------------------------
# i->j
#--------------------------
#COLUMN SECTION PROPERTIES
#resistances
set numstirrups 4
set stirrupsize 3
set stirrups [bar2area $numstirrups $stirrupsize]

puts $outputdir "$numstirrups stirrups of n.$stirrupsize on every element"
# creates columns
foreach col $cols {
    lappend mycols [ConcRect $fc [lindex $col 0] [lindex $col 1] $Cover $Cover]
}
#-original design
# set mybars [list [bar2area 12 12] [bar2area 10 12] [bar2area 10 8] [bar2area 8 8] [bar2area 8 8] [addbars {{2 8} {8 6}}] [addbars {{2 8} {8 6}}] [bar2area 12 5] [bar2area 12 5]]
#- design for 300tm base moment
set mybars [list [bar2area 10 10] [bar2area  3 5] [bar2area 3 5] [bar2area 3 4] [bar2area 2 4] [bar2area 3 4] [bar2area 3 4] [bar2area 3 4] [bar2area 2 4]]
# TODO: error message when below the minimum design steel
#list of bar in columns (symmetrical)
set NCRP_ht [llength $mycols]
set NCRAst {1 1 1 1 1 1 1 1 1}

set mycolfile [open _output/SPO_output/col-strs.out w]
puts $mycolfile "Column yield moments per storey (t-m) \n"
set index 1
foreach bar $mybars {
# assumes uniform yield str in whole column
    set CRT_array($index) [list $index $index $index $index]
    set Myc_in($index) [dict get [ConcFardis [lindex $mycols [expr $index - 1]] $steelStr $bar $bar $stirrups] yieldMoment]
    set Myc_ip($index) $Myc_in($index)
    set Myc_jn($index) $Myc_in($index)
    set Myc_jp($index) $Myc_in($index)
    puts $mycolfile "Col $index [format "%.1f" [expr $Myc_in($index)/9810/1000]]"
    incr index
}

unset index
# moments of inertia
set colIratio  1
set colinertia effInertia
set colarea bruteArea
puts $outputdir "--- Inertias columns ---"
puts $outputdir "column inertia factor is $colIratio"
puts $outputdir "inertia considered in analyses:  $colinertia"
puts $outputdir "area considered in analyses:  $colarea"
# TODO: general inertia and area building procedures, this used /4/ because we had 4 distinct groups
set myindexes [list 0 2 4 7]
set NCSP_ht [llength $myindexes];            #number of stiffness properties arrays per height
set NCSAst {2 2 3 2};
for {set ix 1} {$ix <= 4} {incr ix} {
    set CST_array($ix) [list $ix $ix $ix $ix]
    set Ic($ix) [expr {$colIratio*[dict get [lindex $mycols [lindex $myindexes [expr $ix - 1]]] $colinertia]}]
    puts $mycolfile "Inertia col type $ix is [expr {$Ic($ix)/1e12}] m2"
    set Ac($ix) [dict get [lindex $mycols [lindex $myindexes [expr $ix - 1]]] $colarea]
}
close $mycolfile
#-------------------------------------------------------
# this should yield the same as:
# set Ic(1) [expr 0.5*[dict get [lindex $mycols 0] bruteInertia]]
# set Ic(2) [expr 0.5*[dict get [lindex $mycols 2] bruteInertia]]
# set Ic(3) [expr 0.5*[dict get [lindex $mycols 4] bruteInertia]]
# set Ic(4) [expr 0.5*[dict get [lindex $mycols 7] bruteInertia]]
# set Ac(1) [dict get [lindex $mycols 0] bruteArea]
# set Ac(2) [dict get [lindex $mycols 2] bruteArea]
#

set mycol [Colbuilder 600 300 30 30]
set mysteel [Steelbuilder 100 10 6 0 0 0 0 0 0]





proc storey_nodal_coords {storeys storey_height bays bay_width first_storey_height} {
    set time_zero [clock clicks -millisec]
    for {set storey 0} {$storey <= $storeys} {incr storey} {
        set coords {}
        for {set bay 0} {$bay <= $bays } {incr bay} {
            if {$storey == 0} {
                lappend coords [list [expr $bay*$bay_width] 0]
            } else {
            lappend coords [list [expr $bay*$bay_width] [expr ($storey - 1) * $storey_height + $first_storey_height]]
            }
        }
        dict append Storeycoords $storey $coords
    }
    puts stderr "creating nodal coords took : [expr {([clock clicks -millisec]-$time_zero)/1000.}] (s)" ;# RS
    return $Storeycoords
}
# set Storeycoords [storey_nodal_coords 100 2.8 30 5 6] huge buildings 100+ floors

proc Nodes {storey_heights bay_widths} {
    # storeys and bays format; [list 0 3 6 ..] in meters
    # builds opensees nodes for inelastic 2d rigid-diaphragm frames
    # returns dict containing arrays with nodal tags:
    # nodes{storey0 {joint0 {center:$tag1, up:$tag2 ...} joint-n {center:$tag-n,...}} storey-n {...}}
    # tag's exact numbers are unimportant, and contained within dict keys for easy access
    # counter starts from storey0 (base) and increases upwards, joint0 corresponds to leftmost joint, increasing
    # to the right
    # storey_n node_n(center) is master node
    #           up
    # left    (center)    right
    #           down
    # up,left,right,down are nodes created for zerolength springs
    # they are all constrained to (center), and (center-i) is constrained to (center-j) of the same storey for diaphragm modeling
    set nodes [dict create]
    set counter 0
    for {set storey 0} {$storey < [llength $storey_heights]} {incr storey} {
        for {set bay 0} {$bay < [llength $bay_widths]} {incr bay} {
            set x [lindex $storey_heights $storey]
            set y [lindex $bay_widths $bay]
            set joint$bay [dict create center [incr counter]]
            node $counter $x $y
            if {$storey == 0} {
                dict append joint$bay up [incr counter]
                node $counter $x $y
                equalDOF [dict get [set joint$bay] center] [dict get [set joint$bay] up] 1 2
            } elseif {$storey == [expr [llength $storey_heights] - 1]} {;# last storey, only down,left,right, weird indexing (similar to python)
                if {$bay == 0}
                    foreach direction {right down} {
                        dict append joint$bay $direction [incr counter]
                        node $counter $x $y
                        equalDOF [dict get [set joint$bay] center] [dict get [set joint$bay] $direction] 1 2
                    }
                } elseif {$bay == [expr [llength $bay_widths] - 1]} {
                    foreach direction {left down} {
                        dict append joint$bay $direction [incr counter]
                        node $counter $x $y
                        equalDOF [dict get [set joint$bay] center] [dict get [set joint$bay] $direction] 1 2
                    }
                } else {
                    foreach direction {right left down} {
                        dict append joint$bay $direction [incr counter]
                        node $counter $x $y
                        equalDOF [dict get [set joint$bay] center] [dict get [set joint$bay] $direction] 1 2
                    }
                }
            } else {
                if {$bay == 0} {
                    # center up,down,right node only
                    foreach direction {right up down} {
                        dict append joint$bay $direction [incr counter]
                        node $counter $x $y
                        equalDOF [dict get [set joint$bay] center] [dict get [set joint$bay] $direction] 1 2
                    }
                } elseif {$bay == [expr [llength $bay_widths] - 1]} {
                    foreach direction {up left down} {
                        dict append joint$bay $direction [incr counter]
                        node $counter $x $y
                        equalDOF [dict get [set joint$bay] center] [dict get [set joint$bay] $direction] 1 2
                    }
                } else {
                    foreach direction {right up left down} {
                        dict append joint$bay $direction [incr counter]
                        node $counter $x $y
                        equalDOF [dict get [set joint$bay] center] [dict get [set joint$bay] $direction] 1 2
                    }
                }
            }
                if {$bay > 0} {
                    equalDOF [dict get $joint0 center] [dict get [set joint$bay] center] 1 ; # translational constraint on diaphragm
                }
            dict append nodes storey$storey joint$bay
        }
    }
    return $nodes
}


model BasicBuilder -ndm 2 -ndf 3
set mystoreys {0 3 6 9}
set mybays {0 6 12}
set time_zero [clock clicks -millisec]
set frame1nodes [Nodes $mystoreys $mybays]
puts [dict keys $frame1nodes]
puts [dict values [dict get $frame1nodes]]
# puts [dict values [dict get $frame1nodes storey2]]
# puts [dict get [dict get $frame1nodes storey2] joint0]
# puts [dict get [dict get $frame1nodes storey2] joint1]
puts [dict for {storey joint} $frame1nodes {puts "$storey $joint"}]
puts stderr "creating nodal coords took : [expr {([clock clicks -millisec]-$time_zero)/1000.}] (s)" ; # RS
