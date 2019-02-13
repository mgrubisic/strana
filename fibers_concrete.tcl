proc Concrete_beam {sec_tag x1 x2 cover skin coreID coverID steelID
numBarsTop barAreaTop numBarsBot barAreaBot numBarsIntTot barAreaInt
nfCoreY nfCoreZ nfCoverY nfCoverZ} {
set coverY [expr $x1/2.0];
set coverZ [expr $x2/2.0];
set coreY [expr $coverY-$cover];
set coreZ [expr $coverZ-$skin];
set numBarsInt [expr $numBarsIntTot/2];
section fiberSec $sec_tag {

    patch quadr $coreID $nfCoreZ $nfCoreY -$coreY $coreZ -$coreY -$coreZ $coreY
    -$coreZ $coreY $coreZ

    patch quadr $coverID 2 $nfCoverY -$coverY $coverZ -$coreY $coreZ $coreY $coreZ
    $coverY $coverZ
    patch quadr $coverID 2 $nfCoverY -$coreY -$coreZ -$coverY -$coverZ $coverY
    -$coverZ $coreY -$coreZ
    patch quadr $coverID $nfCoverZ 2 -$coverY $coverZ -$coverY -$coverZ -$coreY
    -$coreZ -$coreY $coreZ
    patch quadr $coverID $nfCoverZ 2 $coreY $coreZ $coreY -$coreZ $coverY -$coverZ
    $coverY $coverZ

    layer straight $steelID $numBarsInt $barAreaInt -$coreY $coreZ $coreY $coreZ;
    layer straight $steelID $numBarsInt $barAreaInt -$coreY -$coreZ $coreY -$coreZ;
    layer straight $steelID $numBarsTop $barAreaTop $coreY $coreZ $coreY -$coreZ;
    layer straight $steelID $numBarsBot $barAreaBot -$coreY $coreZ -$coreY -$coreZ;
}
}
