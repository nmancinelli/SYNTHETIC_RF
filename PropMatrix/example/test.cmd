#!/bin/csh
set xdir=../bin
#ANGL is the incidence angle from vertical
set ANGL=30
#AZM is the clockwise horizontal azimuth from the 1-axis of the elastic
# coefficients
foreach AZM ( 45 )
#
$xdir/synth <<!
aniso.data
$ANGL
1
0.
$AZM
!
#
echo "Done propagation, now convolving source"
#
$xdir/sourc1 <<!
test.$AZM.$ANGL
1 0 0
5.
0.
2
!
#
echo "Done source convolution, now writing to sac"
./saccpt <<!
test.$AZM.$ANGL
0
!

end
echo OK
