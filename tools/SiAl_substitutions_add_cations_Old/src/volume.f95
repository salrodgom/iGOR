real function volume(rv)
 implicit none
 real :: rv(3,3)
 real :: r1x,r1y,r1z,r2x,r2y,r2z,r3x,r3y,r3z
 r1x = rv(1,1)
 r1y = rv(2,1)
 r1z = rv(3,1)
 r2x = rv(1,2)
 r2y = rv(2,2)
 r2z = rv(3,2)
 r3x = rv(1,3)
 r3y = rv(2,3)
 r3z = rv(3,3)
 volume = abs( r1x*(r2y*r3z - r2z*r3y) + r1y*(r3x*r2z - r3z*r2x) + r1z*(r2x*r3y - r2y*r3x) )
 return
end function volume
!
real function area(rv)
 implicit none
 real :: rv(3,3)
 real :: r1x,r1y,r2x,r2y
!
 r1x = rv(1,1)
 r1y = rv(2,1)
 r2x = rv(1,2)
 r2y = rv(2,2)
 area = abs( r1x*r2y - r1y*r2x )
 return
end function area
