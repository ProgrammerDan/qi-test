Computational analysis of the theory of Rindler Horizon hiding mass by Michael McCulloch based on a rotating object.
Proposed here: http://www.ptep-online.com/2019/PP-57-05.PDF

High level summary:
1) Acceleration can reduce the physical distance at which information can reach an object.
2) This "Rindler Horizon" gets closer to an object under acceleration, including rotation related acceleration
3) This implies that even centripetal accelerations of rotating objects will produce Rindler Horizons in the 
   direction opposite of this acceleration.
4) Given these assumptions, an object rotating very, very fast (in the 30k to 400k RPM range) will see
   assymetric hiding of stellar masses.
5) It is possible to model these effects, by assuming an object under uniform rotational acceleration, computing the
   resulting Rindler horizon effects across the object, and determining the reduction in stellar mass that will
   impact the rotation object.
6) The results are in effect the analytic estimation of gravitation impact reduction; if sufficiently large, it should be
   experimentally confirmable.
   
This specific simulation uses estimations, and will be less accurate then a purely mathematic solution. Such a solution is no longer
accessible to me, but the math of a computation model is accessible, so I've done that instead.
Specifically, we do the following:
1) Assume a rotating sphere of uniform density, on the surface of the Earth.
2) Subdivide the cuboid region around the rotating sphere into cuboids of uniform volume and density.
3) Based on each cuboid's centroid distance from the axis of rotation, we can determine its acceleration and corresponding Rindler horizon.
4) Using basic vector math and systems of equations for planes and spheres, we can determine what if anything is predicted to be hidden
   by the Rindler horizon and the plane normal to the acceleration vector.
5) All stellar mass (even fractional) "above" the plane (in the direction of acceleration) are not hidden. Portions of stellar mass below the
   acceleration plane and outside the Rindler horizon are hidden, but only the portions actually hidden.

6) Ultimately I'd like to avoid all point masses! Density times Volume computations where the volumes are computed via sphere-sphere, or sphere-plane, 
   or both are used to determine remaining fractional masses that will interact.
   This is an important point! In the paper above M. E. McCulloch estimates the various rates of rotation that will hide stellar objects but
   makes an important over-simplification, namely, the solution is for a point mass, not a volumetric mass. So, he estimates that the earth will be hidden
   at 3589000 rpm; in fact, this will hide ~50% of the earth's mass from some fraction of the sphere, as the Rindler horizon will be just 
   inside the line formed from the center of the object to the center of the earth. Obviously the earth is not a point mass; while such a small
   error in math can be forgiven it could lead to incorrect predictions about the usefulness of "horizon drives" based on rotation. Unfortunately I am
   fighting an uphill battle to regain rusty math skill so instead of using Gauss's law for gravity which would allow me to properly derive
   directly gravitation impact from the shapes I've derived, I'm resorting internally to "slightly more accurate" point mass / center of mass computations.
   These being inaccurate, however, there are some odd effects on the outcomes as a result.

Errors may exist. I've attempted to account for every kind of intersection possible under the contraints of solutions between sphere-sphere, 
sphere-plane, etc.

June 20, 2020

By ProgrammerDan <programmerdan@gmail.com>
 
Many thanks to the author of Apfloat; this high-precision library is VERY fast and amazing, doing the same with BigDecimal would have been
nearly impossible but was insanely easy with Apfloat. Great work.
