package com.programmerdan.physics.qi;

import java.util.concurrent.atomic.AtomicLong;

import org.apfloat.*;

/**
 * Computational analysis of the theory of Rindler Horizon hiding mass by Michael McCulloch based on a rotating object.
 * Proposed here: http://www.ptep-online.com/2019/PP-57-05.PDF
 * 
 * High level summary:
 * 1) Acceleration can reduce the physical distance at which information can reach an object.
 * 2) This "Rindler Horizon" gets closer to an object under acceleration, including rotation related acceleration
 * 3) This implies that even centripetal accelerations of rotating objects will produce Rindler Horizons in the 
 *    direction opposite of this acceleration.
 * 4) Given these assumptions, an object rotating very, very fast (in the 30k to 400k RPM range) will see
 *    assymetric hiding of stellar masses.
 * 5) It is possible to model these effects, by assuming an object under uniform rotational acceleration, computing the
 *    resulting Rindler horizon effects across the object, and determining the reduction in stellar mass that will
 *    impact the rotation object.
 * 6) The results are in effect the analytic estimation of gravitation impact reduction; if sufficiently large, it should be
 *    experimentally confirmable.
 *    
 *    TODO: Consider Lorentz contraction at high edge velocity. Not sure how to do these computations in Non-Euclidean space, though.
 *     But, this may be something I need to consider anyway, as I'd like to imagine all this in the context of Relativity, anyway.
 *    
 * This specific simulation uses estimations, and will be less accurate then a purely mathematic solution. Such a solution is no longer
 * accessible to me, but the math of a computation model is accessible, so I've done that instead.
 * Specifically, we do the following:
 * 1) Assume a rotating sphere of uniform density, on the surface of the Earth.
 * 2) Subdivide the cuboid region around the rotating sphere into cuboids of uniform volume and density.
 * 3) Based on each cuboid's centroid distance from the axis of rotation, we can determine its acceleration and corresponding Rindler horizon.
 * 4) Using basic vector math and systems of equations for planes and spheres, we can determine what if anything is predicted to be hidden
 *    by the Rindler horizon and the plane normal to the acceleration vector.
 * 5) All stellar mass (even fractional) "above" the plane (in the direction of acceleration) are not hidden. Portions of stellar mass below the
 *    acceleration plane and outside the Rindler horizon are hidden, but only the portions actually hidden.
 * 
 * 6) No point masses are used! Density*Volume computations where the volumes are computed via sphere-sphere, or sphere-plane, or both are used
 *    to determine remaining fractional masses that will interact.
 *    This is an important point! In the paper above M. E. McCulloch estimates the various rates of rotation that will hide stellar objects but
 *    makes an important mistake, namely, the solution is for a point mass, not a volumetric mass. So, he estimates that the earth will be hidden
 *    at 3589000 rpm; in fact, this will hide ~50% of the earth's mass from some fraction of the sphere, as the Rindler horizon will be just 
 *    inside the line formed from the center of the object to the center of the earth. Obviously the earth is not a point mass; while such a small
 *    error in math can be forgiven it could lead to incorrect predictions about the usefulness of "horizon drives" based on rotation.
 * 
 * Errors may exist. I've attempted to account for every kind of intersection possible under the contraints of solutions between sphere-sphere, 
 * sphere-plane, etc.
 * 
 * June 20, 2020
 * 
 * @author ProgrammerDan <programmerdan@gmail.com>
 * 
 * Many thanks to the author of Apfloat; this high-precision library is VERY fast and amazing, doing the same with BigDecimal would have been
 * nearly impossible but was insanely easy with Apfloat. Great work.
 *
 */
public class RindlerHorizonHidingViaAngularMotion implements Runnable {

	/**
	 * This is the radius of the rotating object. Larger values will increase the relative acceleration of points along the edge of the
	 * sphere. I do not look for relativistic effects, e.g. you can construct a rate of rotation and radius that will move faster then light
	 * and I will not stop you.
	 */
	private Apfloat objectRadius =       new Apfloat(".0002", PRECISION); // .0002km or .2 m (20 cm)
	//private Apfloat objectRotationRate = new Apfloat(30000d/60d, PRECISION); // rotations per second -- according to MM, 23krpm should "hide" some of sun from rotating object. Bumping to 30k. 
	//private Apfloat objectRotationRate = new Apfloat(30d/60d, PRECISION); // terrestrial rotations of 30rpm -- should result in no stellar hiding.
	private Apfloat objectRotationRate = new Apfloat(3589000d/60d, PRECISION); // He proposes 3589krpm as enough to "hide" the earth

	private Apfloat objectRotationRateSq = objectRotationRate.multiply(objectRotationRate);
	private Apfloat precomputeRotationFactor = PISQ.multiply(objectRotationRateSq).multiply(new Apfloat(4d, PRECISION));
	
	/**
	 * Change the objectResolution as suits you. Smaller values will take longer to simulate, as the cuboid space occupied by the sphere under
	 * acceleration will be subdivided into more cubes; however the resulting output will be more precise.
	 */
	//private Apfloat objectResolution =   new Apfloat(".000001", PRECISION); // .000001km or .001 m (.1cm, 1mm)
	//private Apfloat objectResolution =   new Apfloat(".00001", PRECISION); // .00001km or .01 m (1cm, 10mm)
	private Apfloat objectResolution =   new Apfloat(".00002", PRECISION); // .00002km or .02 m (2cm, 20mm)

	// At present, object Density has no bearing, as we are determining the impact of stellar masses on us, not the mutual gravitation.
	//private Apfloat objectDensity =      new Apfloat("1800000000000", PRECISION); // carbon fiber kg/km3

	/**
	 * Don't modify this unless you intend to make this a moving simulation, e.g. with system effects. Right now
	 * it's designed to simulate a rotating sphere sitting on the average radius of the earth's surface, while the earth,
	 * moon, and sun are in a line (to maximize the rindler horizon impact).
	 * Other arrangements are possible. Note that this formulation also assumes that the object's rotation is orthoganal to the plane
	 * on which the object and stellar masses are aligned.
	 */
	private Location objectLocation =    new Location(Apfloat.ZERO, Apfloat.ZERO, Apfloat.ZERO);
	
	private static long PRECISION = 1024l;
	private static final FixedPrecisionApfloatHelper MATH = new FixedPrecisionApfloatHelper(PRECISION);
	private static final Apfloat PI = ApfloatMath.pi(PRECISION);
	private static final Apfloat PISQ = PI.multiply(PI);
	private static Apfloat SPEEDOFLIGHT = new Apfloat("299792.458", PRECISION); //km/s
	private static Apfloat SPEEDOFLIGHTSQ = SPEEDOFLIGHT.multiply(SPEEDOFLIGHT);
	private static final Apfloat TWO = new Apfloat("2", PRECISION);
	private static final Apfloat THREE = new Apfloat("3", PRECISION);
	private static final Apfloat FOUR = new Apfloat("4", PRECISION);
	private static final Apfloat SIX = new Apfloat("6", PRECISION);
	private static final Apfloat TWELVE = new Apfloat("12", PRECISION);
	//private static final Apfloat HALF = new Apfloat(".5", PRECISION);
	private static final Apfloat THREEEIGHTS = new Apfloat(".375", PRECISION);
	// Note to self: Apfloat.ZERO is always infinite precision. be careful
	
	private Apfloat sunDensity =         new Apfloat("1408000000000", PRECISION); // kg/km3
	private Apfloat sunRadius =          new Apfloat("695700", PRECISION); // km -- Volumetric mean radius
	private Apfloat sunMass =            FOUR.divide(THREE).multiply(PI).multiply(MATH.pow(sunRadius, THREE)).multiply(sunDensity);
	private Location sunLocation =       new Location(Apfloat.ZERO, Apfloat.ZERO, Apfloat.ZERO);
	
	private Apfloat moonDensity =        new Apfloat("3344000000000", PRECISION); // kg/km3
	private Apfloat moonRadius =         new Apfloat("1737.4", PRECISION); // km -- Volumetric mean radius.
	private Apfloat moonMass =           FOUR.divide(THREE).multiply(PI).multiply(MATH.pow(moonRadius, THREE)).multiply(moonDensity);
	//private Apfloat moonSemimajorAxis =  new Apfloat("384400", PRECISION); // km -- average semimajor axis
	private Apfloat moonPerigee =        new Apfloat("363300", PRECISION); // km -- average perigee axis
	private Apfloat moonApogee =         new Apfloat("405500", PRECISION); // km -- average apogee axis
	//private Apfloat moonInclination =    new Apfloat((5.145 / 360d) * 2d * Math.PI, PRECISION); // in radians, inclination to ecliptic.
	private Location moonLocation =      new Location(Apfloat.ZERO, Apfloat.ZERO, Apfloat.ZERO);
	
	private Apfloat earthDensity =       new Apfloat("5514000000000", PRECISION); // kg/km3
	private Apfloat earthRadius =        new Apfloat("6371", PRECISION); // km -- Volumetric mean radius.
	private Apfloat earthMass =          FOUR.divide(THREE).multiply(PI).multiply(MATH.pow(earthRadius, THREE)).multiply(earthDensity);
	//private Apfloat earthSemimajorAxis = new Apfloat("149600000", PRECISION); //km
	private Apfloat earthPerihelion =    new Apfloat("147090000", PRECISION); // km
	private Apfloat earthAphelion =      new Apfloat("152100000", PRECISION); // km
	//private Apfloat earthTilt =          new Apfloat((23.44 / 360d) * 2d * Math.PI, PRECISION); // equatorial inclination.
	private Location earthLocation =     new Location(Apfloat.ZERO, Apfloat.ZERO, Apfloat.ZERO);
	
	private Apfloat constantG =          new Apfloat("6.67430e-20", PRECISION); // km^3*kg^-1*s^-2
	
	public RindlerHorizonHidingViaAngularMotion(Apfloat apRadiusInKM, Apfloat apRotationRate, Apfloat apResolution,
			Apfloat apLightSpeed, Apfloat apSunRadius, Apfloat apSunDensity, Apfloat apMoonRadius,
			Apfloat apMoonDensity, Apfloat apMoonOrbitalPosition, Apfloat apMoonOffsetAngle, Apfloat apEarthRadius,
			Apfloat apEarthDensity, Apfloat apEarthOrbitalPosition, Apfloat apEarthOffsetAngle, Apfloat apConstantsG, long PRECISION) {

		RindlerHorizonHidingViaAngularMotion.PRECISION = PRECISION;
		
		if (apRadiusInKM != null) {
			objectRadius = apRadiusInKM;
		}
		if (apResolution != null) {
			objectResolution = objectRadius.divide(apResolution);
		}
		if (apRotationRate != null) {
			objectRotationRate = apRotationRate.divide(new Apfloat(60, PRECISION));
			objectRotationRateSq = objectRotationRate.multiply(objectRotationRate);
			precomputeRotationFactor = PISQ.multiply(objectRotationRateSq).multiply(new Apfloat(4d, PRECISION)); 
		}

		if (apLightSpeed != null) {
			RindlerHorizonHidingViaAngularMotion.SPEEDOFLIGHT = apLightSpeed;
			RindlerHorizonHidingViaAngularMotion.SPEEDOFLIGHTSQ = apLightSpeed.multiply(apLightSpeed);
		}
		
		if (apSunRadius != null && apSunDensity != null) {
			sunDensity = apSunDensity;
			sunRadius = apSunRadius;
			sunMass = FOUR.divide(THREE).multiply(PI).multiply(MATH.pow(sunRadius, THREE)).multiply(sunDensity);
		}
		
		if (apMoonRadius != null && apMoonDensity != null) {
			moonDensity = apMoonDensity;
			moonRadius = apMoonRadius;
			moonMass = FOUR.divide(THREE).multiply(PI).multiply(MATH.pow(moonRadius, THREE)).multiply(moonDensity);
		}
		
		if (apEarthRadius != null && apEarthDensity != null) {
			earthDensity = apEarthDensity;
			earthRadius = apEarthRadius;
			earthMass = FOUR.divide(THREE).multiply(PI).multiply(MATH.pow(earthRadius, THREE)).multiply(earthDensity);
		}
		
		if (apConstantsG != null) {
			constantG = apConstantsG;
		}
		
		if (apMoonOffsetAngle != null || apEarthOffsetAngle != null) {
			System.err.println("Warning: moon and earth offset angles are current ignored, ordering is fixed at obj on surface of earth, moon inbetween earth and sun.");
		}
		Apfloat moonDistance = moonApogee;
		Apfloat earthDistance = earthAphelion;
		if (apMoonOrbitalPosition != null) {
			moonDistance = moonPerigee.add( moonApogee.subtract(moonPerigee).multiply(apMoonOrbitalPosition) );
		}
		if (apEarthOrbitalPosition != null) {
			earthDistance = earthPerihelion.add( earthAphelion.subtract(earthPerihelion).multiply(apEarthOrbitalPosition) );
		}
		
		earthLocation = objectLocation.offset(earthRadius.add(objectRadius), Apfloat.ZERO, Apfloat.ZERO);
		moonLocation = earthLocation.offset(moonDistance, Apfloat.ZERO, Apfloat.ZERO);
		sunLocation = earthLocation.offset(earthDistance, Apfloat.ZERO, Apfloat.ZERO);

	}
	
	Location[] gravAcceleration = new Location[] { new Location(Apfloat.ZERO, Apfloat.ZERO, Apfloat.ZERO), new Location(Apfloat.ZERO, Apfloat.ZERO, Apfloat.ZERO), new Location(Apfloat.ZERO, Apfloat.ZERO, Apfloat.ZERO) };
	Location aggregateGravityAcceleration = new Location(Apfloat.ZERO, Apfloat.ZERO, Apfloat.ZERO);
	AtomicLong[] pathsUsed = new AtomicLong[17];
	String[] pathsNamed = new String[] { "Fully Inside Horizon (inc)", "Fully Above Accel Plane (inc)", "Fully Below Accel Plane (exc)", 
										"Half Above Accel Plane (bth)", "Mostly Above Accel Plane (bth)", "Mostly Below Accel Plan (bth)",
										"Lens Intersect Horizon (bth)", "Excl Lens above Accel Plane (exc)", "Incl Lens below Accel Plan (Inc)", 
										"Incl Half Lens (Inc)", "Incl Minority Lens below Accel (bth)", "Incl Majority Lens below Accel (bth)", 
										"Also Half Cap (bth)", "Also Majority Cap (bth)", "Also Minority Cap (bth)", "Only Cap (Inc)", "Cap and Lens (Inc)"};
	AtomicLong currentTakenSteps = new AtomicLong(0l);
	AtomicLong currentDoneSteps = new AtomicLong(0l);
	long stepsOnEdge = 0l;
	long totalSteps = 0l;

	String[] sphereName = new String[]{"Earth", "Moon", "Sun"};
	Location[] sphereCenter = new Location[]{ earthLocation, moonLocation, sunLocation };
	Apfloat[] sphereRadius = new Apfloat[]{ earthRadius, moonRadius, sunRadius };
	Apfloat[] sphereDensity = new Apfloat[]{ earthDensity, moonDensity, sunDensity };
	Apfloat[] sphereMass = new Apfloat[]{ earthMass, moonMass, sunMass };
	
	Apfloat[][] aggregateQIAcceleration = null;
	String[][] aggregateQIPaths = null;

	Location newtonAggregateGravityAcceleration = new Location(Apfloat.ZERO, Apfloat.ZERO, Apfloat.ZERO);
	Apfloat newtonAggregateGravityAccelerationScale = null;
	Location[] newtonGravAcceleration = new Location[] { new Location(Apfloat.ZERO, Apfloat.ZERO, Apfloat.ZERO), new Location(Apfloat.ZERO, Apfloat.ZERO, Apfloat.ZERO), new Location(Apfloat.ZERO, Apfloat.ZERO, Apfloat.ZERO) };
	Apfloat[] newtonGravAccelerationScale = new Apfloat[3];
	
	public void run() {
		for (int i = 0; i < 17; i++) {
			pathsUsed[i] = new AtomicLong(0l);
		}
		/**
		 * Let us consider for a moment a pure math or logical solution.
		 * 
		 * You have a spherical mass "object" spinning along a single axis. Since we assume the spin is uniform subject to a=vr
		 * Consequently, for any particular spherical object "sphere" you'd like to consider, you can ask three questions:
		 * 1) Is the "sphere" inside the Rindler horizon? (fully visible). If it is, you can include its entire mass as a gravitational impact on the object
		 * 2) Is the "sphere" outside the Rindler horizon? If so, some portion of the "object" will see its mass and some portion will not. Specifically,
		 *    the arc formed by rays from the center of the "object" to the edges of the "sphere" form the horizon corridor, this is the area where
		 *    some portion of the "sphere" will be hidden, mass-wise. Given the way the horizon will "sweep" across the arc, 
		 * 
		 */
		/*
		 * Let's arrange our experiment in a number of ways. Basically, we need to position the sun and moon
		 * and our test object.
		 * 
		 * For ease of computation, we will determine locations for earth and moon and test object, then
		 * alter the reference frame to be centered on the rotating test object.
		 * 
		 * Intent here is to determine a reasonable expectation for a rotating uniform sphere in reference to 
		 * QI's proposed Rindler Horizon formation.
		 * 
		 * Since the sphere is uniform, it makes sense that we can compute the acceleration profile of the sphere
		 * once per configuration of rotational speed and sun/earth/moon system, and form an analytic expectation
		 * of resulting non-uniform acceleration of the test sphere, if any. Once we have this computed profile, 
		 * if the math checks out it should be possible to generate, experimentally, tests of this effect.
		 * 
		 * Also, to simplify math initially, our first analysis will assume the object/earth/moon/sun are perfectly aligned with 
		 * the moon at apogee from earth, and earth at Apogee from sun. We can "run" an orbital simulation from there and 
		 * retest impact of local objects based on QI Predictions.
		 * 
		 * We assume solar system orbital plane is X, Y. "up/down" relative to the plane of the orbital plane is Z.
		 */
		
		// we're set up, let's do analytics.
		stepsOnEdge = objectRadius.divide(objectResolution).longValue();
		totalSteps = stepsOnEdge*2*stepsOnEdge*2*stepsOnEdge*2;
		
		//Apfloat resolutionMass = ApfloatMath.pow(objectRadius,3l).multiply(objectDensity); // results in kg, but we don't need it right now.
		final Apfloat objectRadiusSquared = objectRadius.multiply(objectRadius);
		
		sphereName = new String[]{"Earth", "Moon", "Sun"};
		sphereCenter = new Location[]{ earthLocation, moonLocation, sunLocation };
		sphereRadius = new Apfloat[]{ earthRadius, moonRadius, sunRadius };
		sphereDensity = new Apfloat[]{ earthDensity, moonDensity, sunDensity };
		sphereMass = new Apfloat[]{ earthMass, moonMass, sunMass };

		if (  ((long) ((int)totalSteps)) == totalSteps) {
			aggregateQIAcceleration = new Apfloat[(int)totalSteps][sphereName.length];
			aggregateQIPaths = new String[(int)totalSteps][sphereName.length];
		} else {
			System.out.println("Warning: Output space is too dense to allow point tracking.");
		}
		
		// Precompute Newtonian for visualization purposes.
		for (int i = 0; i < sphereName.length; i++) {
			Vector accelLine = Vector.betweenPoints(objectLocation, sphereCenter[i]);
			Apfloat a = constantG.multiply(sphereMass[i]).divide(accelLine.magnitude.multiply(accelLine.magnitude));
			Location contribution = new Location(
					accelLine.unitVector.x.multiply(a),
					accelLine.unitVector.y.multiply(a),
					accelLine.unitVector.z.multiply(a));
			Apfloat gx = newtonAggregateGravityAcceleration.x.add(contribution.x);
			Apfloat gy = newtonAggregateGravityAcceleration.y.add(contribution.y);
			Apfloat gz = newtonAggregateGravityAcceleration.z.add(contribution.z);
			
			newtonAggregateGravityAcceleration = new Location(gx, gy, gz);
			
			gx = newtonGravAcceleration[i].x.add(contribution.x);
			gy = newtonGravAcceleration[i].y.add(contribution.y);
			gz = newtonGravAcceleration[i].z.add(contribution.z);
			
			newtonGravAcceleration[i] = new Location(gx, gy, gz);
			newtonGravAccelerationScale[i] = a;
		}
		
		newtonAggregateGravityAccelerationScale = MATH.sqrt(
				newtonAggregateGravityAcceleration.x.multiply(newtonAggregateGravityAcceleration.x).add(
				newtonAggregateGravityAcceleration.y.multiply(newtonAggregateGravityAcceleration.y)).add(
				newtonAggregateGravityAcceleration.z.multiply(newtonAggregateGravityAcceleration.z))
			);
		
		System.out.println("Object is at " + objectLocation);
		System.out.println("Earth is at " + earthLocation);
		System.out.println("Moon is at " + moonLocation);
		System.out.println("Sun is at " + sunLocation);
		final Apfloat tmpRindlerHorizonDistance = SPEEDOFLIGHTSQ.divide(objectRadius.multiply(precomputeRotationFactor));
		System.out.println(String.format("Tightest Rindler Horizon will be at %#s from Object", tmpRindlerHorizonDistance.precision(20l)));
		
		System.out.println("Prepping simulation, using reference object with " + stepsOnEdge + " steps per edge, resulting in " + totalSteps + " to compute across " + sphereName.length + " grav objects.");
		
		for (double stepX = -stepsOnEdge; stepX < stepsOnEdge; stepX++) {
			Apfloat midX = new Apfloat(stepX + 0.5d, PRECISION);
			Apfloat pointX = midX.multiply(objectResolution);
			for (double stepY = -stepsOnEdge; stepY < stepsOnEdge; stepY++) {
				Apfloat midY = new Apfloat(stepY + 0.5d, PRECISION);
				Apfloat pointY = midY.multiply(objectResolution);
				for (double stepZ = -stepsOnEdge; stepZ < stepsOnEdge; stepZ++) {
					Apfloat midZ = new Apfloat(stepZ + 0.5d, PRECISION);
					Apfloat pointZ = midZ.multiply(objectResolution);
					/* we ignore anything where the midpoint of the resolution block falls outside the sphere.
					 * we imagine each step as such:
					 *    ______
					 *   / _   /
					 *  /____ /|
					 * |     |.|
					 * |  x  | |
					 * |_____|/
					 * 
					 * where the face marks represent the "center" of a small block, with each edge of the box representing the
					 * resolution. The mass of this microbox is edge^3*objectdensity (gives in kg).
					 * 
					 * We compute the box based on
					 *  (stepX,stepY,stepZ)*objectResolution -> (stepX+1,stepY+1,stepZ+1)*objectResolution)
					 * with midpoint:
					 *  (stepX+.5, stepY+.5,stepZ+.5)*objectResolution
					 */
					
					final Location midPoint = new Location(pointX, pointY, pointZ);
					final int xX = (int) stepX;
					final int yY = (int) stepY;
					final int zZ = (int) stepZ;
					if (midPoint.isInsideSphere(objectLocation, objectRadiusSquared)) {
						QIPlatform.processHandler.execute(new Runnable() {
							public void run() {
								StringBuilder paths = new StringBuilder();
								currentTakenSteps.incrementAndGet();
								// it's inside!
								/*
								 * Now we compute the 3D acceleration vector for this box based on rotation.
								 * We assume uniform, stable acceleration of all points.
								 * We assume uniform density.
								 * magnitude of centripetal (center-pointing) acceleration is determined by a=4pi^2rR^2 where r = radius and R = rot/sec
								 * Note that this centripetal vector is based on circular motion around the Z axis of rotation, which we assume for simplicity of math.
								 */
								Vector centerOfRotationPointingVector = Vector.pointDirection(midPoint, new Location(midPoint.x.negate(), midPoint.y.negate(), Apfloat.ZERO) );
								Location centerOfRotationPointingUnitVector = centerOfRotationPointingVector.unitVector;
								Apfloat centerOfRotationPointingRadius = centerOfRotationPointingVector.magnitude;
								Apfloat magnitudeOfCentripetalAcceleration = centerOfRotationPointingRadius.multiply(precomputeRotationFactor); // will be km/s^2
								
								/*
								 * Once we have this acceleration vector, determine the Rindler horizon distance based on acceleration magnitude.
								 * TODO: how do we interpret Rindler horizon for the vector at large? Can we do the math using just matrix math?
								 * 
								 * FOR NOW: we will assume the rindler horizon is "sharp" -- effectively no horizon in direction of accel, and a half-sphere
								 * "behind" our acceleration vector with magnitude at all points using MM's equation: dR = c^2 / a; c is the speed of light,
								 * a is our magnitude of acceleration from above.
								 * 
								 * Note all math units are in km or kg.
								 * 
								 * So, this gives us a point ("midPoint"), and a unit vector (centerOfRotationPointingUnitVector), from that we can determine
								 * a plane. 
								 * 
								 */
								Apfloat rindlerHorizonDistance;
								if (!magnitudeOfCentripetalAcceleration.equals(Apfloat.ZERO)) {
									rindlerHorizonDistance = SPEEDOFLIGHTSQ.divide(magnitudeOfCentripetalAcceleration);
								} else {
									// rindlerHorizon is "unlimited" in this formulation. Although, this would violate QI's minimum speedlimit postulate
									// so perhaps that is the answer? compute using minimum plank length based acceleration. TODO
									rindlerHorizonDistance = earthAphelion.multiply(earthAphelion); // unit in km
								}
								
								/*
								 * Equation of a sphere:
								 * (x-xs)^2+(y-ys)^2+(z-zs)^2-R^2=0
								 * xs, ys, zs is the center location of the object considered.
								 * R is the radius.
								 * 
								 * Equation of a plane:
								 * D=-A*x0-B*y0-C*z0
								 * A, B, C are the unit vector elements normal to the plane.
								 * D is derived as above.
								 * 
								 * 
								 * all distance units in km
								 */
								Apfloat plane_A = centerOfRotationPointingUnitVector.x;
								Apfloat plane_B = centerOfRotationPointingUnitVector.y;
								Apfloat plane_C = centerOfRotationPointingUnitVector.z;
								Apfloat plane_D = midPoint.x.negate().multiply(plane_A)
												.subtract(midPoint.y.multiply(plane_B))
												.subtract(midPoint.z.multiply(plane_C));
								
								Apfloat AsqBsqCsq = plane_A.multiply(plane_A).add(
													plane_B.multiply(plane_B)).add(
													plane_C.multiply(plane_C));
								Apfloat AsqBsqCsqRoot = MATH.sqrt(AsqBsqCsq);
								
								/*
								 * Now per sphere we're considering:
								 * square the radius; 
								 * compute A*xs + B*ys + C*zs + D
								 * 
								 * d = (A*xs + B*ys + C*zs + D) / AsqBsqCsqRoot
								 * if R > |d| then, intersection!
								 * if R = |d| then, tangency (no intersection for us)
								 * if R < |d| then, no intersection at all.
								 * 
								 * d is distance from intersection circle to sphere center (h in other formulas for SphericalCap)
								 * 
								 */
								if (!AsqBsqCsqRoot.equals(Apfloat.ZERO)) { 
									for (int i = 0; i < sphereName.length; i++) {
										paths = new StringBuilder();										
										/*
										 * Consider based on the plane and unit vector, using the point-sphere intersection point, which might be outside
										 * the sphere, is this sphere on the "acceleration" side of the plane? or, the opposite side.
										 * 
										 * if the sphere has no intersection, but is on the acceleration side of the plane, that means its mass applies fully.
										 * How to determine if on the acceleration side of the plane? the unit vector of the plane from sphere center TO plane is
										 * inverse of unit vector of acceleration.
										 * Otherwise; it's on the "rear" and vulnerable to "hiding" by QI. E.g. we continue evaluating this.
										 * 
										 * Note: from https://en.wikipedia.org/wiki/Center_of_mass
										 * R = (1 / (m1+m2)) * (m1*r1+m2*r2) where R is vector point; r1 is vector center of mass 1, r2 is vector center of mass 2. 
										 */
										
										
										/*
										 * Then, ideally using parametric equations if possible, determine what mass fraction impacts this microobject
										 * based on midpoint of microobject & microobject mass and bariometric midpoint of "rindler visible" portion of volume / mass
										 * of each considered object (moon, earth, sun).
										 * We assume all considered objects are perfect spheres for determining visibility; 
										 * equation should be; 
										 * find the plane determine by the acceleration unit vector; 
										 * determine intersection if any of plane and sphere;
										 * determine volume of spherical cap "below" the plane subject to intersection.
										 * determine mass of spherical cap (using density & volume)
										 * should be able to determine midpoint of intersection in sphere; 
										 * generate vector from "true" centroid of sphere through midpoint of intersection.
										 * 
										 * Ok, https://mathworld.wolfram.com/SphericalCap.html has done all the derivations, huzzah.!
										 * So we just need to compute, h, where h is the vector magnitude of the spherical cap's "size".
										 * Then we can compute the distance along the sphere's centroid vector through intersection midpoint that is geometric centroid.
										 * 
										 */
										
										// Here's the outline:
										/*
										 * Direction of acceleration is considered "above plane". Direction against acceleration and subject to
										 * hiding is "below plane". Hiding begins at the Rindler Horizon.
										 * 
										 * High level:
										 * 1) if sphere is fully inside rindler horizon, include mass.
										 * 2) if sphere is fully outside rindler horizon, look for sphere-plane intersection.
										 *   a) If sphere is fully above plane, include mass.
										 *   b) If sphere is fully below plane, do not include mass.
										 *   c) Compute the portion of the sphere that is above the plane, and include that fraction of mass.
										 * 3) If sphere is partially within rindler horizon:
										 *   a) if sphere is fully above plane, include mass.
										 *   b) If sphere is fully below plane, include fraction of mass based on sphere-sphere intersection
										 *   c) Compute sphere-sphere intersection, then compute the portion of the lens that is above the plane, and include that fraction of mass.
										 *
										 * 1) if sphere is fully inside rindler horizon, include mass.
										 * 2) If sphere is fully above plane, include mass.
										 * 3) if sphere is fully outside rindler horizon, look for sphere-plane intersection.
										 *   a) If sphere is fully below plane, do not include mass.
										 *   b) Compute the portion of the sphere that is above the plane, and include that fraction of mass.
										 * 4) If sphere is partially within rindler horizon:
										 *   a) If sphere is fully below plane, include fraction of mass based on sphere-sphere intersection
										 *   b) Compute sphere-sphere intersection, then compute the portion of the lens that is above the plane, and include that fraction of mass.
										 */
										Apfloat sphereRSq = sphereRadius[i].multiply(sphereRadius[i]);
										Location sphereS = sphereCenter[i];
										
										Apfloat impactMass = null;
										Location impactCenter = null;
										
										Vector fromObjectToSphere = Vector.betweenPoints(midPoint, sphereS);
										Apfloat sphereDistance = fromObjectToSphere.magnitude;
										
										Apfloat sphereInnerEdge = sphereDistance.subtract(sphereRadius[i]);
										Apfloat sphereOuterEdge = sphereDistance.add(sphereRadius[i]);
										if (rindlerHorizonDistance.compareTo(sphereOuterEdge) >= 0) {
											pathsUsed[0].incrementAndGet(); paths.append(",0");
											/* 1) if sphere is fully inside rindler horizon, include mass. */
											impactMass = sphereMass[i];
											impactCenter = sphereS;
										} else {
											Apfloat dNumerator = plane_A.multiply(sphereS.x).add(
													plane_B.multiply(sphereS.y)).add(
													plane_C.multiply(sphereS.z)).add(
													plane_D);
											Apfloat compute_d = dNumerator.divide(AsqBsqCsqRoot);
											Apfloat absCompute_d = MATH.abs(compute_d);
											
											/* note here that 2) through 4) all leverage notion of "above" and "below" plane, so we compute that. */
											Apfloat intersectX = sphereS.x.subtract(plane_A.multiply(dNumerator).divide(AsqBsqCsq));
											Apfloat intersectY = sphereS.y.subtract(plane_B.multiply(dNumerator).divide(AsqBsqCsq));
											Apfloat intersectZ = sphereS.z.subtract(plane_C.multiply(dNumerator).divide(AsqBsqCsq));
											Location intersectCenter = new Location(intersectX, intersectY, intersectZ);
											
											Vector fromPlaneToSphere = Vector.betweenPoints(intersectCenter, sphereS);
											Location fromPlaneToSphereUnitVector = fromPlaneToSphere.unitVector;
											
											if (sphereRadius[i].compareTo(absCompute_d) <= 0 && centerOfRotationPointingUnitVector.isEqual(fromPlaneToSphereUnitVector, PRECISION / 2)) {
												pathsUsed[1].incrementAndGet(); paths.append(",1");
												/*2) If sphere is fully above plane, include mass.*/ 
												/* Or, sphere does not intersect plane and is above the plane -- the pointing vector of plane and vector from plane to sphere center are in same direction. */
												impactMass = sphereMass[i];
												impactCenter = sphereS;
											} else if (rindlerHorizonDistance.compareTo(sphereInnerEdge) <= 0) {
												/* 3) if sphere is fully outside rindler horizon, look for sphere-plane intersection. */
												if (sphereRadius[i].compareTo(absCompute_d) <= 0 && !centerOfRotationPointingUnitVector.isEqual(fromPlaneToSphereUnitVector, PRECISION / 2)) {
													pathsUsed[2].incrementAndGet(); paths.append(",2");
													/* a) If sphere is fully below plane, do not include mass. */
													/* Or, sphere does not intersect plane and is below the plane -- the pointing vector of plane and vector from plane to sphere center are inverted (not equal). */
													impactMass = null;
													impactCenter = null;
												} else {
													/* b) Compute the portion of the sphere that is above the plane, and include that fraction of mass. */
													/* Or, we know the sphere is fully outside the rindler horizon and intersects the plane, so find the sphere-plane cap that is above the plane */
													Apfloat intersectRadius = MATH.sqrt(sphereRSq.subtract(compute_d.multiply(compute_d)));
													if (intersectRadius.equalDigits(sphereRadius[i]) > PRECISION/2) {
														pathsUsed[3].incrementAndGet(); paths.append(",3");
														// Special case, "bysection" so we can just do 50% of volume;
														/*
														 * 3 * ( 2R - R )^2 / 4 * ( 3R - R ) = 3 * ( R ^2 ) / 4 * (2 R ) = 3 / 8 * R
														 * this is geometric centroid for this special case.
														 * 
														 * 
														 * This math is adjust from CENTER of SPHERE along intersection vector.
														 */
														
														impactMass = sphereMass[i].divide(TWO);
														Apfloat centroid = THREEEIGHTS.multiply(sphereRadius[i]);
														impactCenter = new Location(
																sphereS.x.add(centerOfRotationPointingUnitVector.x.multiply(centroid)), // distance from sphere center; we use the unit vector. Go ALONG unit as we want "above" plane.
																sphereS.y.add(centerOfRotationPointingUnitVector.y.multiply(centroid)),
																sphereS.z.add(centerOfRotationPointingUnitVector.z.multiply(centroid)));
														// note above we have to use centerOfRotationPointingUnitVector because the vector of the intersection will be effectively 0.
														
													} else {
														/*
														 * OK: we can use the angle between fromPlaneToSphereUnitVector and centerOfRotationPointingUnitVector
														 * to determine to decide if we use the cap or the remainder.
														 * They should either be identical or inverted, based on the properties of intersection of
														 * sphere and plane.  
														 * 
														 */
														Apfloat height = null;
														boolean vectorEquiv = false;
														if (fromPlaneToSphereUnitVector.isEqual(centerOfRotationPointingUnitVector, PRECISION/2)) {
															// if they are equal unit vectors, then we are the "large" bit of the cap. So, h will be bigger then radius.
															vectorEquiv = true;
															height = fromPlaneToSphere.magnitude.add(sphereRadius[i]);
														} else {
															// if they are inverted, then we know we're using the "small" bit of the cap. H will be smaller than radius.
															vectorEquiv = false;
															height = sphereRadius[i].subtract(fromPlaneToSphere.magnitude);
														}
														// now volume; multiply with density; gives mass.
														// Now, volume is V=((PI*h^2)/3) *(3*R-h)
														Apfloat volume = PI.multiply(height.multiply(height)).divide(THREE).multiply(sphereRadius[i].multiply(THREE).subtract(height));
														// centroid will be distance of z = 3*(2*R - h) ^2 / 4 * (3*R - h) along vector from sphere center to intersection center
														Apfloat centroid = 
																THREE.multiply( MATH.pow( sphereRadius[i].multiply(TWO).subtract(height), TWO ) )
															.divide(
																FOUR.multiply( sphereRadius[i].multiply(THREE).subtract(height) )
															);
														
														impactMass = volume.multiply(sphereDensity[i]);
														if (vectorEquiv) {
															pathsUsed[4].incrementAndGet(); paths.append(",4");
															// move along unit vector
															impactCenter = new Location(
																	sphereS.x.add(fromPlaneToSphereUnitVector.x.multiply(centroid)),
																	sphereS.y.add(fromPlaneToSphereUnitVector.y.multiply(centroid)),
																	sphereS.z.add(fromPlaneToSphereUnitVector.z.multiply(centroid)));
														} else {
															pathsUsed[5].incrementAndGet(); paths.append(",5");
															// move against unit vector.
															impactCenter = new Location(
																	sphereS.x.subtract(fromPlaneToSphereUnitVector.x.multiply(centroid)), // distance from sphere center; we use the inverse of the unit vector but reduce to subtraction.
																	sphereS.y.subtract(fromPlaneToSphereUnitVector.y.multiply(centroid)),
																	sphereS.z.subtract(fromPlaneToSphereUnitVector.z.multiply(centroid)));
														}
													}
		
												}
											} else {
												/*4) If sphere is partially within rindler horizon: */
												// in all cases we need to compute sphere-sphere intersection, probably.
												/*
												 * lensVolume = PI*(sphereRadius+rindlerHorizonDistance-sphereDistance)^2*
												 * 					(sphereDistance^2+2*sphereDistance*rindlerHorizonDistance-3*rindlerHorizon^2+2*sphereDistance*sphereRadius+6*rindlerHorizon*sphereRadius-3*sphereRadius^2)/(12*sphereDistance)
												 * 
												 * lensMass = lensVolume*sphereDensity
												 * 
												 */
												Apfloat sphereDistanceSq = sphereDistance.multiply(sphereDistance);
												Apfloat smallSphereRadius = rindlerHorizonDistance.compareTo(sphereRadius[i]) > 0 ? sphereRadius[i] : rindlerHorizonDistance;
												Apfloat smallSphereRSq = smallSphereRadius.multiply(smallSphereRadius);
												Apfloat bigSphereRadius = rindlerHorizonDistance.compareTo(sphereRadius[i]) > 0 ? rindlerHorizonDistance : sphereRadius[i];
												Apfloat bigSphereRSq = bigSphereRadius.multiply(bigSphereRadius);
												Apfloat lensVolume = PI
														.multiply( MATH.pow( bigSphereRadius.add(smallSphereRadius).subtract(sphereDistance),TWO) )
														.multiply(
															sphereDistanceSq
																.add(TWO.multiply(sphereDistance).multiply(smallSphereRadius))
																.subtract(THREE.multiply(smallSphereRSq))
																.add(TWO.multiply(sphereDistance).multiply(bigSphereRadius))
																.add(SIX.multiply(smallSphereRadius).multiply(bigSphereRadius))
																.subtract(THREE.multiply(bigSphereRSq))
														).divide(
															TWELVE.multiply(sphereDistance)
														);
												Apfloat lensMass = lensVolume.multiply(sphereDensity[i]);
												
												// NOTE: here, we will compute an effective lensRadius, that describes a sphere with the same volume
												//       as the lens. They share the same center, but the math to do sphere-plane intersection is 
												//       way more accessible.
												Apfloat lensCircleRadius = MATH.root( THREE.divide(FOUR).multiply(lensVolume).divide(PI), 3l);
												
												
												/*
												 * equivCenterMagnitude = (sphereDistance^2-rindlerHorizonDistance^2+sphereRadius^2) / (2 * sphereDistance)
												 * 
												 * equivCenter = <fromObjectToSphere.unitVector> * equivCenterMagnitude;
												 */
												
												Apfloat lensCenterOffset = (sphereDistanceSq.subtract(smallSphereRSq).add(bigSphereRSq)).divide(TWO.multiply(sphereDistance));
												Location lensCenter = new Location(
															midPoint.x.add(fromObjectToSphere.unitVector.x.multiply(lensCenterOffset)),
															midPoint.y.add(fromObjectToSphere.unitVector.y.multiply(lensCenterOffset)),
															midPoint.z.add(fromObjectToSphere.unitVector.z.multiply(lensCenterOffset))
														);
												
												if (sphereRadius[i].compareTo(absCompute_d) <= 0 && !centerOfRotationPointingUnitVector.isEqual(fromPlaneToSphereUnitVector, PRECISION / 2)) {
													pathsUsed[6].incrementAndGet(); paths.append(",6");
													/* a) If sphere is fully below plane, include fraction of mass based on sphere-sphere intersection */
													// Or, sphere is fully below the plane, but intersects the Rindler horizon, so compute the lens that contributes mass.
													/*
													 * https://mathworld.wolfram.com/Sphere-SphereIntersection.html
													 */
													impactMass = lensMass;
													impactCenter = lensCenter;
												} else {
													// just to reinforce, if we're here, sphere is intersecting with plane (definitely not above, definitely not fully below, so definitely "intersected with")
													// AND, not fully outside, and not fully within, rindler horizon, but inbetween.
													// So the question becomes, we have a full sphere cap, and a partial lens, to include, and what are the contributing volumes.
													
													/* Compute sphere-sphere intersection, then compute the portion of the sphere that is above the plane, and portion of the lens below the plane, and include those fractions of mass. */
													// Or, sphere is partially below the plane and intersects with Rindler horizon, so we have a partial lens and a spherical cap above the plane, and must include both mass fractions.
													// We should do them separately; e.g. use 3)b)'s math in entirety and include its gravity impact; then, compute the lens; cut the lens at the plane. final mass impact is spherecap, and "non-hidden"
													// lens segment.
													
													/*
													 * So after much thought, it occurred to me that the lens formed by the intersection of the rindler horizon and the planet sphere would have a computable volume; due to it being a
													 * composition of spherical caps, and the uniformity of density assumption from above, it is easy to consider instead a sphere with the same volume, centered on the centroid of the
													 * lens; this new sphere-from-lens will have the same mass, but be topologically equivalent and geometrically equivalent; the lens-plane intersection will be mass-equivalent to the
													 * sphere-plane intersection of a sphere of identical mass located at the same point and having the same central radius as the lens.
													 * 
													 * 
													 * Now let's do sphere-plane intersection for the sphere at <equivCenter> radius lensCircleRadius. 
													 * The resulting spherical cap's volume multiplied by the density will then be the desired mass of the lens, if we'd been able to bysect it directly, due to topological equivalence.
													 * 
													 * Note: from https://en.wikipedia.org/wiki/Center_of_mass
													 * R = (1 / (m1+m2)) * (m1*r1+m2*r2) where R is vector point; r1 is vector center of mass 1, r2 is vector center of mass 2.
													 * 
													 *  This new center of mass computation will allow us to combine the mass & center of the major spherical cap, with the mass & center of the equivSphere spherical cap.
													 *  
													 */
																								
													// Note there are three cases here; the plane can intersect with the equivsphere, or it could be fully below, or fully above.
													
													Apfloat equivDNumerator = plane_A.multiply(lensCenter.x).add(
															plane_B.multiply(lensCenter.y)).add(
															plane_C.multiply(lensCenter.z)).add(
															plane_D);
													Apfloat equivCompute_d = equivDNumerator.divide(AsqBsqCsqRoot);
													Apfloat equivAbsCompute_d = MATH.abs(equivCompute_d);
													
													/* note here that we leverage notion of "above" and "below" plane, so we compute that. */
													Apfloat equivIntersectX = lensCenter.x.subtract(plane_A.multiply(equivDNumerator).divide(AsqBsqCsq));
													Apfloat equivIntersectY = lensCenter.y.subtract(plane_B.multiply(equivDNumerator).divide(AsqBsqCsq));
													Apfloat equivIntersectZ = lensCenter.z.subtract(plane_C.multiply(equivDNumerator).divide(AsqBsqCsq));
													Location equivIntersectCenter = new Location(equivIntersectX, equivIntersectY, equivIntersectZ);
													
													Vector fromPlaneToEquiv = Vector.betweenPoints(equivIntersectCenter, lensCenter);
													Location fromPlaneToEquivUnitVector = fromPlaneToEquiv.unitVector;
													
													Apfloat addMass = null;
													Location addCenter = null;
													
													// So contrasting with prior use of this spherical cap computations, we IGNORE the bit that's "visible" due to sphere-plane interaction. We INCLUDE the bit that is "below" the accel plane.
													// This is because we've determined above this mass is inside the Rindler horizon; but the bit above the plane is already accounted for by standard cap computation, which we need to include.
													if (lensCircleRadius.compareTo(equivAbsCompute_d) <= 0 && centerOfRotationPointingUnitVector.isEqual(fromPlaneToEquivUnitVector, PRECISION / 2)) {
														pathsUsed[7].incrementAndGet(); paths.append(",7");
														/* If lens-sphere is fully above plane, do not include mass. */
														/* Or, lens-sphere does not intersect plane and is above the plane -- the pointing vector of plane and vector from plane to equivsphere center are equal. */
														
														// then lens equiv sphere is fully inside the sphere's cap, so we don't contribute from the lens sphere at all.
														
														// here we just use the cap of the sphere that is above the plane.
														addMass = null;
														addCenter = null;
													} else if (lensCircleRadius.compareTo(equivAbsCompute_d) <= 0 && !centerOfRotationPointingUnitVector.isEqual(fromPlaneToEquivUnitVector, PRECISION / 2)) {
														pathsUsed[8].incrementAndGet(); paths.append(",8");
														/* If lensesphere is fully below plane, include all its mass. */
														/* Or, lensesphere does not intersect plane and is below the plane -- the pointing vector of plane and vector from plane to equivsphere center are inverted. */
														
														// lens equiv sphere is full outside the sphere's cap and should be included in full.
														// so do lens mass and cap mass, with center of mass being the centroid of them.
														addMass = lensMass;
														addCenter = lensCenter;
													} else {
														/* lens sphere intersects with plane; we include the portion of mass below the plane.
														/* Compute the portion of the equivsphere that is below the plane, and include that fraction of mass. */
														Apfloat intersectEquivRadius = MATH.sqrt(lensCircleRadius.multiply(lensCircleRadius).subtract(equivCompute_d.multiply(equivCompute_d)));
		
														if (intersectEquivRadius.equalDigits(lensCircleRadius) > PRECISION/2) {
															pathsUsed[9].incrementAndGet(); paths.append(",9");
															// Special case, "bysection" so we can just do 50% of volume;
															/*
															 * 3 * ( 2R - R )^2 / 4 * ( 3R - R ) = 3 * ( R ^2 ) / 4 * (2 R ) = 3 / 8 * R
															 * this is geometric centroid for this special case.
															 * 
															 * 
															 * This math is adjust from CENTER of SPHERE along intersection vector.
															 */
															
															addMass = lensMass.divide(TWO);
															Apfloat centroid = THREEEIGHTS.multiply(lensCircleRadius);
															addCenter = new Location(
																	lensCenter.x.subtract(centerOfRotationPointingUnitVector.x.multiply(centroid)), // distance from sphere center; we use the inverse of the unit vector but reduce to subtraction
																	lensCenter.y.subtract(centerOfRotationPointingUnitVector.y.multiply(centroid)), // going for "below" so inverse accel vector
																	lensCenter.z.subtract(centerOfRotationPointingUnitVector.z.multiply(centroid)));
															// we use centerOfRotationPointingUnitVector because the fromPlaneToEquivUnitVector will be effectively zero.
															
														} else {
															/*
															 * OK: we can use the angle between fromPlaneToSphereUnitVector and centerOfRotationPointingUnitVector
															 * to determine to decide if we use the cap or the remainder.
															 * They should either be identical or inverted, based on the properties of intersection of
															 * sphere and plane.  
															 * 
															 */
															if (equivAbsCompute_d.equalDigits(fromPlaneToEquiv.magnitude) < PRECISION/2) {
																System.err.println(String.format("mismatch in vector magnitutde from lens intersect with plane vs. computed distance: \n%#s vs\n%#s",  fromPlaneToEquiv.magnitude, equivAbsCompute_d));
															}
															Apfloat height = null;
															boolean vectorEquiv = false;
															if (fromPlaneToEquivUnitVector.isEqual(centerOfRotationPointingUnitVector, PRECISION/2)) {
																// if they are equal unit vectors, then the "large" bit of the cap is above the plane. We want below, so h will be smaller than radius
																height = lensCircleRadius.subtract(equivAbsCompute_d); // fromPlaneToEquiv.magnitude);
																vectorEquiv = true;
															} else {
																// if they are inverted, then we know the "small" bit of the cap is above the plane. We want below, so H will be larger than the radius.
																height = lensCircleRadius.add(equivAbsCompute_d); //fromPlaneToEquiv.magnitude.add(lensCircleRadius);
																vectorEquiv = false;
															}
															// now volume; multiply with density; gives mass.
															// Now, volume is V=((PI*h^2)/3) *(3*R-h)
															Apfloat volume = PI.multiply(height.multiply(height)).divide(THREE).multiply(lensCircleRadius.multiply(THREE).subtract(height));
															// centroid will be distance of z = 3*(2*R - h) ^2 / 4 * (3*R - h) along vector from sphere center to intersection center
															Apfloat centroid = 
																	THREE.multiply( MATH.pow( (lensCircleRadius.multiply(TWO)).subtract(height), TWO ) )
																.divide(
																	FOUR.multiply( lensCircleRadius.multiply(THREE).subtract(height) )
																);
															
															addMass = volume.multiply(sphereDensity[i]); // density of lens, and density of equivalent sphere, is the same as sphere's volume is same as lens' volume.
															if (vectorEquiv) {
																pathsUsed[10].incrementAndGet(); paths.append(",10");
																// moving against the planetoequivunitvector 
																addCenter = new Location(
																		lensCenter.x.subtract(fromPlaneToEquivUnitVector.x.multiply(centroid)), // distance from sphere center; we use the inverse of the unit vector but reduce to subtraction.
																		lensCenter.y.subtract(fromPlaneToEquivUnitVector.y.multiply(centroid)),
																		lensCenter.z.subtract(fromPlaneToEquivUnitVector.z.multiply(centroid)));
															} else {
																pathsUsed[11].incrementAndGet(); paths.append(",11");
																// moving with the planetoequivunitvector 
																addCenter = new Location(
																		lensCenter.x.add(fromPlaneToEquivUnitVector.x.multiply(centroid)), // distance from sphere center; we use the unit vector.
																		lensCenter.y.add(fromPlaneToEquivUnitVector.y.multiply(centroid)),
																		lensCenter.z.add(fromPlaneToEquivUnitVector.z.multiply(centroid)));														
															}
														}
		
													}
		
													// ok lets now figure out the sphere's cap. we can then join the sphere's mass with the lens fragment mass, and find their shared center of mass.
													
													/* b) Compute the portion of the sphere that is above the plane, and include that fraction of mass. */
													/* Or, find the sphere-plane cap that is above the plane */
													Apfloat capMass = null;
													Location capCenter = null;
													
													Apfloat intersectRadius = MATH.sqrt(sphereRSq.subtract(compute_d.multiply(compute_d)));
													if (intersectRadius.equalDigits(sphereRadius[i]) > PRECISION/2) {
														pathsUsed[12].incrementAndGet(); paths.append(",12");
														// Special case, "bysection" so we can just do 50% of volume;
														/*
														 * 3 * ( 2R - R )^2 / 4 * ( 3R - R ) = 3 * ( R ^2 ) / 4 * (2 R ) = 3 / 8 * R
														 * this is geometric centroid for this special case.
														 * 
														 * 
														 * This math is adjust from CENTER of SPHERE along intersection vector.
														 */
															
														capMass = sphereMass[i].divide(TWO);
														Apfloat centroid = THREEEIGHTS.multiply(sphereRadius[i]);
														capCenter = new Location(
																sphereS.x.add(centerOfRotationPointingUnitVector.x.multiply(centroid)), // distance from sphere center; we use the unit vector. Go ALONG unit as we want "above" plane.
																sphereS.y.add(centerOfRotationPointingUnitVector.y.multiply(centroid)),
																sphereS.z.add(centerOfRotationPointingUnitVector.z.multiply(centroid)));
														// note above we have to use centerOfRotationPointingUnitVector because the vector of the intersection will be effectively 0.
														
													} else {
														/*
														 * OK: we can use the angle between fromPlaneToSphereUnitVector and centerOfRotationPointingUnitVector
														 * to determine to decide if we use the cap or the remainder.
														 * They should either be identical or inverted, based on the properties of intersection of
														 * sphere and plane.  
														 * 
														 */
														Apfloat height = null;
														boolean vectorEquiv = false;
														if (fromPlaneToSphereUnitVector.isEqual(centerOfRotationPointingUnitVector, PRECISION/2)) {
															// if they are equal unit vectors, then we are the "large" bit of the cap. So, h will be bigger then radius.
															vectorEquiv = true;
															height = fromPlaneToSphere.magnitude.add(sphereRadius[i]);
														} else {
															// if they are inverted, then we know we're using the "small" bit of the cap. H will be smaller than radius.
															vectorEquiv = false;
															height = sphereRadius[i].subtract(fromPlaneToSphere.magnitude);
														}
														// now volume; multiply with density; gives mass.
														// Now, volume is V=((PI*h^2)/3) *(3*R-h)
														Apfloat volume = PI.multiply(height.multiply(height)).divide(THREE).multiply(sphereRadius[i].multiply(THREE).subtract(height));
														// centroid will be distance of z = 3*(2*R - h) ^2 / 4 * (3*R - h) along vector from sphere center to intersection center
														Apfloat centroid = 
																THREE.multiply( MATH.pow( sphereRadius[i].multiply(TWO).subtract(height), TWO ) )
															.divide(
																FOUR.multiply( sphereRadius[i].multiply(THREE).subtract(height) )
															);
														
														capMass = volume.multiply(sphereDensity[i]);
														if (vectorEquiv) {
															pathsUsed[13].incrementAndGet(); paths.append(",13");
															// move along unit vector
															capCenter = new Location(
																	sphereS.x.add(fromPlaneToSphereUnitVector.x.multiply(centroid)),
																	sphereS.y.add(fromPlaneToSphereUnitVector.y.multiply(centroid)),
																	sphereS.z.add(fromPlaneToSphereUnitVector.z.multiply(centroid)));
														} else {
															pathsUsed[14].incrementAndGet(); paths.append(",14");
															// move against unit vector.
															capCenter = new Location(
																	sphereS.x.subtract(fromPlaneToSphereUnitVector.x.multiply(centroid)), // distance from sphere center; we use the inverse of the unit vector but reduce to subtraction.
																	sphereS.y.subtract(fromPlaneToSphereUnitVector.y.multiply(centroid)),
																	sphereS.z.subtract(fromPlaneToSphereUnitVector.z.multiply(centroid)));
														}
													}
													
													// now we have an capmass/center from the real sphere cap; and an addmass/center from the lens equiv sphere cap. let's join tem.
													if (addMass == null) {
														pathsUsed[15].incrementAndGet(); paths.append(",15");
														impactMass = capMass;
														impactCenter = capCenter;
													} else {
														pathsUsed[16].incrementAndGet(); paths.append(",16");
														/*
														 * Note: from https://en.wikipedia.org/wiki/Center_of_mass
														 * R = (1 / (m1+m2)) * (m1*r1+m2*r2) where R is vector point; r1 is vector center of mass 1, r2 is vector center of mass 2.
														 */
														impactMass = capMass.add(addMass);
														impactCenter = new Location(
																(capMass.multiply(capCenter.x).add(addMass.multiply(addCenter.x))).divide(impactMass),
																(capMass.multiply(capCenter.y).add(addMass.multiply(addCenter.y))).divide(impactMass),
																(capMass.multiply(capCenter.z).add(addMass.multiply(addCenter.z))).divide(impactMass)
															);
													}
												}
											}
										}
										
										// OK: now we have impactMass / impactCenter. Let's generate the grav impact.
										if (impactMass != null) {
											contributeAcceleration(xX, yY, zZ, paths.substring(1), midPoint, i, impactMass, impactCenter);
										} else {
											contributeAcceleration(xX, yY, zZ, paths.substring(1), midPoint, i, null, null);
										}
									}
								}
								currentDoneSteps.incrementAndGet();
								if (currentDoneSteps.get()%100 == 0) {
									showPartialProgress();
								}
							}
						});
					} else {
						currentDoneSteps.incrementAndGet();						
					}
				}
			}			
		}
	}
	
	public void showFinalProgress() {
		int z = 0;
		/*for (int x = (int) -stepsOnEdge; x < stepsOnEdge; x++) {
			System.out.print("[");
			for (int y = (int) -stepsOnEdge; y < stepsOnEdge; y++) {
				int index = (int) ((x + stepsOnEdge) + ( (stepsOnEdge * 2) * (y + stepsOnEdge) ) + ( (stepsOnEdge * 2) * (stepsOnEdge * 2) * (z + stepsOnEdge)));
				System.out.print(String.format("%1d", aggregateQIAccelerationUpdCounts[index]));
				if (y < stepsOnEdge -1) {
					System.out.print(",");
				}
			}
			System.out.println("]");
		}*/
		
		Apfloat tmpTakenSteps = new Apfloat(currentTakenSteps.get(), PRECISION);
		Location tempGravAccel = new Location(aggregateGravityAcceleration.x.divide(tmpTakenSteps),
				aggregateGravityAcceleration.y.divide(tmpTakenSteps),
				aggregateGravityAcceleration.z.divide(tmpTakenSteps));		
		System.out.println("Final: Overall grav impact: \nQI: " + tempGravAccel.toString() + " Newton: " + newtonAggregateGravityAcceleration.toString());
		for (int i = 0; i < sphereName.length; i++) {
			tempGravAccel = new Location(gravAcceleration[i].x.divide(tmpTakenSteps),
					gravAcceleration[i].y.divide(tmpTakenSteps),
					gravAcceleration[i].z.divide(tmpTakenSteps));		
			System.out.println("Final: " + sphereName[i] + " grav impact: \nQI: " + tempGravAccel.toString() + " Newton: " + newtonGravAcceleration[i].toString());
		}
		
		System.out.println("Paths:");
		for (int i = 0; i < pathsUsed.length; i++) {
			System.out.print(String.format("%40s:%5d ", pathsNamed[i], pathsUsed[i].get()));
			if (i>0 && (i + 1) % 4 == 0) System.out.println();
		}
	}
	
	private synchronized void contributeAcceleration(final int x, final int y, final int z, String paths, final Location midPoint, int i, Apfloat impactMass, Location impactCenter) {
		// sum up the mass!
		// So given that we are in a high acceleration environment, I'm going to ignore the 2c^2/|a|O term from QI that
		// adjusts inertial mass.
		// a = GM/r^2
		Apfloat a = null;
		if (impactMass != null) {
			Vector accelLine = Vector.betweenPoints(midPoint, impactCenter);
			a = constantG.multiply(impactMass).divide(accelLine.magnitude.multiply(accelLine.magnitude));
			Location contribution = new Location(
					accelLine.unitVector.x.multiply(a),
					accelLine.unitVector.y.multiply(a),
					accelLine.unitVector.z.multiply(a));
			Apfloat gx = aggregateGravityAcceleration.x.add(contribution.x);
			Apfloat gy = aggregateGravityAcceleration.y.add(contribution.y);
			Apfloat gz = aggregateGravityAcceleration.z.add(contribution.z);
	
			aggregateGravityAcceleration = new Location(gx, gy, gz);
			
			
			gx = gravAcceleration[i].x.add(contribution.x);
			gy = gravAcceleration[i].y.add(contribution.y);
			gz = gravAcceleration[i].z.add(contribution.z);
			
			gravAcceleration[i] = new Location(gx, gy, gz);
		} 
		
		if (aggregateQIAcceleration != null) {
			int index = (int) ((x + stepsOnEdge) + ( (stepsOnEdge * 2) * (y + stepsOnEdge) ) + ( (stepsOnEdge * 2) * (stepsOnEdge * 2) * (z + stepsOnEdge)));
			aggregateQIAcceleration[index][i] = a;
			aggregateQIPaths[index][i] = paths;
			int count = 0;
			for (int c = 0; c < aggregateQIPaths[index].length; c++) { count += aggregateQIPaths[index][c] != null ? 1 : 0; }
			if (count % aggregateQIAcceleration[index].length == 0) { // basically, this one is done now.
				float r = 0.5f;
				float g = 0.5f;
				float b = 0.5f;
				if (aggregateQIAcceleration[index][0] != null) {
					if (aggregateQIAcceleration[index][0].compareTo(newtonGravAccelerationScale[0]) > 0) { // our current agg is above newton.
						Apfloat rel = newtonGravAccelerationScale[0].divide(aggregateQIAcceleration[index][0]); // newton / current; value [1, 0] where as current grows larger then newton, value trends to 0. 
						r -= 0.5f - (rel.floatValue() / 2f);
					} else {
						Apfloat rel = aggregateQIAcceleration[index][0].divide(newtonGravAccelerationScale[0]); // current / newton; value [1, 0] where as current is closer to newton, value trends to 1; more distance from newton, value trends to 0.
						r += 0.5f - (rel.floatValue() / 2f);
					}
				}
				if (aggregateQIAcceleration[index][1] != null) {
					if (aggregateQIAcceleration[index][1].compareTo(newtonGravAccelerationScale[1]) > 0) { // our current agg is above newton.
						Apfloat rel = newtonGravAccelerationScale[1].divide(aggregateQIAcceleration[index][1]); // newton / current; value [1, 0] where as current grows larger then newton, value trends to 0. 
						g -= 0.5f - (rel.floatValue() / 2f);
					} else {
						Apfloat rel = aggregateQIAcceleration[index][1].divide(newtonGravAccelerationScale[1]); // current / newton; value [1, 0] where as current is closer to newton, value trends to 1; more distance from newton, value trends to 0.
						g += 0.5f - (rel.floatValue() / 2f);
					}
				}
				if (aggregateQIAcceleration[index][2] != null) {
					if (aggregateQIAcceleration[index][2].compareTo(newtonGravAccelerationScale[2]) > 0) { // our current agg is above newton.
						Apfloat rel = newtonGravAccelerationScale[2].divide(aggregateQIAcceleration[index][2]); // newton / current; value [1, 0] where as current grows larger then newton, value trends to 0. 
						b -= 0.5f - (rel.floatValue() / 2f);
					} else {
						Apfloat rel = aggregateQIAcceleration[index][2].divide(newtonGravAccelerationScale[2]); // current / newton; value [1, 0] where as current is closer to newton, value trends to 1; more distance from newton, value trends to 0.
						b += 0.5f - (rel.floatValue() / 2f);
					}
				}
				QIPlatform.uiQueue.offer(new QIPlatform.QIVisualizationKernel((double) x + stepsOnEdge, (double) y + stepsOnEdge, (double) z + stepsOnEdge, r, g, b, 1.0f / (float) (stepsOnEdge * 2), aggregateQIPaths[index]));
			}
		}
	}

	private synchronized void showPartialProgress() {
		double completePercent = ((double) currentDoneSteps.get() / (double) totalSteps) * 100d;
		System.out.println("Completed " + currentDoneSteps + " with " + currentTakenSteps + " taken, or " + String.format("%6.2f", completePercent) + "% complete. " + QIPlatform.processHandler.getActiveCount() + " threads on the job.");
		if (currentDoneSteps.get()%500 == 0 && currentTakenSteps.get() > 0) {
			Apfloat tmpTakenSteps = new Apfloat(currentTakenSteps.get(), PRECISION);
			Location tempGravAccel = new Location(aggregateGravityAcceleration.x.divide(tmpTakenSteps),
					aggregateGravityAcceleration.y.divide(tmpTakenSteps),
					aggregateGravityAcceleration.z.divide(tmpTakenSteps));
			System.out.println("Overall grav impact: \n" + tempGravAccel.toString());
		}
	}

	public static class Location {
		public Apfloat x;
		public Apfloat y;
		public Apfloat z;
		public Location(Apfloat x, Apfloat y, Apfloat z) {
			this.x = x;
			this.y = y;
			this.z = z;
		}
		
		public Location offset(Apfloat x, Apfloat y, Apfloat z) {
			return new Location(this.x.add(x), this.y.add(y), this.z.add(z));
		}
		public Location offset(Location l) {
			return new Location(this.x.add(l.x), this.y.add(l.y), this.z.add(l.z));
		}
		
		public boolean isInsideSphere(Location sphere, Apfloat sphereRadiusSquare) {
			Apfloat squaredDistance = MATH.pow(x.subtract(sphere.x), 2l)
										.add(
									  MATH.pow(y.subtract(sphere.y), 2l))
									    .add(
									  MATH.pow(z.subtract(sphere.z), 2l));
			return squaredDistance.compareTo(sphereRadiusSquare) < 1; // absolute distance^2 <= sphereRadius^2
		}
		
		public boolean isEqual(Location o, long precision) {
			long xD = x.equalDigits(o.x);
			long yD = y.equalDigits(o.y);
			long zD = z.equalDigits(o.z);
			if (xD >= precision && yD >= precision && zD >= precision) {
				return true;
			}
			return false;
		}
		
		@Override
		public String toString() {
			StringBuilder sb = new StringBuilder();
			sb.append(String.format("<%#s,%#s,%#s>", x.precision(28l), y.precision(28l), z.precision(28l)));
			return sb.toString();
		}
	}
	
	public static class Vector {
		public Location point;
		public Location direction;
		public Apfloat magnitude;
		public Location unitVector;
		public Vector() {	
		}
		
		public static Vector pointDirection(Location point, Location direction) {
			Vector tVec = new Vector();
			tVec.point = point;
			tVec.direction = direction;
			tVec.magnitude = MATH.sqrt(
					direction.x.multiply(direction.x).add(
					direction.y.multiply(direction.y)).add(
					direction.z.multiply(direction.z))
				);
			if (tVec.magnitude.equals(Apfloat.ZERO)) {
				tVec.unitVector = new Location(Apfloat.ZERO, Apfloat.ZERO, Apfloat.ZERO);
			} else {
				tVec.unitVector = new Location(
						tVec.direction.x.divide(tVec.magnitude),
						tVec.direction.y.divide(tVec.magnitude),
						tVec.direction.z.divide(tVec.magnitude));
			}
			return tVec;
		}
		public static Vector betweenPoints(Location point1, Location point2) {
			Location direction = new Location(point2.x.subtract(point1.x), point2.y.subtract(point1.y), point2.z.subtract(point1.z));
			return Vector.pointDirection(point1, direction);
		}
	}
}
