package com.programmerdan.physics.qi;

import java.io.File;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.util.Properties;
import java.util.concurrent.Executors;
import java.util.concurrent.ThreadPoolExecutor;
import java.util.concurrent.TimeUnit;

import org.apfloat.Apfloat;

/** Platform for execution of QI Tests.
 * 
 * @author ProgrammerDan <programmerdan@gmail.com>
 */
public class QIPlatform {

	public static final long PRECISION = 1024l;
	
	private static int multiProcess = Runtime.getRuntime().availableProcessors() * 2 < 1 ? 2 : Runtime.getRuntime().availableProcessors() * 2;
	protected static ThreadPoolExecutor processHandler = (ThreadPoolExecutor) Executors.newFixedThreadPool(multiProcess);
	
	public static void main(String[] args) {
		if (args.length > 0 && "angular".equalsIgnoreCase(args[0])) {
			System.out.println("QI-Test simulation:");
			System.out.println(" running test \"angular\"");
			(new QIPlatform()).runAngularMotionHorizonSimulation(args);
		} else {
			System.out.println("QI-Test simulations:");
			System.out.println("  angular - Rindler Horizon Hiding Via Angular Motion simulation.");
		}
	}

	public void runAngularMotionHorizonSimulation(String[] args) {
		Properties props = new Properties();		
		if (args.length > 1) {
			System.out.println("Assuming argument passed is a configuration file: " + args[1]);
			try {
				File config = new File(args[1]);
				props.load(new FileReader(config));
			} catch (IOException ioe) {
				System.err.println("Failed to load configuration from " + args[1] + ": " + ioe.getMessage());
			}
		}
		
		if (props.isEmpty()) {
			System.out.println("Using all default values. Will result in experiment sphere sitting on surface of earth, with earth-moon-sun in linear");
			System.out.println(" order and each at apoapsis.");
		}
		
		// define knobs
		String radiusInKM = props.getProperty("sphere.Radius", ".0002"); // value in KM.
		System.out.println("Rotating Sphere Radius: " + radiusInKM + "km");
		Apfloat apRadiusInKM = safeSet(props, "sphere.Radius", radiusInKM, ".0002");
		
		String rotationRate = props.getProperty("sphere.RPM", "3589000"); // value in RPM
		System.out.println("Rotating Sphere Rotation: " + rotationRate + "rpm");
		Apfloat apRotationRate = safeSet(props, "sphere.RPM", rotationRate, "3589000");
		
		String resolution = props.getProperty("sphere.Resolution", "2"); // subdivisions of radius to "cuboidize" into. Higher values are slower to compute but more accurate.
		System.out.println("Rotating Sphere Subdivision Per Radius: " + resolution);
		Apfloat apResolution = safeSet(props, "sphere.Resolution", resolution, "2");
		
		String lightSpeed = props.getProperty("lightSpeed", "299792.458");
		System.out.println("Speed of Light in a Vacuum: " + lightSpeed + "km/s");
		Apfloat apLightSpeed = safeSet(props, "lightSpeed", lightSpeed, "299792.458");
		
		String sunRadius = props.getProperty("sun.MeanRadius", "695700");
		System.out.println("Sun Mean Volumetric Radius: " + sunRadius + "km");
		Apfloat apSunRadius = safeSet(props, "sun.MeanRadius", sunRadius, "695700");
		
		String sunDensity = props.getProperty("sun.Density", "1408000000000");
		System.out.println("Sun Mean Density: " + sunDensity + "kg/km^3");
		Apfloat apSunDensity = safeSet(props, "sun.Density",  sunDensity, "1408000000000");
		
		String moonRadius = props.getProperty("moon.MeanRadius", "1737.4");
		System.out.println("Moon Mean Volumetric Radius: " + moonRadius + "km");
		Apfloat apMoonRadius = safeSet(props, "moon.MeanRadius", moonRadius, "1737.4");
		
		String moonDensity = props.getProperty("moon.Density", "3344000000000");
		System.out.println("Moon Mean Density: " + moonDensity + "kg/km^3");
		Apfloat apMoonDensity = safeSet(props, "moon.Density", moonDensity, "3344000000000");
		
		String moonOrbitalPosition = props.getProperty("moon.Position", "1");
		System.out.println("Moon Distance from Earth: 0 Perigee to 1 Apogee: " + moonOrbitalPosition);
		Apfloat apMoonOrbitalPosition = safeSet(props, "moon.Position", moonOrbitalPosition, "1");
		
		String moonOffsetAngle = props.getProperty("moon.Rotation", "0");
		System.out.println("Moon Theta Rotation Around the Earth: " + moonOffsetAngle + " degrees");
		Apfloat apMoonOffsetAngle = safeSet(props, "moon.Rotation", moonOffsetAngle, "0");
		
		String earthRadius = props.getProperty("earth.MeanRadius", "6371");
		System.out.println("Earth Mean Volumetric Radius: " + earthRadius + "km");
		Apfloat apEarthRadius = safeSet(props, "earth.MeanRadius", earthRadius, "6371");
		
		String earthDensity = props.getProperty("earth.Density", "5514000000000");
		System.out.println("Earth Mean Density: " + earthDensity + "kg/km^3");
		Apfloat apEarthDensity = safeSet(props, "earth.Density", earthDensity, "5514000000000");
		
		String earthOrbitalPosition = props.getProperty("earth.Position", "1");
		System.out.println("Earth Distance from Sun: 0 Perihelion to 1 Aphelion: " + earthOrbitalPosition);
		Apfloat apEarthOrbitalPosition = safeSet(props, "earth.Position", earthOrbitalPosition, "1");
		
		String earthOffsetAngle = props.getProperty("earth.Rotation", "0");
		System.out.println("Earth Theta Rotation Around the Sun: " + earthOffsetAngle + " degrees");
		Apfloat apEarthOffsetAngle = safeSet(props, "earth.Rotation", earthOffsetAngle, "0");
		
		String constantsG = props.getProperty("constants.G", "6.67430e-20");
		System.out.println("Constants: G (gravitational constant): " + constantsG + "km^3*kg^-1*s^-2");
		Apfloat apConstantsG = safeSet(props, "constants.G", constantsG, "6.67430e-20");
		
		
		if (args.length > 1) {
			try {
				File config = new File(args[1]);
				props.store(new FileOutputStream(config), "Updated Config with Defaults");
			} catch (IOException ioe) {
				System.err.println("Failed to save configuration to " + args[1] + ": " + ioe.getMessage());
			}
		}
		
		// TODO: offsetangle isn't used yet
		RindlerHorizonHidingViaAngularMotion rHHVAM = new RindlerHorizonHidingViaAngularMotion(
				apRadiusInKM, apRotationRate, apResolution, apLightSpeed,
				apSunRadius, apSunDensity,
				apMoonRadius, apMoonDensity, apMoonOrbitalPosition, apMoonOffsetAngle,
				apEarthRadius, apEarthDensity, apEarthOrbitalPosition, apEarthOffsetAngle,
				apConstantsG, PRECISION);
		
		processHandler.execute(rHHVAM); // ok let's go.
		while (processHandler.getActiveCount() > 0) {
			try {
				Thread.sleep(0);
			} catch (Exception e){}
		}
		processHandler.shutdown(); // ok, it's all on queue, now let's wait.
		while (true) {
			try {
				if (!processHandler.awaitTermination(60l, TimeUnit.MINUTES)) {
					continue;
				} else {
					break;
				}
			} catch (InterruptedException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
		
		rHHVAM.showFinalProgress();
	}
	
	public Apfloat safeSet(Properties props, String propKey, String propValue, String safeValue) {
		Apfloat value = null;
		try {
			value = new Apfloat(propValue, PRECISION);
			props.setProperty(propKey, propValue);
		} catch (Exception e) {
			value = new Apfloat(safeValue, PRECISION);
			props.setProperty(propKey, safeValue);
		}
		return value;
	}
}
