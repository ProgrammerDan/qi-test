package com.programmerdan.physics.qi;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.FontMetrics;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.geom.Point2D;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintStream;
import java.util.LinkedList;
import java.util.Properties;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.Executors;
import java.util.concurrent.LinkedBlockingQueue;
import java.util.concurrent.ThreadPoolExecutor;
import java.util.concurrent.TimeUnit;
import java.util.List;
import java.util.Map;

import javax.swing.JFrame;
import javax.swing.JPanel;

import org.apfloat.Apfloat;

/** Platform for execution of QI Tests.
 * 
 * @author ProgrammerDan <programmerdan@gmail.com>
 */
public class QIPlatform {

	public static final long PRECISION = 1024l;
	
	private static int multiProcess = Runtime.getRuntime().availableProcessors() * 2 < 1 ? 2 : Runtime.getRuntime().availableProcessors() * 2;
	protected static ThreadPoolExecutor processHandler = (ThreadPoolExecutor) Executors.newFixedThreadPool(multiProcess);
	protected static BlockingQueue<QIVisualizationKernel> uiQueue = new LinkedBlockingQueue<>();
	protected static ConcurrentHashMap<Double, List<QIVisualizationKernel>> uiStore = new ConcurrentHashMap<>();
	protected static ConcurrentHashMap<Double, Point2D> uiFrames = new ConcurrentHashMap<>();
	
	protected static PrintStream debugFile = null;
	
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
		System.out.println("Constants:\n G (gravitational constant): " + constantsG + "km^3*kg^-1*s^-2");
		Apfloat apConstantsG = safeSet(props, "constants.G", constantsG, "6.67430e-20");
		
		String showVisualization = props.getProperty("visualize", "true");
		System.out.println("Simulation Constraints:\n Visualization: " + showVisualization);
		Boolean apShowVisualization = safeSetBoolean(props, "visualize", showVisualization, "true");

		String sdotSize = props.getProperty("vis.dotSize", "5");
		System.out.println(" Dot Size: " + sdotSize);
		Integer apDotSize = safeSetInteger(props, "vis.dotSize", sdotSize, "5");

		String sborder = props.getProperty("vis.border", "30");
		System.out.println(" Border Size: " + sborder);
		Integer apBorder = safeSetInteger(props, "vis.border", sborder, "30");

		String storeDebug = props.getProperty("debug", "false");
		System.out.println("Development Aides:\n Debug: " + storeDebug);
		Boolean apStoreDebug = safeSetBoolean(props, "debug", storeDebug, "false");

		String debugFile = props.getProperty("debug.file", "debug.log");
		System.out.println(" Debug File: " + debugFile);
		String apDebugFile = safeSetString(props, "debug.file", debugFile, "debug.log");

		if (apStoreDebug) {
			try {
				QIPlatform.debugFile = new PrintStream(new File(apDebugFile));
			} catch (Exception e) {
				System.err.println(" Failed to initialize Debug File!");
			}
		}
		
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
		
		if (Boolean.TRUE.equals(apShowVisualization)) {
			int dotSize = apDotSize;
			int border = apBorder;
			int displaysOnAnEdge = 2;
			
			int size = apResolution.intValue() * 2 * displaysOnAnEdge * dotSize + border * 3; // four displays, 3x3 pixel per display and 60 pixels border
			
			int cnt = (size / (int) (apResolution.intValue() * 2 * dotSize));
			
			int expansion = (apResolution.intValue() * 2 / cnt) * apResolution.intValue() * 2 * dotSize;
			
			JFrame frame = new JFrame("QIPlatform: Rindler Horizon Hiding Via Angular Motion visualization");
			QIVisualization viz = new QIVisualization(border, dotSize, apResolution.intValue() * 2 * dotSize);
			frame.add(viz, BorderLayout.CENTER);
			viz.setSize(size + expansion, size * 2);
			frame.setSize(size + expansion, size * 2);
			frame.setVisible(true);
			frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
			viz.setBackground(Color.BLACK);
			
			viz.setupBuffer();
			frame.setBackground(Color.BLACK);
			viz.addMouseListener(viz);

			processHandler.execute(rHHVAM); // ok let's go.
			
			while (processHandler.getActiveCount() > 0) {
				try {
					viz.repaint();
					Thread.sleep(15l);
				} catch (Exception e){}
			}
		} else {
			processHandler.execute(rHHVAM); // ok let's go.
			
			while (processHandler.getActiveCount() > 0) {
				try {
					Thread.sleep(0l);
				} catch (Exception e){}
			}
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
	
	public Boolean safeSetBoolean(Properties props, String propKey, String propValue, String safeValue) {
		Boolean value = null;
		try {
			value = Boolean.parseBoolean(propValue);
			props.setProperty(propKey, propValue);
		} catch (Exception e) {
			value = Boolean.parseBoolean(safeValue);
			props.setProperty(propKey,  safeValue);
		}
		return value;
	}

	public Integer safeSetInteger(Properties props, String propKey, String propValue, String safeValue) {
		Integer value = null;
		try {
			value = Integer.parseInt(propValue);
			props.setProperty(propKey, propValue);
		} catch (Exception e) {
			value = Integer.parseInt(safeValue);
			props.setProperty(propKey,  safeValue);
		}
		return value;
	}

	public String safeSetString(Properties props, String propKey, String propValue, String safeValue) {
		String value = null;
		try {
			value = propValue;
			props.setProperty(propKey, propValue);
		} catch (Exception e) {
			value = safeValue;
			props.setProperty(propKey, safeValue);
		}
		return value;
	}
	
	public static class QIVisualizationKernel {
		double x;
		double y;
		double z;
		
		float r;
		float g;
		float b;
		float a;
		
		String[] paths;
		
		public QIVisualizationKernel(double x, double y, double z, float r, float g, float b, float a, String[] paths) {
			this.x = x;
			this.y = y;
			this.z = z;
			
			this.r = r;
			this.g = g;
			this.b = b;
			this.a = a;
			
			this.paths = paths;
		}
	}
	
	public class QIVisualization extends JPanel implements MouseListener {
		
		private int border;
		private int expansion;
		private int panel;
		private BufferedImage buffer;
		
		public QIVisualization(int border, int expansion, int panel) {
			this.border = border;
			this.expansion = expansion;
			this.panel = panel;
		}
		
		public void setupBuffer() {
			buffer = new BufferedImage(this.getWidth(), this.getHeight(), BufferedImage.TYPE_INT_ARGB);
			buffer.createGraphics().setBackground(this.getBackground());
		}
		
		@Override
		protected void paintComponent(Graphics g) {
			super.paintComponent(g);
			
			if (buffer != null) {
				Graphics2D g2 = buffer.createGraphics();
				long start = System.currentTimeMillis();
				while ( !uiQueue.isEmpty() && (System.currentTimeMillis() - start < 15l) ) {
					QIVisualizationKernel kernel = uiQueue.poll();
					if (kernel == null) break;
					int xOff = border;
					int yOff = border;
					Color kernelColor = new Color(kernel.r, kernel.g, kernel.b, kernel.a);
					FontMetrics met = g2.getFontMetrics();
					int row = met.getHeight();
					int base = met.getMaxDescent();
					
					g2.setColor(Color.WHITE);
					g2.drawString("x - y", xOff + (panel - met.stringWidth("x - y")) / 2, yOff - base);					
					g2.setColor(kernelColor);
					// x - y 
					g2.fillRect(xOff + (int) kernel.x * this.expansion, yOff + (int) kernel.y * this.expansion, this.expansion, this.expansion);
					
					xOff += border + panel;
					
					g2.setColor(Color.WHITE);
					g2.drawString("x - z", xOff + (panel - met.stringWidth("x - z")) / 2, yOff - base);
					g2.setColor(kernelColor);
					// x - z
					g2.fillRect(xOff + (int) kernel.x * this.expansion, yOff + (int) kernel.z * this.expansion, this.expansion, this.expansion);
					
					yOff += border + panel;
					
					g2.setColor(Color.WHITE);
					g2.drawString("y - z", xOff + (panel - met.stringWidth("y - z")) / 2, yOff - base);
					g2.setColor(kernelColor);
					// y - z
					g2.fillRect(xOff + (int) kernel.y * this.expansion, yOff + (int) kernel.z * this.expansion, this.expansion, this.expansion);
					
					yOff -= border;
					xOff -= border + panel + border;
					g2.clearRect(xOff, yOff, panel+border+border, row * (4 + kernel.paths.length));
					
					g2.setColor(Color.WHITE);
					g2.drawString(String.format("%3.0f", kernel.x), xOff, yOff + row);
					g2.drawString(String.format("%3.0f", kernel.y), xOff, yOff + row*2);
					g2.drawString(String.format("%3.0f", kernel.z), xOff, yOff + row*3);
					for (int c = 0; c < kernel.paths.length; c++) {
						g2.drawString(String.format("%s", kernel.paths[c]), xOff, yOff + row * (4+c));
					}
					
					yOff -= panel;
					xOff += border + panel + border + panel + border;
					
					int mCount = (panel*4+border*3) / panel;
					
					int yJog = (((int) kernel.z) % mCount) * panel;
					int xJog = (((int) kernel.z) / mCount) * panel;
					
					g2.setColor(Color.WHITE);
					g2.drawString("x - y - layers", xOff + (panel - met.stringWidth("x - y - layers")) / 2, yOff - base);
					float[] kernelRaw = kernelColor.getRGBColorComponents(new float[3]);
					g2.setColor(new Color(kernelRaw[0], kernelRaw[1], kernelRaw[2]));
					// x - y 
					g2.fillRect(xOff + xJog + (int) kernel.x * this.expansion, yOff + yJog + (int) kernel.y * this.expansion, this.expansion, this.expansion);
					
					List<QIVisualizationKernel> uiList = null;
					Point2D uiFrame = null;
					if (uiStore.containsKey(kernel.z)) {
						uiList = uiStore.get(kernel.z);
					} else {
						uiList = new LinkedList<QIVisualizationKernel>();
						uiFrame = new Point2D.Float(xOff + xJog, yOff + yJog);
						//System.out.println("New frame: " + uiFrame.getX() + ", " + uiFrame.getY());
						uiFrames.put(kernel.z, uiFrame);
						uiStore.put(kernel.z, uiList);
					}
					
					uiList.add(kernel);
				}
			}
			g.drawImage(buffer, 0, 0, this);
		}

		@Override
		public void mouseClicked(MouseEvent e) {
			int x = e.getX();
			int y = e.getY();
			//System.out.println("Clicked onto " + x + ", " + y);
			for (Map.Entry<Double, Point2D> region : uiFrames.entrySet() ) {
				Point2D topLeft = region.getValue();
				//System.out.print("<" + topLeft.getX() + ", " + topLeft.getY() + ">-<" + (topLeft.getX() + this.panel) + ", " + (topLeft.getY() + this.panel) + ">"); 
				if (x >= topLeft.getX() && x <= (topLeft.getX() + this.panel) &&
					y >= topLeft.getY() && y <= (topLeft.getY() + this.panel)) {
					//System.out.println("Clicked into frame " + region.getKey());
					List<QIVisualizationKernel> store = uiStore.get(region.getKey());
					if (store != null) {
						int localX = (x - (int) topLeft.getX()) / this.expansion;
						int localY = (y - (int) topLeft.getY()) / this.expansion;
						for (QIVisualizationKernel kernel : store) {
							if (localX == (int) kernel.x && localY == (int) kernel.y) {
								System.out.println(String.join(",", kernel.paths));
							}
						}
					}
				}
			}
			
		}

		@Override
		public void mousePressed(MouseEvent e) {
			
		}

		@Override
		public void mouseReleased(MouseEvent e) {
			
		}

		@Override
		public void mouseEntered(MouseEvent e) {
			
		}

		@Override
		public void mouseExited(MouseEvent e) {
			
		}
	}
}
