package de.jotschi.geoconvert;

import de.jotschi.geoconvert.GeoConvert;
import static org.junit.Assert.assertEquals;

public class GeoConvertTest {

	/**
	 * Tiny test for the toUtm function
	 * 
	 * @throws Exception
	 */
	public void testConvertToUTM() throws Exception {

		double[] utm = GeoConvert.toUtm(48.188767, 16.418349);
		assertEquals(199753.6610440401, utm[0], 0);
		assertEquals(1817294.8635275129, utm[1], 0);
		System.out.println("X: " + utm[0]);
		System.out.println("X: " + utm[1]);

		// double xy1[] = { 48.201566, 16.353718 };
		// double xy2[] = { 48.205227, 16.363395 };
		// int[] cords = GeoConvert.getPos(48.2039039, 16.3574814, xy1, xy2,
		// 300, 300);
		// System.out.println("Cords: " + cords[0] + " " + cords[1] );

	}

}
