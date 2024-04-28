package de.labathome;

import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;

public class TestStaticMethods {

	@Test
	void testStripTrailingZeros() {

		Assertions.assertArrayEquals(null, FortranNamelist.stripTrailingZeros((int[])null));
		Assertions.assertArrayEquals(null, FortranNamelist.stripTrailingZeros((double[])null));

		Assertions.assertArrayEquals(new int[] {}, FortranNamelist.stripTrailingZeros(new int[] {}));
		Assertions.assertArrayEquals(new double[] {}, FortranNamelist.stripTrailingZeros(new double[] {}));

		Assertions.assertArrayEquals(new int[] {42}, FortranNamelist.stripTrailingZeros(new int[] {42}));
		Assertions.assertArrayEquals(new double[] {42.0}, FortranNamelist.stripTrailingZeros(new double[] {42.0}));

		Assertions.assertArrayEquals(new int[] {1, 2, 3}, FortranNamelist.stripTrailingZeros(new int[] {1, 2, 3}));
		Assertions.assertArrayEquals(new double[] {1.0, 2.0, 3.0}, FortranNamelist.stripTrailingZeros(new double[] {1.0, 2.0, 3.0}));

		Assertions.assertArrayEquals(new int[] {0, 0, 1, 2, 3}, FortranNamelist.stripTrailingZeros(new int[] {0, 0, 1, 2, 3}));
		Assertions.assertArrayEquals(new double[] {0.0, 0.0, 1.0, 2.0, 3.0}, FortranNamelist.stripTrailingZeros(new double[] {0.0, 0.0, 1.0, 2.0, 3.0}));

		Assertions.assertArrayEquals(new int[] {}, FortranNamelist.stripTrailingZeros(new int[] {0}));
		Assertions.assertArrayEquals(new double[] {}, FortranNamelist.stripTrailingZeros(new double[] {0.0}));

		Assertions.assertArrayEquals(new int[] {}, FortranNamelist.stripTrailingZeros(new int[] {0, 0, 0}));
		Assertions.assertArrayEquals(new double[] {}, FortranNamelist.stripTrailingZeros(new double[] {0.0, 0.0, 0.0}));

		Assertions.assertArrayEquals(new int[] {1, 2, 3}, FortranNamelist.stripTrailingZeros(new int[] {1, 2, 3, 0}));
		Assertions.assertArrayEquals(new double[] {1.0, 2.0, 3.0}, FortranNamelist.stripTrailingZeros(new double[] {1.0, 2.0, 3.0, 0.0}));

		Assertions.assertArrayEquals(new int[] {0, 0, 1, 2, 3}, FortranNamelist.stripTrailingZeros(new int[] {0, 0, 1, 2, 3, 0}));
		Assertions.assertArrayEquals(new double[] {0.0, 0.0, 1.0, 2.0, 3.0}, FortranNamelist.stripTrailingZeros(new double[] {0.0, 0.0, 1.0, 2.0, 3.0, 0.0}));

		Assertions.assertArrayEquals(new int[] {1, 2, 3}, FortranNamelist.stripTrailingZeros(new int[] {1, 2, 3, 0, 0, 0}));
		Assertions.assertArrayEquals(new double[] {1.0, 2.0, 3.0}, FortranNamelist.stripTrailingZeros(new double[] {1.0, 2.0, 3.0, 0.0, 0.0, 0.0}));

		Assertions.assertArrayEquals(new int[] {0, 0, 1, 2, 3}, FortranNamelist.stripTrailingZeros(new int[] {0, 0, 1, 2, 3, 0, 0, 0}));
		Assertions.assertArrayEquals(new double[] {0.0, 0.0, 1.0, 2.0, 3.0}, FortranNamelist.stripTrailingZeros(new double[] {0.0, 0.0, 1.0, 2.0, 3.0, 0.0, 0.0, 0.0}));

		// --------

		Assertions.assertArrayEquals(new double[] {}, FortranNamelist.stripTrailingZeros(new double[] {0.0}, 0));
		Assertions.assertArrayEquals(new double[] {0.0}, FortranNamelist.stripTrailingZeros(new double[] {0.0}, 1));
		Assertions.assertArrayEquals(new double[] {0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, FortranNamelist.stripTrailingZeros(new double[] {0.0, 0.0}, 6));
	}

	@Test
	void testToFortranAsString() {

		Assertions.assertEquals(".true.", FortranNamelist.toFortranAsString(true));
		Assertions.assertEquals(".false.", FortranNamelist.toFortranAsString(false));

		// --------

		double[] dblArray = { 0.0, 0.0, 1.0, 2.0, 3.0, 0.0, 0.0 };

		String refDoubleArrayFull = "";
		refDoubleArrayFull += String.format("% .20e, ", dblArray[0]);
		refDoubleArrayFull += String.format("% .20e, ", dblArray[1]);
		refDoubleArrayFull += String.format("% .20e, ", dblArray[2]);
		refDoubleArrayFull += String.format("% .20e, ", dblArray[3]);
		refDoubleArrayFull += String.format("% .20e, ", dblArray[4]);
		refDoubleArrayFull += String.format("% .20e, ", dblArray[5]);
		refDoubleArrayFull += String.format("% .20e", dblArray[6]);

		String refDoubleArrayNoTrailingZeros = "";
		refDoubleArrayNoTrailingZeros += String.format("% .20e, ", dblArray[0]);
		refDoubleArrayNoTrailingZeros += String.format("% .20e, ", dblArray[1]);
		refDoubleArrayNoTrailingZeros += String.format("% .20e, ", dblArray[2]);
		refDoubleArrayNoTrailingZeros += String.format("% .20e, ", dblArray[3]);
		refDoubleArrayNoTrailingZeros += String.format("% .20e", dblArray[4]);

		Assertions.assertEquals(refDoubleArrayNoTrailingZeros, FortranNamelist.toFortranAsString(dblArray));
		Assertions.assertEquals(refDoubleArrayFull, FortranNamelist.toFortranAsString(dblArray, false));

		// --------

		int[] intArray = { 0, 0, 1, 2, 3, 0, 0 };

		String refIntArrayFull = "";
		refIntArrayFull += String.format("%d, ", intArray[0]);
		refIntArrayFull += String.format("%d, ", intArray[1]);
		refIntArrayFull += String.format("%d, ", intArray[2]);
		refIntArrayFull += String.format("%d, ", intArray[3]);
		refIntArrayFull += String.format("%d, ", intArray[4]);
		refIntArrayFull += String.format("%d, ", intArray[5]);
		refIntArrayFull += String.format("%d", intArray[6]);

		String refIntArrayNoTrailingZeros = "";
		refIntArrayNoTrailingZeros += String.format("%d, ", intArray[0]);
		refIntArrayNoTrailingZeros += String.format("%d, ", intArray[1]);
		refIntArrayNoTrailingZeros += String.format("%d, ", intArray[2]);
		refIntArrayNoTrailingZeros += String.format("%d, ", intArray[3]);
		refIntArrayNoTrailingZeros += String.format("%d", intArray[4]);

		Assertions.assertEquals(refIntArrayNoTrailingZeros, FortranNamelist.toFortranAsString(intArray));
		Assertions.assertEquals(refIntArrayFull, FortranNamelist.toFortranAsString(intArray, false));
	}

}
