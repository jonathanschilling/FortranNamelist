package de.labathome;

import static org.junit.jupiter.api.Assertions.assertNotNull;

import java.io.IOException;
import java.io.InputStream;
import java.nio.charset.StandardCharsets;
import java.util.Arrays;

import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;

/**
 * Tests for a class to parse Fortran namelists into Java classes. This uses the
 * information given in the class definition to constrain the parser.
 */
public class FortranNamelistTest {

	@Test
	public final void basicTest() {

		MGridParameters testClass = new MGridParameters();

		String namelist = """
				&MGRID_NLI
				  MGRID_EXT = 'w7x_conf17_rev'
				  MGRID_MODE = 'R'
				  LSTELL_SYM = .TRUE.
				  RMIN = 4.30
				  RMAX = 6.30
				  ZMIN = -1.20
				  ZMAX = 1.20
				  IR = 211
				  JZ = 241
				  KP = 36
				/
				""";

		// parse the namelist into a Java object
		FortranNamelist parser = new FortranNamelist(namelist, "mgrid_nli", testClass);
		testClass = (MGridParameters) parser.getParsed();
		assertNotNull(testClass);

		// test for values
		Assertions.assertEquals(testClass.mgrid_ext, "w7x_conf17_rev");
		Assertions.assertEquals(testClass.mgrid_mode, "R");
		Assertions.assertEquals(testClass.stellaratorSymmetric, true);
		Assertions.assertEquals(testClass.Rmin, 4.3);
		Assertions.assertEquals(testClass.Rmax, 6.3);
		Assertions.assertEquals(testClass.Zmin, -1.2);
		Assertions.assertEquals(testClass.Zmax, 1.2);
		Assertions.assertEquals(testClass.numR, 211);
		Assertions.assertEquals(testClass.numZ, 241);
		Assertions.assertEquals(testClass.numPhi, 36);
	}

	@Test
	void vmecInputTest_AsMessyAsPossible() throws IOException {
		String inputFilename = "/de/labathome/FortranNamelist/input/input.AsMessyAsPossible";
		InputStream inputFileStream = getClass().getResourceAsStream(inputFilename);

		String inputFile = new String(inputFileStream.readAllBytes(), StandardCharsets.UTF_8);

		VmecIndataNamelist vmecInput = new VmecIndataNamelist();

		// parse the namelist into a Java object
		FortranNamelist parser = new FortranNamelist(inputFile, "indata", vmecInput);
		vmecInput = (VmecIndataNamelist) parser.getParsed();
		assertNotNull(vmecInput);
		vmecInput.sanitize();

		// Numerical Resolution and Symmetry Assumptions
		Assertions.assertEquals(false, vmecInput.lasym);
		Assertions.assertEquals(5, vmecInput.nfp);
		Assertions.assertEquals(8, vmecInput.mpol);
		Assertions.assertEquals(8, vmecInput.ntor);
		Assertions.assertEquals(22, vmecInput.ntheta);
		Assertions.assertEquals(36, vmecInput.nzeta);

		// Multi-Grid Steps
		int[] ns_array = new int[100];
		ns_array[0] = 51;
		Assertions.assertArrayEquals(ns_array, vmecInput.ns_array);

		double[] ftol_array = new double[100];
		ftol_array[0] = 1.0e-19;
		Assertions.assertArrayEquals(ftol_array, vmecInput.ftol_array);

		int[] niter_array = new int[100];
		Arrays.fill(niter_array, 100000);
		Assertions.assertArrayEquals(niter_array, vmecInput.niter_array);

		// Global Physics Parameters
		Assertions.assertEquals(1.86086390731511, vmecInput.phiedge);
		Assertions.assertEquals(1, vmecInput.ncurr);

		// Profile of Mass or Pressure
		Assertions.assertEquals("two_power", vmecInput.pmass_type.trim());

		double[] am = new double[21];
		am[0] = 66978.0973204654;
		am[1] = 1.90820658424883;
		am[2] = 4.46347452430851;
		Assertions.assertArrayEquals(am, vmecInput.am);

		Assertions.assertArrayEquals(new double[VmecIndataNamelist.nsd], vmecInput.am_aux_s);
		Assertions.assertArrayEquals(new double[VmecIndataNamelist.nsd], vmecInput.am_aux_f);

		Assertions.assertEquals(1.0, vmecInput.pres_scale);
		Assertions.assertEquals(0.0, vmecInput.gamma);
		Assertions.assertEquals(1.0, vmecInput.spres_ped);

		// (Initial Guess for) Rotational Transform Profile
		Assertions.assertEquals("power_series", vmecInput.piota_type.trim());

		double[] ai = new double[21];
		ai[3] = 1.0;
		ai[4] = 1.0;
		Assertions.assertArrayEquals(ai, vmecInput.ai);

		Assertions.assertArrayEquals(new double[VmecIndataNamelist.nsd], vmecInput.ai_aux_s);
		Assertions.assertArrayEquals(new double[VmecIndataNamelist.nsd], vmecInput.ai_aux_f);

		// (Initial Guess for) Toroidal Current Profile
		Assertions.assertEquals("sum_atan", vmecInput.pcurr_type.trim());

		double[] ac = new double[21];
		ac[1] = 1.0;
		ac[2] = 4.65367221956161;
		ac[3] = 2.27775992394378;
		ac[4] = 1.0;
		ac[5] = -0.943974780239535;
		ac[6] = 4.97719715801010;
		ac[7] = 1.0;
		ac[8] = 1.0;
		Assertions.assertArrayEquals(ac, vmecInput.ac);

		Assertions.assertArrayEquals(new double[VmecIndataNamelist.nsd], vmecInput.ac_aux_s);
		Assertions.assertArrayEquals(new double[VmecIndataNamelist.nsd], vmecInput.ac_aux_f);

		Assertions.assertEquals(-913.791680825795, vmecInput.curtor);
		Assertions.assertEquals(1.0, vmecInput.bloat);

		// Free-Boundary Parameters
		Assertions.assertEquals(true, vmecInput.lfreeb);
		Assertions.assertEquals("/u/jsch/v3fit_runs/MGRIDS/mgrid_w7x_v2_10mm_grid.nc", vmecInput.mgrid_file.trim());

		double[] extcur = new double[VmecIndataNamelist.nigroup];
		extcur[0] = -12989.0;
		extcur[1] = -12989.0;
		extcur[2] = -12989.0;
		extcur[3] = -12989.0;
		extcur[4] = -12989.0;
		Assertions.assertArrayEquals(extcur, vmecInput.extcur);

		Assertions.assertEquals(10, vmecInput.nvacskip);

		// Tweaking Parameters
		Assertions.assertEquals(200, vmecInput.nstep);

		double[] aphi = new double[20];
		aphi[0] = 1.0;
		aphi[14] = 1.4;
		Assertions.assertArrayEquals(aphi, vmecInput.aphi);

		Assertions.assertEquals(0.9, vmecInput.delt);
		Assertions.assertEquals(1.0, vmecInput.tcon0);
		Assertions.assertEquals(false, vmecInput.lforbal);

		// Initial Guess for Magnetic Axis Geometry
		double[] raxis_cc = new double[VmecIndataNamelist.ntord + 1];
		raxis_cc[0] = 5.60519782116720;
		raxis_cc[1] = 0.359696693960237;
		raxis_cc[2] = 1.335514497262088E-002;
		raxis_cc[3] = 7.259351998922190E-004;
		raxis_cc[4] = -6.817192552674108E-005;
		raxis_cc[5] = -4.866604255356447E-005;
		raxis_cc[6] = -9.461962764212337E-005;
		raxis_cc[7] = -2.828730415733434E-005;
		raxis_cc[8] = -1.554463199957404E-004;
		Assertions.assertArrayEquals(raxis_cc, vmecInput.raxis_cc);

		double[] zaxis_cs = new double[VmecIndataNamelist.ntord + 1];
		zaxis_cs[1] = -0.301269608014063;
		zaxis_cs[2] = -1.618262861498088E-002;
		zaxis_cs[3] = 1.319614987094394E-004;
		zaxis_cs[4] = 7.237442086550272E-005;
		zaxis_cs[5] = -7.887293984398275E-005;
		zaxis_cs[6] = -5.489994367001566E-006;
		zaxis_cs[7] = 2.463216529452668E-005;
		zaxis_cs[8] = -1.340415732215003E-004;
		Assertions.assertArrayEquals(zaxis_cs, vmecInput.zaxis_cs);

		double[] raxis_cs = new double[VmecIndataNamelist.ntord + 1];
		raxis_cs[1] = -2.093456611578366E-016;
		raxis_cs[2] = -4.361367940788264E-016;
		raxis_cs[3] = -2.267911329209897E-015;
		raxis_cs[4] = -1.500310571631163E-015;
		raxis_cs[5] = -1.395637741052244E-016;
		raxis_cs[6] = 1.290964910473326E-015;
		raxis_cs[7] = -1.116510192841796E-015;
		raxis_cs[8] = 2.233020385683591E-015;
		Assertions.assertArrayEquals(raxis_cs, vmecInput.raxis_cs);

		double[] zaxis_cc = new double[VmecIndataNamelist.ntord + 1];
		zaxis_cc[1] = 3.679904200040097E-017;
		zaxis_cc[2] = -4.906538933386797E-018;
		zaxis_cc[3] = 6.814637407481662E-018;
		zaxis_cc[4] = -2.807630611882445E-017;
		zaxis_cc[5] = -6.187690765993350E-017;
		zaxis_cc[6] = 8.450150385277261E-018;
		zaxis_cc[7] = -5.642519773394817E-017;
		zaxis_cc[8] = 4.497660688937897E-018;
		Assertions.assertArrayEquals(zaxis_cc, vmecInput.zaxis_cc);

		// (Initial Guess for) Boundary Geometry

		double[][] rbc = new double[2 * VmecIndataNamelist.ntord + 1][VmecIndataNamelist.mpold + 1];
		rbc[VmecIndataNamelist.ntord + 0][0] = 5.53767206572562;

		// set to zero again in line below
//		rbc[VmecIndataNamelist.ntord + 1][0] = 0.274510550203623;

		rbc[VmecIndataNamelist.ntord + 2][0] = -6.981826227607880E-003;
		rbc[VmecIndataNamelist.ntord + 3][0] = 1.813036409387450E-003;
		rbc[VmecIndataNamelist.ntord + 4][0] = -1.816182024820829E-003;
		rbc[VmecIndataNamelist.ntord + 5][0] = -9.457063618548155E-007;
		rbc[VmecIndataNamelist.ntord + 6][0] = -3.598476797219176E-004;
		rbc[VmecIndataNamelist.ntord + 7][0] = 5.201212576617179E-005;
		rbc[VmecIndataNamelist.ntord + 8][0] = -3.252552333940336E-004;
		//
		rbc[VmecIndataNamelist.ntord - 8][1] = 1.896612909948546E-004;
		rbc[VmecIndataNamelist.ntord - 7][1] = -2.456907136966591E-004;
		rbc[VmecIndataNamelist.ntord - 6][1] = 2.171005552378400E-004;
		rbc[VmecIndataNamelist.ntord - 5][1] = -3.651594436424253E-004;
		rbc[VmecIndataNamelist.ntord - 4][1] = -5.338074727927220E-004;
		rbc[VmecIndataNamelist.ntord - 3][1] = 2.506984615739400E-004;
		rbc[VmecIndataNamelist.ntord - 2][1] = -8.846005854412054E-004;

		// set to zero again in line below
//		rbc[VmecIndataNamelist.ntord - 1][1] = 1.678418664958850E-002;

		rbc[VmecIndataNamelist.ntord + 0][1] = 0.495662386868802;
		rbc[VmecIndataNamelist.ntord + 1][1] = -0.206910682038032;
		rbc[VmecIndataNamelist.ntord + 2][1] = -1.680826411787574E-002;
		rbc[VmecIndataNamelist.ntord + 3][1] = 7.103573636538324E-004;
		rbc[VmecIndataNamelist.ntord + 4][1] = 4.245820612023922E-004;
		rbc[VmecIndataNamelist.ntord + 5][1] = 6.969490097295169E-004;
		rbc[VmecIndataNamelist.ntord + 6][1] = -7.039468898728505E-004;
		rbc[VmecIndataNamelist.ntord + 7][1] = 5.333996925958589E-004;
		rbc[VmecIndataNamelist.ntord + 8][1] = 6.689061049649439E-005;
		//
		rbc[VmecIndataNamelist.ntord - 8][2] = 1.151620243875582E-004;
		rbc[VmecIndataNamelist.ntord - 7][2] = -3.184094551899957E-004;
		rbc[VmecIndataNamelist.ntord - 6][2] = 4.300662767725242E-004;
		rbc[VmecIndataNamelist.ntord - 5][2] = -3.210600174477836E-004;
		rbc[VmecIndataNamelist.ntord - 4][2] = 2.624531379976409E-004;
		rbc[VmecIndataNamelist.ntord - 3][2] = 5.677789104917515E-005;
		rbc[VmecIndataNamelist.ntord - 2][2] = 2.317403085550396E-003;
		rbc[VmecIndataNamelist.ntord - 1][2] = 9.514153491031787E-003;
		rbc[VmecIndataNamelist.ntord + 0][2] = 3.302699802057374E-002;
		rbc[VmecIndataNamelist.ntord + 1][2] = 4.021410132068875E-002;
		rbc[VmecIndataNamelist.ntord + 2][2] = 6.948281839736793E-002;
		rbc[VmecIndataNamelist.ntord + 3][2] = -8.442177055766123E-004;
		rbc[VmecIndataNamelist.ntord + 4][2] = -1.591791378279287E-004;
		rbc[VmecIndataNamelist.ntord + 5][2] = 2.486332537387373E-005;
		rbc[VmecIndataNamelist.ntord + 6][2] = -7.528422072335360E-005;
		rbc[VmecIndataNamelist.ntord + 7][2] = -8.021064997689406E-005;
		rbc[VmecIndataNamelist.ntord + 8][2] = -5.237201084004603E-005;
		//
		rbc[VmecIndataNamelist.ntord - 8][3] = 1.714830561905995E-004;
		rbc[VmecIndataNamelist.ntord - 7][3] = 7.083574967044192E-005;
		rbc[VmecIndataNamelist.ntord - 6][3] = -1.030190624061461E-004;
		rbc[VmecIndataNamelist.ntord - 5][3] = 7.464723877566652E-005;
		rbc[VmecIndataNamelist.ntord - 4][3] = 1.594997910406184E-004;
		rbc[VmecIndataNamelist.ntord - 3][3] = -1.315594838027009E-004;
		rbc[VmecIndataNamelist.ntord - 2][3] = 1.896710074727297E-004;
		rbc[VmecIndataNamelist.ntord - 1][3] = -4.186302609439121E-004;
		rbc[VmecIndataNamelist.ntord + 0][3] = 3.311501954692046E-004;
		rbc[VmecIndataNamelist.ntord + 1][3] = -8.268632528457621E-003;
		rbc[VmecIndataNamelist.ntord + 2][3] = -2.035417064634164E-002;
		rbc[VmecIndataNamelist.ntord + 3][3] = -1.448421989297009E-002;
		rbc[VmecIndataNamelist.ntord + 4][3] = 7.509697962325974E-004;
		rbc[VmecIndataNamelist.ntord + 5][3] = -1.328389193189970E-004;
		rbc[VmecIndataNamelist.ntord + 6][3] = 1.474276373936834E-004;
		rbc[VmecIndataNamelist.ntord + 7][3] = 1.042690595884966E-004;
		rbc[VmecIndataNamelist.ntord + 8][3] = 6.708871806270065E-005;
		//
		rbc[VmecIndataNamelist.ntord - 8][4] = -1.440375584040593E-004;
		rbc[VmecIndataNamelist.ntord - 7][4] = -1.704752360936465E-004;
		rbc[VmecIndataNamelist.ntord - 6][4] = 9.011920036675072E-005;
		rbc[VmecIndataNamelist.ntord - 5][4] = -1.499685240577500E-004;
		rbc[VmecIndataNamelist.ntord - 4][4] = 5.756550703305372E-005;
		rbc[VmecIndataNamelist.ntord - 3][4] = 1.900144565741371E-004;
		rbc[VmecIndataNamelist.ntord - 2][4] = 4.923760064924362E-005;
		rbc[VmecIndataNamelist.ntord - 1][4] = -5.103652305253303E-004;
		rbc[VmecIndataNamelist.ntord + 0][4] = 2.486409772891553E-003;
		rbc[VmecIndataNamelist.ntord + 1][4] = 3.763697963319822E-003;
		rbc[VmecIndataNamelist.ntord + 2][4] = 9.220272047894581E-003;
		rbc[VmecIndataNamelist.ntord + 3][4] = 4.017321543601945E-003;
		rbc[VmecIndataNamelist.ntord + 4][4] = 9.476230338947471E-004;
		rbc[VmecIndataNamelist.ntord + 5][4] = -7.056521343060718E-004;
		rbc[VmecIndataNamelist.ntord + 6][4] = -6.013036923002932E-005;
		rbc[VmecIndataNamelist.ntord + 7][4] = -8.827308310929046E-005;
		rbc[VmecIndataNamelist.ntord + 8][4] = -7.602682766245118E-005;
		//
		rbc[VmecIndataNamelist.ntord - 8][5] = 2.153814980670956E-005;
		rbc[VmecIndataNamelist.ntord - 7][5] = 1.183952494315462E-005;
		rbc[VmecIndataNamelist.ntord - 6][5] = -7.154739133251302E-005;
		rbc[VmecIndataNamelist.ntord - 5][5] = -1.003063968491245E-004;
		rbc[VmecIndataNamelist.ntord - 4][5] =  1.646390554518521E-004;
		rbc[VmecIndataNamelist.ntord - 3][5] = -5.492928459569241E-005;
		rbc[VmecIndataNamelist.ntord - 2][5] = -1.008690713410809E-004;
		rbc[VmecIndataNamelist.ntord - 1][5] = 5.052245395750897E-004;
		rbc[VmecIndataNamelist.ntord + 0][5] = 6.070892885457450E-004;
		rbc[VmecIndataNamelist.ntord + 1][5] = 7.467412314074342E-004;
		rbc[VmecIndataNamelist.ntord + 2][5] = -1.665044640784653E-003;
		rbc[VmecIndataNamelist.ntord + 3][5] = -1.332894131049668E-003;
		rbc[VmecIndataNamelist.ntord + 4][5] = 8.621325208882533E-004;
		rbc[VmecIndataNamelist.ntord + 5][5] = 2.267011310380025E-004;
		rbc[VmecIndataNamelist.ntord + 6][5] = 1.753152728586062E-004;
		rbc[VmecIndataNamelist.ntord + 7][5] = 2.258512907309003E-005;
		rbc[VmecIndataNamelist.ntord + 8][5] = 4.622663427651913E-005;
		//
		rbc[VmecIndataNamelist.ntord - 8][6] = 2.261179964816800E-005;
		rbc[VmecIndataNamelist.ntord - 7][6] = -1.811117336074824E-005;
		rbc[VmecIndataNamelist.ntord - 6][6] = 7.536505110908006E-005;
		rbc[VmecIndataNamelist.ntord - 5][6] = -3.456196403903005E-006;
		rbc[VmecIndataNamelist.ntord - 4][6] = -5.139872953188579E-005;
		rbc[VmecIndataNamelist.ntord - 3][6] = 3.979845546091016E-005;
		rbc[VmecIndataNamelist.ntord - 2][6] = 1.243950000533798E-004;
		rbc[VmecIndataNamelist.ntord - 1][6] = -1.106699933734676E-004;
		rbc[VmecIndataNamelist.ntord + 0][6] = -1.223947247241717E-004;
		rbc[VmecIndataNamelist.ntord + 1][6] = -7.731238353022811E-005;
		rbc[VmecIndataNamelist.ntord + 2][6] = -1.569941375256721E-003;
		rbc[VmecIndataNamelist.ntord + 3][6] = -1.222191707450016E-003;
		rbc[VmecIndataNamelist.ntord + 4][6] = -1.140339566197411E-003;
		rbc[VmecIndataNamelist.ntord + 5][6] = -9.518279470355112E-004;
		rbc[VmecIndataNamelist.ntord + 6][6] = 2.558431811535456E-004;
		rbc[VmecIndataNamelist.ntord + 7][6] = 2.595101629939122E-005;
		rbc[VmecIndataNamelist.ntord + 8][6] = 1.620499278874425E-005;
		//
		rbc[VmecIndataNamelist.ntord - 8][7] = -1.398163883777544E-005;
		rbc[VmecIndataNamelist.ntord - 7][7] = -3.923123999457625E-005;
		rbc[VmecIndataNamelist.ntord - 6][7] = -2.306491449291393E-005;
		rbc[VmecIndataNamelist.ntord - 5][7] = -1.158175514942671E-004;
		rbc[VmecIndataNamelist.ntord - 4][7] = 9.504349582371491E-005;
		rbc[VmecIndataNamelist.ntord - 3][7] = -7.789561325605551E-005;
		rbc[VmecIndataNamelist.ntord - 2][7] = 3.045433428931733E-005;
		rbc[VmecIndataNamelist.ntord - 1][7] = -1.527742119656307E-005;
		rbc[VmecIndataNamelist.ntord + 0][7] = -1.998021858917175E-005;
		rbc[VmecIndataNamelist.ntord + 1][7] = 4.953044093333177E-005;
		rbc[VmecIndataNamelist.ntord + 2][7] = 4.714962252031772E-004;
		rbc[VmecIndataNamelist.ntord + 3][7] = 5.784116518161847E-004;
		rbc[VmecIndataNamelist.ntord + 4][7] = 7.049240981285555E-004;
		rbc[VmecIndataNamelist.ntord + 5][7] = -1.675739529392045E-004;
		rbc[VmecIndataNamelist.ntord + 6][7] = -5.586250774777524E-005;
		rbc[VmecIndataNamelist.ntord + 7][7] = -2.029459227327149E-004;
		rbc[VmecIndataNamelist.ntord + 8][7] = 5.906682990243203E-006;
		//
		Assertions.assertArrayEquals(rbc, vmecInput.rbc);

		double[][] zbs = new double[2 * VmecIndataNamelist.ntord + 1][VmecIndataNamelist.mpold + 1];
		zbs[VmecIndataNamelist.ntord + 0][0] = 0.0;
		zbs[VmecIndataNamelist.ntord + 1][0] = -0.231215947095178;
		zbs[VmecIndataNamelist.ntord + 2][0] = 3.000448318252753E-004;
		zbs[VmecIndataNamelist.ntord + 3][0] = 2.451759201610496E-004;
		zbs[VmecIndataNamelist.ntord + 4][0] = 2.002740558876558E-003;
		zbs[VmecIndataNamelist.ntord + 5][0] = 2.042500448056660E-004;
		zbs[VmecIndataNamelist.ntord + 6][0] = 1.110653375943014E-004;
		zbs[VmecIndataNamelist.ntord + 7][0] = 1.233633688773586E-004;
		zbs[VmecIndataNamelist.ntord + 8][0] = -5.454884672029991E-004;
		//
		zbs[VmecIndataNamelist.ntord - 8][1] = -6.512120712505940E-004;
		zbs[VmecIndataNamelist.ntord - 7][1] = -2.473300842408858E-004;
		zbs[VmecIndataNamelist.ntord - 6][1] = 1.226721407577620E-003;
		zbs[VmecIndataNamelist.ntord - 5][1] = 1.157564964099680E-005;
		zbs[VmecIndataNamelist.ntord - 4][1] = 1.948945874541883E-004;
		zbs[VmecIndataNamelist.ntord - 3][1] = 1.584625523874076E-003;
		zbs[VmecIndataNamelist.ntord - 2][1] = 5.725438140714447E-003;
		zbs[VmecIndataNamelist.ntord - 1][1] = 6.183781833751985E-003;
		zbs[VmecIndataNamelist.ntord + 0][1] = 0.628787962252855;
		zbs[VmecIndataNamelist.ntord + 1][1] = 0.223524141570276;
		zbs[VmecIndataNamelist.ntord + 2][1] = 1.323272993170185E-002;
		zbs[VmecIndataNamelist.ntord + 3][1] = -6.316519355198045E-004;
		zbs[VmecIndataNamelist.ntord + 4][1] = -4.160839711727137E-004;
		zbs[VmecIndataNamelist.ntord + 5][1] = -1.431427743599354E-004;
		zbs[VmecIndataNamelist.ntord + 6][1] = -3.965264842806391E-004;
		zbs[VmecIndataNamelist.ntord + 7][1] = -2.485908172021094E-005;
		zbs[VmecIndataNamelist.ntord + 8][1] = 3.549072371917396E-004;
		//
		zbs[VmecIndataNamelist.ntord - 8][2] = -3.300731983353361E-004;
		zbs[VmecIndataNamelist.ntord - 7][2] = 4.473218222210373E-004;
		zbs[VmecIndataNamelist.ntord - 6][2] = -6.065256605704031E-004;
		zbs[VmecIndataNamelist.ntord - 5][2] = 3.241308593613608E-004;
		zbs[VmecIndataNamelist.ntord - 4][2] = 6.840233808173711E-005;
		zbs[VmecIndataNamelist.ntord - 3][2] = -9.056185242554235E-005;
		zbs[VmecIndataNamelist.ntord - 2][2] = 1.955219044278804E-004;
		zbs[VmecIndataNamelist.ntord - 1][2] = 7.523213388325822E-003;
		zbs[VmecIndataNamelist.ntord + 0][2] = -5.215130122068161E-003;
		zbs[VmecIndataNamelist.ntord + 1][2] = 2.382604615534300E-002;
		zbs[VmecIndataNamelist.ntord + 2][2] = -5.184435337334851E-002;
		zbs[VmecIndataNamelist.ntord + 3][2] = 2.258335290509644E-003;
		zbs[VmecIndataNamelist.ntord + 4][2] = 6.590696228657478E-004;
		zbs[VmecIndataNamelist.ntord + 5][2] = 9.616300338029773E-005;
		zbs[VmecIndataNamelist.ntord + 6][2] = 1.018064657002039E-004;
		zbs[VmecIndataNamelist.ntord + 7][2] = -6.953157311656576E-005;
		zbs[VmecIndataNamelist.ntord + 8][2] = 9.909633701842818E-005;
		//
		zbs[VmecIndataNamelist.ntord - 8][3] = -4.519578534113128E-005;
		zbs[VmecIndataNamelist.ntord - 7][3] = -8.527405237368145E-005;
		zbs[VmecIndataNamelist.ntord - 6][3] = 1.753730817444603E-004;
		zbs[VmecIndataNamelist.ntord - 5][3] = -2.326348438307375E-005;
		zbs[VmecIndataNamelist.ntord - 4][3] = -2.624975041422071E-004;
		zbs[VmecIndataNamelist.ntord - 3][3] = -1.390802867156250E-005;
		zbs[VmecIndataNamelist.ntord - 2][3] = 8.074618587295024E-004;
		zbs[VmecIndataNamelist.ntord - 1][3] = -1.531636429051088E-003;
		zbs[VmecIndataNamelist.ntord + 0][3] = -1.645890760289522E-003;
		zbs[VmecIndataNamelist.ntord + 1][3] = -8.199761340258446E-003;
		zbs[VmecIndataNamelist.ntord + 2][3] = 5.766395955724331E-003;
		zbs[VmecIndataNamelist.ntord + 3][3] = 1.163122094676383E-002;
		zbs[VmecIndataNamelist.ntord + 4][3] = -5.968673861601821E-004;
		zbs[VmecIndataNamelist.ntord + 5][3] = 1.622756876768911E-004;
		zbs[VmecIndataNamelist.ntord + 6][3] = -3.916688872863889E-005;
		zbs[VmecIndataNamelist.ntord + 7][3] = -6.313722893293626E-005;
		zbs[VmecIndataNamelist.ntord + 8][3] = -3.394782233208609E-005;
		//
		zbs[VmecIndataNamelist.ntord - 8][4] = -2.405330978649966E-005;
		zbs[VmecIndataNamelist.ntord - 7][4] = 2.084090191207986E-004;
		zbs[VmecIndataNamelist.ntord - 6][4] = -6.937710034163128E-005;
		zbs[VmecIndataNamelist.ntord - 5][4] = 3.866322821755758E-005;
		zbs[VmecIndataNamelist.ntord - 4][4] = 1.296929167587253E-004;
		zbs[VmecIndataNamelist.ntord - 3][4] = -1.292084644773919E-004;
		zbs[VmecIndataNamelist.ntord - 2][4] = -9.511070777900972E-005;
		zbs[VmecIndataNamelist.ntord - 1][4] = -2.460373856685125E-004;
		zbs[VmecIndataNamelist.ntord + 0][4] = 4.198477990498741E-003;
		zbs[VmecIndataNamelist.ntord + 1][4] = -2.836264233928369E-003;
		zbs[VmecIndataNamelist.ntord + 2][4] = 9.239409675099113E-003;
		zbs[VmecIndataNamelist.ntord + 3][4] = -4.334921025609129E-003;
		zbs[VmecIndataNamelist.ntord + 4][4] = -1.784089978532210E-003;
		zbs[VmecIndataNamelist.ntord + 5][4] = 4.675626087533873E-004;
		zbs[VmecIndataNamelist.ntord + 6][4] = -3.540448326425835E-005;
		zbs[VmecIndataNamelist.ntord + 7][4] = 5.212974445405975E-005;
		zbs[VmecIndataNamelist.ntord + 8][4] = 2.676132432009841E-005;
		//
		zbs[VmecIndataNamelist.ntord - 8][5] = 3.279177584182842E-005;
		zbs[VmecIndataNamelist.ntord - 7][5] = 2.724539949718756E-005;
		zbs[VmecIndataNamelist.ntord - 6][5] = -1.593085585738485E-004;
		zbs[VmecIndataNamelist.ntord - 5][5] = 2.410957443918596E-004;
		zbs[VmecIndataNamelist.ntord - 4][5] = -1.044020785453555E-004;
		zbs[VmecIndataNamelist.ntord - 3][5] = 2.056825728184085E-005;
		zbs[VmecIndataNamelist.ntord - 2][5] = -3.741888572862706E-005;
		zbs[VmecIndataNamelist.ntord - 1][5] = 3.578584050879604E-004;
		zbs[VmecIndataNamelist.ntord + 0][5] = 6.752128022657002E-004;
		zbs[VmecIndataNamelist.ntord + 1][5] = 1.854756872143505E-003;
		zbs[VmecIndataNamelist.ntord + 2][5] = -1.507318481326176E-003;
		zbs[VmecIndataNamelist.ntord + 3][5] = -9.597505773106182E-004;
		zbs[VmecIndataNamelist.ntord + 4][5] = 9.885927576869936E-004;
		zbs[VmecIndataNamelist.ntord + 5][5] = 1.133619037566474E-005;
		zbs[VmecIndataNamelist.ntord + 6][5] = -1.325292640508833E-004;
		zbs[VmecIndataNamelist.ntord + 7][5] = 2.726109026309649E-005;
		zbs[VmecIndataNamelist.ntord + 8][5] = -1.850542092438669E-005;
		//
		zbs[VmecIndataNamelist.ntord - 8][6] = 7.386763623731291E-007;
		zbs[VmecIndataNamelist.ntord - 7][6] = -2.508309149435569E-005;
		zbs[VmecIndataNamelist.ntord - 6][6] = 2.117983834955785E-005;
		zbs[VmecIndataNamelist.ntord - 5][6] = -5.753278329735706E-005;
		zbs[VmecIndataNamelist.ntord - 4][6] = -1.155249029943377E-005;
		zbs[VmecIndataNamelist.ntord - 3][6] = 1.153115758377326E-004;
		zbs[VmecIndataNamelist.ntord - 2][6] = -2.830337070611086E-005;
		zbs[VmecIndataNamelist.ntord - 1][6] = 1.324540090482543E-005;
		zbs[VmecIndataNamelist.ntord + 0][6] = -7.800374644189332E-004;
		zbs[VmecIndataNamelist.ntord + 1][6] = 3.838347842484833E-004;
		zbs[VmecIndataNamelist.ntord + 2][6] = -6.137053321851298E-004;
		zbs[VmecIndataNamelist.ntord + 3][6] = -1.301407081337896E-003;
		zbs[VmecIndataNamelist.ntord + 4][6] = -4.955356310734652E-004;
		zbs[VmecIndataNamelist.ntord + 5][6] = 3.462627709485832E-004;
		zbs[VmecIndataNamelist.ntord + 6][6] = -9.919741284083595E-005;
		zbs[VmecIndataNamelist.ntord + 7][6] = -1.784091935229791E-005;
		zbs[VmecIndataNamelist.ntord + 8][6] = -2.302543207412257E-005;
		//
		zbs[VmecIndataNamelist.ntord - 8][7] = -1.484141014465394E-005;
		zbs[VmecIndataNamelist.ntord - 7][7] = -1.816080256001114E-005;
		zbs[VmecIndataNamelist.ntord - 6][7] = -3.833112693126000E-006;
		zbs[VmecIndataNamelist.ntord - 5][7] = -1.865844947458311E-006;
		zbs[VmecIndataNamelist.ntord - 4][7] = -9.801486449695234E-007;
		zbs[VmecIndataNamelist.ntord - 3][7] = 8.061675676629167E-006;
		zbs[VmecIndataNamelist.ntord - 2][7] = -1.094928404347744E-006;
		zbs[VmecIndataNamelist.ntord - 1][7] = -3.601316979928921E-005;
		zbs[VmecIndataNamelist.ntord + 0][7] = 2.940524142344392E-004;
		zbs[VmecIndataNamelist.ntord + 1][7] = -8.496033274994724E-004;
		zbs[VmecIndataNamelist.ntord + 2][7] = 6.255838033712031E-004;
		zbs[VmecIndataNamelist.ntord + 3][7] = 3.900613621829464E-004;
		zbs[VmecIndataNamelist.ntord + 4][7] = 1.845971251755716E-003;
		zbs[VmecIndataNamelist.ntord + 5][7] = 2.917639239820137E-004;
		zbs[VmecIndataNamelist.ntord + 6][7] = -5.665554773695777E-005;
		zbs[VmecIndataNamelist.ntord + 7][7] = 1.038283071268155E-004;
		zbs[VmecIndataNamelist.ntord + 8][7] = -2.979906810937261E-006;
		//
		Assertions.assertArrayEquals(zbs, vmecInput.zbs);

		double[][] rbs = new double[2 * VmecIndataNamelist.ntord + 1][VmecIndataNamelist.mpold + 1];
		Assertions.assertArrayEquals(rbs, vmecInput.rbs);

		double[][] zbc = new double[2 * VmecIndataNamelist.ntord + 1][VmecIndataNamelist.mpold + 1];
		Assertions.assertArrayEquals(zbc, vmecInput.zbc);
	}

	// EXTCUR has spaces in indices
	@Test
	void vmecInputTest_BETA_5_ICUR_5K() throws IOException {
		String inputFilename = "/de/labathome/FortranNamelist/input/input.BETA_5_ICUR_5K";
		InputStream inputFileStream = getClass().getResourceAsStream(inputFilename);

		String inputFile = new String(inputFileStream.readAllBytes(), StandardCharsets.UTF_8);

		VmecIndataNamelist vmecInput = new VmecIndataNamelist();

		// parse the namelist into a Java object
		FortranNamelist parser = new FortranNamelist(inputFile, "indata", vmecInput);
		vmecInput = (VmecIndataNamelist) parser.getParsed();
		assertNotNull(vmecInput);
		vmecInput.sanitize();

		// --------------------
		// Numerical Resolution and Symmetry Assumptions
		Assertions.assertEquals(vmecInput.lasym, false);
		Assertions.assertEquals(vmecInput.nfp, 5);
		Assertions.assertEquals(vmecInput.mpol, 8);
		Assertions.assertEquals(vmecInput.ntor, 8);
		Assertions.assertEquals(vmecInput.ntheta, 0);
		Assertions.assertEquals(vmecInput.nzeta, 36);

		// --------------------
		// Multi-Grid Steps
		int[] ns_array = new int[100];
		ns_array[0] = 51;
		Assertions.assertArrayEquals(ns_array, vmecInput.ns_array);

		double[] ftol_array = new double[100];
		ftol_array[0] = 1.0e-19;
		Assertions.assertArrayEquals(ftol_array, vmecInput.ftol_array);

		int[] niter_array = new int[100];
		Arrays.fill(niter_array, 100000);
		Assertions.assertArrayEquals(niter_array, vmecInput.niter_array);

		// --------------------
		// Global Physics Parameters
		Assertions.assertEquals(-1.7104, vmecInput.phiedge);
		Assertions.assertEquals(1, vmecInput.ncurr);

		// --------------------
		// Profile of Mass or Pressure
		Assertions.assertEquals("two_lorentz", vmecInput.pmass_type.toLowerCase());
		double[] am = new double[21];
		am[0] = 1.0;
		am[1] = 1.0;
		am[2] = 1.0;
		am[3] = 1.0;
		am[4] = 1.0;
		am[5] = 1.0;
		am[6] = 1.0;
		am[7] = 1.0;
		Assertions.assertArrayEquals(am, vmecInput.am);
		Assertions.assertArrayEquals(new double[VmecIndataNamelist.nsd], vmecInput.am_aux_s);
		Assertions.assertArrayEquals(new double[VmecIndataNamelist.nsd], vmecInput.am_aux_f);
		Assertions.assertEquals(50.0e3, vmecInput.pres_scale);
		Assertions.assertEquals(0.0, vmecInput.gamma);
		Assertions.assertEquals(1.0, vmecInput.spres_ped);

		// --------------------
		// (Initial Guess for) Rotational Transform Profile
		Assertions.assertEquals("power_series", vmecInput.piota_type.toLowerCase());
		Assertions.assertArrayEquals(new double[21], vmecInput.ai);
		Assertions.assertArrayEquals(new double[VmecIndataNamelist.nsd], vmecInput.ai_aux_s);
		Assertions.assertArrayEquals(new double[VmecIndataNamelist.nsd], vmecInput.ai_aux_f);

		// --------------------
		// (Initial Guess for) Toroidal Current Profile
		Assertions.assertEquals("sum_atan", vmecInput.pcurr_type.toLowerCase());
		double[] ac = new double[21];
		ac[0] = 0.0;
		ac[1] = 1.0;
		ac[2] = 1.0;
		ac[3] = 1.5;
		ac[4] = 1.0;
		Assertions.assertArrayEquals(ac, vmecInput.ac);
		Assertions.assertArrayEquals(new double[VmecIndataNamelist.nsd], vmecInput.ac_aux_s);
		Assertions.assertArrayEquals(new double[VmecIndataNamelist.nsd], vmecInput.ac_aux_f);
		Assertions.assertEquals(5000.0, vmecInput.curtor);
		Assertions.assertEquals(1.0, vmecInput.bloat);

		// --------------------
		// Free-Boundary Parameters
		Assertions.assertEquals(true, vmecInput.lfreeb);
		Assertions.assertEquals("/u/jsch/v3fit_runs/MGRIDS/mgrid_w7x_v2_10mm_grid.nc", vmecInput.mgrid_file);
		double[] extcur = new double[VmecIndataNamelist.nigroup];
		extcur[0] = 1.28700e+04;
		extcur[1] = 1.31220e+04;
		extcur[2] = 1.38960e+04;
		extcur[3] = 1.19640e+04;
		extcur[4] = 1.08440e+04;
		extcur[5] = -1.0;
		extcur[6] = -1.0;
		Assertions.assertArrayEquals(extcur, vmecInput.extcur);
		Assertions.assertEquals(10, vmecInput.nvacskip);

		// --------------------
		// Tweaking Parameters
		Assertions.assertEquals(200, vmecInput.nstep);
		double[] aphi = new double[20];
		aphi[0] = 1.0;
		Assertions.assertArrayEquals(aphi, vmecInput.aphi);
		Assertions.assertEquals(1.0, vmecInput.delt);
		Assertions.assertEquals(1.0, vmecInput.tcon0);
		Assertions.assertEquals(false, vmecInput.lforbal);

		// --------------------
		// Initial Guess for Magnetic Axis Geometry
		double[] raxis_cc = new double[VmecIndataNamelist.ntord + 1];
		raxis_cc[0] = 5.5607E+00;
		Assertions.assertArrayEquals(raxis_cc, vmecInput.raxis_cc);

		double[] zaxis_cs = new double[VmecIndataNamelist.ntord + 1];
		Assertions.assertArrayEquals(zaxis_cs, vmecInput.zaxis_cs);

		double[] raxis_cs = new double[VmecIndataNamelist.ntord + 1];
		Assertions.assertArrayEquals(raxis_cs, vmecInput.raxis_cs);

		double[] zaxis_cc = new double[VmecIndataNamelist.ntord + 1];
		Assertions.assertArrayEquals(zaxis_cc, vmecInput.zaxis_cc);

		// --------------------
		// (Initial Guess for) Boundary Geometry

		double[][] rbc = new double[2 * VmecIndataNamelist.ntord + 1][VmecIndataNamelist.mpold + 1];
		rbc[VmecIndataNamelist.ntord + 0][0] = 5.5210E+00;
		rbc[VmecIndataNamelist.ntord + 1][0] = 2.7849e-01;
		rbc[VmecIndataNamelist.ntord + 0][1] = 4.8900E-01;
		Assertions.assertArrayEquals(rbc, vmecInput.rbc);

		double[][] zbs = new double[2 * VmecIndataNamelist.ntord + 1][VmecIndataNamelist.mpold + 1];
		zbs[VmecIndataNamelist.ntord + 0][0] = 0.0;
		zbs[VmecIndataNamelist.ntord + 1][0] = -2.3504e-01;
		zbs[VmecIndataNamelist.ntord + 0][1] = 6.2496E-01;
		Assertions.assertArrayEquals(zbs, vmecInput.zbs);

		double[][] rbs = new double[2 * VmecIndataNamelist.ntord + 1][VmecIndataNamelist.mpold + 1];
		Assertions.assertArrayEquals(rbs, vmecInput.rbs);

		double[][] zbc = new double[2 * VmecIndataNamelist.ntord + 1][VmecIndataNamelist.mpold + 1];
		Assertions.assertArrayEquals(zbc, vmecInput.zbc);
	}

	// => forgot value of am
	@Test
	void vmecInputTest_dboe_id_1000_1000_1000_1000_() throws IOException {
		String inputFilename = "/de/labathome/FortranNamelist/input/input.dboe_id_1000_1000_1000_1000_+0000_+0000_v_00_pres_00_it_6";
		InputStream inputFileStream = getClass().getResourceAsStream(inputFilename);

		String inputFile = new String(inputFileStream.readAllBytes(), StandardCharsets.UTF_8);

		VmecIndataNamelist vmecInput = new VmecIndataNamelist();

		// parse the namelist into a Java object
		FortranNamelist parser = new FortranNamelist(inputFile, "indata", vmecInput);
		vmecInput = (VmecIndataNamelist) parser.getParsed();
		assertNotNull(vmecInput);
		vmecInput.sanitize();

		// --------------------
		// Numerical Resolution and Symmetry Assumptions
		Assertions.assertEquals(vmecInput.lasym, false);
		Assertions.assertEquals(vmecInput.nfp, 5);
		Assertions.assertEquals(vmecInput.mpol, 12);
		Assertions.assertEquals(vmecInput.ntor, 12);
		Assertions.assertEquals(vmecInput.ntheta, 32);
		Assertions.assertEquals(vmecInput.nzeta, 36);

		// --------------------
		// Multi-Grid Steps
		int[] ns_array = new int[100];
		ns_array[0] = 4;
		ns_array[1] = 9;
		ns_array[2] = 28;
		ns_array[3] = 99;
		Assertions.assertArrayEquals(ns_array, vmecInput.ns_array);

		double[] ftol_array = new double[100];
		ftol_array[0] = 1.0e-3;
		ftol_array[1] = 1.0e-5;
		ftol_array[2] = 1.0e-9;
		ftol_array[3] = 1.0e-14;
		Assertions.assertArrayEquals(ftol_array, vmecInput.ftol_array);

		int[] niter_array = new int[100];
		Arrays.fill(niter_array, 100000);
		Assertions.assertArrayEquals(niter_array, vmecInput.niter_array);

		// --------------------
		// Global Physics Parameters
		Assertions.assertEquals(-1.8, vmecInput.phiedge);
		Assertions.assertEquals(1, vmecInput.ncurr);

		// --------------------
		// Profile of Mass or Pressure
		Assertions.assertEquals("power_series", vmecInput.pmass_type.toLowerCase());
		Assertions.assertArrayEquals(new double[21], vmecInput.am);
		Assertions.assertArrayEquals(new double[VmecIndataNamelist.nsd], vmecInput.am_aux_s);
		Assertions.assertArrayEquals(new double[VmecIndataNamelist.nsd], vmecInput.am_aux_f);
		Assertions.assertEquals(1.0, vmecInput.pres_scale);
		Assertions.assertEquals(0.0, vmecInput.gamma);
		Assertions.assertEquals(1.0, vmecInput.spres_ped);

		// --------------------
		// (Initial Guess for) Rotational Transform Profile
		Assertions.assertEquals("power_series", vmecInput.piota_type.toLowerCase());
		Assertions.assertArrayEquals(new double[21], vmecInput.ai);
		Assertions.assertArrayEquals(new double[VmecIndataNamelist.nsd], vmecInput.ai_aux_s);
		Assertions.assertArrayEquals(new double[VmecIndataNamelist.nsd], vmecInput.ai_aux_f);

		// --------------------
		// (Initial Guess for) Toroidal Current Profile
		Assertions.assertEquals("power_series", vmecInput.pcurr_type.toLowerCase());
		double[] ac = new double[21];
		ac[0] = 0.0;
		ac[1] = 1.0;
		ac[2] = -1.0;
		Assertions.assertArrayEquals(ac, vmecInput.ac);
		Assertions.assertArrayEquals(new double[VmecIndataNamelist.nsd], vmecInput.ac_aux_s);
		Assertions.assertArrayEquals(new double[VmecIndataNamelist.nsd], vmecInput.ac_aux_f);
		Assertions.assertEquals(44000.0, vmecInput.curtor);
		Assertions.assertEquals(1.0, vmecInput.bloat);

		// --------------------
		// Free-Boundary Parameters
		Assertions.assertEquals(true, vmecInput.lfreeb);
		Assertions.assertEquals("mgrid_w7x_nv36_hires.nc", vmecInput.mgrid_file);
		double[] extcur = new double[VmecIndataNamelist.nigroup];
		extcur[0] = 12882.480489739242;
		extcur[1] = 12882.480489739242;
		extcur[2] = 12882.480489739242;
		extcur[3] = 12882.480489739242;
		extcur[4] = 12882.480489739242;
		Assertions.assertArrayEquals(extcur, vmecInput.extcur);
		Assertions.assertEquals(6, vmecInput.nvacskip);

		// --------------------
		// Tweaking Parameters
		Assertions.assertEquals(100, vmecInput.nstep);
		double[] aphi = new double[20];
		aphi[0] = 1.0;
		Assertions.assertArrayEquals(aphi, vmecInput.aphi);
		Assertions.assertEquals(0.6, vmecInput.delt);
		Assertions.assertEquals(1.0, vmecInput.tcon0);
		Assertions.assertEquals(false, vmecInput.lforbal);

		// --------------------
		// Initial Guess for Magnetic Axis Geometry
		double[] raxis_cc = new double[VmecIndataNamelist.ntord + 1];
		raxis_cc[0] = 5.5607;
		raxis_cc[1] = 0.37075;
		Assertions.assertArrayEquals(raxis_cc, vmecInput.raxis_cc);

		double[] zaxis_cs = new double[VmecIndataNamelist.ntord + 1];
		zaxis_cs[1] = -0.30815;
		Assertions.assertArrayEquals(zaxis_cs, vmecInput.zaxis_cs);

		double[] raxis_cs = new double[VmecIndataNamelist.ntord + 1];
		Assertions.assertArrayEquals(raxis_cs, vmecInput.raxis_cs);

		double[] zaxis_cc = new double[VmecIndataNamelist.ntord + 1];
		Assertions.assertArrayEquals(zaxis_cc, vmecInput.zaxis_cc);

		// --------------------
		// (Initial Guess for) Boundary Geometry

		double[][] rbc = new double[2 * VmecIndataNamelist.ntord + 1][VmecIndataNamelist.mpold + 1];
		rbc[VmecIndataNamelist.ntord + 0][0] = 5.5208;
		rbc[VmecIndataNamelist.ntord + 1][0] = 0.27787;
		rbc[VmecIndataNamelist.ntord + 2][0] = -0.0068943;
		rbc[VmecIndataNamelist.ntord + 3][0] = -1.0741E-4;
		rbc[VmecIndataNamelist.ntord + 4][0] = -0.0014489;
		rbc[VmecIndataNamelist.ntord + 5][0] = -7.7056E-5;
		rbc[VmecIndataNamelist.ntord + 6][0] = -3.2743E-4;
		//
		rbc[VmecIndataNamelist.ntord - 6][1] = -5.7628E-4;
		rbc[VmecIndataNamelist.ntord - 5][1] = -1.3638E-4;
		rbc[VmecIndataNamelist.ntord - 4][1] = -0.0011712;
		rbc[VmecIndataNamelist.ntord - 3][1] = -2.206E-4;
		rbc[VmecIndataNamelist.ntord - 2][1] = 0.0023049;
		rbc[VmecIndataNamelist.ntord - 1][1] = 0.024621;
		rbc[VmecIndataNamelist.ntord + 0][1] = 0.48875;
		rbc[VmecIndataNamelist.ntord + 1][1] = -0.21307;
		rbc[VmecIndataNamelist.ntord + 2][1] = -0.018187;
		rbc[VmecIndataNamelist.ntord + 3][1] = 0.0022951;
		rbc[VmecIndataNamelist.ntord + 4][1] = 0.0017087;
		rbc[VmecIndataNamelist.ntord + 5][1] = 2.908E-4;
		rbc[VmecIndataNamelist.ntord + 6][1] = -6.7525E-5;
		//
		rbc[VmecIndataNamelist.ntord - 6][2] = 5.2205E-4;
		rbc[VmecIndataNamelist.ntord - 5][2] = -3.685E-4;
		rbc[VmecIndataNamelist.ntord - 4][2] = -4.4685E-5;
		rbc[VmecIndataNamelist.ntord - 3][2] = 1.3933E-4;
		rbc[VmecIndataNamelist.ntord - 2][2] = 0.0021791;
		rbc[VmecIndataNamelist.ntord - 1][2] = 0.011469;
		rbc[VmecIndataNamelist.ntord + 0][2] = 0.038021;
		rbc[VmecIndataNamelist.ntord + 1][2] = 0.044568;
		rbc[VmecIndataNamelist.ntord + 2][2] = 0.067691;
		rbc[VmecIndataNamelist.ntord + 3][2] = -0.0015827;
		rbc[VmecIndataNamelist.ntord + 4][2] = -7.0661E-4;
		rbc[VmecIndataNamelist.ntord + 5][2] = 1.1135E-4;
		rbc[VmecIndataNamelist.ntord + 6][2] = 1.8763E-4;
		//
		rbc[VmecIndataNamelist.ntord - 6][3] = -1.7461E-4;
		rbc[VmecIndataNamelist.ntord - 5][3] = 6.2895E-5;
		rbc[VmecIndataNamelist.ntord - 4][3] = 2.4326E-4;
		rbc[VmecIndataNamelist.ntord - 3][3] = -2.8097E-5;
		rbc[VmecIndataNamelist.ntord - 2][3] = -2.6826E-5;
		rbc[VmecIndataNamelist.ntord - 1][3] = 0.0015951;
		rbc[VmecIndataNamelist.ntord + 0][3] = -0.0027434;
		rbc[VmecIndataNamelist.ntord + 1][3] = -0.012242;
		rbc[VmecIndataNamelist.ntord + 2][3] = -0.02109;
		rbc[VmecIndataNamelist.ntord + 3][3] = -0.013199;
		rbc[VmecIndataNamelist.ntord + 4][3] = 0.0013953;
		rbc[VmecIndataNamelist.ntord + 5][3] = 7.1688E-5;
		rbc[VmecIndataNamelist.ntord + 6][3] = -5.6587E-6;
		//
		rbc[VmecIndataNamelist.ntord - 6][4] = 1.3336E-4;
		rbc[VmecIndataNamelist.ntord - 5][4] = -1.0498E-4;
		rbc[VmecIndataNamelist.ntord - 4][4] = -4.768E-6;
		rbc[VmecIndataNamelist.ntord - 3][4] = 1.8088E-4;
		rbc[VmecIndataNamelist.ntord - 2][4] = 7.3539E-5;
		rbc[VmecIndataNamelist.ntord - 1][4] = -4.8731E-4;
		rbc[VmecIndataNamelist.ntord + 0][4] = 0.0022625;
		rbc[VmecIndataNamelist.ntord + 1][4] = -0.0011282;
		rbc[VmecIndataNamelist.ntord + 2][4] = 0.0086662;
		rbc[VmecIndataNamelist.ntord + 3][4] = 0.0031656;
		rbc[VmecIndataNamelist.ntord + 4][4] = 4.3233E-4;
		rbc[VmecIndataNamelist.ntord + 5][4] = -6.7706E-4;
		rbc[VmecIndataNamelist.ntord + 6][4] = 3.381E-5;
		//
		rbc[VmecIndataNamelist.ntord - 6][5] = -1.0377E-4;
		rbc[VmecIndataNamelist.ntord - 5][5] = -9.0017E-5;
		rbc[VmecIndataNamelist.ntord - 4][5] = 1.6073E-4;
		rbc[VmecIndataNamelist.ntord - 3][5] = -6.5786E-5;
		rbc[VmecIndataNamelist.ntord - 2][5] = -8.2925E-5;
		rbc[VmecIndataNamelist.ntord - 1][5] = 6.1793E-4;
		rbc[VmecIndataNamelist.ntord + 0][5] = 4.5952E-4;
		rbc[VmecIndataNamelist.ntord + 1][5] = 0.0019919;
		rbc[VmecIndataNamelist.ntord + 2][5] = -5.325E-4;
		rbc[VmecIndataNamelist.ntord + 3][5] = -0.0011767;
		rbc[VmecIndataNamelist.ntord + 4][5] = 0.0010514;
		rbc[VmecIndataNamelist.ntord + 5][5] = 1.5892E-4;
		rbc[VmecIndataNamelist.ntord + 6][5] = 8.5567E-5;
		//
		rbc[VmecIndataNamelist.ntord - 6][6] = 4.7116E-5;
		rbc[VmecIndataNamelist.ntord - 5][6] = -2.8274E-5;
		rbc[VmecIndataNamelist.ntord - 4][6] = -2.3142E-5;
		rbc[VmecIndataNamelist.ntord - 3][6] = 6.8603E-5;
		rbc[VmecIndataNamelist.ntord - 2][6] = 1.0263E-4;
		rbc[VmecIndataNamelist.ntord - 1][6] = -4.8358E-5;
		rbc[VmecIndataNamelist.ntord + 0][6] = 1.4547E-4;
		rbc[VmecIndataNamelist.ntord + 1][6] = -0.0016075;
		rbc[VmecIndataNamelist.ntord + 2][6] = -0.0015746;
		rbc[VmecIndataNamelist.ntord + 3][6] = -0.0013242;
		rbc[VmecIndataNamelist.ntord + 4][6] = -7.94E-4;
		rbc[VmecIndataNamelist.ntord + 5][6] = -6.1283E-4;
		rbc[VmecIndataNamelist.ntord + 6][6] = 1.218E-4;
		//
		Assertions.assertArrayEquals(rbc, vmecInput.rbc);

		double[][] zbs = new double[2 * VmecIndataNamelist.ntord + 1][VmecIndataNamelist.mpold + 1];
		zbs[VmecIndataNamelist.ntord + 0][0] = 0.0;
		zbs[VmecIndataNamelist.ntord + 1][0] = -0.23539;
		zbs[VmecIndataNamelist.ntord + 2][0] = 0.0026848;
		zbs[VmecIndataNamelist.ntord + 3][0] = 0.0020733;
		zbs[VmecIndataNamelist.ntord + 4][0] = 0.0017992;
		zbs[VmecIndataNamelist.ntord + 5][0] = 2.1828E-4;
		zbs[VmecIndataNamelist.ntord + 6][0] = 1.5737E-4;
		//
		zbs[VmecIndataNamelist.ntord - 6][1] = 6.5157E-4;
		zbs[VmecIndataNamelist.ntord - 5][1] = -3.7205E-4;
		zbs[VmecIndataNamelist.ntord - 4][1] = -5.7712E-4;
		zbs[VmecIndataNamelist.ntord - 3][1] = 8.91E-4;
		zbs[VmecIndataNamelist.ntord - 2][1] = 0.0079723;
		zbs[VmecIndataNamelist.ntord - 1][1] = 0.022911;
		zbs[VmecIndataNamelist.ntord + 0][1] = 0.62516;
		zbs[VmecIndataNamelist.ntord + 1][1] = 0.20804;
		zbs[VmecIndataNamelist.ntord + 2][1] = 0.013018;
		zbs[VmecIndataNamelist.ntord + 3][1] = -0.0013049;
		zbs[VmecIndataNamelist.ntord + 4][1] = -1.5574E-4;
		zbs[VmecIndataNamelist.ntord + 5][1] = 3.0286E-4;
		zbs[VmecIndataNamelist.ntord + 6][1] = 2.7985E-4;
		//
		zbs[VmecIndataNamelist.ntord - 6][2] = -6.9642E-4;
		zbs[VmecIndataNamelist.ntord - 5][2] = 2.1661E-4;
		zbs[VmecIndataNamelist.ntord - 4][2] = 1.7501E-4;
		zbs[VmecIndataNamelist.ntord - 3][2] = -5.6433E-5;
		zbs[VmecIndataNamelist.ntord - 2][2] = 2.44E-5;
		zbs[VmecIndataNamelist.ntord - 1][2] = 0.0081389;
		zbs[VmecIndataNamelist.ntord + 0][2] = -0.0042271;
		zbs[VmecIndataNamelist.ntord + 1][2] = 0.019683;
		zbs[VmecIndataNamelist.ntord + 2][2] = -0.050482;
		zbs[VmecIndataNamelist.ntord + 3][2] = 0.0028309;
		zbs[VmecIndataNamelist.ntord + 4][2] = 0.0010155;
		zbs[VmecIndataNamelist.ntord + 5][2] = 2.3052E-4;
		zbs[VmecIndataNamelist.ntord + 6][2] = -5.9997E-5;
		//
		zbs[VmecIndataNamelist.ntord - 6][3] = 2.2405E-4;
		zbs[VmecIndataNamelist.ntord - 5][3] = -1.2172E-5;
		zbs[VmecIndataNamelist.ntord - 4][3] = -2.2395E-4;
		zbs[VmecIndataNamelist.ntord - 3][3] = -2.7206E-5;
		zbs[VmecIndataNamelist.ntord - 2][3] = 6.1237E-4;
		zbs[VmecIndataNamelist.ntord - 1][3] = -4.8551E-4;
		zbs[VmecIndataNamelist.ntord + 0][3] = -0.0014569;
		zbs[VmecIndataNamelist.ntord + 1][3] = -0.0045355;
		zbs[VmecIndataNamelist.ntord + 2][3] = 0.007419;
		zbs[VmecIndataNamelist.ntord + 3][3] = 0.011026;
		zbs[VmecIndataNamelist.ntord + 4][3] = -9.8338E-4;
		zbs[VmecIndataNamelist.ntord + 5][3] = -7.5541E-5;
		zbs[VmecIndataNamelist.ntord + 6][3] = 1.7458E-6;
		//
		zbs[VmecIndataNamelist.ntord - 6][4] = -3.9856E-5;
		zbs[VmecIndataNamelist.ntord - 5][4] = 2.0247E-5;
		zbs[VmecIndataNamelist.ntord - 4][4] = 7.576E-5;
		zbs[VmecIndataNamelist.ntord - 3][4] = -1.7144E-4;
		zbs[VmecIndataNamelist.ntord - 2][4] = 4.2436E-5;
		zbs[VmecIndataNamelist.ntord - 1][4] = 6.1602E-5;
		zbs[VmecIndataNamelist.ntord + 0][4] = 7.808E-4;
		zbs[VmecIndataNamelist.ntord + 1][4] = 0.0019232;
		zbs[VmecIndataNamelist.ntord + 2][4] = 0.0087014;
		zbs[VmecIndataNamelist.ntord + 3][4] = -0.0048804;
		zbs[VmecIndataNamelist.ntord + 4][4] = -0.0014794;
		zbs[VmecIndataNamelist.ntord + 5][4] = 5.5004E-4;
		zbs[VmecIndataNamelist.ntord + 6][4] = 2.5113E-5;
		//
		zbs[VmecIndataNamelist.ntord - 6][5] = -9.4615E-5;
		zbs[VmecIndataNamelist.ntord - 5][5] = 1.9564E-4;
		zbs[VmecIndataNamelist.ntord - 4][5] = -9.7698E-5;
		zbs[VmecIndataNamelist.ntord - 3][5] = 4.1842E-5;
		zbs[VmecIndataNamelist.ntord - 2][5] = -4.1543E-5;
		zbs[VmecIndataNamelist.ntord - 1][5] = 4.2101E-4;
		zbs[VmecIndataNamelist.ntord + 0][5] = 9.837E-4;
		zbs[VmecIndataNamelist.ntord + 1][5] = 1.1787E-4;
		zbs[VmecIndataNamelist.ntord + 2][5] = -0.002834;
		zbs[VmecIndataNamelist.ntord + 3][5] = -4.3649E-4;
		zbs[VmecIndataNamelist.ntord + 4][5] = 0.0011682;
		zbs[VmecIndataNamelist.ntord + 5][5] = 6.0077E-6;
		zbs[VmecIndataNamelist.ntord + 6][5] = -1.5061E-4;
		//
		zbs[VmecIndataNamelist.ntord - 6][6] = 1.933E-5;
		zbs[VmecIndataNamelist.ntord - 5][6] = -6.0921E-5;
		zbs[VmecIndataNamelist.ntord - 4][6] = 1.3194E-5;
		zbs[VmecIndataNamelist.ntord - 3][6] = 1.1157E-4;
		zbs[VmecIndataNamelist.ntord - 2][6] = -2.6581E-5;
		zbs[VmecIndataNamelist.ntord - 1][6] = 7.4858E-7;
		zbs[VmecIndataNamelist.ntord + 0][6] = -9.458E-6;
		zbs[VmecIndataNamelist.ntord + 1][6] = -8.6688E-4;
		zbs[VmecIndataNamelist.ntord + 2][6] = -0.0011037;
		zbs[VmecIndataNamelist.ntord + 3][6] = -9.5527E-4;
		zbs[VmecIndataNamelist.ntord + 4][6] = -5.5342E-4;
		zbs[VmecIndataNamelist.ntord + 5][6] = 9.9096E-5;
		zbs[VmecIndataNamelist.ntord + 6][6] = 3.3856E-5;
		//
		Assertions.assertArrayEquals(zbs, vmecInput.zbs);

		double[][] rbs = new double[2 * VmecIndataNamelist.ntord + 1][VmecIndataNamelist.mpold + 1];
		Assertions.assertArrayEquals(rbs, vmecInput.rbs);

		double[][] zbc = new double[2 * VmecIndataNamelist.ntord + 1][VmecIndataNamelist.mpold + 1];
		Assertions.assertArrayEquals(zbc, vmecInput.zbc);

	}


	@Test
	void vmecInputTest_demo() throws IOException {
		String inputFilename = "/de/labathome/FortranNamelist/input/input.demo";
		InputStream inputFileStream = getClass().getResourceAsStream(inputFilename);

		String inputFile = new String(inputFileStream.readAllBytes(), StandardCharsets.UTF_8);

		VmecIndataNamelist vmecInput = new VmecIndataNamelist();

		// parse the namelist into a Java object
		FortranNamelist parser = new FortranNamelist(inputFile, "indata", vmecInput);
		vmecInput = (VmecIndataNamelist) parser.getParsed();
		assertNotNull(vmecInput);
		vmecInput.sanitize();


		// --------------------
		// Numerical Resolution and Symmetry Assumptions
		Assertions.assertEquals(vmecInput.lasym, false);
		// nfp
		// mpol
		// ntor
		// ntheta
		// nzeta

		// --------------------
		// Multi-Grid Steps
		// ns_array
		// ftol_array
		// niter_array

		// --------------------
		// Global Physics Parameters
		// phiedge
		// ncurr

		// --------------------
		// Profile of Mass or Pressure
		// pmass_type
		// am
		// am_aux_s
		// am_aux_f
		// pres_scale
		// gamma
		// spres_ped

		// --------------------
		// (Initial Guess for) Rotational Transform Profile
		// piota_type
		// ai
		// ai_aux_s
		// ai_aux_f


		// --------------------
		// (Initial Guess for) Toroidal Current Profile
		// pcurr_type
		// ac
		// ac_aux_s
		// ac_aux_f
		// curtor
		// bloat




		// --------------------
		// Free-Boundary Parameters
		// lfreeb
		// mgrid_file
		// extcur
		// nvacskip



		// --------------------
		// Tweaking Parameters
		// nstep
		// aphi
		// delt
		// tcon0
		// lforbal


		// --------------------
		// Initial Guess for Magnetic Axis Geometry
		// raxis_cc
		// zaxis_cs
		// raxis_cs
		// zaxis_cc


		// --------------------
		// (Initial Guess for) Boundary Geometry
		// rbc
		// zbs
		// rbs
		// zbc

	}


	@Test
	void vmecInputTest_demo_vmec6_90() throws IOException {
		String inputFilename = "/de/labathome/FortranNamelist/input/input.demo_vmec6_90";
		InputStream inputFileStream = getClass().getResourceAsStream(inputFilename);

		String inputFile = new String(inputFileStream.readAllBytes(), StandardCharsets.UTF_8);

		VmecIndataNamelist vmecInput = new VmecIndataNamelist();

		// parse the namelist into a Java object
		FortranNamelist parser = new FortranNamelist(inputFile, "indata", vmecInput);
		vmecInput = (VmecIndataNamelist) parser.getParsed();
		assertNotNull(vmecInput);
		vmecInput.sanitize();


		// --------------------
		// Numerical Resolution and Symmetry Assumptions
		Assertions.assertEquals(vmecInput.lasym, false);
		// nfp
		// mpol
		// ntor
		// ntheta
		// nzeta

		// --------------------
		// Multi-Grid Steps
		// ns_array
		// ftol_array
		// niter_array

		// --------------------
		// Global Physics Parameters
		// phiedge
		// ncurr

		// --------------------
		// Profile of Mass or Pressure
		// pmass_type
		// am
		// am_aux_s
		// am_aux_f
		// pres_scale
		// gamma
		// spres_ped

		// --------------------
		// (Initial Guess for) Rotational Transform Profile
		// piota_type
		// ai
		// ai_aux_s
		// ai_aux_f


		// --------------------
		// (Initial Guess for) Toroidal Current Profile
		// pcurr_type
		// ac
		// ac_aux_s
		// ac_aux_f
		// curtor
		// bloat




		// --------------------
		// Free-Boundary Parameters
		// lfreeb
		// mgrid_file
		// extcur
		// nvacskip



		// --------------------
		// Tweaking Parameters
		// nstep
		// aphi
		// delt
		// tcon0
		// lforbal


		// --------------------
		// Initial Guess for Magnetic Axis Geometry
		// raxis_cc
		// zaxis_cs
		// raxis_cs
		// zaxis_cc


		// --------------------
		// (Initial Guess for) Boundary Geometry
		// rbc
		// zbs
		// rbs
		// zbc

	}

	@Test
	void vmecInputTest_demo3d_vmec6_90() throws IOException {
		String inputFilename = "/de/labathome/FortranNamelist/input/input.demo3d_vmec6_90";
		InputStream inputFileStream = getClass().getResourceAsStream(inputFilename);

		String inputFile = new String(inputFileStream.readAllBytes(), StandardCharsets.UTF_8);

		VmecIndataNamelist vmecInput = new VmecIndataNamelist();

		// parse the namelist into a Java object
		FortranNamelist parser = new FortranNamelist(inputFile, "indata", vmecInput);
		vmecInput = (VmecIndataNamelist) parser.getParsed();
		assertNotNull(vmecInput);
		vmecInput.sanitize();


		// --------------------
		// Numerical Resolution and Symmetry Assumptions
		Assertions.assertEquals(vmecInput.lasym, false);
		// nfp
		// mpol
		// ntor
		// ntheta
		// nzeta

		// --------------------
		// Multi-Grid Steps
		// ns_array
		// ftol_array
		// niter_array

		// --------------------
		// Global Physics Parameters
		// phiedge
		// ncurr

		// --------------------
		// Profile of Mass or Pressure
		// pmass_type
		// am
		// am_aux_s
		// am_aux_f
		// pres_scale
		// gamma
		// spres_ped

		// --------------------
		// (Initial Guess for) Rotational Transform Profile
		// piota_type
		// ai
		// ai_aux_s
		// ai_aux_f


		// --------------------
		// (Initial Guess for) Toroidal Current Profile
		// pcurr_type
		// ac
		// ac_aux_s
		// ac_aux_f
		// curtor
		// bloat




		// --------------------
		// Free-Boundary Parameters
		// lfreeb
		// mgrid_file
		// extcur
		// nvacskip



		// --------------------
		// Tweaking Parameters
		// nstep
		// aphi
		// delt
		// tcon0
		// lforbal


		// --------------------
		// Initial Guess for Magnetic Axis Geometry
		// raxis_cc
		// zaxis_cs
		// raxis_cs
		// zaxis_cc


		// --------------------
		// (Initial Guess for) Boundary Geometry
		// rbc
		// zbs
		// rbs
		// zbc

	}

	@Test
	void vmecInputTest_ffhr_d1_nflux_v1() throws IOException {
		String inputFilename = "/de/labathome/FortranNamelist/input/input.ffhr_d1_nflux_v1.1";
		InputStream inputFileStream = getClass().getResourceAsStream(inputFilename);

		String inputFile = new String(inputFileStream.readAllBytes(), StandardCharsets.UTF_8);

		VmecIndataNamelist vmecInput = new VmecIndataNamelist();

		// parse the namelist into a Java object
		FortranNamelist parser = new FortranNamelist(inputFile, "indata", vmecInput);
		vmecInput = (VmecIndataNamelist) parser.getParsed();
		assertNotNull(vmecInput);
		vmecInput.sanitize();


		// --------------------
		// Numerical Resolution and Symmetry Assumptions
		Assertions.assertEquals(vmecInput.lasym, false);
		// nfp
		// mpol
		// ntor
		// ntheta
		// nzeta

		// --------------------
		// Multi-Grid Steps
		// ns_array
		// ftol_array
		// niter_array

		// --------------------
		// Global Physics Parameters
		// phiedge
		// ncurr

		// --------------------
		// Profile of Mass or Pressure
		// pmass_type
		// am
		// am_aux_s
		// am_aux_f
		// pres_scale
		// gamma
		// spres_ped

		// --------------------
		// (Initial Guess for) Rotational Transform Profile
		// piota_type
		// ai
		// ai_aux_s
		// ai_aux_f


		// --------------------
		// (Initial Guess for) Toroidal Current Profile
		// pcurr_type
		// ac
		// ac_aux_s
		// ac_aux_f
		// curtor
		// bloat




		// --------------------
		// Free-Boundary Parameters
		// lfreeb
		// mgrid_file
		// extcur
		// nvacskip



		// --------------------
		// Tweaking Parameters
		// nstep
		// aphi
		// delt
		// tcon0
		// lforbal


		// --------------------
		// Initial Guess for Magnetic Axis Geometry
		// raxis_cc
		// zaxis_cs
		// raxis_cs
		// zaxis_cc


		// --------------------
		// (Initial Guess for) Boundary Geometry
		// rbc
		// zbs
		// rbs
		// zbc

	}


	@Test
	void vmecInputTest_ITER_nflux_verification_0001() throws IOException {
		String inputFilename = "/de/labathome/FortranNamelist/input/input.ITER_nflux_verification_0001.0001";
		InputStream inputFileStream = getClass().getResourceAsStream(inputFilename);

		String inputFile = new String(inputFileStream.readAllBytes(), StandardCharsets.UTF_8);

		VmecIndataNamelist vmecInput = new VmecIndataNamelist();

		// parse the namelist into a Java object
		FortranNamelist parser = new FortranNamelist(inputFile, "indata", vmecInput);
		vmecInput = (VmecIndataNamelist) parser.getParsed();
		assertNotNull(vmecInput);
		vmecInput.sanitize();


		// --------------------
		// Numerical Resolution and Symmetry Assumptions
		Assertions.assertEquals(vmecInput.lasym, false);
		// nfp
		// mpol
		// ntor
		// ntheta
		// nzeta

		// --------------------
		// Multi-Grid Steps
		// ns_array
		// ftol_array
		// niter_array

		// --------------------
		// Global Physics Parameters
		// phiedge
		// ncurr

		// --------------------
		// Profile of Mass or Pressure
		// pmass_type
		// am
		// am_aux_s
		// am_aux_f
		// pres_scale
		// gamma
		// spres_ped

		// --------------------
		// (Initial Guess for) Rotational Transform Profile
		// piota_type
		// ai
		// ai_aux_s
		// ai_aux_f


		// --------------------
		// (Initial Guess for) Toroidal Current Profile
		// pcurr_type
		// ac
		// ac_aux_s
		// ac_aux_f
		// curtor
		// bloat




		// --------------------
		// Free-Boundary Parameters
		// lfreeb
		// mgrid_file
		// extcur
		// nvacskip



		// --------------------
		// Tweaking Parameters
		// nstep
		// aphi
		// delt
		// tcon0
		// lforbal


		// --------------------
		// Initial Guess for Magnetic Axis Geometry
		// raxis_cc
		// zaxis_cs
		// raxis_cs
		// zaxis_cc


		// --------------------
		// (Initial Guess for) Boundary Geometry
		// rbc
		// zbs
		// rbs
		// zbc

	}

	@Test
	void vmecInputTest_JDHtest7() throws IOException {
		String inputFilename = "/de/labathome/FortranNamelist/input/input.JDHtest7";
		InputStream inputFileStream = getClass().getResourceAsStream(inputFilename);

		String inputFile = new String(inputFileStream.readAllBytes(), StandardCharsets.UTF_8);

		VmecIndataNamelist vmecInput = new VmecIndataNamelist();

		// parse the namelist into a Java object
		FortranNamelist parser = new FortranNamelist(inputFile, "indata", vmecInput);
		vmecInput = (VmecIndataNamelist) parser.getParsed();
		assertNotNull(vmecInput);
		vmecInput.sanitize();


		// --------------------
		// Numerical Resolution and Symmetry Assumptions
		Assertions.assertEquals(vmecInput.lasym, false);
		// nfp
		// mpol
		// ntor
		// ntheta
		// nzeta

		// --------------------
		// Multi-Grid Steps
		// ns_array
		// ftol_array
		// niter_array

		// --------------------
		// Global Physics Parameters
		// phiedge
		// ncurr

		// --------------------
		// Profile of Mass or Pressure
		// pmass_type
		// am
		// am_aux_s
		// am_aux_f
		// pres_scale
		// gamma
		// spres_ped

		// --------------------
		// (Initial Guess for) Rotational Transform Profile
		// piota_type
		// ai
		// ai_aux_s
		// ai_aux_f


		// --------------------
		// (Initial Guess for) Toroidal Current Profile
		// pcurr_type
		// ac
		// ac_aux_s
		// ac_aux_f
		// curtor
		// bloat




		// --------------------
		// Free-Boundary Parameters
		// lfreeb
		// mgrid_file
		// extcur
		// nvacskip



		// --------------------
		// Tweaking Parameters
		// nstep
		// aphi
		// delt
		// tcon0
		// lforbal


		// --------------------
		// Initial Guess for Magnetic Axis Geometry
		// raxis_cc
		// zaxis_cs
		// raxis_cs
		// zaxis_cc


		// --------------------
		// (Initial Guess for) Boundary Geometry
		// rbc
		// zbs
		// rbs
		// zbc

	}

	@Test
	void vmecInputTest_minerva() throws IOException {
		String inputFilename = "/de/labathome/FortranNamelist/input/input.minerva-vmec-9aea06f0c9b27276d9a55c2817eda280";
		InputStream inputFileStream = getClass().getResourceAsStream(inputFilename);

		String inputFile = new String(inputFileStream.readAllBytes(), StandardCharsets.UTF_8);

		VmecIndataNamelist vmecInput = new VmecIndataNamelist();

		// parse the namelist into a Java object
		FortranNamelist parser = new FortranNamelist(inputFile, "indata", vmecInput);
		vmecInput = (VmecIndataNamelist) parser.getParsed();
		assertNotNull(vmecInput);
		vmecInput.sanitize();


		// --------------------
		// Numerical Resolution and Symmetry Assumptions
		Assertions.assertEquals(vmecInput.lasym, false);
		// nfp
		// mpol
		// ntor
		// ntheta
		// nzeta

		// --------------------
		// Multi-Grid Steps
		// ns_array
		// ftol_array
		// niter_array

		// --------------------
		// Global Physics Parameters
		// phiedge
		// ncurr

		// --------------------
		// Profile of Mass or Pressure
		// pmass_type
		// am
		// am_aux_s
		// am_aux_f
		// pres_scale
		// gamma
		// spres_ped

		// --------------------
		// (Initial Guess for) Rotational Transform Profile
		// piota_type
		// ai
		// ai_aux_s
		// ai_aux_f


		// --------------------
		// (Initial Guess for) Toroidal Current Profile
		// pcurr_type
		// ac
		// ac_aux_s
		// ac_aux_f
		// curtor
		// bloat




		// --------------------
		// Free-Boundary Parameters
		// lfreeb
		// mgrid_file
		// extcur
		// nvacskip



		// --------------------
		// Tweaking Parameters
		// nstep
		// aphi
		// delt
		// tcon0
		// lforbal


		// --------------------
		// Initial Guess for Magnetic Axis Geometry
		// raxis_cc
		// zaxis_cs
		// raxis_cs
		// zaxis_cc


		// --------------------
		// (Initial Guess for) Boundary Geometry
		// rbc
		// zbs
		// rbs
		// zbc

	}


	@Test
	void vmecInputTest_vmec_asym() throws IOException {
		String inputFilename = "/de/labathome/FortranNamelist/input/input.test.vmec_asym";
		InputStream inputFileStream = getClass().getResourceAsStream(inputFilename);

		String inputFile = new String(inputFileStream.readAllBytes(), StandardCharsets.UTF_8);

		VmecIndataNamelist vmecInput = new VmecIndataNamelist();

		// parse the namelist into a Java object
		FortranNamelist parser = new FortranNamelist(inputFile, "indata", vmecInput);
		vmecInput = (VmecIndataNamelist) parser.getParsed();
		assertNotNull(vmecInput);
		vmecInput.sanitize();


		// --------------------
		// Numerical Resolution and Symmetry Assumptions
		Assertions.assertEquals(vmecInput.lasym, true);
		// nfp
		// mpol
		// ntor
		// ntheta
		// nzeta

		// --------------------
		// Multi-Grid Steps
		// ns_array
		// ftol_array
		// niter_array

		// --------------------
		// Global Physics Parameters
		// phiedge
		// ncurr

		// --------------------
		// Profile of Mass or Pressure
		// pmass_type
		// am
		// am_aux_s
		// am_aux_f
		// pres_scale
		// gamma
		// spres_ped

		// --------------------
		// (Initial Guess for) Rotational Transform Profile
		// piota_type
		// ai
		// ai_aux_s
		// ai_aux_f


		// --------------------
		// (Initial Guess for) Toroidal Current Profile
		// pcurr_type
		// ac
		// ac_aux_s
		// ac_aux_f
		// curtor
		// bloat




		// --------------------
		// Free-Boundary Parameters
		// lfreeb
		// mgrid_file
		// extcur
		// nvacskip



		// --------------------
		// Tweaking Parameters
		// nstep
		// aphi
		// delt
		// tcon0
		// lforbal


		// --------------------
		// Initial Guess for Magnetic Axis Geometry
		// raxis_cc
		// zaxis_cs
		// raxis_cs
		// zaxis_cc


		// --------------------
		// (Initial Guess for) Boundary Geometry
		// rbc
		// zbs
		// rbs
		// zbc

	}





	@Test
	void vmecInputTest_v3fit_jsch_20171207_006_2200ms_pm_10ms_v0_Bprecheck() throws IOException {
		String inputFilename = "/de/labathome/FortranNamelist/input/input.v3fit_jsch_20171207_006_2200ms_pm_10ms_v0_Bprecheck";
		InputStream inputFileStream = getClass().getResourceAsStream(inputFilename);

		String inputFile = new String(inputFileStream.readAllBytes(), StandardCharsets.UTF_8);

		VmecIndataNamelist vmecInput = new VmecIndataNamelist();

		// parse the namelist into a Java object
		FortranNamelist parser = new FortranNamelist(inputFile, "indata", vmecInput);
		vmecInput = (VmecIndataNamelist) parser.getParsed();
		assertNotNull(vmecInput);
		vmecInput.sanitize();

		// Numerical Resolution and Symmetry Assumptions
		Assertions.assertEquals(false, vmecInput.lasym);
		Assertions.assertEquals(5, vmecInput.nfp);
		Assertions.assertEquals(8, vmecInput.mpol);
		Assertions.assertEquals(8, vmecInput.ntor);
		Assertions.assertEquals(22, vmecInput.ntheta);
		Assertions.assertEquals(36, vmecInput.nzeta);

		// Multi-Grid Steps
		int[] ns_array = new int[100];
		ns_array[0] = 51;
		Assertions.assertArrayEquals(ns_array, vmecInput.ns_array);

		double[] ftol_array = new double[100];
		ftol_array[0] = 1.0e-19;
		Assertions.assertArrayEquals(ftol_array, vmecInput.ftol_array);

		int[] niter_array = new int[100];
		Arrays.fill(niter_array, 100000);
		Assertions.assertArrayEquals(niter_array, vmecInput.niter_array);

		// Global Physics Parameters
		Assertions.assertEquals(1.86086390731511, vmecInput.phiedge);
		Assertions.assertEquals(1, vmecInput.ncurr);

		// Profile of Mass or Pressure
		Assertions.assertEquals("two_power", vmecInput.pmass_type.trim());

		double[] am = new double[21];
		am[0] = 66978.0973204654;
		am[1] = 1.90820658424883;
		am[2] = 4.46347452430851;
		Assertions.assertArrayEquals(am, vmecInput.am);

		Assertions.assertArrayEquals(new double[VmecIndataNamelist.nsd], vmecInput.am_aux_s);
		Assertions.assertArrayEquals(new double[VmecIndataNamelist.nsd], vmecInput.am_aux_f);

		Assertions.assertEquals(1.0, vmecInput.pres_scale);
		Assertions.assertEquals(0.0, vmecInput.gamma);
		Assertions.assertEquals(1.0, vmecInput.spres_ped);

		// (Initial Guess for) Rotational Transform Profile
		Assertions.assertEquals("power_series", vmecInput.piota_type.trim());

		Assertions.assertArrayEquals(new double[21], vmecInput.ai);

		Assertions.assertArrayEquals(new double[VmecIndataNamelist.nsd], vmecInput.ai_aux_s);
		Assertions.assertArrayEquals(new double[VmecIndataNamelist.nsd], vmecInput.ai_aux_f);

		// (Initial Guess for) Toroidal Current Profile
		Assertions.assertEquals("sum_atan", vmecInput.pcurr_type.trim());

		double[] ac = new double[21];
		ac[1] = 1.0;
		ac[2] = 4.65367221956161;
		ac[3] = 2.27775992394378;
		ac[4] = 1.0;
		ac[5] = -0.943974780239535;
		ac[6] = 4.97719715801010;
		ac[7] = 1.0;
		ac[8] = 1.0;
		Assertions.assertArrayEquals(ac, vmecInput.ac);

		Assertions.assertArrayEquals(new double[VmecIndataNamelist.nsd], vmecInput.ac_aux_s);
		Assertions.assertArrayEquals(new double[VmecIndataNamelist.nsd], vmecInput.ac_aux_f);

		Assertions.assertEquals(-913.791680825795, vmecInput.curtor);
		Assertions.assertEquals(1.0, vmecInput.bloat);

		// Free-Boundary Parameters
		Assertions.assertEquals(true, vmecInput.lfreeb);
		Assertions.assertEquals("/u/jsch/v3fit_runs/MGRIDS/mgrid_w7x_v2_10mm_grid.nc", vmecInput.mgrid_file.trim());

		double[] extcur = new double[VmecIndataNamelist.nigroup];
		extcur[0] = -12989.0;
		extcur[1] = -12989.0;
		extcur[2] = -12989.0;
		extcur[3] = -12989.0;
		extcur[4] = -12989.0;
		Assertions.assertArrayEquals(extcur, vmecInput.extcur);

		Assertions.assertEquals(10, vmecInput.nvacskip);

		// Tweaking Parameters
		Assertions.assertEquals(200, vmecInput.nstep);

		double[] aphi = new double[20];
		aphi[0] = 1.0;
		Assertions.assertArrayEquals(aphi, vmecInput.aphi);

		Assertions.assertEquals(0.9, vmecInput.delt);
		Assertions.assertEquals(1.0, vmecInput.tcon0);
		Assertions.assertEquals(false, vmecInput.lforbal);

		// Initial Guess for Magnetic Axis Geometry
		double[] raxis_cc = new double[VmecIndataNamelist.ntord + 1];
		raxis_cc[0] = 5.60519782116720;
		raxis_cc[1] = 0.359696693960237;
		raxis_cc[2] = 1.335514497262088E-002;
		raxis_cc[3] = 7.259351998922190E-004;
		raxis_cc[4] = -6.817192552674108E-005;
		raxis_cc[5] = -4.866604255356447E-005;
		raxis_cc[6] = -9.461962764212337E-005;
		raxis_cc[7] = -2.828730415733434E-005;
		raxis_cc[8] = -1.554463199957404E-004;
		Assertions.assertArrayEquals(raxis_cc, vmecInput.raxis_cc);

		double[] zaxis_cs = new double[VmecIndataNamelist.ntord + 1];
		zaxis_cs[1] = -0.301269608014063;
		zaxis_cs[2] = -1.618262861498088E-002;
		zaxis_cs[3] = 1.319614987094394E-004;
		zaxis_cs[4] = 7.237442086550272E-005;
		zaxis_cs[5] = -7.887293984398275E-005;
		zaxis_cs[6] = -5.489994367001566E-006;
		zaxis_cs[7] = 2.463216529452668E-005;
		zaxis_cs[8] = -1.340415732215003E-004;
		Assertions.assertArrayEquals(zaxis_cs, vmecInput.zaxis_cs);

		double[] raxis_cs = new double[VmecIndataNamelist.ntord + 1];
		raxis_cs[1] = -2.093456611578366E-016;
		raxis_cs[2] = -4.361367940788264E-016;
		raxis_cs[3] = -2.267911329209897E-015;
		raxis_cs[4] = -1.500310571631163E-015;
		raxis_cs[5] = -1.395637741052244E-016;
		raxis_cs[6] = 1.290964910473326E-015;
		raxis_cs[7] = -1.116510192841796E-015;
		raxis_cs[8] = 2.233020385683591E-015;
		Assertions.assertArrayEquals(raxis_cs, vmecInput.raxis_cs);

		double[] zaxis_cc = new double[VmecIndataNamelist.ntord + 1];
		zaxis_cc[1] = 3.679904200040097E-017;
		zaxis_cc[2] = -4.906538933386797E-018;
		zaxis_cc[3] = 6.814637407481662E-018;
		zaxis_cc[4] = -2.807630611882445E-017;
		zaxis_cc[5] = -6.187690765993350E-017;
		zaxis_cc[6] = 8.450150385277261E-018;
		zaxis_cc[7] = -5.642519773394817E-017;
		zaxis_cc[8] = 4.497660688937897E-018;
		Assertions.assertArrayEquals(zaxis_cc, vmecInput.zaxis_cc);

		// (Initial Guess for) Boundary Geometry

		double[][] rbc = new double[2 * VmecIndataNamelist.ntord + 1][VmecIndataNamelist.mpold + 1];
		rbc[VmecIndataNamelist.ntord + 0][0] = 5.53767206572562;
		rbc[VmecIndataNamelist.ntord + 1][0] = 0.274510550203623;
		rbc[VmecIndataNamelist.ntord + 2][0] = -6.981826227607880E-003;
		rbc[VmecIndataNamelist.ntord + 3][0] = 1.813036409387450E-003;
		rbc[VmecIndataNamelist.ntord + 4][0] = -1.816182024820829E-003;
		rbc[VmecIndataNamelist.ntord + 5][0] = -9.457063618548155E-007;
		rbc[VmecIndataNamelist.ntord + 6][0] = -3.598476797219176E-004;
		rbc[VmecIndataNamelist.ntord + 7][0] = 5.201212576617179E-005;
		rbc[VmecIndataNamelist.ntord + 8][0] = -3.252552333940336E-004;
		//
		rbc[VmecIndataNamelist.ntord - 8][1] = 1.896612909948546E-004;
		rbc[VmecIndataNamelist.ntord - 7][1] = -2.456907136966591E-004;
		rbc[VmecIndataNamelist.ntord - 6][1] = 2.171005552378400E-004;
		rbc[VmecIndataNamelist.ntord - 5][1] = -3.651594436424253E-004;
		rbc[VmecIndataNamelist.ntord - 4][1] = -5.338074727927220E-004;
		rbc[VmecIndataNamelist.ntord - 3][1] = 2.506984615739400E-004;
		rbc[VmecIndataNamelist.ntord - 2][1] = -8.846005854412054E-004;
		rbc[VmecIndataNamelist.ntord - 1][1] = 1.678418664958850E-002;
		rbc[VmecIndataNamelist.ntord + 0][1] = 0.495662386868802;
		rbc[VmecIndataNamelist.ntord + 1][1] = -0.206910682038032;
		rbc[VmecIndataNamelist.ntord + 2][1] = -1.680826411787574E-002;
		rbc[VmecIndataNamelist.ntord + 3][1] = 7.103573636538324E-004;
		rbc[VmecIndataNamelist.ntord + 4][1] = 4.245820612023922E-004;
		rbc[VmecIndataNamelist.ntord + 5][1] = 6.969490097295169E-004;
		rbc[VmecIndataNamelist.ntord + 6][1] = -7.039468898728505E-004;
		rbc[VmecIndataNamelist.ntord + 7][1] = 5.333996925958589E-004;
		rbc[VmecIndataNamelist.ntord + 8][1] = 6.689061049649439E-005;
		//
		rbc[VmecIndataNamelist.ntord - 8][2] = 1.151620243875582E-004;
		rbc[VmecIndataNamelist.ntord - 7][2] = -3.184094551899957E-004;
		rbc[VmecIndataNamelist.ntord - 6][2] = 4.300662767725242E-004;
		rbc[VmecIndataNamelist.ntord - 5][2] = -3.210600174477836E-004;
		rbc[VmecIndataNamelist.ntord - 4][2] = 2.624531379976409E-004;
		rbc[VmecIndataNamelist.ntord - 3][2] = 5.677789104917515E-005;
		rbc[VmecIndataNamelist.ntord - 2][2] = 2.317403085550396E-003;
		rbc[VmecIndataNamelist.ntord - 1][2] = 9.514153491031787E-003;
		rbc[VmecIndataNamelist.ntord + 0][2] = 3.302699802057374E-002;
		rbc[VmecIndataNamelist.ntord + 1][2] = 4.021410132068875E-002;
		rbc[VmecIndataNamelist.ntord + 2][2] = 6.948281839736793E-002;
		rbc[VmecIndataNamelist.ntord + 3][2] = -8.442177055766123E-004;
		rbc[VmecIndataNamelist.ntord + 4][2] = -1.591791378279287E-004;
		rbc[VmecIndataNamelist.ntord + 5][2] = 2.486332537387373E-005;
		rbc[VmecIndataNamelist.ntord + 6][2] = -7.528422072335360E-005;
		rbc[VmecIndataNamelist.ntord + 7][2] = -8.021064997689406E-005;
		rbc[VmecIndataNamelist.ntord + 8][2] = -5.237201084004603E-005;
		//
		rbc[VmecIndataNamelist.ntord - 8][3] = 1.714830561905995E-004;
		rbc[VmecIndataNamelist.ntord - 7][3] = 7.083574967044192E-005;
		rbc[VmecIndataNamelist.ntord - 6][3] = -1.030190624061461E-004;
		rbc[VmecIndataNamelist.ntord - 5][3] = 7.464723877566652E-005;
		rbc[VmecIndataNamelist.ntord - 4][3] = 1.594997910406184E-004;
		rbc[VmecIndataNamelist.ntord - 3][3] = -1.315594838027009E-004;
		rbc[VmecIndataNamelist.ntord - 2][3] = 1.896710074727297E-004;
		rbc[VmecIndataNamelist.ntord - 1][3] = -4.186302609439121E-004;
		rbc[VmecIndataNamelist.ntord + 0][3] = 3.311501954692046E-004;
		rbc[VmecIndataNamelist.ntord + 1][3] = -8.268632528457621E-003;
		rbc[VmecIndataNamelist.ntord + 2][3] = -2.035417064634164E-002;
		rbc[VmecIndataNamelist.ntord + 3][3] = -1.448421989297009E-002;
		rbc[VmecIndataNamelist.ntord + 4][3] = 7.509697962325974E-004;
		rbc[VmecIndataNamelist.ntord + 5][3] = -1.328389193189970E-004;
		rbc[VmecIndataNamelist.ntord + 6][3] = 1.474276373936834E-004;
		rbc[VmecIndataNamelist.ntord + 7][3] = 1.042690595884966E-004;
		rbc[VmecIndataNamelist.ntord + 8][3] = 6.708871806270065E-005;
		//
		rbc[VmecIndataNamelist.ntord - 8][4] = -1.440375584040593E-004;
		rbc[VmecIndataNamelist.ntord - 7][4] = -1.704752360936465E-004;
		rbc[VmecIndataNamelist.ntord - 6][4] = 9.011920036675072E-005;
		rbc[VmecIndataNamelist.ntord - 5][4] = -1.499685240577500E-004;
		rbc[VmecIndataNamelist.ntord - 4][4] = 5.756550703305372E-005;
		rbc[VmecIndataNamelist.ntord - 3][4] = 1.900144565741371E-004;
		rbc[VmecIndataNamelist.ntord - 2][4] = 4.923760064924362E-005;
		rbc[VmecIndataNamelist.ntord - 1][4] = -5.103652305253303E-004;
		rbc[VmecIndataNamelist.ntord + 0][4] = 2.486409772891553E-003;
		rbc[VmecIndataNamelist.ntord + 1][4] = 3.763697963319822E-003;
		rbc[VmecIndataNamelist.ntord + 2][4] = 9.220272047894581E-003;
		rbc[VmecIndataNamelist.ntord + 3][4] = 4.017321543601945E-003;
		rbc[VmecIndataNamelist.ntord + 4][4] = 9.476230338947471E-004;
		rbc[VmecIndataNamelist.ntord + 5][4] = -7.056521343060718E-004;
		rbc[VmecIndataNamelist.ntord + 6][4] = -6.013036923002932E-005;
		rbc[VmecIndataNamelist.ntord + 7][4] = -8.827308310929046E-005;
		rbc[VmecIndataNamelist.ntord + 8][4] = -7.602682766245118E-005;
		//
		rbc[VmecIndataNamelist.ntord - 8][5] = 2.153814980670956E-005;
		rbc[VmecIndataNamelist.ntord - 7][5] = 1.183952494315462E-005;
		rbc[VmecIndataNamelist.ntord - 6][5] = -7.154739133251302E-005;
		rbc[VmecIndataNamelist.ntord - 5][5] = -1.003063968491245E-004;
		rbc[VmecIndataNamelist.ntord - 4][5] =  1.646390554518521E-004;
		rbc[VmecIndataNamelist.ntord - 3][5] = -5.492928459569241E-005;
		rbc[VmecIndataNamelist.ntord - 2][5] = -1.008690713410809E-004;
		rbc[VmecIndataNamelist.ntord - 1][5] = 5.052245395750897E-004;
		rbc[VmecIndataNamelist.ntord + 0][5] = 6.070892885457450E-004;
		rbc[VmecIndataNamelist.ntord + 1][5] = 7.467412314074342E-004;
		rbc[VmecIndataNamelist.ntord + 2][5] = -1.665044640784653E-003;
		rbc[VmecIndataNamelist.ntord + 3][5] = -1.332894131049668E-003;
		rbc[VmecIndataNamelist.ntord + 4][5] = 8.621325208882533E-004;
		rbc[VmecIndataNamelist.ntord + 5][5] = 2.267011310380025E-004;
		rbc[VmecIndataNamelist.ntord + 6][5] = 1.753152728586062E-004;
		rbc[VmecIndataNamelist.ntord + 7][5] = 2.258512907309003E-005;
		rbc[VmecIndataNamelist.ntord + 8][5] = 4.622663427651913E-005;
		//
		rbc[VmecIndataNamelist.ntord - 8][6] = 2.261179964816800E-005;
		rbc[VmecIndataNamelist.ntord - 7][6] = -1.811117336074824E-005;
		rbc[VmecIndataNamelist.ntord - 6][6] = 7.536505110908006E-005;
		rbc[VmecIndataNamelist.ntord - 5][6] = -3.456196403903005E-006;
		rbc[VmecIndataNamelist.ntord - 4][6] = -5.139872953188579E-005;
		rbc[VmecIndataNamelist.ntord - 3][6] = 3.979845546091016E-005;
		rbc[VmecIndataNamelist.ntord - 2][6] = 1.243950000533798E-004;
		rbc[VmecIndataNamelist.ntord - 1][6] = -1.106699933734676E-004;
		rbc[VmecIndataNamelist.ntord + 0][6] = -1.223947247241717E-004;
		rbc[VmecIndataNamelist.ntord + 1][6] = -7.731238353022811E-005;
		rbc[VmecIndataNamelist.ntord + 2][6] = -1.569941375256721E-003;
		rbc[VmecIndataNamelist.ntord + 3][6] = -1.222191707450016E-003;
		rbc[VmecIndataNamelist.ntord + 4][6] = -1.140339566197411E-003;
		rbc[VmecIndataNamelist.ntord + 5][6] = -9.518279470355112E-004;
		rbc[VmecIndataNamelist.ntord + 6][6] = 2.558431811535456E-004;
		rbc[VmecIndataNamelist.ntord + 7][6] = 2.595101629939122E-005;
		rbc[VmecIndataNamelist.ntord + 8][6] = 1.620499278874425E-005;
		//
		rbc[VmecIndataNamelist.ntord - 8][7] = -1.398163883777544E-005;
		rbc[VmecIndataNamelist.ntord - 7][7] = -3.923123999457625E-005;
		rbc[VmecIndataNamelist.ntord - 6][7] = -2.306491449291393E-005;
		rbc[VmecIndataNamelist.ntord - 5][7] = -1.158175514942671E-004;
		rbc[VmecIndataNamelist.ntord - 4][7] = 9.504349582371491E-005;
		rbc[VmecIndataNamelist.ntord - 3][7] = -7.789561325605551E-005;
		rbc[VmecIndataNamelist.ntord - 2][7] = 3.045433428931733E-005;
		rbc[VmecIndataNamelist.ntord - 1][7] = -1.527742119656307E-005;
		rbc[VmecIndataNamelist.ntord + 0][7] = -1.998021858917175E-005;
		rbc[VmecIndataNamelist.ntord + 1][7] = 4.953044093333177E-005;
		rbc[VmecIndataNamelist.ntord + 2][7] = 4.714962252031772E-004;
		rbc[VmecIndataNamelist.ntord + 3][7] = 5.784116518161847E-004;
		rbc[VmecIndataNamelist.ntord + 4][7] = 7.049240981285555E-004;
		rbc[VmecIndataNamelist.ntord + 5][7] = -1.675739529392045E-004;
		rbc[VmecIndataNamelist.ntord + 6][7] = -5.586250774777524E-005;
		rbc[VmecIndataNamelist.ntord + 7][7] = -2.029459227327149E-004;
		rbc[VmecIndataNamelist.ntord + 8][7] = 5.906682990243203E-006;
		//
		Assertions.assertArrayEquals(rbc, vmecInput.rbc);

		double[][] zbs = new double[2 * VmecIndataNamelist.ntord + 1][VmecIndataNamelist.mpold + 1];
		zbs[VmecIndataNamelist.ntord + 0][0] = 0.0;
		zbs[VmecIndataNamelist.ntord + 1][0] = -0.231215947095178;
		zbs[VmecIndataNamelist.ntord + 2][0] = 3.000448318252753E-004;
		zbs[VmecIndataNamelist.ntord + 3][0] = 2.451759201610496E-004;
		zbs[VmecIndataNamelist.ntord + 4][0] = 2.002740558876558E-003;
		zbs[VmecIndataNamelist.ntord + 5][0] = 2.042500448056660E-004;
		zbs[VmecIndataNamelist.ntord + 6][0] = 1.110653375943014E-004;
		zbs[VmecIndataNamelist.ntord + 7][0] = 1.233633688773586E-004;
		zbs[VmecIndataNamelist.ntord + 8][0] = -5.454884672029991E-004;
		//
		zbs[VmecIndataNamelist.ntord - 8][1] = -6.512120712505940E-004;
		zbs[VmecIndataNamelist.ntord - 7][1] = -2.473300842408858E-004;
		zbs[VmecIndataNamelist.ntord - 6][1] = 1.226721407577620E-003;
		zbs[VmecIndataNamelist.ntord - 5][1] = 1.157564964099680E-005;
		zbs[VmecIndataNamelist.ntord - 4][1] = 1.948945874541883E-004;
		zbs[VmecIndataNamelist.ntord - 3][1] = 1.584625523874076E-003;
		zbs[VmecIndataNamelist.ntord - 2][1] = 5.725438140714447E-003;
		zbs[VmecIndataNamelist.ntord - 1][1] = 6.183781833751985E-003;
		zbs[VmecIndataNamelist.ntord + 0][1] = 0.628787962252855;
		zbs[VmecIndataNamelist.ntord + 1][1] = 0.223524141570276;
		zbs[VmecIndataNamelist.ntord + 2][1] = 1.323272993170185E-002;
		zbs[VmecIndataNamelist.ntord + 3][1] = -6.316519355198045E-004;
		zbs[VmecIndataNamelist.ntord + 4][1] = -4.160839711727137E-004;
		zbs[VmecIndataNamelist.ntord + 5][1] = -1.431427743599354E-004;
		zbs[VmecIndataNamelist.ntord + 6][1] = -3.965264842806391E-004;
		zbs[VmecIndataNamelist.ntord + 7][1] = -2.485908172021094E-005;
		zbs[VmecIndataNamelist.ntord + 8][1] = 3.549072371917396E-004;
		//
		zbs[VmecIndataNamelist.ntord - 8][2] = -3.300731983353361E-004;
		zbs[VmecIndataNamelist.ntord - 7][2] = 4.473218222210373E-004;
		zbs[VmecIndataNamelist.ntord - 6][2] = -6.065256605704031E-004;
		zbs[VmecIndataNamelist.ntord - 5][2] = 3.241308593613608E-004;
		zbs[VmecIndataNamelist.ntord - 4][2] = 6.840233808173711E-005;
		zbs[VmecIndataNamelist.ntord - 3][2] = -9.056185242554235E-005;
		zbs[VmecIndataNamelist.ntord - 2][2] = 1.955219044278804E-004;
		zbs[VmecIndataNamelist.ntord - 1][2] = 7.523213388325822E-003;
		zbs[VmecIndataNamelist.ntord + 0][2] = -5.215130122068161E-003;
		zbs[VmecIndataNamelist.ntord + 1][2] = 2.382604615534300E-002;
		zbs[VmecIndataNamelist.ntord + 2][2] = -5.184435337334851E-002;
		zbs[VmecIndataNamelist.ntord + 3][2] = 2.258335290509644E-003;
		zbs[VmecIndataNamelist.ntord + 4][2] = 6.590696228657478E-004;
		zbs[VmecIndataNamelist.ntord + 5][2] = 9.616300338029773E-005;
		zbs[VmecIndataNamelist.ntord + 6][2] = 1.018064657002039E-004;
		zbs[VmecIndataNamelist.ntord + 7][2] = -6.953157311656576E-005;
		zbs[VmecIndataNamelist.ntord + 8][2] = 9.909633701842818E-005;
		//
		zbs[VmecIndataNamelist.ntord - 8][3] = -4.519578534113128E-005;
		zbs[VmecIndataNamelist.ntord - 7][3] = -8.527405237368145E-005;
		zbs[VmecIndataNamelist.ntord - 6][3] = 1.753730817444603E-004;
		zbs[VmecIndataNamelist.ntord - 5][3] = -2.326348438307375E-005;
		zbs[VmecIndataNamelist.ntord - 4][3] = -2.624975041422071E-004;
		zbs[VmecIndataNamelist.ntord - 3][3] = -1.390802867156250E-005;
		zbs[VmecIndataNamelist.ntord - 2][3] = 8.074618587295024E-004;
		zbs[VmecIndataNamelist.ntord - 1][3] = -1.531636429051088E-003;
		zbs[VmecIndataNamelist.ntord + 0][3] = -1.645890760289522E-003;
		zbs[VmecIndataNamelist.ntord + 1][3] = -8.199761340258446E-003;
		zbs[VmecIndataNamelist.ntord + 2][3] = 5.766395955724331E-003;
		zbs[VmecIndataNamelist.ntord + 3][3] = 1.163122094676383E-002;
		zbs[VmecIndataNamelist.ntord + 4][3] = -5.968673861601821E-004;
		zbs[VmecIndataNamelist.ntord + 5][3] = 1.622756876768911E-004;
		zbs[VmecIndataNamelist.ntord + 6][3] = -3.916688872863889E-005;
		zbs[VmecIndataNamelist.ntord + 7][3] = -6.313722893293626E-005;
		zbs[VmecIndataNamelist.ntord + 8][3] = -3.394782233208609E-005;
		//
		zbs[VmecIndataNamelist.ntord - 8][4] = -2.405330978649966E-005;
		zbs[VmecIndataNamelist.ntord - 7][4] = 2.084090191207986E-004;
		zbs[VmecIndataNamelist.ntord - 6][4] = -6.937710034163128E-005;
		zbs[VmecIndataNamelist.ntord - 5][4] = 3.866322821755758E-005;
		zbs[VmecIndataNamelist.ntord - 4][4] = 1.296929167587253E-004;
		zbs[VmecIndataNamelist.ntord - 3][4] = -1.292084644773919E-004;
		zbs[VmecIndataNamelist.ntord - 2][4] = -9.511070777900972E-005;
		zbs[VmecIndataNamelist.ntord - 1][4] = -2.460373856685125E-004;
		zbs[VmecIndataNamelist.ntord + 0][4] = 4.198477990498741E-003;
		zbs[VmecIndataNamelist.ntord + 1][4] = -2.836264233928369E-003;
		zbs[VmecIndataNamelist.ntord + 2][4] = 9.239409675099113E-003;
		zbs[VmecIndataNamelist.ntord + 3][4] = -4.334921025609129E-003;
		zbs[VmecIndataNamelist.ntord + 4][4] = -1.784089978532210E-003;
		zbs[VmecIndataNamelist.ntord + 5][4] = 4.675626087533873E-004;
		zbs[VmecIndataNamelist.ntord + 6][4] = -3.540448326425835E-005;
		zbs[VmecIndataNamelist.ntord + 7][4] = 5.212974445405975E-005;
		zbs[VmecIndataNamelist.ntord + 8][4] = 2.676132432009841E-005;
		//
		zbs[VmecIndataNamelist.ntord - 8][5] = 3.279177584182842E-005;
		zbs[VmecIndataNamelist.ntord - 7][5] = 2.724539949718756E-005;
		zbs[VmecIndataNamelist.ntord - 6][5] = -1.593085585738485E-004;
		zbs[VmecIndataNamelist.ntord - 5][5] = 2.410957443918596E-004;
		zbs[VmecIndataNamelist.ntord - 4][5] = -1.044020785453555E-004;
		zbs[VmecIndataNamelist.ntord - 3][5] = 2.056825728184085E-005;
		zbs[VmecIndataNamelist.ntord - 2][5] = -3.741888572862706E-005;
		zbs[VmecIndataNamelist.ntord - 1][5] = 3.578584050879604E-004;
		zbs[VmecIndataNamelist.ntord + 0][5] = 6.752128022657002E-004;
		zbs[VmecIndataNamelist.ntord + 1][5] = 1.854756872143505E-003;
		zbs[VmecIndataNamelist.ntord + 2][5] = -1.507318481326176E-003;
		zbs[VmecIndataNamelist.ntord + 3][5] = -9.597505773106182E-004;
		zbs[VmecIndataNamelist.ntord + 4][5] = 9.885927576869936E-004;
		zbs[VmecIndataNamelist.ntord + 5][5] = 1.133619037566474E-005;
		zbs[VmecIndataNamelist.ntord + 6][5] = -1.325292640508833E-004;
		zbs[VmecIndataNamelist.ntord + 7][5] = 2.726109026309649E-005;
		zbs[VmecIndataNamelist.ntord + 8][5] = -1.850542092438669E-005;
		//
		zbs[VmecIndataNamelist.ntord - 8][6] = 7.386763623731291E-007;
		zbs[VmecIndataNamelist.ntord - 7][6] = -2.508309149435569E-005;
		zbs[VmecIndataNamelist.ntord - 6][6] = 2.117983834955785E-005;
		zbs[VmecIndataNamelist.ntord - 5][6] = -5.753278329735706E-005;
		zbs[VmecIndataNamelist.ntord - 4][6] = -1.155249029943377E-005;
		zbs[VmecIndataNamelist.ntord - 3][6] = 1.153115758377326E-004;
		zbs[VmecIndataNamelist.ntord - 2][6] = -2.830337070611086E-005;
		zbs[VmecIndataNamelist.ntord - 1][6] = 1.324540090482543E-005;
		zbs[VmecIndataNamelist.ntord + 0][6] = -7.800374644189332E-004;
		zbs[VmecIndataNamelist.ntord + 1][6] = 3.838347842484833E-004;
		zbs[VmecIndataNamelist.ntord + 2][6] = -6.137053321851298E-004;
		zbs[VmecIndataNamelist.ntord + 3][6] = -1.301407081337896E-003;
		zbs[VmecIndataNamelist.ntord + 4][6] = -4.955356310734652E-004;
		zbs[VmecIndataNamelist.ntord + 5][6] = 3.462627709485832E-004;
		zbs[VmecIndataNamelist.ntord + 6][6] = -9.919741284083595E-005;
		zbs[VmecIndataNamelist.ntord + 7][6] = -1.784091935229791E-005;
		zbs[VmecIndataNamelist.ntord + 8][6] = -2.302543207412257E-005;
		//
		zbs[VmecIndataNamelist.ntord - 8][7] = -1.484141014465394E-005;
		zbs[VmecIndataNamelist.ntord - 7][7] = -1.816080256001114E-005;
		zbs[VmecIndataNamelist.ntord - 6][7] = -3.833112693126000E-006;
		zbs[VmecIndataNamelist.ntord - 5][7] = -1.865844947458311E-006;
		zbs[VmecIndataNamelist.ntord - 4][7] = -9.801486449695234E-007;
		zbs[VmecIndataNamelist.ntord - 3][7] = 8.061675676629167E-006;
		zbs[VmecIndataNamelist.ntord - 2][7] = -1.094928404347744E-006;
		zbs[VmecIndataNamelist.ntord - 1][7] = -3.601316979928921E-005;
		zbs[VmecIndataNamelist.ntord + 0][7] = 2.940524142344392E-004;
		zbs[VmecIndataNamelist.ntord + 1][7] = -8.496033274994724E-004;
		zbs[VmecIndataNamelist.ntord + 2][7] = 6.255838033712031E-004;
		zbs[VmecIndataNamelist.ntord + 3][7] = 3.900613621829464E-004;
		zbs[VmecIndataNamelist.ntord + 4][7] = 1.845971251755716E-003;
		zbs[VmecIndataNamelist.ntord + 5][7] = 2.917639239820137E-004;
		zbs[VmecIndataNamelist.ntord + 6][7] = -5.665554773695777E-005;
		zbs[VmecIndataNamelist.ntord + 7][7] = 1.038283071268155E-004;
		zbs[VmecIndataNamelist.ntord + 8][7] = -2.979906810937261E-006;
		//
		Assertions.assertArrayEquals(zbs, vmecInput.zbs);

		double[][] rbs = new double[2 * VmecIndataNamelist.ntord + 1][VmecIndataNamelist.mpold + 1];
		Assertions.assertArrayEquals(rbs, vmecInput.rbs);

		double[][] zbc = new double[2 * VmecIndataNamelist.ntord + 1][VmecIndataNamelist.mpold + 1];
		Assertions.assertArrayEquals(zbc, vmecInput.zbc);
	}


	@Test
	void vmecInputTest_w7x() throws IOException {
		String inputFilename = "/de/labathome/FortranNamelist/input/input-w7x.txt";
		InputStream inputFileStream = getClass().getResourceAsStream(inputFilename);

		String inputFile = new String(inputFileStream.readAllBytes(), StandardCharsets.UTF_8);

		VmecIndataNamelist vmecInput = new VmecIndataNamelist();

		// parse the namelist into a Java object
		FortranNamelist parser = new FortranNamelist(inputFile, "indata", vmecInput);
		vmecInput = (VmecIndataNamelist) parser.getParsed();
		assertNotNull(vmecInput);
		vmecInput.sanitize();


		// --------------------
		// Numerical Resolution and Symmetry Assumptions
		Assertions.assertEquals(vmecInput.lasym, false);
		// nfp
		// mpol
		// ntor
		// ntheta
		// nzeta

		// --------------------
		// Multi-Grid Steps
		// ns_array
		// ftol_array
		// niter_array

		// --------------------
		// Global Physics Parameters
		// phiedge
		// ncurr

		// --------------------
		// Profile of Mass or Pressure
		// pmass_type
		// am
		// am_aux_s
		// am_aux_f
		// pres_scale
		// gamma
		// spres_ped

		// --------------------
		// (Initial Guess for) Rotational Transform Profile
		// piota_type
		// ai
		// ai_aux_s
		// ai_aux_f


		// --------------------
		// (Initial Guess for) Toroidal Current Profile
		// pcurr_type
		// ac
		// ac_aux_s
		// ac_aux_f
		// curtor
		// bloat




		// --------------------
		// Free-Boundary Parameters
		// lfreeb
		// mgrid_file
		// extcur
		// nvacskip



		// --------------------
		// Tweaking Parameters
		// nstep
		// aphi
		// delt
		// tcon0
		// lforbal


		// --------------------
		// Initial Guess for Magnetic Axis Geometry
		// raxis_cc
		// zaxis_cs
		// raxis_cs
		// zaxis_cc


		// --------------------
		// (Initial Guess for) Boundary Geometry
		// rbc
		// zbs
		// rbs
		// zbc

	}

}
