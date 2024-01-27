package de.labathome;

import static org.junit.jupiter.api.Assertions.assertNotNull;

import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.List;

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
		Assertions.assertEquals(false, vmecInput.lasym);
		Assertions.assertEquals(5, vmecInput.nfp);
		Assertions.assertEquals(8, vmecInput.mpol);
		Assertions.assertEquals(8, vmecInput.ntor);
		Assertions.assertEquals(0, vmecInput.ntheta);
		Assertions.assertEquals(36, vmecInput.nzeta);

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
		Assertions.assertEquals(false, vmecInput.lasym);
		Assertions.assertEquals(5, vmecInput.nfp);
		Assertions.assertEquals(12, vmecInput.mpol);
		Assertions.assertEquals(12, vmecInput.ntor);
		Assertions.assertEquals(32, vmecInput.ntheta);
		Assertions.assertEquals(36, vmecInput.nzeta);

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
		Assertions.assertEquals(false, vmecInput.lasym);
		Assertions.assertEquals(1, vmecInput.nfp);
		Assertions.assertEquals(6, vmecInput.mpol);
		Assertions.assertEquals(0, vmecInput.ntor);
		Assertions.assertEquals(0, vmecInput.ntheta);
		Assertions.assertEquals(0, vmecInput.nzeta);

		// --------------------
		// Multi-Grid Steps
		int[] ns_array = new int[100];
		ns_array[0] = 16;
		ns_array[1] = 31;
		Assertions.assertArrayEquals(ns_array, vmecInput.ns_array);

		double[] ftol_array = new double[100];
		ftol_array[0] = 1.0e-6;
		ftol_array[1] = 5.0e-11;
		Assertions.assertArrayEquals(ftol_array, vmecInput.ftol_array);

		int[] niter_array = new int[100];
		Arrays.fill(niter_array, 600);
		Assertions.assertArrayEquals(niter_array, vmecInput.niter_array);

		// --------------------
		// Global Physics Parameters
		Assertions.assertEquals(-1.3534E+01, vmecInput.phiedge);
		Assertions.assertEquals(0, vmecInput.ncurr);

		// --------------------
		// Profile of Mass or Pressure
		Assertions.assertEquals("power_series", vmecInput.pmass_type.toLowerCase());
		double[] am = new double[21];
		am[0] = 0.2762993401771816E6;
		am[1] = -1.398637555311701E6;
		am[2] = 6.435929178901299E6;
		am[3] = -19.13351221515564E6;
		am[4] = 30.46522063279099E6;
		am[5] = -23.99696824024823E6;
		am[6] = 7.363942466650127E6;
		Assertions.assertArrayEquals(am, vmecInput.am);
		Assertions.assertArrayEquals(new double[VmecIndataNamelist.nsd], vmecInput.am_aux_s);
		Assertions.assertArrayEquals(new double[VmecIndataNamelist.nsd], vmecInput.am_aux_f);
		Assertions.assertEquals(1.0, vmecInput.pres_scale);
		Assertions.assertEquals(0.0, vmecInput.gamma);
		Assertions.assertEquals(1.0, vmecInput.spres_ped);

		// --------------------
		// (Initial Guess for) Rotational Transform Profile
		Assertions.assertEquals("power_series", vmecInput.piota_type.toLowerCase());
		double[] ai = new double[21];
		ai[0] = 2.354498895432447;
		ai[1] = -9.907517672615517;
		ai[2] = 42.86570171033211;
		ai[3] = -122.8693254810872;
		ai[4] = 188.6838680668763;
		ai[5] = -143.1773314731515;
		ai[6] = 42.19838185779675;
		Assertions.assertArrayEquals(ai, vmecInput.ai);
		Assertions.assertArrayEquals(new double[VmecIndataNamelist.nsd], vmecInput.ai_aux_s);
		Assertions.assertArrayEquals(new double[VmecIndataNamelist.nsd], vmecInput.ai_aux_f);

		// --------------------
		// (Initial Guess for) Toroidal Current Profile
		Assertions.assertEquals("power_series", vmecInput.pcurr_type.toLowerCase());
		Assertions.assertArrayEquals(new double[21], vmecInput.ac);
		Assertions.assertArrayEquals(new double[VmecIndataNamelist.nsd], vmecInput.ac_aux_s);
		Assertions.assertArrayEquals(new double[VmecIndataNamelist.nsd], vmecInput.ac_aux_f);
		Assertions.assertEquals(1.582E+06, vmecInput.curtor);
		Assertions.assertEquals(1.0, vmecInput.bloat);

		// --------------------
		// Free-Boundary Parameters
		Assertions.assertEquals(true, vmecInput.lfreeb);
		Assertions.assertEquals("mgrid.tftr", vmecInput.mgrid_file);
		double[] extcur = new double[VmecIndataNamelist.nigroup];
		extcur[0] = -6.800E+01;
		extcur[1] = -1.070E+01;
		extcur[2] = 1.435E+01;
		extcur[3] = -1.405E+01;
		extcur[4] = 2.498E-02;
		Assertions.assertArrayEquals(extcur, vmecInput.extcur);
		Assertions.assertEquals(12, vmecInput.nvacskip);

		// --------------------
		// Tweaking Parameters
		Assertions.assertEquals(100, vmecInput.nstep);
		double[] aphi = new double[20];
		aphi[0] = 1.0;
		Assertions.assertArrayEquals(aphi, vmecInput.aphi);
		Assertions.assertEquals(1.0, vmecInput.delt);
		Assertions.assertEquals(1.0, vmecInput.tcon0);
		Assertions.assertEquals(false, vmecInput.lforbal);

		// --------------------
		// Initial Guess for Magnetic Axis Geometry
		double[] raxis_cc = new double[VmecIndataNamelist.ntord + 1];
		raxis_cc[0] = 2.772E+00;
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
		rbc[VmecIndataNamelist.ntord + 0][0] = 2.599E+00;
		rbc[VmecIndataNamelist.ntord + 0][1] = 9.370E-01;
		Assertions.assertArrayEquals(rbc, vmecInput.rbc);

		double[][] zbs = new double[2 * VmecIndataNamelist.ntord + 1][VmecIndataNamelist.mpold + 1];
		zbs[VmecIndataNamelist.ntord + 0][1] = 9.590E-01;
		Assertions.assertArrayEquals(zbs, vmecInput.zbs);

		double[][] rbs = new double[2 * VmecIndataNamelist.ntord + 1][VmecIndataNamelist.mpold + 1];
		Assertions.assertArrayEquals(rbs, vmecInput.rbs);

		double[][] zbc = new double[2 * VmecIndataNamelist.ntord + 1][VmecIndataNamelist.mpold + 1];
		Assertions.assertArrayEquals(zbc, vmecInput.zbc);
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
		Assertions.assertEquals(false, vmecInput.lasym);
		Assertions.assertEquals(12, vmecInput.nfp);
		Assertions.assertEquals(6, vmecInput.mpol);
		Assertions.assertEquals(6, vmecInput.ntor);
		Assertions.assertEquals(0, vmecInput.ntheta);
		Assertions.assertEquals(0, vmecInput.nzeta);

		// --------------------
		// Multi-Grid Steps
		int[] ns_array = new int[100];
		ns_array[0] = 26;
		ns_array[1] = 51;
		Assertions.assertArrayEquals(ns_array, vmecInput.ns_array);

		double[] ftol_array = new double[100];
		ftol_array[0] = 1.0e-8;
		ftol_array[1] = 1.0e-11;
		Assertions.assertArrayEquals(ftol_array, vmecInput.ftol_array);

		int[] niter_array = new int[100];
		Arrays.fill(niter_array, 3000);
		Assertions.assertArrayEquals(niter_array, vmecInput.niter_array);

		// --------------------
		// Global Physics Parameters
		Assertions.assertEquals(0.23, vmecInput.phiedge);
		Assertions.assertEquals(0, vmecInput.ncurr);

		// --------------------
		// Profile of Mass or Pressure
		Assertions.assertEquals("power_series", vmecInput.pmass_type.toLowerCase());
		double[] am = new double[21];
		am[0] = 3.6500000E+04;
		am[1] = -7.300000E+04;
		am[2] = 3.650000E+04;
		Assertions.assertArrayEquals(am, vmecInput.am);
		Assertions.assertArrayEquals(new double[VmecIndataNamelist.nsd], vmecInput.am_aux_s);
		Assertions.assertArrayEquals(new double[VmecIndataNamelist.nsd], vmecInput.am_aux_f);
		Assertions.assertEquals(1.0, vmecInput.pres_scale);
		Assertions.assertEquals(1.667, vmecInput.gamma);
		Assertions.assertEquals(1.0, vmecInput.spres_ped);

		// --------------------
		// (Initial Guess for) Rotational Transform Profile
		Assertions.assertEquals("power_series", vmecInput.piota_type.toLowerCase());
		double[] ai = new double[21];
		ai[0] = 0.3000000;
		ai[1] = 0.4200000;
		ai[2] = 0.2500000;
		Assertions.assertArrayEquals(ai, vmecInput.ai);
		Assertions.assertArrayEquals(new double[VmecIndataNamelist.nsd], vmecInput.ai_aux_s);
		Assertions.assertArrayEquals(new double[VmecIndataNamelist.nsd], vmecInput.ai_aux_f);

		// --------------------
		// (Initial Guess for) Toroidal Current Profile
		Assertions.assertEquals("power_series", vmecInput.pcurr_type.toLowerCase());
		Assertions.assertArrayEquals(new double[21], vmecInput.ac);
		Assertions.assertArrayEquals(new double[VmecIndataNamelist.nsd], vmecInput.ac_aux_s);
		Assertions.assertArrayEquals(new double[VmecIndataNamelist.nsd], vmecInput.ac_aux_f);
		Assertions.assertEquals(0.0, vmecInput.curtor);
		Assertions.assertEquals(1.0, vmecInput.bloat);

		// --------------------
		// Free-Boundary Parameters
		Assertions.assertEquals(true, vmecInput.lfreeb);
		Assertions.assertEquals("NONE", vmecInput.mgrid_file);
		double[] extcur = new double[VmecIndataNamelist.nigroup];
		extcur[0] = -3.9385E2;
		extcur[1] = 29.4795;
		extcur[2] = -61.3872;
		extcur[3] = -3.133E-2;
		Assertions.assertArrayEquals(extcur, vmecInput.extcur);
		Assertions.assertEquals(12, vmecInput.nvacskip);

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
		raxis_cc[0] = 1.72;
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
		rbc[VmecIndataNamelist.ntord + 0][0] =  1.7183E+00;
		rbc[VmecIndataNamelist.ntord + 0][1] =  2.1404E-01;
		rbc[VmecIndataNamelist.ntord + 1][1] = -5.5576E-02;
		rbc[VmecIndataNamelist.ntord - 1][2] =  9.8162E-04;
		rbc[VmecIndataNamelist.ntord + 0][2] =  4.0247E-03;
		rbc[VmecIndataNamelist.ntord + 1][2] = -6.1395E-03;
		rbc[VmecIndataNamelist.ntord + 0][3] = -2.2444E-04;
		rbc[VmecIndataNamelist.ntord + 1][3] = -7.2860E-04;
		rbc[VmecIndataNamelist.ntord + 0][4] =  2.0864E-04;
		rbc[VmecIndataNamelist.ntord + 1][4] = -5.1596E-04;
		Assertions.assertArrayEquals(rbc, vmecInput.rbc);

		double[][] zbs = new double[2 * VmecIndataNamelist.ntord + 1][VmecIndataNamelist.mpold + 1];
		zbs[VmecIndataNamelist.ntord + 0][0] =  0.0000E+00;
		zbs[VmecIndataNamelist.ntord + 0][1] =  2.5040E-01;
		zbs[VmecIndataNamelist.ntord + 1][1] =  6.2457E-02;
		zbs[VmecIndataNamelist.ntord - 1][2] =  1.3940E-03;
		zbs[VmecIndataNamelist.ntord + 0][2] = -2.6518E-03;
		zbs[VmecIndataNamelist.ntord + 1][2] =  1.0371E-02;
		zbs[VmecIndataNamelist.ntord + 0][3] = -2.2444E-04;
		zbs[VmecIndataNamelist.ntord + 1][3] = -7.2860E-04;
		zbs[VmecIndataNamelist.ntord + 0][4] =  2.0864E-04;
		zbs[VmecIndataNamelist.ntord + 1][4] = -5.1596E-04;
		Assertions.assertArrayEquals(zbs, vmecInput.zbs);

		double[][] rbs = new double[2 * VmecIndataNamelist.ntord + 1][VmecIndataNamelist.mpold + 1];
		Assertions.assertArrayEquals(rbs, vmecInput.rbs);

		double[][] zbc = new double[2 * VmecIndataNamelist.ntord + 1][VmecIndataNamelist.mpold + 1];
		Assertions.assertArrayEquals(zbc, vmecInput.zbc);
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
		Assertions.assertEquals(false, vmecInput.lasym);
		Assertions.assertEquals(10, vmecInput.nfp);
		Assertions.assertEquals(10, vmecInput.mpol);
		Assertions.assertEquals(8, vmecInput.ntor);
		Assertions.assertEquals(32, vmecInput.ntheta);
		Assertions.assertEquals(32, vmecInput.nzeta);

		// --------------------
		// Multi-Grid Steps
		int[] ns_array = new int[100];
		ns_array[0] = 11;
		ns_array[1] = 31;
		ns_array[2] = 61;
		Assertions.assertArrayEquals(ns_array, vmecInput.ns_array);

		double[] ftol_array = new double[100];
		ftol_array[0] = 1.0e-7;
		ftol_array[1] = 1.0e-11;
		ftol_array[2] = 1.0e-15;
		Assertions.assertArrayEquals(ftol_array, vmecInput.ftol_array);

		int[] niter_array = new int[100];
		Arrays.fill(niter_array, 10000);
		Assertions.assertArrayEquals(niter_array, vmecInput.niter_array);

		// --------------------
		// Global Physics Parameters
		Assertions.assertEquals(-8.2523E+01, vmecInput.phiedge);
		Assertions.assertEquals(1, vmecInput.ncurr);

		// --------------------
		// Profile of Mass or Pressure
		Assertions.assertEquals("power_series", vmecInput.pmass_type.toLowerCase());
		double[] am = new double[21];
		am[0] = 1.1359E+06;
		am[1] = -1.5154E+06;
		am[2] = 4.3242E+05;
		am[3] = 1.0275E+06;
		am[4] = -2.3580E+06;
		am[5] = -4.4279E+04;
		am[6] = 2.8985E+06;
		am[7] = 1.7771E+06;
		am[8] = -2.4165E+06;
		am[9] = -4.2657E+06;
		am[10] = 3.3914E+06;
		Assertions.assertArrayEquals(am, vmecInput.am);
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
		Assertions.assertArrayEquals(new double[21], vmecInput.ac);
		Assertions.assertArrayEquals(new double[VmecIndataNamelist.nsd], vmecInput.ac_aux_s);
		Assertions.assertArrayEquals(new double[VmecIndataNamelist.nsd], vmecInput.ac_aux_f);
		Assertions.assertEquals(0.0, vmecInput.curtor);
		Assertions.assertEquals(1.0, vmecInput.bloat);

		// --------------------
		// Free-Boundary Parameters
		Assertions.assertEquals(false, vmecInput.lfreeb);
		Assertions.assertEquals("dummy", vmecInput.mgrid_file);
		Assertions.assertArrayEquals(new double[VmecIndataNamelist.nigroup], vmecInput.extcur);
		Assertions.assertEquals(10, vmecInput.nvacskip);

		// --------------------
		// Tweaking Parameters
		Assertions.assertEquals(500, vmecInput.nstep);
		double[] aphi = new double[20];
		aphi[0] = 1.0;
		Assertions.assertArrayEquals(aphi, vmecInput.aphi);
		Assertions.assertEquals(0.9, vmecInput.delt);
		Assertions.assertEquals(1.0, vmecInput.tcon0);
		Assertions.assertEquals(false, vmecInput.lforbal);

		// --------------------
		// Initial Guess for Magnetic Axis Geometry
		double[] raxis_cc = new double[VmecIndataNamelist.ntord + 1];
		raxis_cc[0] = 0.1473E+02;
		raxis_cc[1] = -0.1431E-01;
		raxis_cc[2] = -0.2578E-04;
		raxis_cc[3] = -0.2536E-05;
		raxis_cc[4] = 0.6638E-06;
		raxis_cc[5] = -0.9506E-06;
		raxis_cc[6] = 0.3329E-06;
		Assertions.assertArrayEquals(raxis_cc, vmecInput.raxis_cc);

		double[] zaxis_cs = new double[VmecIndataNamelist.ntord + 1];
		zaxis_cs[1] = -0.1396E-01;
		zaxis_cs[2] = 0.4200E-05;
		zaxis_cs[3] = 0.7295E-05;
		zaxis_cs[4] = 0.9739E-06;
		zaxis_cs[5] = 0.4574E-05;
		zaxis_cs[6] = 0.8842E-06;
		Assertions.assertArrayEquals(zaxis_cs, vmecInput.zaxis_cs);

		double[] raxis_cs = new double[VmecIndataNamelist.ntord + 1];
		Assertions.assertArrayEquals(raxis_cs, vmecInput.raxis_cs);

		double[] zaxis_cc = new double[VmecIndataNamelist.ntord + 1];
		Assertions.assertArrayEquals(zaxis_cc, vmecInput.zaxis_cc);

		// --------------------
		// (Initial Guess for) Boundary Geometry

		double[][] rbc = new double[2 * VmecIndataNamelist.ntord + 1][VmecIndataNamelist.mpold + 1];
		rbc[VmecIndataNamelist.ntord + 0][0] =  1.4607E+01;
		rbc[VmecIndataNamelist.ntord + 1][0] =  1.7765E-01;
		rbc[VmecIndataNamelist.ntord + 2][0] = -1.3555E-02;
		rbc[VmecIndataNamelist.ntord + 3][0] =  3.9756E-03;
		rbc[VmecIndataNamelist.ntord + 4][0] = -2.9062E-03;
		rbc[VmecIndataNamelist.ntord + 5][0] =  1.3122E-03;
		rbc[VmecIndataNamelist.ntord + 6][0] = -3.1275E-03;
		//
		rbc[VmecIndataNamelist.ntord - 6][1] =  1.3024E-03;
		rbc[VmecIndataNamelist.ntord - 5][1] =  3.3317E-03;
		rbc[VmecIndataNamelist.ntord - 4][1] = -1.4495E-03;
		rbc[VmecIndataNamelist.ntord - 3][1] =  6.2439E-03;
		rbc[VmecIndataNamelist.ntord - 2][1] = -9.3805E-03;
		rbc[VmecIndataNamelist.ntord - 1][1] = -9.2177E-01;
		rbc[VmecIndataNamelist.ntord + 0][1] =  2.3492E+00;
		rbc[VmecIndataNamelist.ntord + 1][1] =  2.6518E-03;
		rbc[VmecIndataNamelist.ntord + 2][1] = -5.6435E-04;
		rbc[VmecIndataNamelist.ntord + 3][1] =  5.4447E-03;
		rbc[VmecIndataNamelist.ntord + 4][1] =  1.7714E-04;
		rbc[VmecIndataNamelist.ntord + 5][1] =  3.0477E-03;
		rbc[VmecIndataNamelist.ntord + 6][1] =  9.0954E-04;
		//
		rbc[VmecIndataNamelist.ntord - 6][2] =  4.2667E-03;
		rbc[VmecIndataNamelist.ntord - 5][2] = -4.9128E-03;
		rbc[VmecIndataNamelist.ntord - 4][2] =  3.9295E-03;
		rbc[VmecIndataNamelist.ntord - 3][2] = -5.0165E-03;
		rbc[VmecIndataNamelist.ntord - 2][2] =  1.2406E-02;
		rbc[VmecIndataNamelist.ntord - 1][2] =  9.1045E-02;
		rbc[VmecIndataNamelist.ntord + 0][2] = -1.9001E-03;
		rbc[VmecIndataNamelist.ntord + 1][2] = -2.9866E-02;
		rbc[VmecIndataNamelist.ntord + 2][2] =  9.6029E-03;
		rbc[VmecIndataNamelist.ntord + 3][2] =  2.0260E-04;
		rbc[VmecIndataNamelist.ntord + 4][2] = -2.2851E-03;
		rbc[VmecIndataNamelist.ntord + 5][2] =  1.0195E-03;
		rbc[VmecIndataNamelist.ntord + 6][2] =  8.8206E-04;
		//
		rbc[VmecIndataNamelist.ntord - 6][3] = -4.0757E-03;
		rbc[VmecIndataNamelist.ntord - 5][3] = -2.6563E-03;
		rbc[VmecIndataNamelist.ntord - 4][3] = -1.6465E-03;
		rbc[VmecIndataNamelist.ntord - 3][3] = -1.8381E-03;
		rbc[VmecIndataNamelist.ntord - 2][3] =  2.0356E-02;
		rbc[VmecIndataNamelist.ntord - 1][3] =  3.2462E-02;
		rbc[VmecIndataNamelist.ntord + 0][3] =  1.8047E-02;
		rbc[VmecIndataNamelist.ntord + 1][3] = -1.7257E-03;
		rbc[VmecIndataNamelist.ntord + 2][3] = -1.7472E-03;
		rbc[VmecIndataNamelist.ntord + 3][3] = -1.4645E-03;
		rbc[VmecIndataNamelist.ntord + 4][3] = -3.0665E-03;
		rbc[VmecIndataNamelist.ntord + 5][3] = -1.3051E-03;
		rbc[VmecIndataNamelist.ntord + 6][3] = -3.1243E-03;
		//
		rbc[VmecIndataNamelist.ntord - 6][4] =  4.7381E-04;
		rbc[VmecIndataNamelist.ntord - 5][4] = -5.2082E-04;
		rbc[VmecIndataNamelist.ntord - 4][4] =  5.8916E-04;
		rbc[VmecIndataNamelist.ntord - 3][4] = -1.5758E-03;
		rbc[VmecIndataNamelist.ntord - 2][4] = -1.5781E-02;
		rbc[VmecIndataNamelist.ntord - 1][4] = -1.8729E-02;
		rbc[VmecIndataNamelist.ntord + 0][4] =  1.1545E-02;
		rbc[VmecIndataNamelist.ntord + 1][4] =  3.0756E-03;
		rbc[VmecIndataNamelist.ntord + 2][4] = -5.3511E-04;
		rbc[VmecIndataNamelist.ntord + 3][4] =  2.1879E-04;
		rbc[VmecIndataNamelist.ntord + 4][4] =  4.0675E-04;
		rbc[VmecIndataNamelist.ntord + 5][4] =  7.3581E-04;
		rbc[VmecIndataNamelist.ntord + 6][4] =  2.2203E-04;
		//
		rbc[VmecIndataNamelist.ntord - 6][5] =  1.0517E-04;
		rbc[VmecIndataNamelist.ntord - 5][5] = -1.6563E-05;
		rbc[VmecIndataNamelist.ntord - 4][5] =  2.2566E-04;
		rbc[VmecIndataNamelist.ntord - 3][5] =  1.6957E-03;
		rbc[VmecIndataNamelist.ntord - 2][5] = -2.8591E-03;
		rbc[VmecIndataNamelist.ntord - 1][5] = -1.9804E-02;
		rbc[VmecIndataNamelist.ntord + 0][5] = -2.7778E-04;
		rbc[VmecIndataNamelist.ntord + 1][5] =  1.5405E-03;
		rbc[VmecIndataNamelist.ntord + 2][5] =  4.4580E-04;
		rbc[VmecIndataNamelist.ntord + 3][5] =  8.7347E-04;
		rbc[VmecIndataNamelist.ntord + 4][5] = -5.1359E-04;
		rbc[VmecIndataNamelist.ntord + 5][5] =  4.2326E-04;
		rbc[VmecIndataNamelist.ntord + 6][5] = -2.7550E-04;
		//
		rbc[VmecIndataNamelist.ntord - 6][6] =  4.6757E-04;
		rbc[VmecIndataNamelist.ntord - 5][6] =  1.4952E-04;
		rbc[VmecIndataNamelist.ntord - 4][6] = -4.1890E-04;
		rbc[VmecIndataNamelist.ntord - 3][6] =  7.1438E-04;
		rbc[VmecIndataNamelist.ntord - 2][6] =  6.2976E-04;
		rbc[VmecIndataNamelist.ntord - 1][6] = -1.2376E-02;
		rbc[VmecIndataNamelist.ntord + 0][6] = -3.7337E-03;
		rbc[VmecIndataNamelist.ntord + 1][6] =  6.4686E-04;
		rbc[VmecIndataNamelist.ntord + 2][6] =  4.0144E-04;
		rbc[VmecIndataNamelist.ntord + 3][6] = -2.4674E-05;
		rbc[VmecIndataNamelist.ntord + 4][6] = -1.3095E-04;
		rbc[VmecIndataNamelist.ntord + 5][6] =  5.4803E-04;
		rbc[VmecIndataNamelist.ntord + 6][6] =  1.4177E-04;
		//
		rbc[VmecIndataNamelist.ntord - 6][7] = -1.8961E-04;
		rbc[VmecIndataNamelist.ntord - 5][7] = -3.0268E-04;
		rbc[VmecIndataNamelist.ntord - 4][7] =  9.9123E-05;
		rbc[VmecIndataNamelist.ntord - 3][7] =  8.3561E-04;
		rbc[VmecIndataNamelist.ntord - 2][7] =  2.8116E-03;
		rbc[VmecIndataNamelist.ntord - 1][7] = -5.6180E-03;
		rbc[VmecIndataNamelist.ntord + 0][7] = -1.1974E-03;
		rbc[VmecIndataNamelist.ntord + 1][7] = -5.0556E-04;
		rbc[VmecIndataNamelist.ntord + 2][7] = -2.0245E-05;
		rbc[VmecIndataNamelist.ntord + 3][7] =  1.6912E-04;
		rbc[VmecIndataNamelist.ntord + 4][7] =  1.3252E-04;
		rbc[VmecIndataNamelist.ntord + 5][7] =  5.4404E-04;
		rbc[VmecIndataNamelist.ntord + 6][7] =  1.5693E-04;
		//
		Assertions.assertArrayEquals(rbc, vmecInput.rbc);

		double[][] zbs = new double[2 * VmecIndataNamelist.ntord + 1][VmecIndataNamelist.mpold + 1];
		zbs[VmecIndataNamelist.ntord + 0][0] =  0.0000E+00;
		zbs[VmecIndataNamelist.ntord + 1][0] =  1.7520E-01;
		zbs[VmecIndataNamelist.ntord + 2][0] = -5.5886E-03;
		zbs[VmecIndataNamelist.ntord + 3][0] =  2.1715E-03;
		zbs[VmecIndataNamelist.ntord + 4][0] = -2.8229E-03;
		zbs[VmecIndataNamelist.ntord + 5][0] =  4.6609E-03;
		zbs[VmecIndataNamelist.ntord + 6][0] = -1.2143E-03;
		//
		zbs[VmecIndataNamelist.ntord - 6][1] =  1.1669E-04;
		zbs[VmecIndataNamelist.ntord - 5][1] = -3.7551E-03;
		zbs[VmecIndataNamelist.ntord - 4][1] =  2.6225E-03;
		zbs[VmecIndataNamelist.ntord - 3][1] = -4.1894E-03;
		zbs[VmecIndataNamelist.ntord - 2][1] =  1.6500E-02;
		zbs[VmecIndataNamelist.ntord - 1][1] =  9.3818E-01;
		zbs[VmecIndataNamelist.ntord + 0][1] =  2.4639E+00;
		zbs[VmecIndataNamelist.ntord + 1][1] =  1.7039E-02;
		zbs[VmecIndataNamelist.ntord + 2][1] =  7.7155E-03;
		zbs[VmecIndataNamelist.ntord + 3][1] = -3.3930E-03;
		zbs[VmecIndataNamelist.ntord + 4][1] =  1.0016E-03;
		zbs[VmecIndataNamelist.ntord + 5][1] = -3.4721E-03;
		zbs[VmecIndataNamelist.ntord + 6][1] =  5.0819E-04;
		//
		zbs[VmecIndataNamelist.ntord - 6][2] = -1.5798E-03;
		zbs[VmecIndataNamelist.ntord - 5][2] = -2.4941E-03;
		zbs[VmecIndataNamelist.ntord - 4][2] = -3.2112E-04;
		zbs[VmecIndataNamelist.ntord - 3][2] = -1.4767E-03;
		zbs[VmecIndataNamelist.ntord - 2][2] =  5.4646E-03;
		zbs[VmecIndataNamelist.ntord - 1][2] =  1.6500E-02;
		zbs[VmecIndataNamelist.ntord + 0][2] = -6.8627E-02;
		zbs[VmecIndataNamelist.ntord + 1][2] = -2.5114E-03;
		zbs[VmecIndataNamelist.ntord + 2][2] =  6.3576E-03;
		zbs[VmecIndataNamelist.ntord + 3][2] = -2.3277E-03;
		zbs[VmecIndataNamelist.ntord + 4][2] = -2.7521E-03;
		zbs[VmecIndataNamelist.ntord + 5][2] =  2.0969E-03;
		zbs[VmecIndataNamelist.ntord + 6][2] =  3.0422E-04;
		//
		zbs[VmecIndataNamelist.ntord - 6][3] =  2.1192E-03;
		zbs[VmecIndataNamelist.ntord - 5][3] = -5.7772E-04;
		zbs[VmecIndataNamelist.ntord - 4][3] =  2.8470E-03;
		zbs[VmecIndataNamelist.ntord - 3][3] = -5.7875E-03;
		zbs[VmecIndataNamelist.ntord - 2][3] = -2.0740E-02;
		zbs[VmecIndataNamelist.ntord - 1][3] =  2.6512E-02;
		zbs[VmecIndataNamelist.ntord + 0][3] = -2.3445E-02;
		zbs[VmecIndataNamelist.ntord + 1][3] = -2.0347E-03;
		zbs[VmecIndataNamelist.ntord + 2][3] =  6.8606E-03;
		zbs[VmecIndataNamelist.ntord + 3][3] = -1.5784E-03;
		zbs[VmecIndataNamelist.ntord + 4][3] =  3.1090E-03;
		zbs[VmecIndataNamelist.ntord + 5][3] = -5.5238E-04;
		zbs[VmecIndataNamelist.ntord + 6][3] =  2.1421E-03;
		//
		zbs[VmecIndataNamelist.ntord - 6][4] = -3.5303E-04;
		zbs[VmecIndataNamelist.ntord - 5][4] =  6.8739E-04;
		zbs[VmecIndataNamelist.ntord - 4][4] =  2.7879E-04;
		zbs[VmecIndataNamelist.ntord - 3][4] = -3.6348E-04;
		zbs[VmecIndataNamelist.ntord - 2][4] =  5.9639E-03;
		zbs[VmecIndataNamelist.ntord - 1][4] =  2.4917E-02;
		zbs[VmecIndataNamelist.ntord + 0][4] = -1.6743E-02;
		zbs[VmecIndataNamelist.ntord + 1][4] =  2.3455E-03;
		zbs[VmecIndataNamelist.ntord + 2][4] =  9.0423E-04;
		zbs[VmecIndataNamelist.ntord + 3][4] = -1.1928E-03;
		zbs[VmecIndataNamelist.ntord + 4][4] =  1.8489E-04;
		zbs[VmecIndataNamelist.ntord + 5][4] = -8.6419E-04;
		zbs[VmecIndataNamelist.ntord + 6][4] =  3.7669E-04;
		//
		zbs[VmecIndataNamelist.ntord - 6][5] = -5.2479E-04;
		zbs[VmecIndataNamelist.ntord - 5][5] = -3.7521E-04;
		zbs[VmecIndataNamelist.ntord - 4][5] = -7.0973E-05;
		zbs[VmecIndataNamelist.ntord - 3][5] =  5.0243E-05;
		zbs[VmecIndataNamelist.ntord - 2][5] =  2.6809E-03;
		zbs[VmecIndataNamelist.ntord - 1][5] =  1.7578E-02;
		zbs[VmecIndataNamelist.ntord + 0][5] = -3.9804E-03;
		zbs[VmecIndataNamelist.ntord + 1][5] = -2.0682E-03;
		zbs[VmecIndataNamelist.ntord + 2][5] =  2.8655E-04;
		zbs[VmecIndataNamelist.ntord + 3][5] = -1.5860E-03;
		zbs[VmecIndataNamelist.ntord + 4][5] =  5.2804E-04;
		zbs[VmecIndataNamelist.ntord + 5][5] = -4.3421E-04;
		zbs[VmecIndataNamelist.ntord + 6][5] =  1.4665E-04;
		//
		zbs[VmecIndataNamelist.ntord - 6][6] = -2.4682E-04;
		zbs[VmecIndataNamelist.ntord - 5][6] = -5.5887E-05;
		zbs[VmecIndataNamelist.ntord - 4][6] =  6.3329E-04;
		zbs[VmecIndataNamelist.ntord - 3][6] =  7.7064E-05;
		zbs[VmecIndataNamelist.ntord - 2][6] = -4.3674E-03;
		zbs[VmecIndataNamelist.ntord - 1][6] =  7.3160E-03;
		zbs[VmecIndataNamelist.ntord + 0][6] = -3.5003E-03;
		zbs[VmecIndataNamelist.ntord + 1][6] = -4.4304E-04;
		zbs[VmecIndataNamelist.ntord + 2][6] =  2.0925E-04;
		zbs[VmecIndataNamelist.ntord + 3][6] =  2.3550E-04;
		zbs[VmecIndataNamelist.ntord + 4][6] = -8.2796E-05;
		zbs[VmecIndataNamelist.ntord + 5][6] = -1.6727E-04;
		zbs[VmecIndataNamelist.ntord + 6][6] =  4.7081E-04;
		//
		zbs[VmecIndataNamelist.ntord - 6][7] =  1.6118E-05;
		zbs[VmecIndataNamelist.ntord - 5][7] =  2.6743E-04;
		zbs[VmecIndataNamelist.ntord - 4][7] = -7.3447E-04;
		zbs[VmecIndataNamelist.ntord - 3][7] =  1.1699E-04;
		zbs[VmecIndataNamelist.ntord - 2][7] = -2.3165E-03;
		zbs[VmecIndataNamelist.ntord - 1][7] = -1.3629E-03;
		zbs[VmecIndataNamelist.ntord + 0][7] =  9.8401E-04;
		zbs[VmecIndataNamelist.ntord + 1][7] = -3.7117E-04;
		zbs[VmecIndataNamelist.ntord + 2][7] =  2.8225E-04;
		zbs[VmecIndataNamelist.ntord + 3][7] =  1.3625E-04;
		zbs[VmecIndataNamelist.ntord + 4][7] = -3.7989E-04;
		zbs[VmecIndataNamelist.ntord + 5][7] = -3.2861E-05;
		zbs[VmecIndataNamelist.ntord + 6][7] =  4.2161E-05;
		//
		Assertions.assertArrayEquals(zbs, vmecInput.zbs);

		double[][] rbs = new double[2 * VmecIndataNamelist.ntord + 1][VmecIndataNamelist.mpold + 1];
		Assertions.assertArrayEquals(rbs, vmecInput.rbs);

		double[][] zbc = new double[2 * VmecIndataNamelist.ntord + 1][VmecIndataNamelist.mpold + 1];
		Assertions.assertArrayEquals(zbc, vmecInput.zbc);
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
		Assertions.assertEquals(false, vmecInput.lasym);
		Assertions.assertEquals(1, vmecInput.nfp);
		Assertions.assertEquals(5, vmecInput.mpol);
		Assertions.assertEquals(0, vmecInput.ntor);
		Assertions.assertEquals(0, vmecInput.ntheta);
		Assertions.assertEquals(1, vmecInput.nzeta);

		// --------------------
		// Multi-Grid Steps
		int[] ns_array = new int[100];
		ns_array[0] = 100;
		Assertions.assertArrayEquals(ns_array, vmecInput.ns_array);

		double[] ftol_array = new double[100];
		ftol_array[0] = 1.0e-20;
		Assertions.assertArrayEquals(ftol_array, vmecInput.ftol_array);

		int[] niter_array = new int[100];
		Arrays.fill(niter_array, 25000);
		Assertions.assertArrayEquals(niter_array, vmecInput.niter_array);

		// --------------------
		// Global Physics Parameters
		Assertions.assertEquals(1.54721E3, vmecInput.phiedge);
		Assertions.assertEquals(1, vmecInput.ncurr);

		// --------------------
		// Profile of Mass or Pressure
		Assertions.assertEquals("two_power", vmecInput.pmass_type.toLowerCase());
		double[] am = new double[21];
		am[0] = 1.0;
		am[1] = 5.0;
		am[2] = 10.0;
		Assertions.assertArrayEquals(am, vmecInput.am);
		Assertions.assertArrayEquals(new double[VmecIndataNamelist.nsd], vmecInput.am_aux_s);
		Assertions.assertArrayEquals(new double[VmecIndataNamelist.nsd], vmecInput.am_aux_f);
		Assertions.assertEquals(9.21730E5, vmecInput.pres_scale);
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
		Assertions.assertEquals("two_power", vmecInput.pcurr_type.toLowerCase());
		double[] ac = new double[21];
		ac[0] = 1.0;
		ac[1] = 5.0;
		ac[2] = 10.0;
		Assertions.assertArrayEquals(ac, vmecInput.ac);
		Assertions.assertArrayEquals(new double[VmecIndataNamelist.nsd], vmecInput.ac_aux_s);
		Assertions.assertArrayEquals(new double[VmecIndataNamelist.nsd], vmecInput.ac_aux_f);
		Assertions.assertEquals(18.595E6, vmecInput.curtor);
		Assertions.assertEquals(1.0, vmecInput.bloat);

		// --------------------
		// Free-Boundary Parameters
		Assertions.assertEquals(false, vmecInput.lfreeb);
		Assertions.assertEquals("NONE", vmecInput.mgrid_file);
		Assertions.assertArrayEquals(new double[VmecIndataNamelist.nigroup], vmecInput.extcur);
		Assertions.assertEquals(9, vmecInput.nvacskip);

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
		raxis_cc[0] = 6.21;
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
		rbc[VmecIndataNamelist.ntord + 0][0] = 6.21;
		rbc[VmecIndataNamelist.ntord + 0][1] = 2.00;
		rbc[VmecIndataNamelist.ntord + 0][2] = 0.4952;
		Assertions.assertArrayEquals(rbc, vmecInput.rbc);

		double[][] zbs = new double[2 * VmecIndataNamelist.ntord + 1][VmecIndataNamelist.mpold + 1];
		zbs[VmecIndataNamelist.ntord + 0][1] = 3.684;
		Assertions.assertArrayEquals(zbs, vmecInput.zbs);

		double[][] rbs = new double[2 * VmecIndataNamelist.ntord + 1][VmecIndataNamelist.mpold + 1];
		Assertions.assertArrayEquals(rbs, vmecInput.rbs);

		double[][] zbc = new double[2 * VmecIndataNamelist.ntord + 1][VmecIndataNamelist.mpold + 1];
		Assertions.assertArrayEquals(zbc, vmecInput.zbc);
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
		Assertions.assertEquals(false, vmecInput.lasym);
		Assertions.assertEquals(3, vmecInput.nfp);
		Assertions.assertEquals(9, vmecInput.mpol);
		Assertions.assertEquals(6, vmecInput.ntor);
		Assertions.assertEquals(0, vmecInput.ntheta);
		Assertions.assertEquals(32, vmecInput.nzeta);

		// --------------------
		// Multi-Grid Steps
		int[] ns_array = new int[100];
		ns_array[0] = 13;
		ns_array[1] = 25;
		ns_array[2] = 49;
		Assertions.assertArrayEquals(ns_array, vmecInput.ns_array);

		double[] ftol_array = new double[100];
		ftol_array[0] = 1.0e-6;
		ftol_array[1] = 1.0e-6;
		ftol_array[2] = 1.0e-11;
		Assertions.assertArrayEquals(ftol_array, vmecInput.ftol_array);

		int[] niter_array = new int[100];
		Arrays.fill(niter_array, 2500);
		Assertions.assertArrayEquals(niter_array, vmecInput.niter_array);

		// --------------------
		// Global Physics Parameters
		Assertions.assertEquals(5.07375027776403E-01, vmecInput.phiedge);
		Assertions.assertEquals(1, vmecInput.ncurr);

		// --------------------
		// Profile of Mass or Pressure
		Assertions.assertEquals("power_series", vmecInput.pmass_type.toLowerCase());
		double[] am = new double[21];
		am[0] = 5.37278520000000E+04;
		am[1] = -4.01304760000000E+03;
		am[2] = -2.83335960000000E+04;
		am[3] = -3.71706440000000E+05;
		am[4] = 1.40196850000000E+06;
		am[5] = -2.51988380000000E+06;
		am[6] = 2.10862500000000E+06;
		am[7] = -6.40367620000000E+05;
		Assertions.assertArrayEquals(am, vmecInput.am);
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
		ac[0] = 8.18395699999999E+03;
		ac[1] = 1.43603560000000E+06;
		ac[2] = -1.07407140000000E+07;
		ac[3] = 7.44389200000000E+07;
		ac[4] = -3.22215650000000E+08;
		ac[5] = 8.81050800000000E+08;
		ac[6] = -1.49389660000000E+09;
		ac[7] = 1.52746800000000E+09;
		ac[8] = -8.67901590000000E+08;
		ac[9] = 2.10351200000000E+08;
		Assertions.assertArrayEquals(ac, vmecInput.ac);
		Assertions.assertArrayEquals(new double[VmecIndataNamelist.nsd], vmecInput.ac_aux_s);
		Assertions.assertArrayEquals(new double[VmecIndataNamelist.nsd], vmecInput.ac_aux_f);
		Assertions.assertEquals(-1.74250000000000E+05, vmecInput.curtor);
		Assertions.assertEquals(1.0, vmecInput.bloat);

		// --------------------
		// Free-Boundary Parameters
		Assertions.assertEquals(true, vmecInput.lfreeb);
		Assertions.assertEquals("mgrid_c01r00_z16mod.nc", vmecInput.mgrid_file);
		double[] extcur = new double[VmecIndataNamelist.nigroup];
		extcur[0] = 6.96686783152366E+05;
		extcur[1] = 6.53525037003549E+05;
		extcur[2] = 5.66650763486760E+05;
		extcur[5] = 1.50284713000000E+05;
		extcur[6] = 1.29214167000000E+05;
		extcur[7] = 5.70890000000000E+04;
		extcur[8] = -3.87772300000000E+03;
		extcur[9] = 2.17549077362584E+04;
		Assertions.assertArrayEquals(extcur, vmecInput.extcur);
		Assertions.assertEquals(6, vmecInput.nvacskip);

		// --------------------
		// Tweaking Parameters
		Assertions.assertEquals(200, vmecInput.nstep);
		double[] aphi = new double[20];
		aphi[0] = 1.0;
		Assertions.assertArrayEquals(aphi, vmecInput.aphi);
		Assertions.assertEquals(0.9, vmecInput.delt);
		Assertions.assertEquals(1.0, vmecInput.tcon0);
		Assertions.assertEquals(false, vmecInput.lforbal);

		// --------------------
		// Initial Guess for Magnetic Axis Geometry
		double[] raxis_cc = new double[VmecIndataNamelist.ntord + 1];
		raxis_cc[0] = 1.45203311141139E+00;
		raxis_cc[1] = 1.08463357427697E-01;
		raxis_cc[2] = 6.77770122237979E-03;
		raxis_cc[3] = -5.75204471728962E-04;
		raxis_cc[4] = -2.02185338333559E-04;
		raxis_cc[5] = -6.47478490773586E-05;
		raxis_cc[6] = 1.32082168157713E-05;
		Assertions.assertArrayEquals(raxis_cc, vmecInput.raxis_cc);

		double[] zaxis_cs = new double[VmecIndataNamelist.ntord + 1];
		zaxis_cs[1] = -5.65210082402504E-02;
		zaxis_cs[2] = -4.28840244477725E-03;
		zaxis_cs[3] = 3.25741366836217E-04;
		zaxis_cs[4] = 1.33306425930255E-04;
		zaxis_cs[5] = -9.41641794761882E-06;
		zaxis_cs[6] = -1.44303491413310E-05;
		Assertions.assertArrayEquals(zaxis_cs, vmecInput.zaxis_cs);

		double[] raxis_cs = new double[VmecIndataNamelist.ntord + 1];
		Assertions.assertArrayEquals(raxis_cs, vmecInput.raxis_cs);

		double[] zaxis_cc = new double[VmecIndataNamelist.ntord + 1];
		Assertions.assertArrayEquals(zaxis_cc, vmecInput.zaxis_cc);

		// --------------------
		// (Initial Guess for) Boundary Geometry

		double[][] rbc = new double[2 * VmecIndataNamelist.ntord + 1][VmecIndataNamelist.mpold + 1];
		rbc[VmecIndataNamelist.ntord + 0][0] =  1.36857749802181E+00;
		rbc[VmecIndataNamelist.ntord + 1][0] =  2.41152702255788E-02;
		rbc[VmecIndataNamelist.ntord + 2][0] = -2.12317317987761E-03;
		rbc[VmecIndataNamelist.ntord + 3][0] =  2.83432164971144E-03;
		rbc[VmecIndataNamelist.ntord + 4][0] = -2.92618385211359E-04;
		rbc[VmecIndataNamelist.ntord + 5][0] = -2.56149208980095E-04;
		rbc[VmecIndataNamelist.ntord + 6][0] = -2.82914898677351E-05;
		//
		rbc[VmecIndataNamelist.ntord - 6][1] = -1.73121409572297E-04;
		rbc[VmecIndataNamelist.ntord - 5][1] =  4.69971527973080E-04;
		rbc[VmecIndataNamelist.ntord - 4][1] = -7.12763573181719E-04;
		rbc[VmecIndataNamelist.ntord - 3][1] =  1.99390198612186E-03;
		rbc[VmecIndataNamelist.ntord - 2][1] = -1.77451429017279E-02;
		rbc[VmecIndataNamelist.ntord - 1][1] =  6.61831005874843E-03;
		rbc[VmecIndataNamelist.ntord + 0][1] =  2.93167538410578E-01;
		rbc[VmecIndataNamelist.ntord + 1][1] = -1.06483996582863E-01;
		rbc[VmecIndataNamelist.ntord + 2][1] = -5.14109720700927E-03;
		rbc[VmecIndataNamelist.ntord + 3][1] = -2.90970702428460E-03;
		rbc[VmecIndataNamelist.ntord + 4][1] =  2.18054220267298E-04;
		rbc[VmecIndataNamelist.ntord + 5][1] =  2.89454010253765E-04;
		rbc[VmecIndataNamelist.ntord + 6][1] =  3.97531150948874E-04;
		//
		rbc[VmecIndataNamelist.ntord - 6][2] = -1.17125864636197E-04;
		rbc[VmecIndataNamelist.ntord - 5][2] =  9.01104603372894E-05;
		rbc[VmecIndataNamelist.ntord - 4][2] =  4.36688592773377E-04;
		rbc[VmecIndataNamelist.ntord - 3][2] = -1.20652288948884E-03;
		rbc[VmecIndataNamelist.ntord - 2][2] = -2.40419871025780E-03;
		rbc[VmecIndataNamelist.ntord - 1][2] = -3.46251579426286E-03;
		rbc[VmecIndataNamelist.ntord + 0][2] =  5.45513943746032E-02;
		rbc[VmecIndataNamelist.ntord + 1][2] =  7.42363199100095E-02;
		rbc[VmecIndataNamelist.ntord + 2][2] =  5.49637213388619E-02;
		rbc[VmecIndataNamelist.ntord + 3][2] =  5.86024948782749E-03;
		rbc[VmecIndataNamelist.ntord + 4][2] =  1.67226814070710E-03;
		rbc[VmecIndataNamelist.ntord + 5][2] = -2.68462946556597E-04;
		rbc[VmecIndataNamelist.ntord + 6][2] = -2.10121319304193E-04;
		//
		rbc[VmecIndataNamelist.ntord - 6][3] = -1.41998037840183E-05;
		rbc[VmecIndataNamelist.ntord - 5][3] =  6.10298767275015E-05;
		rbc[VmecIndataNamelist.ntord - 4][3] = -2.38014886237573E-04;
		rbc[VmecIndataNamelist.ntord - 3][3] = -5.50177908090369E-04;
		rbc[VmecIndataNamelist.ntord - 2][3] =  2.96079022665501E-04;
		rbc[VmecIndataNamelist.ntord - 1][3] =  8.20049883462086E-03;
		rbc[VmecIndataNamelist.ntord + 0][3] = -1.10296798720414E-02;
		rbc[VmecIndataNamelist.ntord + 1][3] = -1.39050720566585E-02;
		rbc[VmecIndataNamelist.ntord + 2][3] = -2.04914067892290E-02;
		rbc[VmecIndataNamelist.ntord + 3][3] = -1.08551820563319E-02;
		rbc[VmecIndataNamelist.ntord + 4][3] = -8.36533977325674E-05;
		rbc[VmecIndataNamelist.ntord + 5][3] = -2.84363473236722E-04;
		rbc[VmecIndataNamelist.ntord + 6][3] =  1.20525424025838E-04;
		//
		rbc[VmecIndataNamelist.ntord - 6][4] = -1.90101862806390E-05;
		rbc[VmecIndataNamelist.ntord - 5][4] =  1.35300770190253E-04;
		rbc[VmecIndataNamelist.ntord - 4][4] = -1.22743240916024E-04;
		rbc[VmecIndataNamelist.ntord - 3][4] = -3.69583324783397E-04;
		rbc[VmecIndataNamelist.ntord - 2][4] =  2.06425422451682E-03;
		rbc[VmecIndataNamelist.ntord - 1][4] = -3.14239053920561E-03;
		rbc[VmecIndataNamelist.ntord + 0][4] =  1.09761060295478E-03;
		rbc[VmecIndataNamelist.ntord + 1][4] = -7.43826891812780E-04;
		rbc[VmecIndataNamelist.ntord + 2][4] =  5.47006268566223E-03;
		rbc[VmecIndataNamelist.ntord + 3][4] =  1.51743579924205E-03;
		rbc[VmecIndataNamelist.ntord + 4][4] =  1.18078525294347E-03;
		rbc[VmecIndataNamelist.ntord + 5][4] = -5.71943411527878E-04;
		rbc[VmecIndataNamelist.ntord + 6][4] =  4.74694834534509E-05;
		//
		rbc[VmecIndataNamelist.ntord - 6][5] = -7.86166021172168E-07;
		rbc[VmecIndataNamelist.ntord - 5][5] =  3.89759685186833E-05;
		rbc[VmecIndataNamelist.ntord - 4][5] = -1.67460059931849E-04;
		rbc[VmecIndataNamelist.ntord - 3][5] =  5.56426951644379E-04;
		rbc[VmecIndataNamelist.ntord - 2][5] = -1.23972340797831E-03;
		rbc[VmecIndataNamelist.ntord - 1][5] =  1.67592872278822E-03;
		rbc[VmecIndataNamelist.ntord + 0][5] = -2.27531396605288E-03;
		rbc[VmecIndataNamelist.ntord + 1][5] =  3.05085683157440E-03;
		rbc[VmecIndataNamelist.ntord + 2][5] =  3.38544193108488E-03;
		rbc[VmecIndataNamelist.ntord + 3][5] =  1.09305938193071E-03;
		rbc[VmecIndataNamelist.ntord + 4][5] =  1.46840673941703E-03;
		rbc[VmecIndataNamelist.ntord + 5][5] =  1.30008648760872E-04;
		rbc[VmecIndataNamelist.ntord + 6][5] =  3.49473141596692E-04;
		//
		rbc[VmecIndataNamelist.ntord - 6][6] = -3.31048834134425E-05;
		rbc[VmecIndataNamelist.ntord - 5][6] =  4.99961503409912E-05;
		rbc[VmecIndataNamelist.ntord - 4][6] =  1.22368824831979E-04;
		rbc[VmecIndataNamelist.ntord - 3][6] = -4.68737889144547E-04;
		rbc[VmecIndataNamelist.ntord - 2][6] =  3.42867629753431E-04;
		rbc[VmecIndataNamelist.ntord - 1][6] = -4.15211024757546E-04;
		rbc[VmecIndataNamelist.ntord + 0][6] =  2.14693152888931E-03;
		rbc[VmecIndataNamelist.ntord + 1][6] = -1.23425698102826E-03;
		rbc[VmecIndataNamelist.ntord + 2][6] = -2.13050629026877E-03;
		rbc[VmecIndataNamelist.ntord + 3][6] = -2.05963196473025E-03;
		rbc[VmecIndataNamelist.ntord + 4][6] = -5.94445587874021E-04;
		rbc[VmecIndataNamelist.ntord + 5][6] = -4.63047073800230E-04;
		rbc[VmecIndataNamelist.ntord + 6][6] =  7.44887172064422E-05;
		//
		rbc[VmecIndataNamelist.ntord - 6][7] = -4.49438952219980E-05;
		rbc[VmecIndataNamelist.ntord - 5][7] =  1.08004528938234E-04;
		rbc[VmecIndataNamelist.ntord - 4][7] = -1.72791199664687E-04;
		rbc[VmecIndataNamelist.ntord - 3][7] =  5.14331239761795E-05;
		rbc[VmecIndataNamelist.ntord - 2][7] =  2.39510733146159E-05;
		rbc[VmecIndataNamelist.ntord - 1][7] =  4.41181654888970E-04;
		rbc[VmecIndataNamelist.ntord + 0][7] = -5.53912937885497E-04;
		rbc[VmecIndataNamelist.ntord + 1][7] =  1.53097959705758E-05;
		rbc[VmecIndataNamelist.ntord + 2][7] =  2.48441918459947E-04;
		rbc[VmecIndataNamelist.ntord + 3][7] = -1.74382686690892E-03;
		rbc[VmecIndataNamelist.ntord + 4][7] = -1.06587916422560E-03;
		rbc[VmecIndataNamelist.ntord + 5][7] = -9.91257607002373E-04;
		rbc[VmecIndataNamelist.ntord + 6][7] = -2.64094675646115E-04;
		//
		rbc[VmecIndataNamelist.ntord - 6][8] =  1.52352124103436E-05;
		rbc[VmecIndataNamelist.ntord - 5][8] = -2.73037734772156E-05;
		rbc[VmecIndataNamelist.ntord - 4][8] =  5.24404987591455E-06;
		rbc[VmecIndataNamelist.ntord - 3][8] = -1.55839428845508E-05;
		rbc[VmecIndataNamelist.ntord - 2][8] =  3.19126211925698E-04;
		rbc[VmecIndataNamelist.ntord - 1][8] = -4.26047111168129E-04;
		rbc[VmecIndataNamelist.ntord + 0][8] =  5.65716891618795E-05;
		rbc[VmecIndataNamelist.ntord + 1][8] =  4.80256924725329E-04;
		rbc[VmecIndataNamelist.ntord + 2][8] = -8.97562619683718E-04;
		rbc[VmecIndataNamelist.ntord + 3][8] =  3.54166404464476E-04;
		rbc[VmecIndataNamelist.ntord + 4][8] =  1.08210251854160E-03;
		rbc[VmecIndataNamelist.ntord + 5][8] =  4.77644916433333E-04;
		rbc[VmecIndataNamelist.ntord + 6][8] =  4.03892509040954E-04;
		//
		Assertions.assertArrayEquals(rbc, vmecInput.rbc);

		double[][] zbs = new double[2 * VmecIndataNamelist.ntord + 1][VmecIndataNamelist.mpold + 1];
		zbs[VmecIndataNamelist.ntord + 0][0] =  0.00000000000000E+00;
		zbs[VmecIndataNamelist.ntord + 1][0] = -1.95628026443856E-02;
		zbs[VmecIndataNamelist.ntord + 2][0] =  1.28757135233300E-02;
		zbs[VmecIndataNamelist.ntord + 3][0] = -2.44433180433167E-03;
		zbs[VmecIndataNamelist.ntord + 4][0] =  1.39033896396214E-03;
		zbs[VmecIndataNamelist.ntord + 5][0] =  1.43391023894103E-04;
		zbs[VmecIndataNamelist.ntord + 6][0] = -9.11088052887535E-05;
		//
		zbs[VmecIndataNamelist.ntord - 6][1] = -4.20225007239080E-04;
		zbs[VmecIndataNamelist.ntord - 5][1] =  2.46590162823638E-04;
		zbs[VmecIndataNamelist.ntord - 4][1] =  6.73892865947826E-04;
		zbs[VmecIndataNamelist.ntord - 3][1] =  6.15547224946056E-04;
		zbs[VmecIndataNamelist.ntord - 2][1] = -4.69515747168374E-03;
		zbs[VmecIndataNamelist.ntord - 1][1] = -4.35653052804666E-02;
		zbs[VmecIndataNamelist.ntord + 0][1] =  4.36403907496931E-01;
		zbs[VmecIndataNamelist.ntord + 1][1] =  2.04603176926176E-01;
		zbs[VmecIndataNamelist.ntord + 2][1] =  1.99515171615169E-02;
		zbs[VmecIndataNamelist.ntord + 3][1] =  6.12900186275585E-03;
		zbs[VmecIndataNamelist.ntord + 4][1] =  8.08506200033142E-04;
		zbs[VmecIndataNamelist.ntord + 5][1] = -9.02822885583646E-04;
		zbs[VmecIndataNamelist.ntord + 6][1] = -8.57716540177626E-05;
		//
		zbs[VmecIndataNamelist.ntord - 6][2] =  1.35440143932485E-04;
		zbs[VmecIndataNamelist.ntord - 5][2] =  7.37505527491143E-05;
		zbs[VmecIndataNamelist.ntord - 4][2] =  3.00176185750889E-04;
		zbs[VmecIndataNamelist.ntord - 3][2] = -7.84563881185602E-04;
		zbs[VmecIndataNamelist.ntord - 2][2] = -6.96389101314910E-03;
		zbs[VmecIndataNamelist.ntord - 1][2] =  1.90979930539848E-02;
		zbs[VmecIndataNamelist.ntord + 0][2] =  1.13493111759993E-02;
		zbs[VmecIndataNamelist.ntord + 1][2] = -5.30212500614690E-03;
		zbs[VmecIndataNamelist.ntord + 2][2] = -1.94987355684378E-02;
		zbs[VmecIndataNamelist.ntord + 3][2] = -1.03151084236557E-03;
		zbs[VmecIndataNamelist.ntord + 4][2] = -1.12859000704959E-03;
		zbs[VmecIndataNamelist.ntord + 5][2] = -9.18429803433583E-05;
		zbs[VmecIndataNamelist.ntord + 6][2] =  2.52007970188467E-04;
		//
		zbs[VmecIndataNamelist.ntord - 6][3] = -1.58128592678358E-04;
		zbs[VmecIndataNamelist.ntord - 5][3] =  6.76197245854339E-05;
		zbs[VmecIndataNamelist.ntord - 4][3] =  9.26279271028277E-05;
		zbs[VmecIndataNamelist.ntord - 3][3] = -9.86787242780619E-04;
		zbs[VmecIndataNamelist.ntord - 2][3] =  4.04129157771415E-03;
		zbs[VmecIndataNamelist.ntord - 1][3] = -6.65176024252472E-03;
		zbs[VmecIndataNamelist.ntord + 0][3] =  1.76716485147735E-03;
		zbs[VmecIndataNamelist.ntord + 1][3] =  1.30740675896453E-03;
		zbs[VmecIndataNamelist.ntord + 2][3] =  5.37562887686355E-03;
		zbs[VmecIndataNamelist.ntord + 3][3] =  6.74474773606641E-03;
		zbs[VmecIndataNamelist.ntord + 4][3] =  4.83119885960677E-04;
		zbs[VmecIndataNamelist.ntord + 5][3] =  2.83834239709773E-04;
		zbs[VmecIndataNamelist.ntord + 6][3] = -8.93028891847923E-05;
		//
		zbs[VmecIndataNamelist.ntord - 6][4] =  9.72586976270553E-05;
		zbs[VmecIndataNamelist.ntord - 5][4] = -5.43902374508031E-05;
		zbs[VmecIndataNamelist.ntord - 4][4] = -2.75779105622785E-04;
		zbs[VmecIndataNamelist.ntord - 3][4] =  5.15092368519431E-04;
		zbs[VmecIndataNamelist.ntord - 2][4] =  7.19495375166338E-05;
		zbs[VmecIndataNamelist.ntord - 1][4] = -6.82342815450098E-04;
		zbs[VmecIndataNamelist.ntord + 0][4] = -1.59508941399147E-03;
		zbs[VmecIndataNamelist.ntord + 1][4] =  2.92442691854345E-03;
		zbs[VmecIndataNamelist.ntord + 2][4] =  1.37433091928370E-03;
		zbs[VmecIndataNamelist.ntord + 3][4] = -1.11546938717248E-03;
		zbs[VmecIndataNamelist.ntord + 4][4] = -1.34099512042491E-03;
		zbs[VmecIndataNamelist.ntord + 5][4] =  1.11894928004839E-04;
		zbs[VmecIndataNamelist.ntord + 6][4] = -6.79281384985767E-05;
		//
		zbs[VmecIndataNamelist.ntord - 6][5] = -1.97583128167806E-05;
		zbs[VmecIndataNamelist.ntord - 5][5] = -2.35125450627463E-05;
		zbs[VmecIndataNamelist.ntord - 4][5] =  1.20443272625496E-04;
		zbs[VmecIndataNamelist.ntord - 3][5] =  1.89358353470569E-04;
		zbs[VmecIndataNamelist.ntord - 2][5] = -7.13920973881904E-04;
		zbs[VmecIndataNamelist.ntord - 1][5] =  3.09664459049274E-04;
		zbs[VmecIndataNamelist.ntord + 0][5] =  2.89872578851657E-04;
		zbs[VmecIndataNamelist.ntord + 1][5] =  1.70054415516227E-03;
		zbs[VmecIndataNamelist.ntord + 2][5] = -2.78554355922176E-04;
		zbs[VmecIndataNamelist.ntord + 3][5] =  2.17574640716384E-04;
		zbs[VmecIndataNamelist.ntord + 4][5] =  6.55897480907826E-04;
		zbs[VmecIndataNamelist.ntord + 5][5] =  3.28920149251306E-04;
		zbs[VmecIndataNamelist.ntord + 6][5] = -9.10050055565279E-05;
		//
		zbs[VmecIndataNamelist.ntord - 6][6] = -5.39233185118701E-05;
		zbs[VmecIndataNamelist.ntord - 5][6] =  8.18778365653110E-05;
		zbs[VmecIndataNamelist.ntord - 4][6] = -3.18694179250392E-05;
		zbs[VmecIndataNamelist.ntord - 3][6] = -1.19860793341686E-04;
		zbs[VmecIndataNamelist.ntord - 2][6] = -8.41992546728642E-05;
		zbs[VmecIndataNamelist.ntord - 1][6] =  4.62123282814622E-04;
		zbs[VmecIndataNamelist.ntord + 0][6] =  5.83170410220200E-04;
		zbs[VmecIndataNamelist.ntord + 1][6] = -1.24515973144034E-03;
		zbs[VmecIndataNamelist.ntord + 2][6] = -8.62757753151694E-04;
		zbs[VmecIndataNamelist.ntord + 3][6] = -1.08330956107557E-03;
		zbs[VmecIndataNamelist.ntord + 4][6] = -5.66437806990121E-04;
		zbs[VmecIndataNamelist.ntord + 5][6] = -4.68125310802505E-04;
		zbs[VmecIndataNamelist.ntord + 6][6] = -9.31140200352417E-05;
		//
		zbs[VmecIndataNamelist.ntord - 6][7] =  7.74632767666931E-05;
		zbs[VmecIndataNamelist.ntord - 5][7] =  9.03944020758115E-06;
		zbs[VmecIndataNamelist.ntord - 4][7] = -5.54927970109239E-05;
		zbs[VmecIndataNamelist.ntord - 3][7] = -5.02366241482580E-05;
		zbs[VmecIndataNamelist.ntord - 2][7] =  2.49466342908641E-05;
		zbs[VmecIndataNamelist.ntord - 1][7] =  5.18216789203044E-04;
		zbs[VmecIndataNamelist.ntord + 0][7] = -8.93894110075399E-04;
		zbs[VmecIndataNamelist.ntord + 1][7] = -3.63625710061707E-05;
		zbs[VmecIndataNamelist.ntord + 2][7] =  1.01543763330733E-04;
		zbs[VmecIndataNamelist.ntord + 3][7] =  1.11960034378585E-03;
		zbs[VmecIndataNamelist.ntord + 4][7] =  1.22624186244426E-03;
		zbs[VmecIndataNamelist.ntord + 5][7] =  4.65068227747870E-04;
		zbs[VmecIndataNamelist.ntord + 6][7] =  2.66027879388423E-04;
		//
		zbs[VmecIndataNamelist.ntord - 6][8] = -1.48940185824720E-06;
		zbs[VmecIndataNamelist.ntord - 5][8] =  3.00261934904714E-05;
		zbs[VmecIndataNamelist.ntord - 4][8] = -7.87559272120201E-05;
		zbs[VmecIndataNamelist.ntord - 3][8] =  1.60049638912638E-04;
		zbs[VmecIndataNamelist.ntord - 2][8] = -1.33145192383227E-04;
		zbs[VmecIndataNamelist.ntord - 1][8] = -2.11635932867683E-04;
		zbs[VmecIndataNamelist.ntord + 0][8] =  4.87630348200717E-04;
		zbs[VmecIndataNamelist.ntord + 1][8] = -8.01736119211754E-04;
		zbs[VmecIndataNamelist.ntord + 2][8] =  7.73008785638218E-04;
		zbs[VmecIndataNamelist.ntord + 3][8] =  1.05352811801637E-03;
		zbs[VmecIndataNamelist.ntord + 4][8] =  3.27365422328267E-04;
		zbs[VmecIndataNamelist.ntord + 5][8] =  1.87096527574704E-04;
		zbs[VmecIndataNamelist.ntord + 6][8] =  2.81418831862846E-04;
		//
		Assertions.assertArrayEquals(zbs, vmecInput.zbs);

		double[][] rbs = new double[2 * VmecIndataNamelist.ntord + 1][VmecIndataNamelist.mpold + 1];
		Assertions.assertArrayEquals(rbs, vmecInput.rbs);

		double[][] zbc = new double[2 * VmecIndataNamelist.ntord + 1][VmecIndataNamelist.mpold + 1];
		Assertions.assertArrayEquals(zbc, vmecInput.zbc);
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
		Assertions.assertEquals(false, vmecInput.lasym);
		Assertions.assertEquals(5, vmecInput.nfp);
		Assertions.assertEquals(9, vmecInput.mpol);
		Assertions.assertEquals(9, vmecInput.ntor);
		Assertions.assertEquals(31, vmecInput.ntheta);
		Assertions.assertEquals(36, vmecInput.nzeta);

		// --------------------
		// Multi-Grid Steps
		int[] ns_array = new int[100];
		ns_array[0] = 51;
		Assertions.assertArrayEquals(ns_array, vmecInput.ns_array);

		double[] ftol_array = new double[100];
		ftol_array[0] = 1.0e-19;
		Assertions.assertArrayEquals(ftol_array, vmecInput.ftol_array);

		int[] niter_array = new int[100];
		Arrays.fill(niter_array, 40000);
		Assertions.assertArrayEquals(niter_array, vmecInput.niter_array);

		// --------------------
		// Global Physics Parameters
		Assertions.assertEquals(1.7699999999999996, vmecInput.phiedge);
		Assertions.assertEquals(1, vmecInput.ncurr);

		// --------------------
		// Profile of Mass or Pressure
		Assertions.assertEquals("cubic_spline", vmecInput.pmass_type.toLowerCase());
		Assertions.assertArrayEquals(new double[21], vmecInput.am);
		double[] am_aux_s = new double[VmecIndataNamelist.nsd];
		am_aux_s[0] = 0.0;
		am_aux_s[1] = 0.05263157894736842;
		am_aux_s[2] = 0.10526315789473684;
		am_aux_s[3] = 0.15789473684210525;
		am_aux_s[4] = 0.21052631578947367;
		am_aux_s[5] = 0.2631578947368421;
		am_aux_s[6] = 0.3157894736842105;
		am_aux_s[7] = 0.3684210526315789;
		am_aux_s[8] = 0.42105263157894735;
		am_aux_s[9] = 0.47368421052631576;
		am_aux_s[10] = 0.5263157894736842;
		am_aux_s[11] = 0.5789473684210527;
		am_aux_s[12] = 0.631578947368421;
		am_aux_s[13] = 0.6842105263157894;
		am_aux_s[14] = 0.7368421052631579;
		am_aux_s[15] = 0.7894736842105263;
		am_aux_s[16] = 0.8421052631578947;
		am_aux_s[17] = 0.894736842105263;
		am_aux_s[18] = 0.9473684210526315;
		am_aux_s[19] = 1.0;
		Assertions.assertArrayEquals(am_aux_s, vmecInput.am_aux_s);
		double[] am_aux_f = new double[VmecIndataNamelist.nsd];
		am_aux_f[0] = 15000.0;
		am_aux_f[1] = 14909.166131292528;
		am_aux_f[2] = 14741.635846137355;
		am_aux_f[3] = 14521.820600835814;
		am_aux_f[4] = 14257.134278803771;
		am_aux_f[5] = 13950.829398004116;
		am_aux_f[6] = 13604.112173624924;
		am_aux_f[7] = 13216.837862521901;
		am_aux_f[8] = 12787.7527095645;
		am_aux_f[9] = 12314.517356476668;
		am_aux_f[10] = 11793.571927230603;
		am_aux_f[11] = 11219.830592877699;
		am_aux_f[12] = 10586.131895214381;
		am_aux_f[13] = 9882.27330756975;
		am_aux_f[14] = 9093.24651905931;
		am_aux_f[15] = 8195.74281770274;
		am_aux_f[16] = 7150.3162692801325;
		am_aux_f[17] = 5879.993807303139;
		am_aux_f[18] = 4186.566653391287;
		am_aux_f[19] = 0.0;
		Assertions.assertArrayEquals(am_aux_f, vmecInput.am_aux_f);
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
		Assertions.assertEquals("cubic_spline_i", vmecInput.pcurr_type.toLowerCase());
		Assertions.assertArrayEquals(new double[21], vmecInput.ac);
		double[] ac_aux_s = new double[VmecIndataNamelist.nsd];
		ac_aux_s[0] = 0.0;
		ac_aux_s[1] = 0.05263157894736842;
		ac_aux_s[2] = 0.10526315789473684;
		ac_aux_s[3] = 0.15789473684210525;
		ac_aux_s[4] = 0.21052631578947367;
		ac_aux_s[5] = 0.2631578947368421;
		ac_aux_s[6] = 0.3157894736842105;
		ac_aux_s[7] = 0.3684210526315789;
		ac_aux_s[8] = 0.42105263157894735;
		ac_aux_s[9] = 0.47368421052631576;
		ac_aux_s[10] = 0.5263157894736842;
		ac_aux_s[11] = 0.5789473684210527;
		ac_aux_s[12] = 0.631578947368421;
		ac_aux_s[13] = 0.6842105263157894;
		ac_aux_s[14] = 0.7368421052631579;
		ac_aux_s[15] = 0.7894736842105263;
		ac_aux_s[16] = 0.8421052631578947;
		ac_aux_s[17] = 0.894736842105263;
		ac_aux_s[18] = 0.9473684210526315;
		ac_aux_s[19] = 1.0;
		Assertions.assertArrayEquals(ac_aux_s, vmecInput.ac_aux_s);
		double[] ac_aux_f = new double[VmecIndataNamelist.nsd];
		ac_aux_f[0] = 0.0;                   ;
		ac_aux_f[1] = -2.2717840184115134E-17;
		ac_aux_f[2] = -9.762180051166474E-15;
		ac_aux_f[3] = -3.7509888840171997E-13;
		ac_aux_f[4] = -5.413128334700753E-12;
		ac_aux_f[5] = -4.6070250442265955E-11;
		ac_aux_f[6] = -2.8314965468191934E-10;
		ac_aux_f[7] = -1.4024542799398035E-9;
		ac_aux_f[8] = -5.9887936564489915E-9;
		ac_aux_f[9] = -2.3083631767847313E-8;
		ac_aux_f[10] = -8.3120038746978E-8;
		ac_aux_f[11] = -2.8768703205051676E-7;
		ac_aux_f[12] = -9.82864337149123E-7;
		ac_aux_f[13] = -3.4098943349774834E-6;
		ac_aux_f[14] = -1.24408071653565E-5;
		ac_aux_f[15] = -5.018748688406884E-5;
		ac_aux_f[16] = -2.4362756424749988E-4;
		ac_aux_f[17] = -0.001694662989805786;
		ac_aux_f[18] = -0.0278172191370382;
		ac_aux_f[19] = -46.450308147784405;
		Assertions.assertArrayEquals(ac_aux_f, vmecInput.ac_aux_f);
		Assertions.assertEquals(-1500.0, vmecInput.curtor);
		Assertions.assertEquals(1.0, vmecInput.bloat);

		// --------------------
		// Free-Boundary Parameters
		Assertions.assertEquals(true, vmecInput.lfreeb);
		Assertions.assertEquals("mgrid_w7x_nv36_hires.nc", vmecInput.mgrid_file);
		double[] extcur = new double[VmecIndataNamelist.nigroup];
		extcur[0] = -13068.0;
		extcur[1] = -13068.0;
		extcur[2] = -13067.0;
		extcur[3] = -13068.0;
		extcur[4] = -13067.0;
		extcur[5] = 698.0000000000001;
		extcur[6] = 699.0000000000001;
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
		raxis_cc[0] = 5.5617;
		raxis_cc[1] = 0.37114;
		Assertions.assertArrayEquals(raxis_cc, vmecInput.raxis_cc);

		double[] zaxis_cs = new double[VmecIndataNamelist.ntord + 1];
		zaxis_cs[0] = -0.0;
		zaxis_cs[1] = -0.30675;
		Assertions.assertArrayEquals(zaxis_cs, vmecInput.zaxis_cs);

		double[] raxis_cs = new double[VmecIndataNamelist.ntord + 1];
		Assertions.assertArrayEquals(raxis_cs, vmecInput.raxis_cs);

		double[] zaxis_cc = new double[VmecIndataNamelist.ntord + 1];
		Assertions.assertArrayEquals(zaxis_cc, vmecInput.zaxis_cc);

		// --------------------
		// (Initial Guess for) Boundary Geometry

		double[][] rbc = new double[2 * VmecIndataNamelist.ntord + 1][VmecIndataNamelist.mpold + 1];
		rbc[VmecIndataNamelist.ntord + 0][0] = 5.52;
		rbc[VmecIndataNamelist.ntord + 1][0] = 0.28143;
		rbc[VmecIndataNamelist.ntord + 2][0] = -0.0048433;
		rbc[VmecIndataNamelist.ntord + 3][0] = -9.1286E-4;
		rbc[VmecIndataNamelist.ntord + 4][0] = -0.0014998;
		rbc[VmecIndataNamelist.ntord + 5][0] = 4.9311E-5;
		rbc[VmecIndataNamelist.ntord + 6][0] = -2.4697E-4;
		//
		rbc[VmecIndataNamelist.ntord - 6][1] = -4.698E-4;
		rbc[VmecIndataNamelist.ntord - 5][1] = -1.2051E-5;
		rbc[VmecIndataNamelist.ntord - 4][1] = -8.6534E-4;
		rbc[VmecIndataNamelist.ntord - 3][1] = -0.001272;
		rbc[VmecIndataNamelist.ntord - 2][1] = 0.0060171;
		rbc[VmecIndataNamelist.ntord - 1][1] = 0.029677;
		rbc[VmecIndataNamelist.ntord + 0][1] = 0.48552;
		rbc[VmecIndataNamelist.ntord + 1][1] = -0.23402;
		rbc[VmecIndataNamelist.ntord + 2][1] = -0.0091058;
		rbc[VmecIndataNamelist.ntord + 3][1] = 0.003289;
		rbc[VmecIndataNamelist.ntord + 4][1] = 3.4635E-4;
		rbc[VmecIndataNamelist.ntord + 5][1] = 1.1576E-4;
		rbc[VmecIndataNamelist.ntord + 6][1] = -1.6641E-5;
		//
		rbc[VmecIndataNamelist.ntord - 6][2] = 4.6997E-4;
		rbc[VmecIndataNamelist.ntord - 5][2] = -2.8196E-4;
		rbc[VmecIndataNamelist.ntord - 4][2] = 3.1877E-5;
		rbc[VmecIndataNamelist.ntord - 3][2] = 1.0664E-4;
		rbc[VmecIndataNamelist.ntord - 2][2] = 0.0029138;
		rbc[VmecIndataNamelist.ntord - 1][2] = 0.012928;
		rbc[VmecIndataNamelist.ntord + 0][2] = 0.041011;
		rbc[VmecIndataNamelist.ntord + 1][2] = 0.044821;
		rbc[VmecIndataNamelist.ntord + 2][2] = 0.065192;
		rbc[VmecIndataNamelist.ntord + 3][2] = -0.0049404;
		rbc[VmecIndataNamelist.ntord + 4][2] = -1.8067E-4;
		rbc[VmecIndataNamelist.ntord + 5][2] = 3.1057E-4;
		rbc[VmecIndataNamelist.ntord + 6][2] = 6.1451E-6;
		//
		rbc[VmecIndataNamelist.ntord - 6][3] = -1.6953E-4;
		rbc[VmecIndataNamelist.ntord - 5][3] = 1.1373E-4;
		rbc[VmecIndataNamelist.ntord - 4][3] = 1.8704E-4;
		rbc[VmecIndataNamelist.ntord - 3][3] = -9.5654E-5;
		rbc[VmecIndataNamelist.ntord - 2][3] = 1.1619E-4;
		rbc[VmecIndataNamelist.ntord - 1][3] = 0.0010936;
		rbc[VmecIndataNamelist.ntord + 0][3] = -0.0036822;
		rbc[VmecIndataNamelist.ntord + 1][3] = -0.01239;
		rbc[VmecIndataNamelist.ntord + 2][3] = -0.020066;
		rbc[VmecIndataNamelist.ntord + 3][3] = -0.011376;
		rbc[VmecIndataNamelist.ntord + 4][3] = 0.0023417;
		rbc[VmecIndataNamelist.ntord + 5][3] = -2.4847E-4;
		rbc[VmecIndataNamelist.ntord + 6][3] = -1.7681E-5;
		//
		rbc[VmecIndataNamelist.ntord - 6][4] = 1.3002E-4;
		rbc[VmecIndataNamelist.ntord - 5][4] = -1.1019E-4;
		rbc[VmecIndataNamelist.ntord - 4][4] = 1.8622E-5;
		rbc[VmecIndataNamelist.ntord - 3][4] = 1.6933E-4;
		rbc[VmecIndataNamelist.ntord - 2][4] = 3.8313E-5;
		rbc[VmecIndataNamelist.ntord - 1][4] = -3.7628E-4;
		rbc[VmecIndataNamelist.ntord + 0][4] = 0.002604;
		rbc[VmecIndataNamelist.ntord + 1][4] = -1.612E-4;
		rbc[VmecIndataNamelist.ntord + 2][4] = 0.0085646;
		rbc[VmecIndataNamelist.ntord + 3][4] = 0.0020155;
		rbc[VmecIndataNamelist.ntord + 4][4] = 7.6437E-5;
		rbc[VmecIndataNamelist.ntord + 5][4] = -5.8086E-4;
		rbc[VmecIndataNamelist.ntord + 6][4] = 1.051E-4;
		//
		rbc[VmecIndataNamelist.ntord - 6][5] = -1.0675E-4;
		rbc[VmecIndataNamelist.ntord - 5][5] = -4.648E-5;
		rbc[VmecIndataNamelist.ntord - 4][5] = 1.6249E-4;
		rbc[VmecIndataNamelist.ntord - 3][5] = -9.6931E-5;
		rbc[VmecIndataNamelist.ntord - 2][5] = -1.4994E-5;
		rbc[VmecIndataNamelist.ntord - 1][5] = 7.3516E-4;
		rbc[VmecIndataNamelist.ntord + 0][5] = 5.5656E-4;
		rbc[VmecIndataNamelist.ntord + 1][5] = 0.0019807;
		rbc[VmecIndataNamelist.ntord + 2][5] = -6.7779E-4;
		rbc[VmecIndataNamelist.ntord + 3][5] = -7.0688E-4;
		rbc[VmecIndataNamelist.ntord + 4][5] = 0.0012466;
		rbc[VmecIndataNamelist.ntord + 5][5] = -1.2292E-5;
		rbc[VmecIndataNamelist.ntord + 6][5] = 3.4932E-5;
		//
		rbc[VmecIndataNamelist.ntord - 6][6] = 5.2421E-5;
		rbc[VmecIndataNamelist.ntord - 5][6] = -3.9072E-5;
		rbc[VmecIndataNamelist.ntord - 4][6] = -1.7258E-5;
		rbc[VmecIndataNamelist.ntord - 3][6] = 7.39E-5;
		rbc[VmecIndataNamelist.ntord - 2][6] = 9.1887E-5;
		rbc[VmecIndataNamelist.ntord - 1][6] = -5.2882E-5;
		rbc[VmecIndataNamelist.ntord + 0][6] = -6.49E-5;
		rbc[VmecIndataNamelist.ntord + 1][6] = -0.0018093;
		rbc[VmecIndataNamelist.ntord + 2][6] = -0.0016026;
		rbc[VmecIndataNamelist.ntord + 3][6] = -0.0012294;
		rbc[VmecIndataNamelist.ntord + 4][6] = -7.2298E-4;
		rbc[VmecIndataNamelist.ntord + 5][6] = -4.6175E-4;
		rbc[VmecIndataNamelist.ntord + 6][6] = 2.1078E-4;
		//
		Assertions.assertArrayEquals(rbc, vmecInput.rbc);

		double[][] zbs = new double[2 * VmecIndataNamelist.ntord + 1][VmecIndataNamelist.mpold + 1];
		zbs[VmecIndataNamelist.ntord + 0][0] = -0.0;
		zbs[VmecIndataNamelist.ntord + 1][0] = -0.23796;
		zbs[VmecIndataNamelist.ntord + 2][0] = 8.1503E-4;
		zbs[VmecIndataNamelist.ntord + 3][0] = 0.0026029;
		zbs[VmecIndataNamelist.ntord + 4][0] = 0.0015708;
		zbs[VmecIndataNamelist.ntord + 5][0] = 1.8537E-4;
		zbs[VmecIndataNamelist.ntord + 6][0] = 2.1705E-4;
		//
		zbs[VmecIndataNamelist.ntord - 6][1] = 7.9072E-4;
		zbs[VmecIndataNamelist.ntord - 5][1] = -9.8388E-5;
		zbs[VmecIndataNamelist.ntord - 4][1] = -6.0039E-5;
		zbs[VmecIndataNamelist.ntord - 3][1] = -4.9979E-5;
		zbs[VmecIndataNamelist.ntord - 2][1] = 0.010386;
		zbs[VmecIndataNamelist.ntord - 1][1] = 0.032068;
		zbs[VmecIndataNamelist.ntord + 0][1] = 0.62135;
		zbs[VmecIndataNamelist.ntord + 1][1] = 0.1838;
		zbs[VmecIndataNamelist.ntord + 2][1] = 0.013137;
		zbs[VmecIndataNamelist.ntord + 3][1] = 3.109E-4;
		zbs[VmecIndataNamelist.ntord + 4][1] = -7.1302E-4;
		zbs[VmecIndataNamelist.ntord + 5][1] = -1.5327E-4;
		zbs[VmecIndataNamelist.ntord + 6][1] = 1.4661E-4;
		//
		zbs[VmecIndataNamelist.ntord - 6][2] = -6.2128E-4;
		zbs[VmecIndataNamelist.ntord - 5][2] = 2.4333E-4;
		zbs[VmecIndataNamelist.ntord - 4][2] = 6.5569E-5;
		zbs[VmecIndataNamelist.ntord - 3][2] = -1.2046E-5;
		zbs[VmecIndataNamelist.ntord - 2][2] = 6.5241E-4;
		zbs[VmecIndataNamelist.ntord - 1][2] = 0.0089215;
		zbs[VmecIndataNamelist.ntord + 0][2] = -0.0054534;
		zbs[VmecIndataNamelist.ntord + 1][2] = 0.018121;
		zbs[VmecIndataNamelist.ntord + 2][2] = -0.052132;
		zbs[VmecIndataNamelist.ntord + 3][2] = 0.0060893;
		zbs[VmecIndataNamelist.ntord + 4][2] = 5.0937E-4;
		zbs[VmecIndataNamelist.ntord + 5][2] = -6.9752E-5;
		zbs[VmecIndataNamelist.ntord + 6][2] = 3.4368E-5;
		//
		zbs[VmecIndataNamelist.ntord - 6][3] = 2.5083E-4;
		zbs[VmecIndataNamelist.ntord - 5][3] = -6.9072E-5;
		zbs[VmecIndataNamelist.ntord - 4][3] = -2.2969E-4;
		zbs[VmecIndataNamelist.ntord - 3][3] = 7.4446E-5;
		zbs[VmecIndataNamelist.ntord - 2][3] = 5.809E-4;
		zbs[VmecIndataNamelist.ntord - 1][3] = -0.0010221;
		zbs[VmecIndataNamelist.ntord + 0][3] = -0.0011369;
		zbs[VmecIndataNamelist.ntord + 1][3] = -0.0043846;
		zbs[VmecIndataNamelist.ntord + 2][3] = 0.0080932;
		zbs[VmecIndataNamelist.ntord + 3][3] = 0.010319;
		zbs[VmecIndataNamelist.ntord + 4][3] = -0.0019609;
		zbs[VmecIndataNamelist.ntord + 5][3] = 1.9694E-4;
		zbs[VmecIndataNamelist.ntord + 6][3] = 4.9105E-5;
		//
		zbs[VmecIndataNamelist.ntord - 6][4] = -9.394E-5;
		zbs[VmecIndataNamelist.ntord - 5][4] = 4.5925E-5;
		zbs[VmecIndataNamelist.ntord - 4][4] = 7.2163E-5;
		zbs[VmecIndataNamelist.ntord - 3][4] = -1.6596E-4;
		zbs[VmecIndataNamelist.ntord - 2][4] = 4.5852E-5;
		zbs[VmecIndataNamelist.ntord - 1][4] = 4.7079E-5;
		zbs[VmecIndataNamelist.ntord + 0][4] = 0.001406;
		zbs[VmecIndataNamelist.ntord + 1][4] = 0.0018031;
		zbs[VmecIndataNamelist.ntord + 2][4] = 0.0082571;
		zbs[VmecIndataNamelist.ntord + 3][4] = -0.0058909;
		zbs[VmecIndataNamelist.ntord + 4][4] = -7.0749E-4;
		zbs[VmecIndataNamelist.ntord + 5][4] = 7.0162E-4;
		zbs[VmecIndataNamelist.ntord + 6][4] = -1.4223E-4;
		//
		zbs[VmecIndataNamelist.ntord - 6][5] = -7.7646E-5;
		zbs[VmecIndataNamelist.ntord - 5][5] = 2.1039E-4;
		zbs[VmecIndataNamelist.ntord - 4][5] = -1.2919E-4;
		zbs[VmecIndataNamelist.ntord - 3][5] = 3.8293E-5;
		zbs[VmecIndataNamelist.ntord - 2][5] = -1.348E-5;
		zbs[VmecIndataNamelist.ntord - 1][5] = 6.1354E-4;
		zbs[VmecIndataNamelist.ntord + 0][5] = 9.2855E-4;
		zbs[VmecIndataNamelist.ntord + 1][5] = 3.3226E-5;
		zbs[VmecIndataNamelist.ntord + 2][5] = -0.0026372;
		zbs[VmecIndataNamelist.ntord + 3][5] = -7.4071E-5;
		zbs[VmecIndataNamelist.ntord + 4][5] = 0.0011784;
		zbs[VmecIndataNamelist.ntord + 5][5] = -2.2804E-4;
		zbs[VmecIndataNamelist.ntord + 6][5] = -1.053E-4;
		//
		zbs[VmecIndataNamelist.ntord - 6][6] = 1.7109E-5;
		zbs[VmecIndataNamelist.ntord - 5][6] = -5.9836E-5;
		zbs[VmecIndataNamelist.ntord - 4][6] = 3.0786E-5;
		zbs[VmecIndataNamelist.ntord - 3][6] = 9.3985E-5;
		zbs[VmecIndataNamelist.ntord - 2][6] = -3.6056E-5;
		zbs[VmecIndataNamelist.ntord - 1][6] = 1.2491E-5;
		zbs[VmecIndataNamelist.ntord + 0][6] = -2.3626E-4;
		zbs[VmecIndataNamelist.ntord + 1][6] = -8.2965E-4;
		zbs[VmecIndataNamelist.ntord + 2][6] = -0.0011174;
		zbs[VmecIndataNamelist.ntord + 3][6] = -8.229E-4;
		zbs[VmecIndataNamelist.ntord + 4][6] = -3.8131E-4;
		zbs[VmecIndataNamelist.ntord + 5][6] = 1.9311E-4;
		zbs[VmecIndataNamelist.ntord + 6][6] = 9.7126E-7;
		//
		Assertions.assertArrayEquals(zbs, vmecInput.zbs);

		double[][] rbs = new double[2 * VmecIndataNamelist.ntord + 1][VmecIndataNamelist.mpold + 1];
		Assertions.assertArrayEquals(rbs, vmecInput.rbs);

		double[][] zbc = new double[2 * VmecIndataNamelist.ntord + 1][VmecIndataNamelist.mpold + 1];
		Assertions.assertArrayEquals(zbc, vmecInput.zbc);
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
		Assertions.assertEquals(true, vmecInput.lasym);
		Assertions.assertEquals(5, vmecInput.nfp);
		Assertions.assertEquals(5, vmecInput.mpol);
		Assertions.assertEquals(4, vmecInput.ntor);
		Assertions.assertEquals(0, vmecInput.ntheta);
		Assertions.assertEquals(36, vmecInput.nzeta);

		// --------------------
		// Multi-Grid Steps
		int[] ns_array = new int[100];
		ns_array[0] = 15;
		Assertions.assertArrayEquals(ns_array, vmecInput.ns_array);

		double[] ftol_array = new double[100];
		ftol_array[0] = 1.0e-20;
		Assertions.assertArrayEquals(ftol_array, vmecInput.ftol_array);

		int[] niter_array = new int[100];
		Arrays.fill(niter_array, -1);
		niter_array[0] = 52;
		Assertions.assertArrayEquals(niter_array, vmecInput.niter_array);

		// --------------------
		// Global Physics Parameters
		Assertions.assertEquals(-0.035, vmecInput.phiedge);
		Assertions.assertEquals(1, vmecInput.ncurr);

		// --------------------
		// Profile of Mass or Pressure
		Assertions.assertEquals("two_power", vmecInput.pmass_type.toLowerCase());
		double[] am = new double[21];
		am[0] = 1.0;
		am[1] = 5.0;
		am[2] = 10.0;
		Assertions.assertArrayEquals(am, vmecInput.am);
		Assertions.assertArrayEquals(new double[VmecIndataNamelist.nsd], vmecInput.am_aux_s);
		Assertions.assertArrayEquals(new double[VmecIndataNamelist.nsd], vmecInput.am_aux_f);
		Assertions.assertEquals(432.29080924603676, vmecInput.pres_scale);
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
		Assertions.assertEquals("two_power", vmecInput.pcurr_type.toLowerCase());
		double[] ac = new double[21];
		ac[0] = 1.0;
		ac[1] = 5.0;
		ac[2] = 10.0;
		Assertions.assertArrayEquals(ac, vmecInput.ac);
		Assertions.assertArrayEquals(new double[VmecIndataNamelist.nsd], vmecInput.ac_aux_s);
		Assertions.assertArrayEquals(new double[VmecIndataNamelist.nsd], vmecInput.ac_aux_f);
		Assertions.assertEquals(43229.08092460368, vmecInput.curtor);
		Assertions.assertEquals(1.0, vmecInput.bloat);

		// --------------------
		// Free-Boundary Parameters
		Assertions.assertEquals(true, vmecInput.lfreeb);
		Assertions.assertEquals("mgrid_test.nc", vmecInput.mgrid_file);
		double[] extcur = new double[VmecIndataNamelist.nigroup];
		extcur[0] = 4700.0;
		extcur[1] = 1000.0;
		Assertions.assertArrayEquals(extcur, vmecInput.extcur);
		Assertions.assertEquals(9, vmecInput.nvacskip);

		// --------------------
		// Tweaking Parameters
		Assertions.assertEquals(200, vmecInput.nstep);
		double[] aphi = new double[20];
		aphi[0] = 1.0;
		Assertions.assertArrayEquals(aphi, vmecInput.aphi);
		Assertions.assertEquals(0.7, vmecInput.delt);
		Assertions.assertEquals(1.0, vmecInput.tcon0);
		Assertions.assertEquals(false, vmecInput.lforbal);

		// --------------------
		// Initial Guess for Magnetic Axis Geometry
		double[] raxis_cc = new double[VmecIndataNamelist.ntord + 1];
		raxis_cc[0] = 0.786037734951267;
		raxis_cc[1] = -0.0302978726119071;
		raxis_cc[2] = 0.000915258048528683;
		raxis_cc[3] = -1.91959906744039e-05;
		raxis_cc[4] = 1.45930777845745e-06;
		Assertions.assertArrayEquals(raxis_cc, vmecInput.raxis_cc);

		double[] zaxis_cs = new double[VmecIndataNamelist.ntord + 1];
		zaxis_cs[0] = 0.0;
		zaxis_cs[1] = 0.0273158409510113;
		zaxis_cs[2] = -0.000937096584979097;
		zaxis_cs[3] = 1.741833421328e-05;
		zaxis_cs[4] = -5.91222841432118e-07;
		Assertions.assertArrayEquals(zaxis_cs, vmecInput.zaxis_cs);

		double[] raxis_cs = new double[VmecIndataNamelist.ntord + 1];
		Assertions.assertArrayEquals(raxis_cs, vmecInput.raxis_cs);

		double[] zaxis_cc = new double[VmecIndataNamelist.ntord + 1];
		Assertions.assertArrayEquals(zaxis_cc, vmecInput.zaxis_cc);

		// --------------------
		// (Initial Guess for) Boundary Geometry

		double[][] rbc = new double[2 * VmecIndataNamelist.ntord + 1][VmecIndataNamelist.mpold + 1];
		rbc[VmecIndataNamelist.ntord + 0][0] = 0.780906309727434;
		rbc[VmecIndataNamelist.ntord + 1][0] = -0.046151739531816;
		rbc[VmecIndataNamelist.ntord + 2][0] = 0.0025346622908875;
		rbc[VmecIndataNamelist.ntord + 3][0] = -0.000294570027126834;
		rbc[VmecIndataNamelist.ntord + 4][0] = 4.87628905278232e-05;
		//
		rbc[VmecIndataNamelist.ntord - 4][1] = 1.10991985982551e-05;
		rbc[VmecIndataNamelist.ntord - 3][1] = -6.69760263500533e-05;
		rbc[VmecIndataNamelist.ntord - 2][1] = 0.000344811566716634;
		rbc[VmecIndataNamelist.ntord - 1][1] = -0.00566091264003757;
		rbc[VmecIndataNamelist.ntord + 0][1] = 0.131435746235883;
		rbc[VmecIndataNamelist.ntord + 1][1] = -0.039386399008392;
		rbc[VmecIndataNamelist.ntord + 2][1] = 0.00488869740987901;
		rbc[VmecIndataNamelist.ntord + 3][1] = -0.000438124054375238;
		rbc[VmecIndataNamelist.ntord + 4][1] = 4.81645985535357e-05;
		//
		rbc[VmecIndataNamelist.ntord - 4][2] = -1.53011178035308e-06;
		rbc[VmecIndataNamelist.ntord - 3][2] = 2.35413167655254e-05;
		rbc[VmecIndataNamelist.ntord - 2][2] = -8.83026333729552e-05;
		rbc[VmecIndataNamelist.ntord - 1][2] = 0.000146058452530242;
		rbc[VmecIndataNamelist.ntord + 0][2] = 0.00147170066093047;
		rbc[VmecIndataNamelist.ntord + 1][2] = -0.00320611057583166;
		rbc[VmecIndataNamelist.ntord + 2][2] = 0.00228943846580024;
		rbc[VmecIndataNamelist.ntord + 3][2] = -0.000319173467542497;
		rbc[VmecIndataNamelist.ntord + 4][2] = 2.60138638307058e-05;
		//
		rbc[VmecIndataNamelist.ntord - 4][3] = -2.22269133805289e-06;
		rbc[VmecIndataNamelist.ntord - 3][3] = 8.93182437579942e-07;
		rbc[VmecIndataNamelist.ntord - 2][3] = -1.86926512003518e-05;
		rbc[VmecIndataNamelist.ntord - 1][3] = -0.000103222718482529;
		rbc[VmecIndataNamelist.ntord + 0][3] = -0.00010590693425821;
		rbc[VmecIndataNamelist.ntord + 1][3] = -3.70460344545912e-05;
		rbc[VmecIndataNamelist.ntord + 2][3] = 0.000384952481942609;
		rbc[VmecIndataNamelist.ntord + 3][3] = -8.73464736915271e-05;
		rbc[VmecIndataNamelist.ntord + 4][3] = 1.70559212798558e-05;
		//
		rbc[VmecIndataNamelist.ntord - 4][4] = -4.11611838328142e-06;
		rbc[VmecIndataNamelist.ntord - 3][4] = -3.670147571078e-06;
		rbc[VmecIndataNamelist.ntord - 2][4] = 3.27663252086114e-06;
		rbc[VmecIndataNamelist.ntord - 1][4] = -5.49897884990649e-05;
		rbc[VmecIndataNamelist.ntord + 0][4] = 0.000611036997986113;
		rbc[VmecIndataNamelist.ntord + 1][4] = 0.000157752536759372;
		rbc[VmecIndataNamelist.ntord + 2][4] = 0.000238939203660074;
		rbc[VmecIndataNamelist.ntord + 3][4] = -7.45864378517666e-05;
		rbc[VmecIndataNamelist.ntord + 4][4] = 1.0767515627784e-05;
		//
		Assertions.assertArrayEquals(rbc, vmecInput.rbc);

		double[][] zbs = new double[2 * VmecIndataNamelist.ntord + 1][VmecIndataNamelist.mpold + 1];
		zbs[VmecIndataNamelist.ntord + 0][0] = 0.0;
		zbs[VmecIndataNamelist.ntord + 1][0] = 0.0449223151020507;
		zbs[VmecIndataNamelist.ntord + 2][0] = -0.00157959039976886;
		zbs[VmecIndataNamelist.ntord + 3][0] = 0.000110096448338774;
		zbs[VmecIndataNamelist.ntord + 4][0] = -3.41523759959853e-05;
		//
		zbs[VmecIndataNamelist.ntord - 4][1] = 8.86243753757977e-06;
		zbs[VmecIndataNamelist.ntord - 3][1] = -7.61994068511444e-06;
		zbs[VmecIndataNamelist.ntord - 2][1] = -0.000102387369227698;
		zbs[VmecIndataNamelist.ntord - 1][1] = -0.00308184739411174;
		zbs[VmecIndataNamelist.ntord + 0][1] = 0.164601748843378;
		zbs[VmecIndataNamelist.ntord + 1][1] = 0.0303393558516301;
		zbs[VmecIndataNamelist.ntord + 2][1] = -0.00451585182685075;
		zbs[VmecIndataNamelist.ntord + 3][1] = 0.000322392099487897;
		zbs[VmecIndataNamelist.ntord + 4][1] = -3.48737658782808e-05;
		//
		zbs[VmecIndataNamelist.ntord - 4][2] = -1.28442261664991e-06;
		zbs[VmecIndataNamelist.ntord - 3][2] = 2.93065492924715e-05;
		zbs[VmecIndataNamelist.ntord - 2][2] = -0.000193078473966611;
		zbs[VmecIndataNamelist.ntord - 1][2] = 0.000917965755419898;
		zbs[VmecIndataNamelist.ntord + 0][2] = -0.00272248379340565;
		zbs[VmecIndataNamelist.ntord + 1][2] = 0.00355657869194155;
		zbs[VmecIndataNamelist.ntord + 2][2] = -0.00205297288922009;
		zbs[VmecIndataNamelist.ntord + 3][2] = 0.000328755917997435;
		zbs[VmecIndataNamelist.ntord + 4][2] = -3.57925786395997e-05;
		//
		zbs[VmecIndataNamelist.ntord - 4][3] = 3.83110072934457e-06;
		zbs[VmecIndataNamelist.ntord - 3][3] = -3.24731939556698e-06;
		zbs[VmecIndataNamelist.ntord - 2][3] = -3.77850265219013e-06;
		zbs[VmecIndataNamelist.ntord - 1][3] = -7.43052945541678e-07;
		zbs[VmecIndataNamelist.ntord + 0][3] = 0.000120213894931637;
		zbs[VmecIndataNamelist.ntord + 1][3] = 9.70588501764729e-05;
		zbs[VmecIndataNamelist.ntord + 2][3] = -0.000305563766958157;
		zbs[VmecIndataNamelist.ntord + 3][3] = 7.12819152728525e-05;
		zbs[VmecIndataNamelist.ntord + 4][3] = -1.14192148882868e-05;
		//
		zbs[VmecIndataNamelist.ntord - 4][4] = 2.79111955871817e-06;
		zbs[VmecIndataNamelist.ntord - 3][4] = -1.0123457043081e-06;
		zbs[VmecIndataNamelist.ntord - 2][4] = 1.33742394885452e-05;
		zbs[VmecIndataNamelist.ntord - 1][4] = -6.56302557453055e-05;
		zbs[VmecIndataNamelist.ntord + 0][4] = 0.000408527289143525;
		zbs[VmecIndataNamelist.ntord + 1][4] = -0.00113568343044368;
		zbs[VmecIndataNamelist.ntord + 2][4] = 0.000453498195180046;
		zbs[VmecIndataNamelist.ntord + 3][4] = -5.07003662052117e-05;
		zbs[VmecIndataNamelist.ntord + 4][4] = 5.22398328609138e-06;
		//
		Assertions.assertArrayEquals(zbs, vmecInput.zbs);

		double[][] rbs = new double[2 * VmecIndataNamelist.ntord + 1][VmecIndataNamelist.mpold + 1];
		Assertions.assertArrayEquals(rbs, vmecInput.rbs);

		double[][] zbc = new double[2 * VmecIndataNamelist.ntord + 1][VmecIndataNamelist.mpold + 1];
		Assertions.assertArrayEquals(zbc, vmecInput.zbc);
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
		Assertions.assertEquals(false, vmecInput.lasym);
		Assertions.assertEquals(5, vmecInput.nfp);
		Assertions.assertEquals(12, vmecInput.mpol);
		Assertions.assertEquals(12, vmecInput.ntor);
		Assertions.assertEquals(0, vmecInput.ntheta);
		Assertions.assertEquals(36, vmecInput.nzeta);

		// --------------------
		// Multi-Grid Steps
		int[] ns_array = new int[100];
		ns_array[0] = 4;
		ns_array[1] = 9;
		ns_array[2] = 28;
		ns_array[3] = 99;
		Assertions.assertArrayEquals(ns_array, vmecInput.ns_array);

		double[] ftol_array = new double[100];
		ftol_array[0] = 1.0e-5;
		ftol_array[1] = 1.0e-7;
		ftol_array[2] = 1.0e-8;
		ftol_array[3] = 1.0e-14;
		Assertions.assertArrayEquals(ftol_array, vmecInput.ftol_array);

		int[] niter_array = new int[100];
		Arrays.fill(niter_array, 8000);
		Assertions.assertArrayEquals(niter_array, vmecInput.niter_array);

		// --------------------
		// Global Physics Parameters
		Assertions.assertEquals(-1.29, vmecInput.phiedge);
		Assertions.assertEquals(1, vmecInput.ncurr);

		// --------------------
		// Profile of Mass or Pressure
		Assertions.assertEquals("power_series", vmecInput.pmass_type.toLowerCase());
		double[] am = new double[21];
		am[0] = 1.0e-6;
		am[1] = -1.0e-6;
		Assertions.assertArrayEquals(am, vmecInput.am);
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
		Assertions.assertArrayEquals(new double[21], vmecInput.ac);
		Assertions.assertArrayEquals(new double[VmecIndataNamelist.nsd], vmecInput.ac_aux_s);
		Assertions.assertArrayEquals(new double[VmecIndataNamelist.nsd], vmecInput.ac_aux_f);
		Assertions.assertEquals(0.0, vmecInput.curtor);
		Assertions.assertEquals(1.0, vmecInput.bloat);

		// --------------------
		// Free-Boundary Parameters
		Assertions.assertEquals(true, vmecInput.lfreeb);
		Assertions.assertEquals("mgrid_w7x_nv36_hires.nc", vmecInput.mgrid_file);
		double[] extcur = new double[VmecIndataNamelist.nigroup];
		extcur[0] = 14000.0;
		extcur[1] = 14000.0;
		extcur[2] = 13160.0;
		extcur[3] = 12950.0;
		extcur[4] = 12390.0;
		extcur[5] = -9660.0;
		extcur[6] = -9660.0;
		Assertions.assertArrayEquals(extcur, vmecInput.extcur);
		Assertions.assertEquals(6, vmecInput.nvacskip);

		// --------------------
		// Tweaking Parameters
		Assertions.assertEquals(100, vmecInput.nstep);
		double[] aphi = new double[20];
		aphi[0] = 1.0;
		Assertions.assertArrayEquals(aphi, vmecInput.aphi);
		Assertions.assertEquals(0.9, vmecInput.delt);
		Assertions.assertEquals(1.0, vmecInput.tcon0);
		Assertions.assertEquals(false, vmecInput.lforbal);

		// --------------------
		// Initial Guess for Magnetic Axis Geometry
		double[] raxis_cc = new double[VmecIndataNamelist.ntord + 1];
		raxis_cc[0] = 5.5461E+00;
		raxis_cc[1] = 3.8004E-01;
		Assertions.assertArrayEquals(raxis_cc, vmecInput.raxis_cc);

		double[] zaxis_cs = new double[VmecIndataNamelist.ntord + 1];
		zaxis_cs[0] = 0.0;
		zaxis_cs[1] = -3.0390E-01;
		Assertions.assertArrayEquals(zaxis_cs, vmecInput.zaxis_cs);

		double[] raxis_cs = new double[VmecIndataNamelist.ntord + 1];
		Assertions.assertArrayEquals(raxis_cs, vmecInput.raxis_cs);

		double[] zaxis_cc = new double[VmecIndataNamelist.ntord + 1];
		Assertions.assertArrayEquals(zaxis_cc, vmecInput.zaxis_cc);

		// --------------------
		// (Initial Guess for) Boundary Geometry

		double[][] rbc = new double[2 * VmecIndataNamelist.ntord + 1][VmecIndataNamelist.mpold + 1];
		rbc[VmecIndataNamelist.ntord + 0][0] = 5.5289e+00;
		rbc[VmecIndataNamelist.ntord + 1][0] = 2.6718e-01;
		rbc[VmecIndataNamelist.ntord + 2][0] = -2.1739e-03;
		rbc[VmecIndataNamelist.ntord + 3][0] = -2.8507e-04;
		rbc[VmecIndataNamelist.ntord + 4][0] = -1.8275e-03;
		//
		rbc[VmecIndataNamelist.ntord - 4][1] = 1.6341e-03;
		rbc[VmecIndataNamelist.ntord - 3][1] = 1.2163e-03;
		rbc[VmecIndataNamelist.ntord - 2][1] = 3.0809e-03;
		rbc[VmecIndataNamelist.ntord - 1][1] = 3.4528e-02;
		rbc[VmecIndataNamelist.ntord + 0][1] = 4.4872e-01;
		rbc[VmecIndataNamelist.ntord + 1][1] = -2.7085e-01;
		rbc[VmecIndataNamelist.ntord + 2][1] = -4.9284e-03;
		rbc[VmecIndataNamelist.ntord + 3][1] = 4.9911e-04;
		rbc[VmecIndataNamelist.ntord + 4][1] = 6.5184e-05;
		//
		rbc[VmecIndataNamelist.ntord - 4][2] = 4.8734e-04;
		rbc[VmecIndataNamelist.ntord - 3][2] = 1.4054e-03;
		rbc[VmecIndataNamelist.ntord - 2][2] = 4.8047e-03;
		rbc[VmecIndataNamelist.ntord - 1][2] = 1.5134e-02;
		rbc[VmecIndataNamelist.ntord + 0][2] = 3.5687e-02;
		rbc[VmecIndataNamelist.ntord + 1][2] = 4.9779e-02;
		rbc[VmecIndataNamelist.ntord + 2][2] = 6.5200e-02;
		rbc[VmecIndataNamelist.ntord + 3][2] = -1.1350e-02;
		rbc[VmecIndataNamelist.ntord + 4][2] = -1.6119e-03;
		//
		rbc[VmecIndataNamelist.ntord - 4][3] = 1.0339e-04;
		rbc[VmecIndataNamelist.ntord - 3][3] = -3.2332e-04;
		rbc[VmecIndataNamelist.ntord - 2][3] = -3.4468e-04;
		rbc[VmecIndataNamelist.ntord - 1][3] = -1.5729e-03;
		rbc[VmecIndataNamelist.ntord + 0][3] = -2.0611e-03;
		rbc[VmecIndataNamelist.ntord + 1][3] = -1.4756e-02;
		rbc[VmecIndataNamelist.ntord + 2][3] = -1.9949e-02;
		rbc[VmecIndataNamelist.ntord + 3][3] = -8.5802e-03;
		rbc[VmecIndataNamelist.ntord + 4][3] = 3.8516e-03;
		//
		rbc[VmecIndataNamelist.ntord - 4][4] = 1.1352e-04;
		rbc[VmecIndataNamelist.ntord - 3][4] = 3.6285e-04;
		rbc[VmecIndataNamelist.ntord - 2][4] = 2.4647e-04;
		rbc[VmecIndataNamelist.ntord - 1][4] = 6.2828e-04;
		rbc[VmecIndataNamelist.ntord + 0][4] = 2.7421e-03;
		rbc[VmecIndataNamelist.ntord + 1][4] = 4.9943e-03;
		rbc[VmecIndataNamelist.ntord + 2][4] = 7.4223e-03;
		rbc[VmecIndataNamelist.ntord + 3][4] = -5.0041e-04;
		rbc[VmecIndataNamelist.ntord + 4][4] = -5.9196e-04;
		//
		Assertions.assertArrayEquals(rbc, vmecInput.rbc);

		double[][] zbs = new double[2 * VmecIndataNamelist.ntord + 1][VmecIndataNamelist.mpold + 1];
		zbs[VmecIndataNamelist.ntord + 0][0] = -0.0000e+00;
		zbs[VmecIndataNamelist.ntord + 1][0] = -2.0666e-01;
		zbs[VmecIndataNamelist.ntord + 2][0] = -4.1458e-03;
		zbs[VmecIndataNamelist.ntord + 3][0] = -1.9493e-03;
		zbs[VmecIndataNamelist.ntord + 4][0] = 2.0185e-03;
		//
		zbs[VmecIndataNamelist.ntord - 4][1] = 3.4400e-03;
		zbs[VmecIndataNamelist.ntord - 3][1] = 7.7506e-03;
		zbs[VmecIndataNamelist.ntord - 2][1] = 1.4961e-02;
		zbs[VmecIndataNamelist.ntord - 1][1] = 4.0227e-02;
		zbs[VmecIndataNamelist.ntord + 0][1] = 5.6892e-01;
		zbs[VmecIndataNamelist.ntord + 1][1] = 2.0596e-01;
		zbs[VmecIndataNamelist.ntord + 2][1] = -5.7604e-03;
		zbs[VmecIndataNamelist.ntord + 3][1] = -5.6140e-03;
		zbs[VmecIndataNamelist.ntord + 4][1] = -4.2485e-03;
		//
		zbs[VmecIndataNamelist.ntord - 4][2] = 4.5363e-04;
		zbs[VmecIndataNamelist.ntord - 3][2] = 3.1625e-04;
		zbs[VmecIndataNamelist.ntord - 2][2] = 3.5963e-04;
		zbs[VmecIndataNamelist.ntord - 1][2] = 8.3725e-03;
		zbs[VmecIndataNamelist.ntord + 0][2] = 1.1405e-03;
		zbs[VmecIndataNamelist.ntord + 1][2] = 2.3889e-02;
		zbs[VmecIndataNamelist.ntord + 2][2] = -6.0502e-02;
		zbs[VmecIndataNamelist.ntord + 3][2] = 8.9796e-03;
		zbs[VmecIndataNamelist.ntord + 4][2] = 9.1004e-04;
		//
		zbs[VmecIndataNamelist.ntord - 4][3] = -3.7464e-04;
		zbs[VmecIndataNamelist.ntord - 3][3] = -4.2385e-05;
		zbs[VmecIndataNamelist.ntord - 2][3] = 5.3668e-04;
		zbs[VmecIndataNamelist.ntord - 1][3] = -1.7563e-03;
		zbs[VmecIndataNamelist.ntord + 0][3] = -4.2733e-03;
		zbs[VmecIndataNamelist.ntord + 1][3] = -4.4707e-03;
		zbs[VmecIndataNamelist.ntord + 2][3] = 9.5155e-03;
		zbs[VmecIndataNamelist.ntord + 3][3] = 1.0233e-02;
		zbs[VmecIndataNamelist.ntord + 4][3] = -2.8137e-03;
		//
		zbs[VmecIndataNamelist.ntord - 4][4] = 1.2480e-04;
		zbs[VmecIndataNamelist.ntord - 3][4] = -8.7567e-05;
		zbs[VmecIndataNamelist.ntord - 2][4] = 7.6525e-05;
		zbs[VmecIndataNamelist.ntord - 1][4] = 6.1672e-04;
		zbs[VmecIndataNamelist.ntord + 0][4] = 3.6261e-03;
		zbs[VmecIndataNamelist.ntord + 1][4] = -2.8280e-03;
		zbs[VmecIndataNamelist.ntord + 2][4] = 7.3549e-03;
		zbs[VmecIndataNamelist.ntord + 3][4] = -5.6303e-03;
		zbs[VmecIndataNamelist.ntord + 4][4] = -2.8346e-04;
		//
		Assertions.assertArrayEquals(zbs, vmecInput.zbs);

		double[][] rbs = new double[2 * VmecIndataNamelist.ntord + 1][VmecIndataNamelist.mpold + 1];
		Assertions.assertArrayEquals(rbs, vmecInput.rbs);

		double[][] zbc = new double[2 * VmecIndataNamelist.ntord + 1][VmecIndataNamelist.mpold + 1];
		Assertions.assertArrayEquals(zbc, vmecInput.zbc);
	}

	@Test
	void testParsingW7xReferenceRuns() {

		String refFolder = "/data/jonathan/vmecWebservice/reference_runs_vmec_w7x/";
		if (!new File(refFolder).exists()) {
			// only test if folder is accessible
			return;
		}

		String refIdListFile = "ref_id_list_from_webservice.txt";
		try {
			List<String> refIdLines = Files.readAllLines(Paths.get(refFolder, refIdListFile));
			for (String refId: refIdLines) {
				refId = refId.strip();
				String[] shortId_subFolder = refId.split("\\s+");

//				String shortId = shortId_subFolder[0];
				String subFolder = shortId_subFolder[1];

				if (subFolder.endsWith("/")) {
					subFolder = subFolder.substring(0, subFolder.length() - 1);
				}

				String longId = subFolder.replace("/", ".");

				String inputFilename = "input." + longId;

				String inputFile = Files.readString(Paths.get(refFolder, "reference_runs", subFolder, inputFilename));

				// --------------
				// now try parsing the INDATA namelist

				VmecIndataNamelist vmecInput = new VmecIndataNamelist();

				// parse the namelist into a Java object
				FortranNamelist parser = new FortranNamelist(inputFile, "indata", vmecInput);
				vmecInput = (VmecIndataNamelist) parser.getParsed();
				assertNotNull(vmecInput);
				vmecInput.sanitize();

//				System.out.println(shortId);
			}
		} catch (Exception e) {
			throw new RuntimeException(e);
		}
	}

}
