package nov15;




import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;

import java.util.Random;

/**
 * This program simulated brownian dynamics of integrin particles interactin with ligands on a substrate
 * 
 * @author tamarabidone
 *
 * 
 */

	/* An ARRAY of integrin and ligands instances */
	/*
	 * in this simulation the arrays are not varying lenght, but I would like to
	 * keep the arraylist in C++ cause I'll have eventually addition/removal of
	 * integrins and ligands over the course of the simualtions
	 */
public class Main_App {	
	ArrayList<Integrin_3D> sIntegrin = new ArrayList<Integrin_3D>();
	ArrayList<Ligand_3D> sLigand = new ArrayList<Ligand_3D>();

	double eta; // Viscosity
	double timestep;
	int nOfIntegrins, nOfLigands;
	Random rand = new Random();
	double initialtime = 0;
	int norm=1;
	
	
	public static void main(String[] args) {
		Main_App simulation = new Main_App();
		simulation.initialize();
		for (int i = 0; i < (int) (100 / simulation.timestep); i++) {
			simulation.doStep();
		}

		System.exit(0);

	}
	
	
	public void initialize() {

		timestep = 0.0001;
		nOfIntegrins = 100;
		nOfLigands = 300;

		for (int i = 0; i < nOfIntegrins; i++) {
			sIntegrin.add(new Integrin_3D());

			sIntegrin.get(i).DomainSize = 1;
			sIntegrin.get(i).dt = timestep;
			sIntegrin.get(i).temperature = 300;
			sIntegrin.get(i).zIntegrin = 0;
			sIntegrin.get(i).zetaI = 0.0142;// pN*s/um //friction corresponding to free diffusion of beta-3 integrin
											// 0.29 um^2/s sIntegrin.get(i).LigandNum=-1;
			sIntegrin.get(i).BoundLigandYN = 0;
			sIntegrin.get(i).TimesUnbound = 0;
			sIntegrin.get(i).contaAT = 0;
			sIntegrin.get(i).kub = 0;
			sIntegrin.get(i).FinalLife = 0;
			sIntegrin.get(i).V = 0.0;
			double[] initialpos = new double[3];
			initialpos[0] = -sIntegrin.get(i).DomainSize / 2 + Math.random() * (sIntegrin.get(i).DomainSize);
			initialpos[1] = -sIntegrin.get(i).DomainSize / 2 + Math.random() * (sIntegrin.get(i).DomainSize);
			initialpos[2] = sIntegrin.get(i).zIntegrin;
			sIntegrin.get(i).start(initialpos);
		}

		double posx=-0.5;
		double posy=-0.5;
		int totalx = (int) (sIntegrin.get(0).DomainSize/0.05);
		System.out.println(totalx);
		for (int i = 0; i < nOfLigands; i++) {
			sLigand.add(new Ligand_3D());
			sLigand.get(i).DomainSize = 1.0;
			sLigand.get(i).dt = 1E-4;// time step
			sLigand.get(i).zLigand = -0.025;
			sLigand.get(i).integrinBoundNum = -1;
			sLigand.get(i).radius = 0.018;
			sLigand.get(i).zetaL = 500;
			double[] initialpos = new double[3];
			initialpos[0] = -(sLigand.get(i).DomainSize / 2) + Math.random() * sLigand.get(i).DomainSize;
			initialpos[1] = -(sLigand.get(i).DomainSize / 2) + Math.random() * sLigand.get(i).DomainSize;
			initialpos[2] = sLigand.get(i).zLigand;
			sLigand.get(i).start(initialpos);
			sLigand.get(i).Y = 2100*norm;
		}
	}

	int typeCB = 1;
	int rep = 1;




	
	
	
	
	/* The step is done at each time step */
	public void doStep() {

		/* Initiate the class with in interactions between integrins and ligands */
		Force_Interaction forceinter = new Force_Interaction(sIntegrin, sLigand, nOfIntegrins, nOfLigands);

		for (int j = 0; j < nOfLigands; j++)
			sLigand.get(j).step();

		for (int j = 0; j < nOfIntegrins; j++) 
			sIntegrin.get(j).step();

		forceinter.IL_interact();
		forceinter.CatchBond(typeCB);

		/////This is for output files////
		//I would like to have this in another file for readability
		int numLigandBoundIntegrin = 0;
		double t = 0;
		int unboun = 0;

		int numBound = 0;

		for (int i = 0; i < nOfLigands; i++) {
			if (sLigand.get(i).BoundYN == 1) {
				numBound = numBound + 1;
				t = t + sLigand.get(i).FF_interact;

			} else {
				unboun = unboun + 1;

			}

		}
		t = t / numBound;

		if (numBound > nOfIntegrins) 
			System.out.println("error");
	

		int ccc = 0;
		double[] Thresh = { 0.1, 0.2, 0.3, 0.4, 0.5, 1.0, 1.5, 2.0, 2.5, 3 };
		int uno = 0;
		int unocinque = 0;
		int due = 0;
		int duecinque = 0;
		int tre = 0;

		int zerocinque = 0;
		int zerouno = 0, zerodue = 0, zerotre = 0, zeroquattro = 0;
		for (int i = 0; i < nOfIntegrins; i++) {
			if (sIntegrin.get(i).FinalLife * sIntegrin.get(i).dt <= Thresh[0]) {
				if (sIntegrin.get(i).FinalLife * sIntegrin.get(i).dt > 0)
					zerouno = zerouno + 1;
			}
			if (sIntegrin.get(i).FinalLife * sIntegrin.get(i).dt > Thresh[0]
					&& sIntegrin.get(i).FinalLife * sIntegrin.get(i).dt <= Thresh[1]) {
				zerodue = zerodue + 1;
			}
			if (sIntegrin.get(i).FinalLife * sIntegrin.get(i).dt > Thresh[1]
					&& sIntegrin.get(i).FinalLife * sIntegrin.get(i).dt <= Thresh[2]) {
				zerotre = zerotre + 1;
			}
			if (sIntegrin.get(i).FinalLife * sIntegrin.get(i).dt > Thresh[2]
					&& sIntegrin.get(i).FinalLife * sIntegrin.get(i).dt <= Thresh[3]) {
				zeroquattro = zeroquattro + 1;
			}
			if (sIntegrin.get(i).FinalLife * sIntegrin.get(i).dt > Thresh[3]
					&& sIntegrin.get(i).FinalLife * sIntegrin.get(i).dt <= Thresh[4]) {
				zerocinque = zerocinque + 1;
			}
			if (sIntegrin.get(i).FinalLife * sIntegrin.get(i).dt > Thresh[4]
					&& sIntegrin.get(i).FinalLife * sIntegrin.get(i).dt <= Thresh[5]) {
				uno = uno + 1;
			}
			if (sIntegrin.get(i).FinalLife * sIntegrin.get(i).dt > Thresh[5]
					&& sIntegrin.get(i).FinalLife * sIntegrin.get(i).dt <= Thresh[6]) {
				unocinque = unocinque + 1;
			}
			if (sIntegrin.get(i).FinalLife * sIntegrin.get(i).dt > Thresh[6]
					&& sIntegrin.get(i).FinalLife * sIntegrin.get(i).dt <= Thresh[7]) {
				due = due + 1;
			}
			if (sIntegrin.get(i).FinalLife * sIntegrin.get(i).dt > Thresh[7]
					&& sIntegrin.get(i).FinalLife * sIntegrin.get(i).dt <= Thresh[8]) {
				duecinque = duecinque + 1;
			}

			if (sIntegrin.get(i).FinalLife * sIntegrin.get(i).dt > Thresh[8]
					&& sIntegrin.get(i).FinalLife * sIntegrin.get(i).dt <= Thresh[9]) {
				tre = tre + 1;
			}
			if (sIntegrin.get(i).TimesUnbound > 0) {
				ccc = ccc + sIntegrin.get(i).TimesUnbound;

			}
			ccc = ccc / nOfIntegrins;

			if (sIntegrin.get(i).BoundLigandYN == 1) {
				numLigandBoundIntegrin = numLigandBoundIntegrin + 1;
			}
		}

		double[] dd_min = new double[nOfLigands];
		double[] dd_avg = new double[nOfLigands];
		for (int i = 0; i < nOfLigands - 1; i++) {
			if (sLigand.get(i).BoundYN == 1) {
				double minu = 10000;
				int count = 0;

				for (int j = i + 1; j < nOfLigands; j++) {
					if (sLigand.get(j).BoundYN == 1) {
						double dtt = Math.sqrt(
								Math.pow(sLigand.get(i).Bead_PositionL[0] - sLigand.get(j).Bead_PositionL[0], 2) + Math
										.pow(sLigand.get(i).Bead_PositionL[1] - sLigand.get(j).Bead_PositionL[1], 2));
						count = count + 1;
						dd_avg[i] = dd_avg[i] + dtt;
						if (dtt < minu) {
							dd_min[i] = dtt;
							minu = dtt;
						}
					}

				}
				dd_avg[i] = dd_avg[i] / count;
			} else {
				dd_avg[i] = 0;
				dd_min[i] = 0;
			}
		}


		String l2 = "Fig_0_" + Integer.toString(nOfLigands) + "_" + Integer.toString(nOfIntegrins) + "_" + Integer.toString(typeCB) + "_Integrin_"
				+ Double.toString(sIntegrin.get(0).V * 1E3) + "_v_" + Double.toString(sLigand.get(0).Y) + "_"
				+ Double.toString(rep);
		String l3 = "Fig_0_" + Integer.toString(nOfLigands) + "_" + Integer.toString(nOfIntegrins) + "_" + Integer.toString(typeCB) + "_IntegrinFinals_"
				+ Double.toString(sIntegrin.get(0).V * 1E3) + "_v_" + Double.toString(sLigand.get(0).Y) + "_"
				+ Double.toString(rep);
		String ld = "Forces_on_Ligands_" + Integer.toString(nOfLigands) + "_" + Integer.toString(nOfIntegrins) + "_" + Integer.toString(typeCB) + "_L_"
				+ Double.toString(sIntegrin.get(0).V * 1E3) + "_v_" + Integer.toString(sLigand.get(0).Y) + "_"
				+ Integer.toString(rep);
		String l4 = " Fig_0_" + Integer.toString(nOfLigands) + "_" + Integer.toString(nOfIntegrins) + "_General_" + Integer.toString(typeCB) + "_L_"
				+ Double.toString(sIntegrin.get(0).V * 1E3) + "_v_" + Integer.toString(sLigand.get(0).Y) + "_"
				+ Integer.toString(rep);

//// Saving output files every second
		if (sIntegrin.get(0).time % 1 < 0.0001) {
			// System.out.println("Working Directory = " + System.getProperty("user.dir"));
			{
				try {
					FileWriter fw2 = new FileWriter(l2, true);
					FileWriter fw3 = new FileWriter(l3, true);
					FileWriter fw4 = new FileWriter(l4, true);
					fw4.write(String.format("%.3f", sIntegrin.get(0).time) + "    " + "   "
							+ String.format("%d", numLigandBoundIntegrin) + "   "+ String.format("%.3f", t)+"  "+ String.format("%d", zerouno) + "   "
							+ String.format("%d", zerodue) + "   " + String.format("%d", zerotre) + "     "                                                                                                       
							+ String.format("%d", zeroquattro) + "     " + String.format("%d", zerocinque) + "  "
							+ String.format("%d", uno) + "     " + String.format("%d", unocinque) + "     "
							+ String.format("%d", due) + "     " + String.format("%d", duecinque) + "     "
							+ String.format("%d", tre) + " \n");// appends the string to the file
					fw4.close();
					for (int i = 0; i < nOfIntegrins; i++) {
						fw2.write(String.format("%.3f", sIntegrin.get(0).time) + "    " + i + "   "
								+ String.format("%d", sIntegrin.get(i).BoundLigandYN) + "   "
								+ String.format("%.3f", sIntegrin.get(i).Bead_Position[0]) + "   "
								+ String.format("%.3f", sIntegrin.get(i).Bead_Position[1]) + "   "
								+ String.format("%.3f", sIntegrin.get(i).Bead_Position[2]) + "     " + "    "
								+ String.format("%d", sIntegrin.get(i).contaAT) + "     "
								+ String.format("%.6f", sIntegrin.get(i).finalTension) + "  "
								+ String.format("%.6f", sIntegrin.get(i).FinalLife * sIntegrin.get(i).dt) + "     "
								+ String.format("%d", sIntegrin.get(i).TimesUnbound) + "     "
								+ String.format("%.4f", sIntegrin.get(i).contaBoundLifetimeTot) + " \n");// appends the
																											// string to
																											// the file

						if (sIntegrin.get(i).FinalLife > 0)
							fw3.write(String.format("%.3f", sIntegrin.get(0).time) + "    " + i + "   "
									+ String.format("%.6f", sIntegrin.get(i).finalTension) + "  "
									+ String.format("%.6f", sIntegrin.get(i).FinalLife * sIntegrin.get(i).dt) + "     "
									+ " \n");// appends the string to the file
					}
					fw2.close();
					fw3.close();
					FileWriter fw1 = new FileWriter(ld, true);

					for (int i = 0; i < nOfLigands; i++) {
						if (sLigand.get(i).BoundYN == 1)
							fw1.write(String.format("%.3f", sIntegrin.get(0).time) + "    " + String.format("%d", i)
									+ "    " + String.format("%d", sLigand.get(i).BoundYN) + "     "
									+ String.format("%.3f", sLigand.get(i).Bead_PositionL[0]) + "    "
									+ String.format("%.3f", sLigand.get(i).Bead_PositionL[1]) + "    "
									+ String.format("%.3f", sLigand.get(i).Bead_PositionL[2]) + "    "
									+ String.format("%d", sLigand.get(i).integrinBoundNum) + "    "
									+ String.format("%3.4f", sLigand.get(i).FF_interact) + "    "
									+ " \n");// appends the string to the file
					}
					fw1.close();
				} catch (IOException ioe) {
					System.err.println("IOException: " + ioe.getMessage());
				}

			}
		}
		
	}

}
