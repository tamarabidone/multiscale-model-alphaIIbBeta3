package nov15;

import java.util.ArrayList;

public class Force_Interaction{
	
	ArrayList<Integrin_3D> sIntegrin = new ArrayList<Integrin_3D>();
	ArrayList<Ligand_3D> sLigand = new ArrayList<Ligand_3D>();
	int nOfIntegrins,nOfLigands;

	public Force_Interaction(ArrayList<Integrin_3D> sIntegrin,  ArrayList<Ligand_3D> sLigand, int nOfIntegrins, int nOfLigands){
		this.sIntegrin= sIntegrin;
		this.sLigand= sLigand;
		this.nOfIntegrins = nOfIntegrins;	
		this.nOfLigands = nOfLigands;	
	}

	public void IL_interact(){
		for(int i=0;i<nOfIntegrins;i++){		
			double r_min=10000;
			double A=	sLigand.get(0).radius*sLigand.get(0).radius*Math.PI;
			double k_it=sLigand.get(0).Y*(A)/(sIntegrin.get(0).initialp[2]-sLigand.get(0).initialp[2]);//Y*A/L	
//System.out.println("   is it rigth?"+k_it);
			if (sIntegrin.get(i).BoundLigandYN==1){
				int lig1=sIntegrin.get(i).LigandNum;
				double rr=distance(sIntegrin.get(i).Bead_Position,sLigand.get(lig1).Bead_PositionL);
				double force=k_it*( rr-(0.025));
				sLigand.get(lig1).FF_interact=force;
				sIntegrin.get(i).Bead_Position=move2(sIntegrin.get(i).Bead_Position, sLigand.get(lig1).Bead_PositionL, force, rr, sIntegrin.get(i).dt, sIntegrin.get(i).zetaI);
				sLigand.get(lig1).Bead_PositionL=move2(sLigand.get(lig1).Bead_PositionL, sIntegrin.get(i).Bead_Position,  force, rr, sLigand.get(lig1).dt, sLigand.get(lig1).zetaL);	
				sLigand.get(lig1).BoundYN=1;
			}
			
			int lig=-1;
			r_min=1000;
			if (sIntegrin.get(i).BoundLigandYN==0) {
				for(int j=0;j<nOfLigands;j++){	
				double dist=distance(sIntegrin.get(i).Bead_Position,sLigand.get(j).Bead_PositionL);
				if (sIntegrin.get(i).BoundLigandYN==0  && sLigand.get(j).BoundYN==0 && dist<r_min){						
					r_min=dist;
					lig=j;
				}
			}
			if (r_min<=0.027){
				sLigand.get(lig).integrinBoundNum=i;
				sLigand.get(lig).BoundYN=1;
				sIntegrin.get(i).BoundLigandYN=1;
				sIntegrin.get(i).LigandNum=lig;
				double force=k_it*( r_min-(0.025));
				//System.out.println(k_it);
				sLigand.get(lig).FF_interact=force;
				sLigand.get(lig).Bead_PositionL=move2(sLigand.get(lig).Bead_PositionL, sIntegrin.get(i).Bead_Position,  force, r_min, sLigand.get(lig).dt, sLigand.get(lig).zetaL);	
				sIntegrin.get(i).Bead_Position=move2(sIntegrin.get(i).Bead_Position, sLigand.get(lig).Bead_PositionL, force, r_min, sIntegrin.get(i).dt, sIntegrin.get(i).zetaI);
			}
			}
		}
	}
	
	
protected void CatchBond(int typeCB){
	double slope=0, base=0, A=0, slope2=0, B=0, C=0;


		slope =8.2;
	base=2;
	slope2=1;
	A=-0.012;
	B=0.3;
	C=0.5E-5;
	
	for(int i=0;i<nOfLigands;i++){
		if (sLigand.get(i).BoundYN==1){
		 double P;
		 double force=sLigand.get(i).FF_interact;
		// System.out.println(force);
	  sIntegrin.get(sLigand.get(i).integrinBoundNum).kub=base*Math.exp(A*slope*force)+C*Math.exp(B*slope2*force);	
	  sIntegrin.get(sLigand.get(i).integrinBoundNum).FinalLife=0;
		 	 P=(sIntegrin.get(sLigand.get(i).integrinBoundNum).kub*sLigand.get(i).dt);
		 	sIntegrin.get(sLigand.get(i).integrinBoundNum).contaAT=sIntegrin.get(sLigand.get(i).integrinBoundNum).contaAT+1;
		 	 if (Math.random()<P){	
		 		sLigand.get(i).BoundYN=0;
		 		sIntegrin.get(sLigand.get(i).integrinBoundNum).finalTension=sLigand.get(i).FF_interact;
		 		sLigand.get(i).FF_interact=0;
		 		sIntegrin.get(sLigand.get(i).integrinBoundNum).BoundLigandYN=0;
		 		sIntegrin.get(sLigand.get(i).integrinBoundNum).FinalLife=sIntegrin.get(sLigand.get(i).integrinBoundNum).contaAT;
		 		sIntegrin.get(sLigand.get(i).integrinBoundNum).contaBoundLifetimeTot=(sIntegrin.get(sLigand.get(i).integrinBoundNum).contaBoundLifetimeTot+sIntegrin.get(sLigand.get(i).integrinBoundNum).contaAT*sIntegrin.get(sLigand.get(i).integrinBoundNum).dt);
		 		sIntegrin.get(sLigand.get(i).integrinBoundNum).contaAT=0;
		 		sIntegrin.get(sLigand.get(i).integrinBoundNum).kub=0;
		 		sLigand.get(i).integrinBoundNum=-1;
		 		}		
		 	 }
	}
}

	
private double[] move2(double[] beadA, double[] beadB, double force, double r, double dt, double zeta) {
	double temp[] = new double[3];
	temp[0] = beadA[0]-(force)*(beadA[0]-beadB[0])/r*dt/zeta;
	temp[1] = beadA[1]-force*(beadA[1]-beadB[1])/r*dt/zeta;
    temp[2] = beadA[2];//force*(beadA[2]-beadB[2])/r*dt/zeta;
	return temp;
}
	private double distance(double[] vecA, double[] vecB) {
		return Math.sqrt((vecA[0]-vecB[0])*(vecA[0]-vecB[0])+(vecA[1]-vecB[1])*(vecA[1]-vecB[1])+(vecA[2]-vecB[2])*(vecA[2]-vecB[2]));
	}
}