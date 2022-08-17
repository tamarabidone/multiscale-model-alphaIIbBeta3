package nov15;


import java.util.*;

/**
 * The class Integrin_3D includes the properties of a single
 * integrin. Thermal and velocity motions 
 * are included. Main class is in Mechanosensing_3DApp.
 * @author tamarabidone
 *
 */

public class Integrin_3D{
	double V;
	double time = 0, dt;
	double Bead_Position[];  
	double Bead_Force[];	
	double initialp[];	//Initial positions
	double kub=0;
	int LigandNum, TimesUnbound;	
	int contaAT, BoundLigandYN;
	double zetaI,  contaBoundLifetimeTot, FinalLife, finalTension=0;
	double DomainSize;
	double temperature;
    Random rand = new Random();
    double zIntegrin;
    public void start(double[] initialpos){
    	initialp = initialpos;
	    Bead_Position=new double[3];
	    Bead_Position[0]=initialpos[0];
	    Bead_Position[1]=initialpos[1];
	    Bead_Position[2]=initialpos[2];   
    }
    

    /*Update forces and move along the 2D surface*/
	public void step(){		
	    Bead_Force=new double[3];
		time+=dt;
	
		if (BoundLigandYN==0){
			Boundary_Force_Update();
			Thermal_Force_Update();
			kub=0;
			}
		else{
		Thermal_Force_Update();
		Flow_Force_Update();
		}

		move();
    }

	protected void Thermal_Force_Update(){
			double thermal_Force=Math.sqrt(2*0.0000138*temperature*zetaI/dt);
			Bead_Force[0]+=thermal_Force*rand.nextGaussian();
			Bead_Force[1]+=thermal_Force*rand.nextGaussian();
		//	Bead_Force[2]=thermal_Force*rand.nextGaussian();
		
	}	

	protected void Flow_Force_Update(){
		Bead_Force[0]+=-zetaI*V+0;
	}

	protected void Boundary_Force_Update(){
		if (Bead_Position[0]>=DomainSize/2 && BoundLigandYN==0){
			Bead_Position[0]=Bead_Position[0]-DomainSize;
		}
		if (Bead_Position[0]<-DomainSize/2 && BoundLigandYN==0){
			Bead_Position[0]=Bead_Position[0]+DomainSize;
		}
		
		if (Bead_Position[1]>=DomainSize/2  && BoundLigandYN==0){
			Bead_Position[1]=Bead_Position[1]-DomainSize;
		}
		if (Bead_Position[1]<-DomainSize/2  && BoundLigandYN==0){
			Bead_Position[1]=Bead_Position[1]+DomainSize;
		}
	}
	
	protected void move(){
			for(int j=0; j<2;j++){
				Bead_Position[j]=Bead_Position[j]-Bead_Force[j]*dt/zetaI;
		}
	}
	
	//FORCES********************************************************
	protected double Boundary_Force_z(double z){		
		
		double f0 = 100*(z-zIntegrin);	//harmonic restrain for the lipid bilayer. 100pN/um
		return f0;
	}	

}
