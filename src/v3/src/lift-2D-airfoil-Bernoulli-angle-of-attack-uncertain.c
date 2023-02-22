/*v3 - assume that angle of attack is uncertain (0°,5° or 10°) (pressure coefficient distributions will be selected accordingly) */
#include <math.h>
#include <stdint.h>
#include <stdio.h>

#include <string.h>
#include <stdlib.h>
#include <uncertain.h>

/*  Overview: 
 *	Computation of generated lift force for a 2D NACA 2412 airfoil based on Bernoulli s equation (applicable only for inviscid and incompressible dry air flow)
 *  Reference levels h0 (reference height is set at sea level) = 0, P0 (reference pressure) = 1 atm
 *
 *	Inputs v3 [NACA 2412, uncertain angle of attack (0°,5° or 10°)]:
 *	The inputs and their ranges are:
 *	-	`A`:		0.23 m^2 - area of the 2D airfoil
 *	-	`h`:	    0.0 m - elevation (troposphere)
 *  -   `T`:	    15 °C- ambient temperature uncertain 
 *  -   `Rh`:	    0.0 - humidity level (dry air)
 *	-	`V`:		30 m/s - free stream velocity below supersonic speed
 *	-	`Сp1`:		~-2.8 to 1.0 - coefficient for pressurre distribution over an airfoil
 *	-	`Сp2`:		~-0.54 to 1.14 - coefficient for pressurre distribution under an airfoil
 *
 *
 *  Velocities are being calculated based on pressure coefficient distributions 
 *  𝑣x = V * sqrt(|1-Cpx|)
 *   
 *  Air density r(kg/m^3) calculation process:
 *
 *  air pressure in Pa [Pair = P0 × exp(-g × M × (h - h0)/(R × T)] = result in atm (1 atm = 101,325 Pa) 
 *  where M (molar mass of air) = 0.0289644 kg/mol, R (universal gas constant) = 8.31432 N · m/(mol · K), g (acceleration due to gravity) = 9.80665 m/(s^2), T - temperature in Kelvins, reference height h0 = 0 m and reference pressure at sea level P0 = 1 atm
 *
 *  r(humid air) = (Pd/(Rd*T))+(Pv/(Rv*T)) 
 *  where Rd = 287.058 J/(kg·K), Rv = 461.495 J/(kg·K),  T - temperature in Kelvins, Pv (water vapor pressure in P) = p1*Rh,  Pd (pressure of dry air in Pa) = Pair - Pv, Psat (saturation vapor pressure in Pa)= 6.1078*10^(7,5*T/(T+237,3))
 *
 *
 *  Lift force is being calculated by multiplying the net pressure by the wing area.
 *  
 *	Outputs:
 *  - 'Fl' -  Lift force (N)
 *	Fl = 1/2 * 𝜌 * a  * ((𝑣2)^2- (𝑣1)^2)  
 *
 *
 */

enum {
     sampleCount     = 3,
     row     = 140,
	 col     = 7,
     totalLength     = (row-1)*2,
};

void read_csv(int row, int col, char *filename, double **data){
	FILE *file;
	file = fopen(filename, "r");

	int i = 0;
    char line[4098];
	while (fgets(line, 4098, file) && (i < row))
    {
        char* tmp = strdup(line);

	    int j = 0;
	    char* tok;
	    for (tok = strtok(line, ";"); tok && *tok; j++, tok = strtok(NULL, ";\n"))
	    {
            int index=0;
            while(tok[index]!='\0')
            {
                if(tok[index]==',')
                {
                    tok[index]='.';
                }
                index++;
            }
	        data[i][j] = atof(tok);
	    }

        free(tmp);
        i++;
    }
}

static void
loadInputs(double *  A, double *  v1, double * v2, double * r, double *  Cp1, double *  Cp2, double **data)
{
    double V  = 30.0;
    double Rh = 0.0; 
    double h  = 0.0; 
    double T  = 15.0;

    /* air pressure in Pascals Pair = P0 × exp(-g × M × (h - h0)/(R × T) * 101,325 (1 atm );
    *
    * Rd = 287.058 J/(kg·K), Rv = 461.495 J/(kg·K),  
    * Pv = Psat*Rh,  
    * Pd = Pair - Pv, 
    * Psat = 6.1078*10^(7,5*T/(T+237,3))
    *
    */
    double Pair = exp((-9.81 * 0.0289644 * h)/(8.31432 * (T+273.15))) * 101325.0;// atm * sea level pressure 101325 hPa
    double Psat = 6.1078*pow(10.0,7.5*T/(T+237.3));
    double Pv = Psat*Rh;
    double Pd = Pair - Pv;
    // printf("T=%f\n", T);
    // printf("h=%f\n", h);
    // printf("Rh=%f\n", Rh);
    // printf("Pa=%f\n", Pair);
    // printf("P1=%f\n", Psat);  
    // printf("Pv=%f\n", Pv);
    // printf("Pd=%f\n", Pd);

    double empricalPressureOver10AOA[row-1];
    double empricalPressureOver5AOA[row-1];
    double empricalPressureOver0AOA[row-1];
    double empricalPressureUnder0AOA[row-1];
    double empricalPressureUnder5AOA[row-1];
    double empricalPressureUnder10AOA[row-1];
    for(int i = 1; i < 140; i++)
    {
        empricalPressureOver10AOA[i-1] = data[i][1];
        empricalPressureOver5AOA[i-1] = data[i][2];
        empricalPressureOver0AOA[i-1] = data[i][3];
        empricalPressureUnder0AOA[i-1] = data[i][4];
        empricalPressureUnder5AOA[i-1] = data[i][5];
        empricalPressureUnder10AOA[i-1] = data[i][6];
    }
    double empricalPressure10AOA[totalLength];
    double empricalPressure5AOA[totalLength];
    double empricalPressure0AOA[totalLength];

    //create arrays that have over (from index 0 to row-1) and lower (from index row -1 to end of the array ) pressure coefficients for certain angles of attack
    for (int i = 0; i < totalLength/2; i++)
    {
        empricalPressure10AOA[i]  =empricalPressureOver10AOA[i];
        empricalPressure10AOA[i+row-1]=empricalPressureUnder10AOA[i];

        empricalPressure5AOA[i]  =empricalPressureOver5AOA[i];
        empricalPressure5AOA[i+row-1]=empricalPressureUnder5AOA[i];

        empricalPressure0AOA[i]  =empricalPressureOver0AOA[i];
        empricalPressure0AOA[i+row-1]=empricalPressureUnder0AOA[i];
        // printf("empricalPressure10AOA[i]=%f\n", empricalPressure10AOA[i]);   
        // printf("empricalPressure10AOA[i+row-1]=%f\n", empricalPressure10AOA[i+row-1]);   
    }

    double empiricalPressureCoefficientsUncertain[sampleCount][totalLength];
  
    for(int i = 0; i < totalLength; i++){
            empiricalPressureCoefficientsUncertain[0][i] = empricalPressure10AOA[i];
            empiricalPressureCoefficientsUncertain[1][i] = empricalPressure5AOA[i];
            empiricalPressureCoefficientsUncertain[2][i] = empricalPressure0AOA[i];
    }

    double empiricalPressureCoefficients[totalLength];

	libUncertainDoubleDistFromMultidimensionalSamples(
			empiricalPressureCoefficients,
			(void *) empiricalPressureCoefficientsUncertain,
			sampleCount,
			totalLength);

    /*Vx = Vstream * sqrt(|1-Cpx|)*/ 
	for (int i = 0; i < sizeof(empiricalPressureCoefficients)/sizeof(double)/2; i++) //(row-1) * sizeof(double)
	{ 
		*v1 += V * sqrt(fabs(1-empiricalPressureCoefficients[i]));
        *v2 += V * sqrt(fabs(1-empiricalPressureCoefficients[row-1+i]));
	}
    
	*v1 /= sizeof(empiricalPressureCoefficients)/sizeof(double);
    *v2 /= sizeof(empiricalPressureCoefficients)/sizeof(double);
    // printf("v1=%f\n", *v1);
    // printf("v2=%f\n", *v2);
	*A		= 2.3E-1;
    // printf("area=%f\n", *A);

    /*  r = (Pd/(Rd*T))+(Pv/(Rv*T)). */
    *r = (Pd/(287.058*(T+273.15)))+(Pv/(461.495*(T+273.15)));
    // printf("density=%f\n", *r);

}

int main(int argc, char * argv[])
{
	if (argc < 1){
		printf("Please specify the CSV file as an input.\n");
		exit(0);
	}

    char fname[256]; strcpy(fname, argv[1]);

    double ** data;
	data = (double **)malloc(row * sizeof(double *));
	for (int i = 0; i < row; ++i){
		data[i] = (double *)malloc(col * sizeof(double));
	}

	double	A, v1, v2, r, Cp1, Cp2, liftForce;


	read_csv(row, col, fname, data);
	loadInputs(&A, &v1, &v2, &r, &Cp1, &Cp2, data);

    /*	Fl = 1/2 * 𝜌 * a  * ((𝑣1)^2- (𝑣2)^2) */
	liftForce = r*A*(pow(v1, 2)-pow(v2, 2)) / 2.0;

	printf("Lift force = %f\n", liftForce);

	return 0;
}



