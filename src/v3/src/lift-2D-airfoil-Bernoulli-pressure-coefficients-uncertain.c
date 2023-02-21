/*v3 - assume that temperature, elevation, humidity and, therefore, fluid density are uncertain
* V = 30 m/s
* A = 0.23 m^2
* empirical coefficient for pressure distributions at 10° angle of attack
*/
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
 *	-	`A`:		0.1 to 1 m^2 - area of the 2D airfoil
 *	-	`h`:	    0.0 to 11019 m - elevation (troposphere)
 *  -   `T`:	   -50 °C to 50 °C - ambient temperature uncertain 
 *	-	`V`:		10 to 343 m/s - velocity below supersonic speed
 *	-	`Сp1`:		~-2.8 to 1.0 - coefficient for pressurre distribution over an airfoil
 *	-	`Сp2`:		~-0.54 to 1.14 - coefficient for pressurre distribution under an airfoil
 *
 *
 *	Bernouilli law 
 *  P1+1/2*𝜌*(𝑣1)^2+𝜌*𝑔*𝑦1=P2+1/2*𝜌*(𝑣2)^2+𝜌*𝑔*𝑦2.
 *	Assumimg the difference in height between the top and bottom of the airfoil is neglidgible, there hydrostatic term is going to be canceled out and leave us just static and dynamic pressure terms
 *	P1+1/2*𝜌*(𝑣1)^2=P2+1/2*𝜌*(𝑣2)^2.
 *  After moving static pressure terms to the left side and dynamic pressure terms to the right side, we get difference in pressure between the bottom and top of the airfoil
 *  P1-P2 = 1/2 * 𝜌 * ((𝑣2)^2- (𝑣1)^2)
 *  Velocities are being calculated based on pressure coefficient distributions 
 *  V = |Vstream| - sqrt(Cp - 1)
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
 *	Fl = 1/2 * 𝜌 * a  * ((𝑣2)^2- (𝑣1)^2)
 *
 *
 */
 void read_csv(int row, int col, char *filename, double **data){
	FILE *file;
	file = fopen(filename, "r");

	int i = 0;
    char line[2024];
	while (fgets(line, 4098, file) && (i < row))
    {
    	// double row[ssParams->nreal + 1];
        char* tmp = strdup(line);

	    int j = 0;
	    char* tok;
        //printf("%s\t", line);
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
	        //printf("%f\t", data[i][j]);
	    }
	    //printf("\n");

        free(tmp);
        i++;
    }
}

static void
loadInputs(double *  A, double *  v1, double * v2, double * r, double *  Cp1, double *  Cp2)
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
    printf("T=%f\n", T);
    printf("h=%f\n", h);
    printf("Rh=%f\n", Rh);
    printf("Pa=%f\n", Pair);
    printf("P1=%f\n", Psat);  
    printf("Pv=%f\n", Pv);
    printf("Pd=%f\n", Pd);

    int const sampleCount     = 3;
    int const row     = 140;
	int const col     = 7;
    int const totalLength     = (row-1)*2;
	char fname[256] = "all_angles.csv";

	double **data;
	data = (double **)malloc(row * sizeof(double *));
	for (int i = 0; i < row; ++i){
		data[i] = (double *)malloc(col * sizeof(double));
	}

	read_csv(row, col, fname, data);

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

    for (int i = 0; i < totalLength/2; i++)
    {
        empricalPressure10AOA[i]  =empricalPressureOver10AOA[i];
        empricalPressure10AOA[i+row-1]=empricalPressureUnder10AOA[i];

        empricalPressure5AOA[i]  =empricalPressureOver5AOA[i];
        empricalPressure5AOA[i+row-1]=empricalPressureUnder5AOA[i];

        empricalPressure0AOA[i]  =empricalPressureOver0AOA[i];
        empricalPressure0AOA[i+row-1]=empricalPressureUnder0AOA[i];
        //printf("empricalPressure10AOA[i]=%f\n", empricalPressure10AOA[i]);   
        //printf("empricalPressure10AOA[i+row-1]=%f\n", empricalPressure10AOA[i+row-1]);   
    }

    double empiricalPressureCoefficientsUncertain[sampleCount][2*(row-1)]  = { 
        *empricalPressure10AOA ,
        *empricalPressure5AOA ,
        *empricalPressure0AOA
    };

    double empiricalPressureCoefficients[2*(row-1)];

	libUncertainDoubleDistFromMultidimensionalSamples(
			empiricalPressureCoefficients,
			(void *) empiricalPressureCoefficientsUncertain,
			sampleCount,
			2*(row-1));


    /*Vx = Vstream * sqrt(|1-Cpx|)*/ 
	for (int i = 0; i < sizeof(empiricalPressureCoefficients)/sizeof(double)/2; i++) //(row-1) * sizeof(double)
	{  //printf("empricalPressure10AOA[row-1+i]=%f\n", empricalPressure10AOA[row-1+i]);   
		*v1 += V * sqrt(fabs(1-empiricalPressureCoefficients[i]));
        //printf("v1=%f\n", *v1); 
        *v2 += V * sqrt(fabs(1-empiricalPressureCoefficients[row-1+i]));
        //printf("v2=%f\n", *v2); 
	}
    
	*v1 /= sizeof(empiricalPressureCoefficients)/sizeof(double);
    *v2 /= sizeof(empiricalPressureCoefficients)/sizeof(double);
    printf("v1=%f\n", *v1);
    printf("v2=%f\n", *v2);
	*A		= 2.3E-1;
    printf("area=%f\n", *A);

    /*  r = (Pd/(Rd*T))+(Pv/(Rv*T)). */
    *r = (Pd/(287.058*(T+273.15)))+(Pv/(461.495*(T+273.15)));
    printf("density=%f\n", *r);

}

int main(int argc, char *	argv[])
{
	double	A, v1, v2, r, Cp1, Cp2, liftForce;

	loadInputs(&A, &v1, &v2, &r, &Cp1, &Cp2);

    /*	Fl = 1/2 * 𝜌 * a  * ((𝑣1)^2- (𝑣2)^2)*/
	liftForce = r*A*(pow(v1, 2)-pow(v2, 2)) / 2.0;

	printf("Lift force = %f\n", liftForce);

	return 0;
}



