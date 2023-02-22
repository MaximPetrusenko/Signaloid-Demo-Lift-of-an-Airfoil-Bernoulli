/*v1 - no uncertainties */ 
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <uncertain.h>

/*  Overview: 
 *	Computation of generated lift force for a 2D NACA 2412 airfoil based on Bernoulli s equation (applicable only for inviscid and incompressible dry air flow)
 *  Reference levels h0 (reference height is set at sea level) = 0, P0 (reference pressure) = 1 atm
 *
 *	Inputs v1 [NACA 2412, 10° angle of attack]:
 *	The inputs and their ranges are:
 *	-	`A`:		0.23 m^2 - area of the 2D airfoil
 *	-	`h`:	   	0.0 m - elevation (troposphere)
 *  -   `T`:	    15 °C - ambient temperature
 *  -   `Rh`:	    0.0 - humidity level (dry air)
 *	-	`V`:		30 m/s - free stream velocity below supersonic speed
 *	-	`Сp1`:		~-0.54 to 1.14 - coefficient for pressurre distribution under an airfoil (digitized plot) at 10° angle of attack
 *	-	`Сp2`:		~-2.8 to 1.0 - coefficient for pressurre distribution over an airfoil (digitized plot) at 10° angle of attack
 *
 *  Velocities are being calculated based on pressure coefficient distributions 
 *  𝑣x = Vstream * sqrt(|1-Cpx|)
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


static void
loadInputs(double *  A, double *  v1, double * v2, double * r, double *  Cp1, double *  Cp2)
{
    double V  = 30.0;
    double Rh = 0.0;
    double h  = 0.0;
    double T  = 15.0;
    //air pressure
    double Pair = exp((-9.81 * 0.0289644 * h)/(8.31432 * (T+273.15))) * 101325.0; // atm * sea level pressure 101325 hPa
    //saturation vapor pressure
    double Psat = 6.1078*pow(10.0,7.5*T/(T+237.3));
    //water vapor pressure
    double Pv = Psat*Rh;
    //pressure of dry air
    double Pd = Pair - Pv;

	double Cp2[] = {
        -2.3444, -2.4402, -2.5411, -2.577, -2.7322, -2.7316, -2.5977,
        -2.575, -2.5415, -2.3405, -2.3121, -2.2061, -2.1597, -2.0826,
        -1.9988, -1.9037, -1.7997, -1.7692, -1.63, -1.6235, -1.4999
        -1.4769, -1.4098, -1.3809, -1.3528, -1.3367, -1.3181, -1.2695,
        -1.239, -1.1633, -1.1599, -1.0807, -1.0715, -1.0127, -0.9936,
        -0.9336, -0.8987, -0.8544, -0.8222, -0.7642, -0.7355, -0.6851,
        -0.645, -0.6061, -0.5636, -0.538, -0.4927, -0.4825, -0.4468,
        -0.4431, -0.4454, -0.444, -0.4329, -0.4205, -0.4094, -0.3889,
        -0.3636, -0.349, -0.3179, -0.2992, -0.2832, -0.2727, -0.2596,
        -0.2451, -0.2248, -0.2195, -0.2012, -0.1998, -0.1808, -0.1781,
        -0.1831, -0.1885, -0.1837, -0.1769, -0.1889, -0.1865, -0.1799,
        -0.1841, -0.1785, -0.1838, -0.1742, -0.1779, -0.1823, -0.1789
	};

        double Cp1[] = {
        0.8111, 0.9226, 1.0007, 0.9934, 0.8905, 0.8737, 0.7471,
        0.7336, 0.714, 0.6252, 0.6152, 0.5857, 0.5611, 0.4833,
        0.429, 0.403, 0.3861, 0.3781, 0.3431, 0.3423, 0.3439,
        0.3448, 0.3393, 0.3353, 0.3354, 0.3368, 0.3345, 0.3272,
        0.3228, 0.3067, 0.3057, 0.2782, 0.275, 0.2539, 0.2432,
        0.2017, 0.1893, 0.2187, 0.2461, 0.2578, 0.2585, 0.2675,
        0.2711, 0.2632, 0.2392, 0.2206, 0.1912, 0.1853, 0.1643,
        0.1539, 0.1427, 0.1439, 0.1585, 0.1679, 0.1675, 0.1579,
        0.1564, 0.159, 0.164, 0.1606, 0.1438, 0.1286, 0.1278,
        0.1299, 0.1214, 0.1199, 0.1316, 0.1322, 0.1217, 0.1134,
        0.1001, 0.102, 0.1118, 0.1173, 0.1216, 0.1122, 0.1017,
        0.1134, 0.1031, 0.1022, 0.1164, 0.1036, 0.1032, 0.1174
	}; 
    
    /*Vx = V * sqrt(|1-Cpx|)*/
	for (int i = 0; i < sizeof(Cp2)/sizeof(double); i++)
	{
		*v1 += V * sqrt(fabs(1-Cp1[i]));
        *v2 += V * sqrt(fabs(1-Cp2[i]));
	}
	*v1 /= sizeof(Cp1)/sizeof(double);
    *v2 /= sizeof(Cp2)/sizeof(double);
	*A	= 2.3E-1;

    /*  air density kg/m^3
    r = (Pd/(Rd*T))+(Pv/(Rv*T)). */
    *r = (Pd/(287.058*(T+273.15)))+(Pv/(461.495*(T+273.15)));

}

int main(int argc, char *	argv[])
{
	double	A, v1, v2, r, Cp1, Cp2, liftForce;

	loadInputs(&A, &v1, &v2, &r, &Cp1, &Cp2);

    /*	Fl = 1/2 * 𝜌 * a  * ((𝑣1)^2- (𝑣2)^2) */
	liftForce = r*A*(pow(v1, 2)-pow(v2, 2)) / 2.0;

	printf("Lift force = %f N\n", liftForce);

	return 0;
}
