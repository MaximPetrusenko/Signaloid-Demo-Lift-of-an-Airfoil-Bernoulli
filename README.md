# Signaloid-Demo-Lift-of-an-Airfoil-Bernoulli

Model for lift generation of a 2D NACA 2412 airfoil using Bernoulli equation

[<img src="https://assets.signaloid.io/add-to-signaloid-cloud-logo-dark-v6.png#gh-dark-mode-only" alt="[Add to signaloid.io]" height="30">](https://signaloid.io/repositories?connect=https://github.com/MaximPetrusenko/Signaloid-Demo-Lift-of-an-Airfoil-Bernoulli#gh-dark-mode-only)
[<img src="https://assets.signaloid.io/add-to-signaloid-cloud-logo-light-v6.png#gh-light-mode-only" alt="[Add to signaloid.io]" height="30">](https://signaloid.io/repositories?connect=https://github.com/MaximPetrusenko/Signaloid-Demo-Lift-of-an-Airfoil-Bernoulli#gh-light-mode-only)

# Lift generation model based on Bernoulli equation
This example shows how uncertainties in empirical model parameters affect the uncertainty distribution of the model's output, for a model of a physical process. The example calculates lift force over 2-dimensional airfoil. The equation for lift force derived from the Bernoulli equation, therefore it is assumed that fluid(air) flow is inviscid and incompressible, which means that it is not compressible (density changes in the air at different points on an airfoil are suggested to be small and not be taken into the account) and has no viscosity (internal friction). 
Other assumptions are:
- Reference levels for height and pressure are h0 = 0, P0 = 1 atm (101325 Pa)
- Attitude 0°
- Possible angles of attack: 0°, 5°, 10° (empirical data based on those AOA)

## Inputs
The inputs and their ranges are:
  -	`A`:		0.1 to 1 m^2 - area of the 2D airfoil
  -	`h`:	  0.0 to 11019 m - elevation (troposphere)
  - `T`:	  -50 °C to 50 °C - ambient temperature
  -	`V`:		10 to 343 m/s - velocity below supersonic speed
  -	`Сp1`:	~-2.8 to ~1.0 - coefficient for pressurre distribution over an airfoil
  -	`Сp2`:	~-0.54 to ~1.14 - coefficient for pressurre distribution under an airfoil 
  
![Coefficient-of-Pressure-vs-Percent-Chord-Length](https://user-images.githubusercontent.com/72452482/220147400-506b3916-9d18-4a7d-9139-2214d3a66bc1.png)

## Outputs
Lift force is being calculated by multiplying the net pressure by the wing area.
Fl = 1/2 * 𝜌 * a  * ((𝑣2)^2- (𝑣1)^2)
 
## Repository Tree Structure
The repository has all the model parameters as point-valued numbers, with the 𝑣1 and 𝑣2 computed as the mean values from the empirical distributions of pressure coefficients [^0] over and under an airfoil respectively.

```
.
├── README.md
├── v1
│   ├── README.md
│   └── src
│       ├── README.md
│       └── lift-2D-airfoil-Bernoulli-no-distributions.c

```
<br/>

[^0]: Matsson, O. & Voth, John & McCain, Connor & McGraw, Connor. (2016). [Aerodynamic Performance of the NACA 2412 Airfoil at Low Reynolds Number](https://www.researchgate.net/publication/319271205_Aerodynamic_Performance_of_the_NACA_2412_Airfoil_at_Low_Reynolds_Number).
