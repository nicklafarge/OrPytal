# OrPytal

OrPytal is a Python package for orbital mechanics. It began as a pet project to find a way to easily create two-body orbits and states, without needt of parameters. For example, what if you wanted to know what an orbit would look like with a particular radius of apoapsis and specific energy? How about semilatus rectum and eccentricity? OrPytal makes these problems trivial and syntactically easy.

While this has been a fun project to sporadically work on, it is not something that I am actively working on. 


## A quick example

The following examle shows just how easy it is to create orbits and states using OrPytal

```	
import orpytal as op

orbit = op.Orbit(op.bodies.saturn)
orbit.p = 41000 * op.units.miles
orbit.e = 0.9

# Alternatively op.Orbit(op.bodies.saturn, p=41000*op.units.miles, e=0.9)

example_state = orbit.get_state(t_since_rp=orbit.period/2)

print(example_state)
```

The resulting output is
```
-------------------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------------------
                                          Keplarian State Data                                         
-------------------------------------------------------------------------------------------------------
                                            Meta Information                                           
  Orbit Name: .................
  State Name: .................
Central Body: Saturn
        Type: Elliptic
  Equitorial: .................
   Ascending: True
-------------------------------------------------------------------------------------------------------
                                            State Parameters                                           
    Pos. Mag. |r|:  6.5983104000e+05         km |          True Anomaly:  3.1415926536e+00        rad |     
    Vel. Mag. |v|:  2.3976288910e+00       km/s |     Eccentric Anomaly:  3.1415926536e+00        rad |     
Flight Path Angle:  0.0000000000e+00        rad |          Mean Anomaly:  3.1415926536e+00        rad |     
Flight Path Angle:  0.0000000000e+00        rad |       Hyperb. Anomaly: ................. .......... |     
 Time Since Peri.:  1.0439272879e+05          s |      Arg. of Latitude: ................. .......... |     
-------------------------------------------------------------------------------------------------------
                                            Perifocal Frame                                            
         Position: ................. .......... |              Velocity: ................. .......... |     
                e: -6.5983104000e+05         km |                     e: -2.9362485470e-15       km/s |     
                p:  8.0805997111e-11         km |                     p: -2.3976288910e+00       km/s |     
                h:  0.0000000000e+00         km |                     h:  0.0000000000e+00       km/s |     
-------------------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------------------
                                            Orbit Parameters                                           
     Eccentricity:  9.0000000000e-01         nd |           Inclination: ................. .......... |     
   Semimajor Axis:  3.4727949474e+05         km |        Ascending Node: ................. .......... |     
   Semiminor Axis:  1.5137562227e+05         km |      Arg of Periapsis: ................. .......... |     
   Rad. Periapsis:  3.4727949474e+04         km |         Rad. Apoapsis:  6.5983104000e+05         km |     
 Semilatus Rectum:  6.5983104000e+04         km |      Angular Momentum:  1.5820299647e+06     km^2/s |     
      Mean Motion:  3.0093979629e-05      rad/s |        Orbital Period:  2.0878545758e+05          s |     
  Specific Energy: -5.4611930842e+01   km^2/s^2 |           Flyby Angle: ................. .......... |     
       V Infinity: ................. .......... |           TA Infinity: ................. .......... |     
-------------------------------------------------------------------------------------------------------
```
