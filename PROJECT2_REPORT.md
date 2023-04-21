<a name="br1"></a>The Computational Analysis of a Spinning Top

and

The Utilization of Euler’s Equations to Graphically Simulate its Motion.

Samuel Pilon

AE 426, Section 1, Dongeun Seo, Department of Aerospace Engineering, Embry-Riddle Aeronautical University, Daytona Beach, FL.

March 2023

` `The infamous spinning top with friction is a classic example of rotational motion, and understandingits dynamics is crucial for a variety of applications in engineering and physics. To analyze the motion of aspinning top with friction, Euler's equations of motion are utilized. Euler's equations of motion are a set ofdifferential equations that describe the rotational motion of a rigid body in three dimensions. Presented in thispaper is computational analysis of the spinning top with friction about the spin axis using Euler's equations ofmotion. The model was constructed using numerical methods, and various parameters were adjusted tosimulate the motion of the spinning top. The results obtained from the simulations demonstrate the effect offriction on the motion of the spinning top, and the computational model provides insight into the complexmodeling behavior of a spinning top. This paper will be useful for engineers and physicists working in the fieldof rotational dynamics, as it provides a detailed analysis of the spinning top with friction problem using Euler'sequations of motion.

||Introduction||<p>C were determined using an example out of</p><p>“Orbital Mechanics for Engineering Students”</p>|
| :- | :- | :- | :- |
|<p>` `The end goal of this computational</p><p>simulation is to directly identify and understand</p><p>the way friction alters the motion of a spinning top.</p><p>Within the simulation, two different types of</p><p>measurements were recorded in which resulted in</p><p>the angle and angular velocity components.</p><p>Therefore, to fully understand the scope of this</p><p>simulation, further evaluation must take place to</p><p>identify why and where these trends happen and</p><p>what makes them happen the way they do.</p><p>The location of this experiment was</p><p>simulated to be in Daytona Beach Florida where</p><p>the local gravitational acceleration due to gravity</p><p>ꢀꢁ</p><p>(g) is 9.79 ꢂ!ꢀ Using MATLAB a simulation was</p><p>created to simulate friction in a spinning top. The</p><p>the mass of the top was set to be 0.5 kg and the</p><p>body was assumed to be homogenous with a</p><p>constant density which allows the inertia tensor to</p><p>be constant in two directions. The values for A and</p>||<p>` `Using the same inertial tensor equation, a</p><p>parallel axis component is added to accounts for</p><p>the torque generated about the center of mass to</p><p>the pivotal point. To initialize the state vectors of</p><p>the three main Euler angles, assumptions were</p><p>made so that there was an initial spin1, precession2,</p><p>and nutation3 rate. It should be noted that in a non-</p><p>friction case it is assumed that the top will have</p><p>continuous motion forever, whereas a friction case</p><p>as stated in the simulation requirements, will have</p><p>a linear decrease of spin with a torque applied in</p><p>the e3 axis leading to the spin to stop in 60</p><p>seconds. The simulation is set to stop at the 5o</p><p>angle with respect to the axis of symmetry which</p><p>is set to have the top hit the floor. 4</p><p>In the recorded simulation a total of six graphs</p><p>were plotted to demonstration the motion of the</p><p>top. Using this information, it was then used to</p><p>analyze and develop a conclusion regarding how</p>|
1





|<p><a name="br2"></a>the top spins different in an ideal and non-ideal</p><p>case.</p>||<p>corresponding angles and description of motion.</p><p>[1]</p>|
| :- | :- | :- |
||<p>1 Precession: Angular movement of axis due to external</p><p>` `torques over time.</p><p>2 Spin: Rotational motion around an object's axis</p><p>` `characterized by angular velocity and momentum.</p><p>3 Nutation: Periodic variation in axis orientation due to</p><p>` `external gravitational forces acting on non-spherical</p><p>` `shape.</p><p>4 85 degrees relative to top from 90 degrees.</p>||<p>The description of motion in an angular setting is</p><p>represented by ω. Which in a physical form is the</p><p>rate of change in precession and nutation with</p><p>respect to time.</p><p>(1)</p>|
||Simulation Method||<p>Unit vectors can be incorporated to keep the</p><p>convention seen in the diagram consistent which</p><p>ultimately results in equation 2.</p>|
|<p>` `The simulation method was deemed</p><p>successful with the data and insight that was</p><p>obtained. The method that was simulated use laws</p><p>of attitude dynamics that can be used to accurately</p><p>describe the motion and angular velocity of the</p><p>spinning top. To make sure that the simulation</p><p>stays consistent, initial values were set and</p><p>maintained thought-out the simulation.</p>||<p>(2)</p><p>However, since it is in a body following frame, the</p><p>relative angular velocity with respect to the co</p><p>moving system can be described as.</p><p>(3)</p>|
||Variable Value Units||

|<p>g</p><p>m</p><p>d</p><p>A</p><p>C</p><p>Initial</p><p>nutation</p>||<p>9\.79</p><p>.5</p><p>.5</p><p>12e-4</p><p>4\.5e-4</p><p>1e-8</p>||<p>(m/s^2)</p><p>(kg)</p><p>(m)</p><p>(kg-m^2)</p><p>(kg-m^2)</p><p>degrees</p>|
| :- | :- | :- | :- | :- |

||||<p>Ultimately simplifying to equation 4 which</p><p>represent the angular velocity of the spinning top.</p><p>(4)</p><p>To simulate the contributing factors that attribute</p><p>to the slowdown of the top the forces on the</p><p>system were analyzed and factored into</p><p>consideration.</p>|
| :- | :- | :- | :- |
|<p>Table 1: Represented in the table are the variables</p><p>that were kept constant.</p>||(5) F = −mꢃKꢁ = −mgEꢂꢄ|
|<p>To accurately describe the motion of a spinning</p><p>top a convention must be chosen and continuously</p><p>used throughout the simulation. In Image 1, a</p><p>diagram is set in place to allow the convention to</p><p>stay consistent.</p>||<p>(6) ꢃ = −ꢄꢅ sin θ ꢆꢇꢅ − ꢄꢅ cos θ ꢆꢇꢄ</p><p>Using the relationship shown in equation 6 can be</p><p>incorporated into equation 7 which shows the net</p><p>moment that acts on the body resulting in equation</p><p>8\.</p>|
(7) ꢈꢉꢉꢉꢉꢉꢆꢉꢉꢇꢉꢉꢈ⃗ = (ꢊ × ꢃ) + ꢈꢉ

(8) ꢀꢀꢀꢀꢀꢀꢀꢀꢀꢀꢀꢀꢀꢀꢀꢀꢀꢀꢀꢉꢈꢉꢉꢉꢉꢆꢉꢉꢇꢉꢉꢈ⃗ = ꢄꢅꢋ sin θ ꢆꢇꢊ + ꢃ ꢆꢇꢄ

|<p>Image 1: Represented in the image is a diagram of</p><p>the spinning top problem labeled with the</p>||<p>Since it was assumed that the inertial tensor would</p><p>consist of an AAC diagonal vector, it can be noted</p><p>by equation 9, and then the resulting angular</p><p>momentum can be calculated using equation 10.</p>|
| :- | :- | :- |
2





|<p><a name="br3"></a> ꢎ 0 0</p><p>(9) ꢌ = ꢍ0 ꢎ 0ꢐ</p><p>` `0 0 ꢏ</p><p>ꢔ̇</p><p>ꢎ 0 0</p><p>ꢉ⃑ = ꢌꢒ = ꢍ ꢕ̇ꢖꢗꢘꢔ</p><p>(10). ꢑ 0 ꢎ 0ꢐ ꢓ ꢜ</p><p>0 0 ꢏ</p><p>ꢕ̇ꢙꢚꢖꢔ + ꢛ̇</p><p>The resulting matrix of equation 10, is noted in</p><p>equation 11, this was calculated using matrix</p><p>multiplication.</p>||<p>ꢂꢃꢄ sin ꢅ</p><p>(15)ꢀꢀꢀꢀꢀꢁ ꢇ =</p><p>0</p><p>ꢆ</p><p>ꢉꢅ̈</p><p>ꢉꢊ̈ꢋꢌꢍꢅ + ꢉꢊ̇ꢅ̇ꢎꢏꢋꢅ</p><p>ꢀꢀꢀꢀꢀꢀꢀꢀꢀꢀꢀꢀꢀꢀꢀꢀꢈ ꢒꢀ</p><p>` `ꢐꢊ̇ꢎꢏꢋꢅ − ꢐꢊ̇ꢅ̇ꢋꢌꢍꢅ + ꢑ̈</p><p>ꢐꢊ̇ꢋꢌꢍꢅ(ꢐꢊ̇ꢎꢏꢋꢅ + ꢐꢑ̇) − ꢉꢊ̇"ꢋꢌꢍꢅꢎꢏꢋꢅ</p><p>+ ꢈ</p><p>−ꢅ̇(ꢐꢊ̇ꢎꢏꢋꢅ + ꢐꢑ̇) + ꢉꢊ̇ꢅ̇ꢎꢏꢋꢅ ꢒ</p><p>0</p>|
| :- | :- | :- |
|<p>ꢎꢔ̇</p><p>ꢉ⃑ = ꢓ ꢎꢕ̇ꢖꢗꢘꢔ</p><p>(11) ꢑ ꢜ</p><p>ꢏꢕ̇ꢙꢚꢖꢔ + ꢏꢛ̇</p><p>Following the calculation of the angular</p><p>momentum the net torque of the system can be</p><p>calculated by taking the time derivative with</p><p>respect to each variable. This results in the net</p><p>torque matrix represented in equation 12.</p>|||

|(12)||<p>ꢋ</p><p>ꢋꢈ</p>||<p>ꢎꢔ̈</p><p>ꢉ⃑ = ꢓ ꢎꢕ̈ꢖꢗꢘꢔ + ꢎꢕ̇ꢔ̇ꢙꢚꢖꢔ</p><p>ꢑ ꢜ</p><p>ꢏꢕ̇ꢙꢚꢖꢔ − ꢏꢕ̇ꢔ̇ꢖꢗꢘꢔ + ꢛ̈</p>|
| :- | :- | :- | :- | :- |

|<p></p><p>Using Euler equations of motion, a system of</p><p>equations can be created to represent the</p><p>dependencies of motion on each of the Euler</p><p>angles of the spinning top. This system is</p><p>developed using the initial condition for the total</p><p>net moment on the system. This can be represented</p><p>in equation 13 and solved in equation 14.</p>||<p>The system is now set up to be pulled apart into 3</p><p>separate system of differential equations. This is</p><p>noted in equations 16-18.</p><p>(16) ꢅ̈ꢀ=ꢀ ! (ꢂꢃꢄ sin ꢅ + (ꢉ − ꢐ)ꢀꢊ̇"ꢀꢋꢌꢍꢅꢎꢏꢋꢅ −</p><p>$</p><p>ꢀꢐꢊ̇ꢑ̇ꢋꢌꢍꢅ)</p><p>̇ ꢐ̇ꢑꢒꢂꢐꢓꢌꢏ̇ꢐ̇ꢍꢎꢏ̇ꢐ̇ꢑꢒꢂꢐ</p><p>` `ꢎꢂꢔꢆꢐ</p><p>(17) ꢕ̈=ꢀ (ꢌꢍꢎ)ꢏ</p><p>(18)ꢀꢀꢀꢀꢀꢀꢀꢀꢀꢀꢀꢀꢛ̈ꢀ=ꢀ ꢊ (ꢃ − ꢏꢀꢕ̈ꢙꢚꢖꢔ − ꢀꢏꢕ̇ꢔ̇ꢖꢗꢘꢔ)</p><p>ꢌ</p><p>It is crucial to this type of problem that the DCM</p><p>be found prior to this calculation, the reason for</p><p>this is due to the nature of the system. The</p><p>conversion must take place to explain the motion</p><p>of the spinning top fully and accurately by</p><p>changings the relative attitude with respect to the</p><p>reference frame.</p><p>This DCM is known as a typical 3-1-3 rotation</p><p>with precession, nutation, and spin as the rotation</p><p>angles. This can be visualized in equation 19.</p>|
| :- | :- | :- |
|<p>(13) ꢀꢀꢀꢀꢀꢀꢀꢀꢀꢀꢀꢀꢀꢀꢀꢀꢀꢀꢀꢈꢉꢉꢉꢉꢉꢆꢉꢉꢇꢉꢉꢈ⃗ =</p><p>ꢋ</p><p>ꢑ</p><p>ꢉ⃑ + Ωꢀ ×ꢀꢑꢉ⃑</p><p>ꢋꢈ</p><p>ꢂꢃꢄ sin ꢅ</p><p>(14)ꢀꢀꢀꢀꢀꢁ ꢇ =</p><p>0</p><p>ꢆ</p><p>ꢉꢅ̈</p><p>ꢉꢊ̈ꢋꢌꢍꢅ + ꢉꢊ̇ꢅ̇ꢎꢏꢋꢅ</p><p>ꢈ ꢒ</p><p>ꢐꢊ̇ꢎꢏꢋꢅ − ꢐꢊ̇ꢅ̇ꢋꢌꢍꢅ + ꢑ̈</p><p>` `ꢓ! ꢓ" ꢓ#</p><p>ꢅ̇ꢊ̇ꢋꢌꢍꢅ ꢊ̇ꢎꢏꢋꢅ</p><p>+ ꢈ ꢒ</p><p>ꢉꢅ̇ꢉꢊ̇ꢋꢌꢍꢅ ꢐꢊ̇ꢎꢏꢋꢅ + ꢐꢑ̇</p>||<p>(19)</p><p>ꢎꢏꢋ(ꢅ) −ꢋꢌꢍ(ꢅ) 0 ꢎꢏꢋ(ꢑ) −ꢋꢌꢍ(ꢑ) 0</p><p>1 0 0</p><p>ꢁꢋꢌꢍ(ꢅ) ꢎꢏꢋ(ꢅ) 0ꢇ ꢁ ꢇ ꢁꢋꢌꢍ(ꢑ) ꢎꢏꢋ(ꢑ) 0ꢇ</p><p>0 ꢎꢏꢋ(ꢊ) −ꢋꢌꢍ(ꢊ)</p><p>0 0 1 0 0 1</p><p>0 ꢋꢌꢍ(ꢊ) ꢎꢏꢋ(ꢊ)</p><p>These matrices are multiplied in the order of the</p><p>sequence and then multiplied with the angular</p><p>velocity of the body frame which was noted in</p><p>equation 4. The result of this gives the angular</p><p>velocity of the body frame with respect to that of</p><p>the inertial frame.</p>|
|<p>Equation 14 can be simplified and re written to be</p><p>more easily understood which is noted in equation</p><p>15\.</p>||<p>Using ode45() the system of three differential</p><p>equations was solved simultaneously and then</p><p>plotted against time until the stopping condition</p><p>was met. To use the ode45() function in</p><p>MATLAB, the system was split into multiple first</p><p>order differential equations. The “spinningtop.m”</p><p>function was created and called using</p><p>spinningtop(k) this contained the angle and the rate</p>|
3





|<p><a name="br4"></a>at each given time instance. Then the ode was</p><p>solved and passed to the original script where it</p><p>was then plotted.</p><p>Discussion On Results</p><p>The output of the graphs for the spinning</p><p>top shows some interesting trends that occurred</p><p>over a span of about 60 seconds. The precession</p><p>graph exhibits a sine wave trend and shows an</p><p>increase in precession over time. This sine wave</p><p>trend is caused by the top rotating around a vertical</p><p>axis while being pulled by a force that causes the</p><p>axis to rotate around a horizontal axis. The top</p><p>moves in a circular motion around the vertical axis</p><p>because of this force, which is perpendicular to the</p><p>axis of rotation. The periodic nature of this motion</p><p>results in a sine wave trend, with the amplitude and</p><p>frequency depending on the force and inertial</p><p>characteristics of the top.</p>||<p>nutation motion, the angular velocity experiences a</p><p>sudden change, resulting in a spike in the graph.</p><p>Then, after more time the angular velocity would</p><p>eventually decreases linearly due to the</p><p>conservation of angular momentum, causing the</p><p>graph to have a repeating pattern of a spike and</p><p>then a linear decrease over time.</p><p>On the other hand, the angular velocity of</p><p>the precession graph also shows a repeating pattern</p><p>of a spike and then a linear decrease with respect</p><p>to time. This pattern is because the top experiences</p><p>a sudden change in angular velocity when it</p><p>reaches the maximum of its precession motion,</p><p>resulting in a spike in the graph. Then, the angular</p><p>velocity decreases linearly due to the conservation</p><p>of angular momentum, causing the graph to follow</p><p>a similar pattern as the nutation graph. The</p><p>amplitude and frequency of the spikes and the</p><p>overall trend of the graph depend on the force</p><p>acting on the top and its inertial properties.</p>|
| :- | :- | :- |
|<p>` `The nutation graph started small and</p><p>exponentially got bigger over time. This happens</p><p>because the top's axis of rotation starts to wobble</p><p>due to an external force acting on it, which is the</p><p>reason why outside the simulation tops tend to</p><p>deviate from their original axis of rotation. Within</p><p>the bounds of the simulation, the graphs show a</p><p>trend where it starts small and exponentially gets</p><p>bigger over time because the top's axis of rotation</p><p>starts to wobble due to the precession motion. This</p><p>wobbling is known as nutation, and it occurs when</p><p>the axis of rotation changes its orientation in space,</p><p>while still maintaining its precession motion. The</p><p>exponential increase is since nutation is a self-</p><p>reinforcing motion, and the amplitude builds up</p><p>over time.</p><p>The spin graph had a very similar trend to that of</p><p>the precession in that it was a sine wave. This</p><p>happens because the top's angular momentum is</p><p>conserved, so as the precession increases, the spin</p><p>decreases, and vice versa. This relationship</p><p>between the precession and spin motion causes the</p><p>sine wave trend in the spin graph.</p><p>` `The angular velocity of the nutation graph</p><p>followed a very similar trend to that of the nutation</p><p>graph itself. As the top's axis of rotation starts to</p><p>wobble, the angular velocity of nutation starts</p><p>small and exponentially increases due to the self-</p><p>reinforcing nature of nutation. At the peak of the</p>||<p>Lastly, the angular velocity of the spin held a</p><p>linear decrease over time as that is due to friction.</p><p>The angular velocity of the spin graph holds a</p><p>linear decrease over time due to friction. The</p><p>friction between the top and the surface it is</p><p>spinning on causes a loss of energy, which results</p><p>in a decrease in the spin velocity over time.</p><p>Concluding Statements</p><p>` `In the experiment, the rotational motion of</p><p>a spinning top with friction was modeled and</p><p>examined using Euler angles. The precession,</p><p>nutation, and spin of the top were successfully</p><p>captured by the simulation. Yet, the presumption</p><p>that the pivotal point stays stable throughout</p><p>rotation resulted in the discovery of a significant</p><p>source of error.</p><p>Three angles known as Euler angles are used to</p><p>describe how a rigid body is oriented in three-</p><p>dimensional space. These were used in the</p><p>experiment to replicate the rotational motion of the</p><p>top. The simulation allowed the researchers to</p><p>learn important things about the top's behavior and</p><p>interactions with its surroundings.</p>|
4



|<p><a name="br5"></a>The assumption that the pivotal point stays</p><p>stationary was found to be the main cause of</p><p>inaccuracy, despite the simulation producing</p><p>accurate plots. In fact, when the top rotates, the</p><p>pivotal point shifts, which changes how the top</p><p>moves. This inaccuracy draws attention to a flaw</p><p>in the model, but it has no impact on the</p><p>experiment's overall importance.</p>||<p>References</p><p>[1] Curtis, H. D. (2021). Orbital Mechanics for</p><p>` `Engineering Students: Revised fourth</p><p>` `edition. Elsevi</p>|
| :- | :- | :- |
The experiment proved useful in providing athorough investigation of the motion and behaviorof the spinning top notwithstanding the discoveredinaccuracy. The parameters that govern the mobility of the top and the effects of friction werebetter understood by the researchers. Futureadjustments to the model might focus on how thepivot moves, which would result in even moreprecise simulations of spinning tops with friction.

5




<a name="br6"></a>Appendix

SAMUELꢀRꢀPILON

AEꢀ426ꢀ­ꢀProject­2ꢀDateꢀ­ꢀ03/29/2023ꢀForꢀpersonalꢀuseꢀofꢀSAMUELꢀPILONꢀonly.ꢀNotꢀtoꢀbeꢀdistributed

**Contents**

QUESTIONꢀSTATEMENT

MAINꢀFUNCTIONꢀ(SPINNINGꢀTOP)

SETꢀGIVENꢀVALUES

EULERSꢀEQUATIONSꢀOFꢀMOTION

DIRECTꢀCOSINEꢀMATRIX

ODEꢀ45ꢀ(ORDINARYꢀDIFFERENTIALꢀEQUATIONꢀSOLVER)

PLOTS

**QUESTIONꢀSTATEMENT**

%{

Derive the equation of motion using Euler's rotational equation in terms ofprecession (φ), nutation (θ), and spin (ψ) angles. The center-of-mass of the topis located at r\_G = lˆe3 in Fig. 1.

Consider that the friction with the oor is adding resistance in such a waythat there is a linear decrease of the spin with a torque applied in the ˆe3 axisleading to the spin stop in 1 minute. Graph the solution to the non-linearequations using MATLAB as a function of precession, nutation and spin andstop the motion when the surface of the top (5º with respect to the axis ofsymmetry) touches the oor. Explain the physical motion of the top with thesimulation.

%}

**MAINꢀFUNCTIONꢀ(SPINNINGꢀTOP)**

%{

This program uses numerical methods to solve the equations of motionknown as Euler's equations, which describe the rotational motion of aspinning top. By integrating these equations over time, the programcan predict the future orientation and angular velocity of the top basedon its initial conditions and the forces acting upon it.

The quaternion allows for efficient and accurate updates to the top'sorientation as it undergoes rotational motion, providing a time history ofits orientation throughout the simulation.

User subfunction required: Spinningtop.m

%}

clear all; close all; clc

**SETꢀGIVENꢀVALUES**

ThisꢀsectionꢀofꢀcodeꢀinitializesꢀvariablesꢀandꢀconstantsꢀrelatedꢀtoꢀaꢀSpiningꢀtopꢀproblem

A = 12.e-4; % moment of inertia about one axis of rotation (kg\*m^2)C = 4.5e-4; % moment of inertia about another axis of rotation (kg\*m^2)m = 0.5; % mass (kg)

g = 9.79; % acceleration due to gravity (m/s^2)

d = 0.05; % distance from origin to center of mass (m)

initial\_nutations = [1e-8]; % initial nutations in degrees

tspan = linspace(0,60,100); %time span of 60 seconds incremented in steps of 100

syms precession spin nutation precessiondot spindot nutationdot thetaddot psiddot phiddot mu

**EULERSꢀEQUATIONSꢀOFꢀMOTION**

Defineꢀleft­handꢀandꢀright­handꢀsidesꢀofꢀEulerꢀequationsꢀofꢀmotion

Momentright\_body = [A\*thetaddot, A\*(phiddot\*sin(nutation) + precessiondot\*nutationdot\*cos(nutation)), ... C\*(psiddot + phiddot\*cos(nutation) - precessiondot\*nutationdot\*sin(nutation))] + ...

cross([nutationdot, precessiondot\*sin(nutation), precessiondot\*cos(nutation)], ...

` `[A\*nutationdot, A\*precessiondot\*sin(nutation), C\*(spindot + precessiondot\*cos(nutation))]); % torques acting on the spinning topMomentleft\_body = [m\*g\*d\*sin(nutation), 0, -mu]; % forces acting on the spinning top

**DIRECTꢀCOSINEꢀMATRIX**

%Relate the coordinates of a vector in one reference

%frame to the coordinates of the same vector in another reference frame.

%Represents a rotation around the z-axis by an angle nutation.

DCM\_1\_nutation = [cos(nutation),-sin(nutation),0;sin(nutation),cos(nutation),0;0,0,1];

%Represents a rotation around the x-axis by an angle precession.

6




<a name="br7"></a>DCM\_2\_precession = [1,0,0;0,cos(precession),-sin(precession);0,sin(precession),cos(precession)];

%Represents a rotation around the z-axis by an angle spin.

DCM\_3\_spin = [cos(spin),-sin(spin),0;sin(spin),cos(spin),0;0,0,1];

%create a single 3x3 matrix DCM that represents the overall%rotation from the inertial frame to the body frame.DCM = DCM\_1\_nutation\*DCM\_2\_precession\*DCM\_3\_spin;

%Represents the angular velocity of the body frame relative%to the inertial frame, expressed in the body frame

w\_body = [nutationdot,precessiondot\*sin(nutation),precessiondot\*cos(nutation)];

%Computes the angular velocity vector by transforming from the body frame to the inertial frame

%using the DCM matrix.w\_inertial = w\_body\*DCM;

% Solve for the angular accelerations

thetaddot = solve(Momentleft\_body(1) == Momentright\_body(1), thetaddot)phiddot = solve(Momentleft\_body(2) == Momentright\_body(2), phiddot)psiddot = solve(Momentleft\_body(3) == Momentright\_body(3), psiddot)

**ODEꢀ45ꢀ(ORDINARYꢀDIFFERENTIALꢀEQUATIONꢀSOLVER)**

figure('Name', 'Euler Angles and Their Rates', 'color', [1 1 1])

opt = odeset('RelTol', 1e-9,'AbsTol',1e-9); % Sets the relative and absolute error tolerances

for p = 1

initial\_nutations = initial\_nutations(p); % Setting the initial nutation angle for this iterationMu\_initial = 0; % Setting the initial friction factor to zero

% Initial conditions in order: precession, nutation, spin, precessiondot, nutationdot, spindot

ini = [0, initial\_nutations\*pi/180, 0, 0, 0, 100\*pi/180]; % Converting initial nutation angle from degrees to radians

% Increasing friction factor until the spin stops

` `defining\_angle = [9,9;9,9]; % Defining an angle matrix to fall in loopwhile defining\_angle(end,end) > 0 % While the spin angular velocity is greater than zero Mu\_initial = Mu\_initial + 0.001; % Increase the friction factor

` `[t,defining\_angle] = ode45(@(t,k) Spiningtop(t,k,Mu\_initial), tspan, ini, opt); % Solve the ODE using the current friction factorend

% This section of the code is used to identify when the nutation angle reaches 85 degrees% and then to extract the corresponding time and angle values up to that point,% effectively stopping the motion.

t\_index = find(defining\_angle(:,2) > 85\*pi/180); % Find the index where the nutation angle reaches 85 degreest\_index = t\_index(1); % Get the first index (i.e., the earliest time) where the nutation angle exceeds 85 degrees

t = t(1:t\_index); % Keep time values up to the point where nutation angle exceeds 85 degrees

defining\_angle = defining\_angle(1:t\_index,:); % Keep angle values up to the point where nutation angle exceeds 85 degrees

**PLOTS**

% Plots

figure(p); % Create a new figure for each initial nutation angle

titles = {'Precession (φ)', 'Nutation (θ)', 'Spin (ψ)', 'Angular Velocity (Precession)', 'Angular Velocity (Nutation)', 'Angular Velocity (Spin)'}

for i = 1:6

% Create subplot for current angle/velocity subplot\_index = subplot(2,3,i); subplot\_handle = subplot(subplot\_index);

% Plot angle/velocity vs time with red color and linewidth of 1.5 plot\_handle = plot(t, defining\_angle(:,i), 'r-', 'LineWidth', 1.5);

% Customize subplot properties

set(subplot\_handle, 'GridLineStyle', ':', 'GridColor', 'k', 'XMinorGrid', 'on');

xlabel('Time (s)');

if i < 4

ylabel('Angle (rad)');

else

ylabel('Angular Velocity (rad/s)');

end

title(titles(i));

end

if p == 1

` `sgtitle('Spinning Top Simulation for Small Initial Nutation Angle (0.01e-5)'); % Add a super title to the first figureend

end

*ꢀ*

*PublishedꢀwithꢀMATLAB®ꢀR2020a*

7




<a name="br8"></a>**Contents**

MOTIONꢀOFꢀTHEꢀSPININGꢀTOP

**MOTIONꢀOFꢀTHEꢀSPININGꢀTOP**

Thisꢀfunctionꢀdefinesꢀtheꢀsystemꢀofꢀordinaryꢀdifferentialꢀequationsꢀ(ODEs)ꢀthatꢀdescribeꢀtheꢀmotionꢀofꢀtheꢀspinningꢀtop.

INPUTS:ꢀ­ꢀt:ꢀcurrentꢀtimeꢀ­ꢀk:ꢀvectorꢀofꢀcurrentꢀvaluesꢀforꢀtheꢀsixꢀstateꢀvariablesꢀ[phi,ꢀtheta,ꢀpsi,ꢀphi\_dot,ꢀtheta\_dot,ꢀpsi\_dot]ꢀ­ꢀfriction:theꢀfrictionꢀfactorꢀofꢀtheꢀspinningꢀtop

OUTPUT:ꢀ­ꢀdifft:ꢀvectorꢀofꢀtheꢀtimeꢀderivativeꢀofꢀtheꢀsixꢀstateꢀvariables

TheꢀfunctionꢀusesꢀtheꢀinputꢀvaluesꢀtoꢀcomputeꢀtheꢀtimeꢀderivativesꢀofꢀtheꢀsixꢀstateꢀvariablesꢀinꢀorderꢀtoꢀsolveꢀtheꢀODEs.ꢀTheequationsꢀforꢀtheseꢀtimeꢀderivativesꢀwereꢀderivedꢀfromꢀtheꢀequationsꢀofꢀmotionꢀforꢀtheꢀspinningꢀtop,ꢀtakingꢀintoꢀaccountꢀtheꢀeffectsꢀofprecession,ꢀnutation,ꢀandꢀspin.ꢀTheꢀfrictionꢀfactorꢀisꢀalsoꢀincludedꢀinꢀtheꢀequations.

function difft = Spiningtop(t,k,friction)

difft(1) = k(4);

difft(2) = k(5);

difft(3) = k(6);

difft(4) = (3\*k(6)\*k(5) - 13\*k(4)\*k(5)\*cos(k(2)))/(8\*sin(k(2)));difft(5) = (5\*cos(k(2))\*sin(k(2))\*k(5)^2)/8 - (3\*k(4)\*sin(k(2))\*((k(6))/8 ... + (4895\*k(4)\*cos(k(2)))/24));

difft(6) = k(4)\*k(5)\*sin(k(2)) - difft(4)\*cos(k(2)) - (20000\*friction)/9;

difft = difft';

end

*ꢀ*

*PublishedꢀwithꢀMATLAB®ꢀR2020a ꢀ*

8




<a name="br9"></a>thetaddot =

(5\*cos(nutation)\*sin(nutation)\*precessiondot^2)/8 - (3\*spindot\*sin(nutation)\*precessiondot)/8 + (4895\*sin(nutation))/24

phiddot =

(3\*nutationdot\*spindot - 13\*nutationdot\*precessiondot\*cos(nutation))/(8\*sin(nutation))

psiddot =

nutationdot\*precessiondot\*sin(nutation) - phiddot\*cos(nutation) - (20000\*mu)/9

9
