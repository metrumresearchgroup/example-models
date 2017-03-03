<b> Torsten </b> is a library of C++ functions that support applications of Stan in Pharmacometrics. The current prototype provides:
* One and Two Compartment Model functions that compute solutions analytically
* General Linear Compartment Model function that computes a matrix exponential solution
* General Compartment Model functions that compute solutions numerically
  * Runge-Kutta 4th/5th order for non-stiff ODE systems
  * Backward Differentiation for stiff ODE systems
  
** This prototype is still under development ** and has been uploaded to facilitate working with the community of Stan developers. The current version was written by Charles Margossian, Bill Gillespie, and Metrum Research Group, LLC. We have recieved extensive help from the Stan development team.

See the user manual (`torstenManual.pdf`) for more information and guidance on the examples. If you have any questions, please raise an issue on GitHub or send me an e-mail at charlesm@metrumrg.com. 

Licensing
---------
The Torsten library is open-source and licensed under the BSD 3-clause license. 


Install
-------
To install cmdStan with torsten, run the shell script script `setupTorsten.sh`.

We are working with Stan's development team to create a system to add and share Stan packages. In the mean time, users can download a forked version of Stan with Torsten from GitHub. The latest version of Torsten (v0.82) is compatible with Stan v2.14.0. Torsten is agnostic to which Stan interface you use. The setupTorsten file installs Torsten with cmdStan.


Examples
---------
For each model, we provide the following files:
* *modelName*.stan
* *modelName*.data.R
* *modelName*.init.R
* *modelName*Simulation.R 

The simulation file can be used to create the data and initial estimate files. 

There are three Stan files for the  two compartment model: `TwoCptModel`, `LinTwoCptModel`, and `GenTwoCptModel`. This is probably a good place to start. `effCptModel` tackles a PKPD model based on a linear ODE system, and the Friberg-Karlsson model a semi-mechanist model described by a nonlinear ODE system. Finally, `TwoCptModelPopulation` extends the two compartment model to the case where we have data from multiple patients with inter-individual variability in the parameters. 

See the manual for more assistance.

Under the R directory, we provide tools to run the examples via cmdStan and look at diagnostic plots for a two compartment model and a two compartment population model.

C++ Code
--------
The C++ code for Torsten can be found on the following repos:

Math Library: https://github.com/charlesm93/math/tree/torsten-master

Stan Grammar: https://github.com/charlesm93/stan/tree/torsten-master

Updates
-------
03/02/2017
* Update examples, user manual, and other files for Torsten v 0.82.

10/31/2016
* Add Linear Compartment model function.
* Revise the user's manual.
* Torsten is now compatible with a development version, post v 2.12. 

08/02/2016
* Update Stan and Math branch to match stan-dev. This is important because of the recent bug fix in stan 2.10. 

