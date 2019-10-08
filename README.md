# Modeling of Tensegrity Structures (MOTES)

### **Welcome to **MOTES** software!**

#### General Information
 
Our research team focuses on integrating structure and control design using Tensegrity structure. The group focuses on designing tensegrity structures to meet the specified objectives. These objectives can vary from minimizing the mass of the structure to controlling the structure to meet specific performance. This software is intended to study the statics and dynamics of tensegrity systems. The authors would like to make this as open-source software to help other researchers who are also interested in this field. 

Contribution of this software is mainly in these two aspects: First, the static analysis provides the minimum mass of the tensegrity structure by optimizing for the tensile force in the strings and compressive force in the bars for no external force (self-equilibrium state) and in the presence of a given external force. The optimization problem is written as a Linear Programming to solve for the minimum mass required under yielding constraints. The software also allows solving for the minimum mass under buckling and yielding failure criteria through a non-linear optimization solver. 
Second, the dynamic analysis uses a second-order matrix differential equation to simulate the dynamics of any complexity of tensegrity structure. This dynamic model assumes the bars to be rigid and strings to show elastic behavior (Hookean).

Undergraduate linear algebra and some basic knowledge of MATLAB is required to understand the codes well. This software is developed based on:

•	64-bit Windows
•	MATLAB 
•	MATLAB Optimization Toolbox
Note: Win7/Win10/Mac OS/Linux/Win XP/Win Vista compatible with a MATLAB version later than 2009a should work fine. However, we encourage the user to run the software with the latest MATLAB release if possible. 

#### LICENSE

    /* This Source Code Form is subject to the terms of the Mozilla Public
     * License, v. 2.0. If a copy of the MPL was not distributed with this
     * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 
---

***<font size=4>Modeling of Tensegrity Structrues (MOTES) folder contains the following parts:</font>***

---

#### setup.m 
If one wants to start using the MOTES software, run the 'setup.m' first.
Open MATLAB and run 'setup.m' file, it will:

- Add all library functions to MATLAB path, 
- Add the Software Verification and Examples,
- Add User Guide,
- Add Videos folder,
- Add JOSS paper.

Note: setup.m must be run every time MATLAB opens, before any other file.

#### JossPaper

This folder contains the journal paper corresponding to the software, including source file documents and reference. The journal paper will provide the background introduction, a summary of our work, applications, and references, etc. 

#### Function Library

This folder contains:

All the library functions for tensegrity statics and dynamics analysis are organized in this folder. By following instructions of statics and dynamics analysis from ''User_Guide.pdf", one can perform the analysis.

#### Software Verification and Examples

This folder contains:

1. Dynamics Examples

Here, we give examples to verify and demonstrate the dynamics of this software.

2. Statics Examples

Here, we give examples to verify and demonstrate the statics of this software.

#### User Guide

This folder contains "User_Guide.pdf". The file provides a detailed description, how to use the software, various applications, and how to become a developer as well.

#### Videos
Some interesting tensegrity animation examples are shown in this folder.

---

### Help Desk:

We are open and willing to answer any question. Please state your problem clearly and use the following emails to contact:<br>
Raman Goyal: <ramaniitrgoyal92@tamu.edu>, Muhao Chen: <muhaochen@tamu.edu>. Thank you!

----

### Join MOTES Community and Contribute

#### How to contribute

Feedbacks and contributions are appreciated. Please use the same nomenclature so that everybody could be on the same page.

1. Fork it
2. Submit a pull request OR send emails to help desk.

We will reply to you ASAP.

#### Coding standards

* MATLAB (>= 2009a)
* Function input and output comments
* Use the same Nomenclature as follows

#### Nomenclature

##### Geometry: 
    N: initial node positions
    C_b: bar connectivity
    C_s: string connectivity
    C_r: Bar center of mass 'connectivity'
    C_sb: string connections to bars
    C_ss: string connections to string nodes
    C_bb: non-zero columns of C_b (old C_b matrix)
    C_nb: satisfies R_b = N*C_nb^T*C_r'
    B: bar matrix
    S: string matrix
    s_0: resting string lengths
    n: number of nodes
    beta: number of bars
    alpha: number of string members
    sigma: number of string point mass nodes
##### Force
    gamma: string member force densities
    lambda: bar member force densities 
    P: constraint matrix
    D: constraint matrix
    M: mass matrix used in dynamics
    Minv: inverse of M
    W: external node forces
    k: string stiffness coefficients
    constants: various constants throughout the simulation
##### Initial Conditions
	Nd0: initial node velocities
    m: bar masses
    ms: string node point masses
    mgyro: gyro wheel masses
    gyro_r: gyro wheel radii
    gyro_h: gyro wheel heights
    gyro_omega: gyro angular rates
    Jt_hat: transaxial bar+gyro moments of inertia
    Ja_hat: axial bar+gyro moments of inertia
    len_hat: initial bar lengths
##### Simulation
    tf: simulation time duration
    dt: simulation time step
    W: external node forces
    k: string stiffness coefficients
    N0: initial node positions
    Nd0: initial node velocities
    sim: integration specification variables
