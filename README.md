# Modeling of Tensegrity Structures (MOTES)

### **Welcome to **MOTES** software!**

#### General Information
 
Our research team focuses on integrating structure and control design using Tensegrity structure. The group focuses on designing tensegrity structures to meet the specified objectives. These objectives can vary from minimizing the mass of the structure to controlling the structure to meet specific performance. This software is intended to study the statics and dynamics of tensegrity systems. The authors would like to make this as open-source software to help other researchers who are also interested in this field. 

Contribution of this software is mainly in these two aspects: First, the static analysis provides the minimum mass of the tensegrity structure by optimizing for the tensile force in the strings and compressive force in the bars for no external force (self-equilibrium state) and in the presence of a given external force. The optimization problem is written as a Linear Programming to solve for the minimum mass required under yielding constraints. The software also allows solving for the minimum mass under buckling and yielding failure criteria through a non-linear optimization solver. 
Second, the dynamic analysis uses a second-order matrix differential equation to simulate the dynamics of any complexity of tensegrity structure. This dynamic model assumes the bars to be rigid and strings to show elastic behavior (Hookean).


#### LICENSE

    /* This Source Code Form is subject to the terms of the Mozilla Public
     * License, v. 2.0. If a copy of the MPL was not distributed with this
     * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 
---

***<font size=4>Modeling of Tensegrity Structrues (MOTES) folder contains the following parts:</font>***

---

#### setup.m 
If one wants to start using the MOTES software, 'setup.m' is to be run only the first time.
Open MATLAB and run 'setup.m' file, it will:

- Add all library functions to MATLAB path, 
- Add the Software Verification and Examples,
- Add User Guide,
- Add Videos folder,
- Add JOSS paper.

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
