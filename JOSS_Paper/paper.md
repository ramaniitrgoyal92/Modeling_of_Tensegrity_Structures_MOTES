---

title: 'TEAM: Tensegrity Engineering Analysis Master'

tags:
  - Tensegrity systems
  - Multibody dynamics
  - Flexible structures
  - Prestressable structures
  - Statics Analysis
  - MATLAB

authors:
  - name: Raman Goyal
    orcid: 0000-0002-8128-3051
    affiliation: 1

  - name: Muhao Chen
    orcid: 0000-0003-1812-6835
    affiliation: 1

  - name: Manoranjan Majji
    orcid: 0000-0002-8218-2631
    affiliation: 1

  - name: Robert E. Skelton
    orcid: 0000-0001-6503-9115
    affiliation: 1
    
affiliations:
 - name: Department of Aerospace Engineering, Texas A&M University, College Station, Texas, USA
   index: 1

date: 15 June 2019
bibliography: paper.bib
---

# Summary

Tensegrity systems [@Skelton_2009_Tensegrity_Book] dynamics is a subset of the class of Multi-body dynamics which includes cylindrical rigid bodies (bars) and elastic members (strings) arranged in a stabilizable topology. The name ‘Tensegrity’ is coined by Buckminster Fuller [@Fuller_1959]. The name comes by adding ‘Tensile +Integrity = Tensegrity’ for the art-form created by Ioganson (1921) and Snelson (1948) [@Snelson_1965]. Tensegrity structures are created by methodically arranging tensile members (strings) and compressive members (bars) to form a stable system. Skelton defines "Class-1" tensegrity system where no compressive members touch each other, and a "Class-k" system, where k compressive members are connected at a node [@Skelton_2009_Tensegrity_Book]. The multiple compressive members are connected through ball joints, causing all the members in a tensegrity structure to be axially loaded, i.e., no moment is present on any individual member. Tensegrity structures can also be prestressed to have uni-directional loading for all the individual members giving the freedom to design tension and compression members separately. This freedom provides good structural efficiency (high strength to mass ratio) to tensegrity structures [@Skelton_2009_Tensegrity_Book]. A particular topology of tensegrity structure can provide various self-equilibrium solutions corresponding to different values of prestress in the structure. The various stable prestress values provide a domain set to minimize the mass of the structure. This minimization has proved tensegrity to be an optimal mass solution for various loading conditions. Moreover, different prestress values correspond to different stiffness, allowing to change the stiffness without changing the shape. These properties of tensegrity structure have led the researchers to use tensegrity concepts in various applications from Civil engineering bridges [@Carpentieri_2015_Bridge] to robotic manipulators [@Karnan_Goyal_2017_IROS] to various space applications; landers [@Vytas_2013,Goyal_2019_Buckling], deployable structures [@Tibert_Pellegrino_2003_Mast], and, space habitat designs [@Goyal_2017_Habitat,Chen_2018_Habitat].

Most approaches in the field of Multi-body dynamics use a minimal coordinate representation. This eliminates redundant variables as the body coordinates are used based on connection, from the first body to the last body in a particular order. One disadvantage of this approach is the use of transcendental functions, which can cause significant computational errors. The other disadvantage is the limited configuration of the rigid body connections to a topological tree. The proposed formulation uses a non-minimal coordinate system to describe the dynamics of the system. The approach uses 6 DOF to define a rigid symmetric cylindrical bar with no rotation about the longitudinal axis. For a "Class-k" tensegrity structure, there would be no moment or torque along the bar axis as the bars are connected by frictionless ball joints. 

Moreover, all the strings are also connected to the node centers of the bars, resulting in no torque (no rotation) along the axis of the bar. The use of non-minimal system results in a second-order matrix differential equation to describe the dynamics of any tensegrity system. One disadvantage of such a method is the violation of the rigid body (constant bar length) constraint in the presence of computational errors. In order to satisfy the constant bar length constraint, we use an analytical correction step to eliminate the approximated errors at each integration step. 

Tensegrity Engineering Analysis Master (TEAM) provides two categories for the analysis of any tensegrity structure. First, the static analysis provides the required tension in the strings and force in the bars for a self-equilibrium state or the equilibrium condition in the presence of an external force. This analysis uses \emph{Linear Programming} to solve for the minimum mass required under yielding constraints. The software also allows solving for minimum mass under buckling and yielding failure criteria using \emph{non-linear optimization}. 

Second, the dynamic analysis uses a second-order matrix differential equation to simulate the dynamics of any complexity of tensegrity structure. The software runs a self-developed modified Runge-Kutta integration package to solve the nonlinear differential equations. A bar length correction scheme is used to correct the dynamics response that might incur random errors because of computational limitations. This analytical correction step also restricts the errors in connection constraints for class-k structures. The class-k, bar-to-bar connections (ball joints) are formulated as linear constraints. These constraints in the motion space give rise to constraint forces in the structure. An analytical solution is provided to solve for the constraint forces, and a reduced order model is developed to simulate the dynamics for this restricted motion space.
The mass in the strings is also included in the dynamics by discretizing the string into several point masses along the length of the string. This complete mathematical formulation for the tensegrity dynamics is developed in our recent dynamics paper [@Goyal_Dynamics_2019].

The software is being used to develop various tensegrity structures like Tensegrity robotic arm, Tensegrity antenna, Tensegrity lander, and Space Habitat, where we integrate structure and control design to get the required performance.
The software was also used as a part of the class curriculum for the course \textit{AERO-489/689 Design Elective: Advanced Statics and Dynamics of Flexible Structures: Tensegrity Systems} which was offered in spring 2019.

# Acknowledgements

The authors would like to thank NASA NIAC Phase II grant, Dr. Edwin Peraza Hernandez (UC Irvine), Dr. Maziar Izadi, and Mr. James Henrickson for their help and suggestions during the development of this software. 

# References

<!--Please contact the authors at ramaniitrgoyal92@tamu.edu, muhaochen@tamu.edu, mmajji@tamu.edu or bobskelton@tamu.edu for a copy of the submitted paper. -->
<!-- # - name: Assistant Professor, Director of LASR Laboratory, Texas A&M University #   index: 2 -->
<!-- # - name: TEES Eminent Professor, Member National Academy of Engineering, Texas A&M University #   index: 3 -->
