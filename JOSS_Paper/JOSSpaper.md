---
title: 'TEAM: Tensegrity Engineering Analysis Master'

tags:
  - tensegrity
  - Lagrangian
  - statics
  - dynamics
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
    orcid: 
    affiliation: 1
    
affiliations:
 - name: Department of Aerospace Engineering, Texas A&M University, College Station, Texas, USA
   index: 1

 # - name: Associate Professor, Texas A&M University
 #   index: 2

date: 13 September 2018
bibliography: JOSSpaper.bib
---

# Summary

A tensegrity system is an arrangement of axially-loaded elements (no element bends, even though the overall structure bends), that we loosely characterize as a network of bars and cables. The bars take compressive axial loads and the cables handle tensile loads. Since failure due to axial stresses happens at higher loads than at bending, a tensegrity structure has a higher strength-to-weight ratio. The famous architect Buckminster Fuller, in the 60's, coined the term tensegrity, combining the words tensile and integrity. Since then, tensegrity principles have found applications in diverse domains like architecture [@ingber1998architecture], biological modeling [@scarr2014biotensegrity] as well as civil engineering design [@tenseBook]. Tensegrity structures, through use of pre-stresses in the bars and cables, can also achieve controlled stiffness in the structure, which makes it attractive in applications such as soft-robotics and prosthetics. In essence, tensegrity principles can be applied in the design of any structure where mass is premium, a high strength-to-weight ratio is critical, and structural stiffness needs to be tailored in both space and time. These include several applications from various engineering sectors such as aerospace (morphing airframes), energy (wind turbine blades, off-shore structures) as well as biomedical engineering (stents, minimally invasive surgical tools) and many more. Clearly, a framework is required that can efficiently model the dynamics of tensegrity structures directly from the topology of bars and cables.

The dynamics of tensegrity systems is governed by multi-body dynamics, given by a set of ordinary differential equations. We have developed a Lagrangian formulation for deriving these differential equations directly from the given topology of members (bars and strings), and their mass and geometric properties. Three key features of classical tensegrity systems are a) actuations only occur via cables, b) bar-to-bar connections are pin joints, and c) the bars do not spin about their longitudinal axis. These properties are exploited to simplify the equations of motion. However, the Lagrangian framework presented here is general enough to allow modeling of general multi-body systems with actuated joints.

`STEDY` is a MATLAB package for conducting numerically accurate tensegrity dynamics. The equations of motion are developed in non-minimum Cartesian coordinates. Although the choice of non-minimum coordinates is fairly controversial in the robotics community on account of a reduction in the system's configuration space in the presence of joints, adopting Cartesian coordinates precludes any singularities developed in the mass matrix [@tenseBook] and facilitates derivations of elegant differential-algebraic equations governing the multi-body dynamics. However, usage of non-minimum coordinates is likely to result in the solution drifting away from the constraint space because of integration errors [@ConstraintViolation]. The software package includes an implementation of the direct correction scheme that minimizes numerically induced constraint violations and energy conservation errors. We have demonstrated superiority of the method in terms of accuracy over the commercially available dynamics simulator, Simscape Multibody [@matlab].

`STEDY` was designed to be used by researchers looking to simulate tensegrity dynamics without having to derive the differential equations for any structure in particular. The inputs that are required from the user are simply the parameters that uniquely define the tensegrity structure (initial nodal configuration, connectivity), the material and geometric properties, and simulation environment properties (inertially fixed nodes, external forces, duration, ODE solver options). Currently, `STEDY` is being used by TAMU researchers to develop [space habitats] and [wind turbine blades] using tensegrity principles.

The software documentation and the manuscript [@LagCMAME] comprehensively cover the functionality and the theory behind the Lagrangian formulation and the novel constraint correction method implemented in the package. We shall seek to incorporate techniques introducing control of dynamic tensegrity systems in future versions of the software.

# Acknowledgements

This work was supported by NSF IUSE/PFE: RED: REvolutionizing Diversity Of Engineering (REDO-E)Award Number:1730693; and NASA NIAC Phase II grant, on Tensegrity Approaches to In-Space Construction of a 1g Growable Habitat.

# References
<!-- 1. Shao-Chen Hsu, Vaishnav Tadiparthi, and Raktim Bhattacharya, "A Lagrangian Formulation for Constrained Multibody Dynamics in Tensegrity Systems", Manuscript submitted for publication.

    Please contact the authors at addyhsu@tamu.edu, vaishnavtv@tamu.edu, or raktim@tamu.edu for a copy of the submitted paper. -->

[space habitats]: https://www.nasa.gov/feature/tensegrity-approaches-to-in-space-construction-of-a-1g-growable-habitat

[wind turbine blades]:https://www.nsf.gov/awardsearch/showAward?AWD_ID=1762825&HistoricalAwards=false
