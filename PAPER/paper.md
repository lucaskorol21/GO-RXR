---
title: 'GO-RXR: Global Optimization of Resonant X-ray Reflectometry'
tags:
  - Python
  - resonant x-ray reflectometry
  - material science
  - thin film heterostructure
authors:
  - name: Lucas Korol
    orchid: 0009-0007-1329-6065
    affiliation: "1"
  - name: Robert J. Green
    orchid: 0000-0003-1849-1576
    affiliation: "1, 3"
  - name: Jesus P. Curbelo
    orchid: 0000-0003-2417-071X
    affiliation: "2"
  - name: Raymond J. Spiteri
    orchid: 0000-0002-3513-6237
    affiliation: "2"
affiliations:
  - name: Department of Physics & Engineering Physics, University of Saskatchewan, Saskatoon, Canada S7N 5E2
    index: 1
  - name: Department of Computer Science, University of Saskatchewan, Saskatoon, Canada S7N 5E2
    index: 2
  - name: Stewart Blusson Quantum Matter Institute, University of British Columbia, Vancouver, Canada V6T 1Z1 
    index: 3
date: 2 April 2024
bibliography: paper.bib
---

# Summary

Resonant x-ray reflectometry (RXR) is a cutting-edge synchrotron technique used to characterize the depth-dependent structure of quantum materials [@keimer_moore_NPh_2017; @green-etal_SRN_2020]. However, the main challenge impeding the success of RXR data analysis lies in its extreme complexity, driven by complicated model construction and the fitting of numerous independent variables. This complexity results in prolonged analysis periods that demand significant engagement from researchers. In response to these challenges, the Global Optimization of Resonant X-ray Reflectometry (GO-RXR) software emerged from rigorous development efforts as a main contribution from the work by [@korol_MSc_2023]. GO-RXR streamlines data analysis, enhances visualization, and reduces the expertise required, offering researchers a more efficient means to analyze RXR data. 

This paper presents an overview of GO-RXR, highlighting its functionality, example use-cases, and impact in materials science research. Through its comprehensive approach and user-friendly design, GO-RXR offers researchers an efficient tool for analyzing RXR data, facilitating breakthroughs in understanding complex material systems. Additionally, publications and ongoing research utilizing GO-RXR underscore its versatility and impact in advancing scientific exploration.

# Statement of Need

RXR offers unique insights into the depth-dependent crystal, electronic, and magnetic structures of quantum materials, enabling the investigation of nanoscale characteristics of new candidate materials with a precision unmatched by any other current experimental technique **[REF - new?]**. Despite its potential, RXR remains significantly underutilized, with far fewer publications compared to techniques such as x-ray absorption spectroscopy. The main challenge hindering the widespread adoption of RXR lies in the extreme difficulty of data analysis, which requires both large-scale computational quantum mechanics simulations and the fitting of reflectivity models with numerous independent variables, such as layer thickness, interfacial roughness, and complex refractive indices that vary with energy. Each of these parameters must be finely adjusted to match the experimental data, making the process intricate and labor-intensive. For example, the oscillations in the Kiessig fringes during theta/two-theta reflectivity scans are directly related to the thickness of the film, while the decay of these fringes indicates the roughness of various interfaces within the material. This makes the analysis process for each sample highly demanding and time-consuming. Consequently, experimental advancements have far outpaced the progress in analytical methods, leaving a substantial amount of collected data unexplored **[REF???]**.

The complexity of RXR analysis extends beyond computational demands, requiring a deep understanding of material properties, such as ferromagnetic order or electronic reconstruction in polar-mismatched heterostructures, and the physics of light-matter interactions, such as energy-dependent absorption edge shifts and dichroism effects. This expertise is pivotal because it provides intuition about parameter adjustments and guides the direction of data analysis to achieve desired outcomes. In addressing this challenge, GO-RXR integrates global optimization algorithms, thereby lowering the expertise threshold necessary for effective data analysis. Through extensive development, diverse global optimization algorithms and unique objective functions were thoroughly explored. The software's capability to capture features in experimental data without exhaustive parameter adjustments significantly reduces the expertise required. It effectively models the oscillations in Kiessig fringes, which are directly linked to film thickness, while preserving the overall shape of the theta/two-theta reflectivity scans using a total variation penalty term. This ensures that both fine details and broader trends in the data are accurately represented. Additionally, GO-RXR offers enhanced flexibility in modeling strain at interfaces by allowing different form factors to be applied to distinct layers within the same element, enabling a more precise depiction of complex interfacial phenomena. GO-RXR serves as a valuable scientific tool for material scientists, offering advanced capabilities to streamline data analysis, reduce the expertise barrier, and ultimately facilitate breakthrough discoveries in the field of materials science.

One notable example is the analysis of LaMnO<sub>3</sub>/SrTiO<sub>3</sub> thin-film heterostructures addressed by [@korol_MSc_2023]. Data collected at the resonant elastic and inelastic x-ray scattering beamline (REIXS) at the Canadian Light Source (CLS) in 2017 had remained unanalyzed despite multiple attempts using available tools in 2021. The data analysis was particularly challenging due to the complex element-specific interactions and the presence of a magnetic dead layer at the surface. The analysis was particularly challenging due to the complex element-specific interactions and the need to model the magnetic dead layer at the interface, which existing tools could not achieve effectively. However, in 2023, the use of GO-RXR enabled a successful analysis, demonstrating the software's efficacy in overcoming longstanding barriers and highlighting its impact on advancing RXR studies.


# Comparison

To the best of our knowledge, no existing software tool comprehensively addresses the general problem of RXR data analysis across a wide range of materials and conditions. Most of the currently available tools are designed for very specific tasks and cannot be extensively applied to the diverse challenges encountered in RXR studies. For instance, tools like GenX [@bjorck2007_genx] and Motofit [@nelson2006_motofit], while effective for neutron and non-resonant x-ray reflectometry, lack the capability to model key aspects of RXR data, such as magnetic scattering contributions and energy-dependent changes in refractive indices, which are critical for interpreting element-specific interactions in complex material. Similarly, more modern tools like Refl1D [@refl1d_github] and BornAgain [@pospelov2020_bornagain] offer powerful features for modeling and fitting specific types of reflectometry data but still lack the flexibility and comprehensive approach needed for the broad range of challenges in RXR analysis.

GO-RXR was developed to fill this gap by providing a comprehensive tool for RXR data analysis. It incorporates advanced global optimization algorithms that efficiently handle the high dimensionality and complexity of RXR datasets, avoiding common issues like local minima. GO-RXR’s ability to customize objective functions and integrate fitting techniques, such as total variation penalty methods, enables accurate modeling of complex electronic and magnetic structures in thin films. This makes it uniquely suited for a broad range of material systems and experimental conditions **[examples here]**.

Additionally, GO-RXR features a user-friendly graphical interface with intuitive plotting options, enhancing data interpretation and visualization. This combination of advanced analytical capabilities and accessibility bridges the gap between experimental data collection and theoretical modeling, offering researchers a powerful tool for studying complex material systems in RXR.

# Functionality

GO-RXR is available for installation through its GitHub repository, providing users with easy access to its functionality. The accompanying documentation offers thorough guidance, including detailed installation instructions, a user guide, and tutorials featuring two example cases solved step-by-step. This comprehensive support ensures that users can quickly grasp the software's features and efficiently apply them to their data analysis tasks.

## Basic concepts

The primary goal of GO-RXR is to simplify and optimize the data analysis process for resonant x-ray reflectometry. Initially conceived as a command-line tool, it evolved into a GUI-based software with enhanced visualization capabilities. Developed using Python with PyQt5 for the interface, this software integrates the Pythonreflectivity [@pythonreflectivity] open-source package to carry out reflectivity calculations efficiently. GO-RXR has been tested extensively on Ubuntu 22.04 with Python 3.10, ensuring reliable performance on this operating system. After running the GUI_GO.py file, the start screen in \autoref{fig:start-screen} will be displayed.

![Start screen of GO-RXR. The interface is divided into three main sections: 1) Toolbar, 2) Workspace Navigator, and 3) Workspace Area. \label{fig:start-screen}](../FIGURES/go-rxr-start.png)

\autoref{fig:flowchart} displays the main steps GO-RXR uses to go from input to output. The inputs of the software are the experimental data and are defined by the reflected intensity, incident grazing angle, photon energy, and polarization. The first step in the data analysis is the data selection step. In this step, the user selects the experimental datasets to include into the data analysis. The next step is the parameter selection. In this step, the user selects the model parameters to vary. The final step is the data fitting. In this step, the user selects the global optimization algorithm and its parameters and fits the data from the selected experimental and simulated datasets. The output of GO-RXR is the depth-dependent density profile, as defined by the model parameters.

![Flowchart of GO-RXR. The flowchart illustrates the path the software uses to take the experimental data and convert them into a depth-dependent density profile.\label{fig:flowchart}](../FIGURES/go-rxr-flowchart.png)

GO-RXR incorporates three global optimization algorithms specifically designed to tackle the inherent complexities of RXR data analysis. These are Differential Evolution (DE) [@storn1997_de], Simplicial Homology Global Optimization (SHGO) [@endres2018_shgo], and Dual Annealing [@xiang1997_dual_annealing], each chosen for their effectiveness in handling the high-dimensional optimization challenges typical of RXR. A significant enhancement provided by GO-RXR is the integration of boundary and weight functions, which enable the selective emphasis of specific data regions, thereby improving the precision and relevance of the optimization results [@korol_MSc_2023]. Additionally, the objective function includes a total variation penalty, ensuring that the optimization not only fits the data but also preserves its physical trends, resulting in more accurate and reliable outcomes [@korol_MSc_2023]. These global optimization algorithms, combined with  local optimization nonlinear least-squares algorithms for post-processing refinement, allow GO-RXR to effectively model the structural and electronic properties of thin-film materials. The algorithms are implemented within GO-RXR using the SciPy library [@virtanen2020_scipy]. \autoref{fig:optimization} displays the Optimization workspace, which includes sections for parameter boundaries, objective function selection, and algorithm parameters, allowing users to fine-tune the optimization process to achieve the most accurate and physically meaningful results.

![Optimization workspace. \label{fig:optimization}](../FIGURES/go-rxr-optimization.png)


## Example use-cases

The GO-RXR software has been instrumental in analyzing RXR data for diverse applications. For instance, it has facilitated the in-depth study of the electrochemical water splitting catalyst La<sub>0.7</sub> Sr<sub>0.3</sub> MnO3/SrTiO<sub>3</sub> (LSMO/STO), revealing insights into the material's structural, electronic, and magnetic depth profiles. Moreover, GO-RXR's utilization in investigating the relationship between film thickness and the presence of ferromagnetism in the LaMnO<sub>3</sub>/SrTiO<sub>3</sub> heterostructure has shed new light on the mechanisms underlying magnetic phase transitions in ultra-thin films.

The adaptability of GO-RXR extends to analyzing a wide array of quantum materials, including semiconductors, superconductors, and magnetic materials, allowing researchers to explore their depth-dependent structures and properties. More generally, GO-RXR could be applied to materials science research, particularly in comprehending the structural and electronic attributes for designing cutting-edge electronic devices and functional materials, especially thin films. Furthermore, its ability to streamline data analysis processes and enhance visualization can benefit researchers across different disciplines, including physics, chemistry, and materials engineering, by providing them with a robust tool for studying complex material systems with greater efficiency and accuracy.


# Publications and ongoing research

The GO-RXR software package was utilized for analyzing the RXR data in the paper titled "The effect of intrinsic magnetic order on electrochemical water splitting" published in the journal Applied Physics Reviews by [@vanderMinne_etal_2023]. Additionally, GO-RXR played a valuable role in a complementary experimental study investigating morphological and chemical disorder in epitaxial La<sub>0.67</sub>Sr<sub>0.33</sub>MnO<sub>3</sub>. This study has been accepted for publication in the journal ACS Applied Materials & Interfaces, with its preprint format available online [@verhage_etal_2023]. These instances underscore the versatility and impact of GO-RXR in advancing scientific exploration, with ongoing research endeavors continuing to harness its capabilities for further discoveries.

# Acknowledgments

Robert J. Green and Raymond J. Spiteri acknowledge the support
from the Natural Sciences and Engineering Research Council of
Canada (NSERC) Discovery Grant program. Lucas Korol
acknowledges the support from the NSERC CREATE to INSPIRE
program. 

# References
