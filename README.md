# Urate_Project2024

The contents of this repository relate to the manuscript 'Evaluating the causal relationships between urate, blood pressure, and kidney function in the general population: a two-sample Mendelian Randomization study' by Haotian Tang, Venexia M Walker and Tom R Gaunt.

# Link to the manuscript

Will copy it here once we submit to arxiv

# Abstract

Associations between blood urate levels, blood pressure (BP), and kidney function have previously been reported in observational studies. However, causal inference between these three traits is challenging due to potentially bidirectional relationships. We applied bidirectional univariable Mendelian randomization (UVMR) to assess the causal relationships between urate levels, BP, and kidney function, proxied by estimated glomerular filtration %rate (eGFR), by using genetic associations from both UK Biobank and CKDGen. We performed multivariable MR (MVMR) to assess the independent effects of urate and BP on eGFR. Effect estimates are presented as standard deviation (SD) change in outcome per SD increase in exposure [95% confidence interval]. UVMR analysis suggested a bidirectional causal effect between urate and eGFR (urate on log(eGFR): beta=-0.10 [-0.22 to 0.02]; log(eGFR) on urate: beta=-0.11 [-0.17 to -0.04]). There was also strong evidence of bidirectional causal effects between urate and SBP (urate on SBP: beta=0.08 [0.04 to 0.11]; SBP on urate: beta=0.13 [0.08 to 0.18]). Similar bidirectional causal effects were identified between urate and DBP (urate on DBP: beta=0.09 [0.05 to 0.14]; DBP on urate: beta=0.13 [0.08 to 0.18]). However, there was weak evidence of a causal effect of eGFR on BP or BP on eGFR. MVMR results suggested the causal effect of urate on eGFR was independent of BP. Our results provide evidence for bidirectional causal effects between urate and both eGFR and BP, suggesting urate control as a potential intervention to reduce BP and decline in kidney function in the general population. We found little evidence of a causal relationship between BP and eGFR.

# Questions
Please send any questions to Haotian Tang (haotian.tang@bristol.ac.uk).

# Script documentation

Each script within the script repo includes a documentation at the beginning that clarifies the purpose of the script. The numbers in each script name indicate the order in which the code should be run.

Version of R: 

4.2.2

Main packages that need to be installed: 

[TwoSampleMR 0.5.11](https://github.com/MRCIEU/TwoSampleMR)

[MVMR 0.4](https://github.com/WSpiller/MVMR)

[forestplot 3.1.3](https://github.com/gforge/forestplot)
