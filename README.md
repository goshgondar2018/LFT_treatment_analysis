# Code for "The impact of repeated rapid test strategies on the effectiveness of at-home antiviral treatments for SARS-CoV-2" by Tigist F. Menkir and Christl A. Donnelly

# estimating_drug_efficacy_curves.R #

Code for estimating risk ratios of hospitalization using the Pfizer EPIC-HR latest summary findings 

# generating_EVs_main_may2.R #

Code for estimating weighted risk ratios (RRs) of hospitalization under our main treatment efficacy scenario, further assuming no delays from testing to treatment and full treatment coverage, for each of the testing strategies. Intermediate output includes weighted RRs across simulations; the median, LQ, and UQ of the simulated weighted RRs; estimated proportions given treatment across simulations; the median, LQ, and UQ of the simulated proportions given treatment; estimated proportions benefiting from treatment across simulations; and the median, LQ, and UQ of the simulated proportions benefiting from treatment. The final synthesized output combines results across testing strategies. 

# generating_EVs_sensitivity_analysis_may2.R #

Code for estimating weighted risk ratios (RRs) of hospitalization under two additional scenarios for treatment efficacy ('fast decline to zero' and 'efficacy preserved'). Intermediate output includes, by scenario, weighted RRs across simulations; the median, LQ, and UQ of the simulated weighted RRs; estimated proportions given treatment across simulations; the median, LQ, and UQ of the simulated proportions given treatment; estimated proportions benefiting from treatment across simulations; and the median, LQ, and UQ of the simulated proportions benefiting from treatment. The final synthesized output combines results across testing strategies, again by scenario.

# generating_EVs_sensitivity_analyses_access_apr8.R #

Code for estimating weighted risk ratios (RRs) of hospitalization under a range of test-to-treatment delays and treatment coverage proportions. Output includes weighted RRs for every combination of assumed treatment coverage and test-to-treatment delays, a contour plot of weighted RRs, treatment coverage, and test-to-treatment delays, and a 3-d point plot of weighted RRs, treatment coverage, and test-to-treatment delays.

# test_coverage_scenarios_may2.R #

Code for estimating proportion given treatment and benefit as a function of proportion tested, across strategies.

# generating_EVs_figures_may10.R #

Code to produce all remaining manuscript and supplement figures 
