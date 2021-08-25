# BPTA_polymerization_sequestration
Matlab scripts for fitting a supramolecular polymerization in competition with a sequestrating event, as described in our manuscript 'Temperature-dependent modulation by biaryl-based monomers of the chain length and morphology of biphenyl-based supramolecular polymers'

# Usage
Either download the scripts or pull them using 'git clone git://github.com/mathijs-m/BPTA_polymerization_sequestration.git
Run Fit_polymerization_sequestration.m. The script will ask for the dataset that you want to fit. 

# Input format
The script requires a plain text file consisting of three tab or comma separated columns as input. The first column should contain the mol% of the additive in the (fixed) total concentration, the second column should contain the temperature and the third column the CD intensity. The script will ask for the total concentration of the system.

# Output
After running successfully, the script will output the fitted CD trace, a speciation plot and txt files of all data that can be used in other programs, such as Origin, Prism or Excel.

# Questions?
If you have problems running the script or you want to change something but you're not sure about what, where and how, don't hesitate to reach out to me. Either here, or by sending an e-mail.
