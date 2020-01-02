from md_analyzer import *

######## USER PROGRAM ##########
FN=2000

plot_structure=False
calc_velocities=True
calc_bond_stats=False

#analyze_md_output("XDATCAR100",100, plot_structure, calc_velocities, calc_bond_stats)
analyze_md_output("XDATCAR300",300, plot_structure, calc_velocities, calc_bond_stats)

#show all plots (if any)
plt.show()
