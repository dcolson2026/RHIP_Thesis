import ROOT
import time
import numpy as np
from particle import Particle
from ROOT_analysis_functions import caldeltaphi, particle_lineage_dfs, calpT, calphi, getParticleName, delta_phi_correlation, caleta, following_charm_dfs

PI = np.pi

# Open the ROOT file
root_file = ROOT.TFile("/home/daniel/LibraFiles/CleanThesis/PythiaData/pythia_1mil_events.root", "READ")
#root_file = ROOT.TFile("pythia_100_events_pT_hard.root", "READ")
tree = root_file.Get("events") # fetches the tree

# Create variables to hold data
pdg = ROOT.std.vector('int')()
status = ROOT.std.vector('int')()
px = ROOT.std.vector('float')()
py = ROOT.std.vector('float')()
pz = ROOT.std.vector('float')()
mother_1 = ROOT.std.vector('int')()
mother_2 = ROOT.std.vector('int')()
daughter_1 = ROOT.std.vector('int')()
daughter_2 = ROOT.std.vector('int')()


# Set branch addresses
# when get entry is called, this fills the variables from above with
# values from "branch" e.g. "pid"
tree.SetBranchAddress("pid", pdg)
tree.SetBranchAddress("status", status)
tree.SetBranchAddress("px", px)
tree.SetBranchAddress("py", py)
tree.SetBranchAddress("pz", pz)
tree.SetBranchAddress("mother_1", mother_1)
tree.SetBranchAddress("mother_2", mother_2)
tree.SetBranchAddress("daughter_1", daughter_1)
tree.SetBranchAddress("daughter_2", daughter_2)

# Define PDG IDs of interest
C_PDG = {4}  # c quark
CBAR_PDG = {-4} # cbar quark

# Create a ROOT histogram
h1_all_ccbar_delta_phi = ROOT.TH1F("h1_all_ccbar_delta_phi", "#Delta#phi Correlation: c - cbar quarks; #Delta#phi (radians); Counts", 30, -PI/2, 3*PI/2)

# Loop over all events
m_events = tree.GetEntries()
print(m_events, "events")
for i in range(50000): # m_entries, edit to look at single event
    ccbar_pair_count = 0
    ccbar_hadron_daughters = 0
    if i % 10000 == 0:
        print(str(i) + " events have completed")
    tree.GetEntry(i)
    
    # Collect phi values of c quarks and cbar quarks in the event
    c_phi_list = []
    cbar_phi_list = []
    all_pion_pT_list = []
    
    # Collect daughter indices
    c_quark_daughter_indices_list = []
    cbar_quark_daughter_indices_list = []

    c_quark_pre_hadronization_index = -1
    cbar_quark_pre_hadronization_index = -1
    c_hadron_daughters = []
    cbar_hadron_daughters = []

    # 09.10.25
    c_descendant_hadron_indices = []
    cbar_descendant_hadron_indices = []
    c_descendant_parton_indices = []
    cbar_descendant_parton_indices = []
    c_descendant_init_showers_indices = []
    cbar_descendant_init_showers_indices = []
    c_descendant_final_showers_indices = []
    cbar_descendant_final_showers_indices = []
    
    n_particles = len(pdg)
    for j in range(n_particles): # particle loop, total number of particles in event "i"
        # Particle loop variables
        particle_px = px[j]
        particle_py = py[j]
        particle_pz = pz[j]
        particle_pT = calpT(particle_px, particle_py)
        particle_pdg = pdg[j]
        particle_status = status[j]
        particle_phi = calphi(y_momentum=particle_py, x_momentum=particle_px)
        particle_eta = caleta(x_momentum=particle_px, y_momentum=particle_py, z_momentum=particle_pz)
        next_particle_pdg = pdg[(j+1)%n_particles]
        next_particle_phi = calphi(y_momentum=py[(j+1)%n_particles], x_momentum=px[(j+1)%n_particles])
        
        mother_1_index = mother_1[j]
        mother_1_pdg = pdg[mother_1_index]
        mother_2_index = mother_2[j]
        mother_2_pdg = pdg[mother_2_index]
        next_particle_mother_1_index = mother_1[(j+1)%n_particles] # mod

        daughter_1_index = daughter_1[j]
        daughter_1_pdg = pdg[daughter_1_index]
        daughter_2_index = daughter_2[j]
        daughter_2_pdg = pdg[daughter_2_index]
        daughter_1_status = status[daughter_1_index]
        daughter_2_status = status[daughter_2_index]
        
        
        # if we have a c or cbar, check the status of its daughter. if its a primary hadron from hadronization, fill
        if pdg[j] in C_PDG:
            if 81 <= np.abs(daughter_1_status) <= 89:
                ccbar_hadron_daughters+=1
                c_hadron_daughters.append(daughter_1_index)
                c_quark_pre_hadronization_index = j
            if 81 <= np.abs(daughter_2_status) <= 89:
                ccbar_hadron_daughters+=1
                c_hadron_daughters.append(daughter_2_index)
        elif pdg[j] in CBAR_PDG:
            if 81 <= np.abs(daughter_1_status) <= 89:
                ccbar_hadron_daughters+=1
                cbar_hadron_daughters.append(daughter_1_index)
                cbar_quark_pre_hadronization_index = j
            if 81 <= np.abs(daughter_2_status) <= 89:
                ccbar_hadron_daughters+=1
                cbar_hadron_daughters.append(daughter_2_index)
        
        # finds ccbar pair and confirms that they came from the same mother particle (so they are paired)
        if particle_pdg == 4 and next_particle_pdg == -4 and mother_1_index == next_particle_mother_1_index:
            ccbar_pair_count += 1
            h1_all_ccbar_delta_phi.Fill(caldeltaphi(particle_phi, next_particle_phi))
            following_charm_dfs(j, pdg, daughter_1, daughter_2, c_quark_daughter_indices_list)
            following_charm_dfs(j+1, pdg, daughter_1, daughter_2, cbar_quark_daughter_indices_list)
        
        # finds ccbar this time if ordering is reversed (doubtful, but just in case)
        # very same to assume that there is only one ccbar pair per event (if not zero) due to low cross section
        if particle_pdg == -4 and next_particle_pdg == 4 and mother_1_index == next_particle_mother_1_index:
            ccbar_pair_count += 1
            h1_all_ccbar_delta_phi.Fill(caldeltaphi(particle_phi, next_particle_phi))
            following_charm_dfs(j, pdg, daughter_1, daughter_2, cbar_quark_daughter_indices_list)
            following_charm_dfs(j+1, pdg, daughter_1, daughter_2, c_quark_daughter_indices_list)
        ### you are now exiting the particle loop. welcome to event loop, population: you
    
    if len(c_quark_daughter_indices_list) != 0 or len(cbar_quark_daughter_indices_list) != 0:
        print(len(c_quark_daughter_indices_list))
        print(len(cbar_quark_daughter_indices_list))
    


# Style for ROOT histogram
# h1_delta_phi.SetLineColor(ROOT.kBlue)
h1_all_ccbar_delta_phi.SetLineWidth(2)

# Save the histogram as a ROOT file
output_file = ROOT.TFile("/home/daniel/LibraFiles/CleanThesis/RootOutputs/09_30_2025_following_charm.root", "RECREATE")

# Write out all histograms
h1_all_ccbar_delta_phi.Write()


output_file.Close()
