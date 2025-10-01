import ROOT
import time
import numpy as np
from particle import Particle
from ROOT_analysis_functions import caldeltaphi, particle_lineage_dfs, calpT, calphi, getParticleName, delta_phi_correlation, caleta

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
PION_PDGS = {111, 211, -211} # pions
D0_PDGS = {421, 411, -411}

# Create a ROOT histogram
h1_all_ccbar_delta_phi = ROOT.TH1F("h1_all_ccbar_delta_phi", "#Delta#phi Correlation: c - cbar quarks; #Delta#phi (radians); Counts", 30, -PI/2, 3*PI/2)
#h2_delta_pT_cquark_chadrons = ROOT.TH1F("h2_delta_pT_cquark_chadrons", "#DeltapT Correlation: cquark - chadrons; #DeltapT (GeV/c); Counts", 30, -PI/2, 3*PI/2)
h3_all_pion_pT = ROOT.TH1F("h3_all_pion_pT", "pT of all #pi^{+}, #pi^{-}, #pi^{0},", 30, 0, 20)
h4_all_ccbar_daughters = ROOT.TH1I("h4_all_ccbar_daughters", "Immediate Daughters of c and cbar", 5000, -2500, 2500)
h5_all_D_pT = ROOT.TH1F("h5_all_D_pT", "pT of all #D^{+}, #D^{-}, #D^{0},", 30, 0, 20)
h6_final_state_hadrons_of_ccbar = ROOT.TH1I("h6_final_state_hadrons_of_ccbar", "Final State Hadrons of c and cbar", 5000, -2500, 2500)
h7_final_state_hadrons_of_ccbar_pT = ROOT.TH1F("h7_final_state_hadrons_of_ccbar_pT", "pT of Final State Hadrons of c and cbar", 30, 0, 20)
h8_all_pion_eta = ROOT.TH1F("h8_all_pion_eta", "#eta of all particles", 1024, -50, 50)

h9_all_primary_hadron_ccbar_daughters = ROOT.TH1I("h9_all_primary_hadron_ccbar_daughters", "Daughters of c and cbar- primary hadrons produced by hadronization process", 5000, -2500, 2500) # status code 81-89
h10_all_primary_hadron_ccbar_daughter_dphi = ROOT.TH1F("h10_all_primary_hadron_ccbar_daughter_dphi", "Dphi of ccbar daughters- primary hadrons produced by hadronization process", 30, -PI/2, 3*PI/2) # status code 81-89
h11_ccbar_pair_before_hadronization_dphi = ROOT.TH1F("h11_ccbar_pair_before_hadronization_dphi", "Dphi of ccbar daughters right before hadronization process", 30, -PI/2, 3*PI/2) # status code 81-89
h12_ccbar_hadron_descendants_dphi = ROOT.TH1F("h12_ccbar_hadron_descendants_dphi", "Dphi of ccbar descendant hadrons", 30, -PI/2, 3*PI/2) # status code 81-89
h13_ccbar_parton_descendants_dphi = ROOT.TH1F("h13_ccbar_parton_descendants_dphi", "Dphi of ccbar descendant partons", 30, -PI/2, 3*PI/2) # status code 71-79
h15_ccbar_final_showers_descendants_dphi = ROOT.TH1F("h15_ccbar_final_showers_descendants_dphi", "Dphi of ccbar descendant final shower particles", 30, -PI/2, 3*PI/2) # status code 51-59
h14_ccbar_init_showers_descendants_dphi = ROOT.TH1F("h14_ccbar_init_showers_descendants_dphi", "Dphi of ccbar descendant initial shower particles", 30, -PI/2, 3*PI/2) # status code 41-49

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
        
        # Collecting for correlation lists
        if pdg[j] in C_PDG:
            c_phi_list.append(particle_phi)
        elif pdg[j] in CBAR_PDG:
            cbar_phi_list.append(particle_phi)
        
        # pT histograms
        if particle_pdg in PION_PDGS:
            h3_all_pion_pT.Fill(particle_pT)
            h8_all_pion_eta.Fill(particle_eta)
        elif (particle_pdg in C_PDG) or (particle_pdg in CBAR_PDG):
            h4_all_ccbar_daughters.Fill(daughter_1_pdg)
            h4_all_ccbar_daughters.Fill(daughter_2_pdg)
        elif particle_pdg in D0_PDGS:
            h5_all_D_pT.Fill(particle_pT)
        
        # if we have a c or cbar, check the status of its daughter. if its a primary hadron from hadronization, fill
        if pdg[j] in C_PDG:
            if 81 <= np.abs(daughter_1_status) <= 89:
                h9_all_primary_hadron_ccbar_daughters.Fill(daughter_1_pdg)
                ccbar_hadron_daughters+=1
                c_hadron_daughters.append(daughter_1_index)
                c_quark_pre_hadronization_index = j
            if 81 <= np.abs(daughter_2_status) <= 89:
                h9_all_primary_hadron_ccbar_daughters.Fill(daughter_2_pdg)
                ccbar_hadron_daughters+=1
                c_hadron_daughters.append(daughter_2_index)
        elif pdg[j] in CBAR_PDG:
            if 81 <= np.abs(daughter_1_status) <= 89:
                h9_all_primary_hadron_ccbar_daughters.Fill(daughter_1_pdg)
                ccbar_hadron_daughters+=1
                cbar_hadron_daughters.append(daughter_1_index)
                cbar_quark_pre_hadronization_index = j
            if 81 <= np.abs(daughter_2_status) <= 89:
                h9_all_primary_hadron_ccbar_daughters.Fill(daughter_2_pdg)
                ccbar_hadron_daughters+=1
                cbar_hadron_daughters.append(daughter_2_index)
        
        # finds ccbar pair and confirms that they came from the same mother particle (so they are paired)
        if particle_pdg == 4 and next_particle_pdg == -4 and mother_1_index == next_particle_mother_1_index:
            ccbar_pair_count += 1
            h1_all_ccbar_delta_phi.Fill(caldeltaphi(particle_phi, next_particle_phi))
            particle_lineage_dfs(j, daughter_1, daughter_2, c_quark_daughter_indices_list)
            particle_lineage_dfs(j+1, daughter_1, daughter_2, cbar_quark_daughter_indices_list)
        
        # finds ccbar this time if ordering is reversed (doubtful, but just in case)
        # very same to assume that there is only one ccbar pair per event (if not zero) due to low cross section
        if particle_pdg == -4 and next_particle_pdg == 4 and mother_1_index == next_particle_mother_1_index:
            ccbar_pair_count += 1
            h1_all_ccbar_delta_phi.Fill(caldeltaphi(particle_phi, next_particle_phi))
            particle_lineage_dfs(j, daughter_1, daughter_2, cbar_quark_daughter_indices_list)
            particle_lineage_dfs(j+1, daughter_1, daughter_2, c_quark_daughter_indices_list)
        ### you are now exiting the particle loop. welcome to event loop, population: you
    
    # these are the correlated hadrons that proceed after the ccbars
    if ccbar_hadron_daughters == 4 and ccbar_pair_count == 1:
        # print(ccbar_hadron_daughters, "matches", len(c_hadron_daughters)+len(cbar_hadron_daughters))
        for o in c_hadron_daughters:
            for p in cbar_hadron_daughters:
                c_phi_temp = calphi(y_momentum=py[o], x_momentum=px[o])
                cbar_phi_temp = calphi(y_momentum=py[p], x_momentum=px[p])
                dphi_temp = caldeltaphi(c_phi_temp, cbar_phi_temp)
                h10_all_primary_hadron_ccbar_daughter_dphi.Fill(dphi_temp)
        # dphi of c and cbar RIGHT before hadronization
        c_phi_pre_hadronization = calphi(y_momentum=py[c_quark_pre_hadronization_index], x_momentum=px[c_quark_pre_hadronization_index])
        cbar_phi_pre_hadronization = calphi(y_momentum=py[cbar_quark_pre_hadronization_index], x_momentum=px[cbar_quark_pre_hadronization_index])
        h11_ccbar_pair_before_hadronization_dphi.Fill(caldeltaphi(c_phi_pre_hadronization, cbar_phi_pre_hadronization))

    # new attempt. i have all of the descendants of the ccbar pair as c_quark_daughter_indices_list and cbar_quark_daughter_indices_list. now
    # try idea from 09.10.25
    # selections based on status codes:
    for i in c_quark_daughter_indices_list:
        c_daughter_pT = calpT(px[i], py[i])
        c_daughter_pdg = pdg[i]
        # if c_daughter_pT < 2: continue # applies a pT cut
        if 81 <= np.abs(status[i]) <= 89:
            c_descendant_hadron_indices.append(i)
        elif 71 <= np.abs(status[i]) <= 79:
            c_descendant_parton_indices.append(i)
        elif 51 <= np.abs(status[i]) <= 59:
            c_descendant_final_showers_indices.append(i)
        elif 41 <= np.abs(status[i]) <= 49:
            c_descendant_init_showers_indices.append(i)
    for i in cbar_quark_daughter_indices_list:
        cbar_daughter_pT = calpT(px[i], py[i])
        cbar_daughter_pdg = pdg[i]
        # if cbar_daughter_pT < 2: continue # applies a pT cut
        if 81 <= np.abs(status[i]) <= 89:
            cbar_descendant_hadron_indices.append(i)
        elif 71 <= np.abs(status[i]) <= 79:
            cbar_descendant_parton_indices.append(i)
        elif 51 <= np.abs(status[i]) <= 59:
            cbar_descendant_final_showers_indices.append(i)
        elif 41 <= np.abs(status[i]) <= 49:
            cbar_descendant_init_showers_indices.append(i)
    # hadron descendants dphi
    for i in c_descendant_hadron_indices:
        for j in cbar_descendant_hadron_indices:
            # print("made it!")
            c_descendant_hadron_phi = calphi(y_momentum=py[i], x_momentum=px[i])
            cbar_descendant_hadron_phi = calphi(y_momentum=py[j], x_momentum=px[j])
            dphi_hadron_descendants = caldeltaphi(c_descendant_hadron_phi, cbar_descendant_hadron_phi)
            h12_ccbar_hadron_descendants_dphi.Fill(dphi_hadron_descendants)
    # parton descendants dphi
    for i in c_descendant_parton_indices:
        for j in cbar_descendant_parton_indices:
            # print("made it! again")
            c_descendant_parton_phi = calphi(y_momentum=py[i], x_momentum=px[i])
            cbar_descendant_parton_phi = calphi(y_momentum=py[j], x_momentum=px[j])
            dphi_parton_descendants = caldeltaphi(c_descendant_parton_phi, cbar_descendant_parton_phi)
            h13_ccbar_parton_descendants_dphi.Fill(dphi_parton_descendants)
    # final shower descendants dphi
    for i in c_descendant_final_showers_indices:
        for j in cbar_descendant_final_showers_indices:
            # print("made it! finally")
            c_descendant_final_showers_phi = calphi(y_momentum=py[i], x_momentum=px[i])
            cbar_descendant_final_showers_phi = calphi(y_momentum=py[j], x_momentum=px[j])
            dphi_final_showers_descendants = caldeltaphi(c_descendant_final_showers_phi, cbar_descendant_final_showers_phi)
            h15_ccbar_final_showers_descendants_dphi.Fill(dphi_final_showers_descendants)
    # init shower descendants dphi
    for i in c_descendant_init_showers_indices:
        for j in cbar_descendant_init_showers_indices:
            # print("made it! again")
            c_descendant_init_showers_phi = calphi(y_momentum=py[i], x_momentum=px[i])
            cbar_descendant_init_showers_phi = calphi(y_momentum=py[j], x_momentum=px[j])
            dphi_init_showers_descendants = caldeltaphi(c_descendant_init_showers_phi, cbar_descendant_init_showers_phi)
            h14_ccbar_init_showers_descendants_dphi.Fill(dphi_init_showers_descendants)
    


    #delta_phi_correlation(c_phi_list, cbar_phi_list, h1_all_ccbar_delta_phi)
    
    for indx in c_quark_daughter_indices_list:
        if status[indx] > 0: # indicates final state hadron
            temp_pT = calpT(px[indx], py[indx])
            h6_final_state_hadrons_of_ccbar.Fill(pdg[indx])
            h7_final_state_hadrons_of_ccbar_pT.Fill(temp_pT)
        #if 

    for indx in cbar_quark_daughter_indices_list:
        if status[indx] > 0: # indicates final state hadron
            temp_pT = calpT(px[indx], py[indx])
            h6_final_state_hadrons_of_ccbar.Fill(pdg[indx])
            h7_final_state_hadrons_of_ccbar_pT.Fill(temp_pT)



# Style for ROOT histogram
# h1_delta_phi.SetLineColor(ROOT.kBlue)
h1_all_ccbar_delta_phi.SetLineWidth(2)
h3_all_pion_pT.SetLineWidth(2)
h4_all_ccbar_daughters.SetLineWidth(2)
h5_all_D_pT.SetLineWidth(2)
h6_final_state_hadrons_of_ccbar.SetLineWidth(2)
h7_final_state_hadrons_of_ccbar_pT.SetLineWidth(2)
h8_all_pion_eta.SetLineWidth(2)
h9_all_primary_hadron_ccbar_daughters.SetLineWidth(2)
h10_all_primary_hadron_ccbar_daughter_dphi.SetLineWidth(2)
h11_ccbar_pair_before_hadronization_dphi.SetLineWidth(2)
h12_ccbar_hadron_descendants_dphi.SetLineWidth(2)
h13_ccbar_parton_descendants_dphi.SetLineWidth(2)
h14_ccbar_init_showers_descendants_dphi.SetLineWidth(2)
h15_ccbar_final_showers_descendants_dphi.SetLineWidth(2)

# Save the histogram as a ROOT file
output_file = ROOT.TFile("/home/daniel/LibraFiles/CleanThesis/RootOutputs/full_09_26_2025.root", "RECREATE")

# Write out all histograms
h1_all_ccbar_delta_phi.Write()
h3_all_pion_pT.Write()
h4_all_ccbar_daughters.Write()
h5_all_D_pT.Write()
h6_final_state_hadrons_of_ccbar.Write()
h7_final_state_hadrons_of_ccbar_pT.Write()
h8_all_pion_eta.Write()
h9_all_primary_hadron_ccbar_daughters.Write()
h10_all_primary_hadron_ccbar_daughter_dphi.Write()
h11_ccbar_pair_before_hadronization_dphi.Write()
h12_ccbar_hadron_descendants_dphi.Write()
h13_ccbar_parton_descendants_dphi.Write()
h14_ccbar_init_showers_descendants_dphi.Write()
h15_ccbar_final_showers_descendants_dphi.Write()

output_file.Close()
