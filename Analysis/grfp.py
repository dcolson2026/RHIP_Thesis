import ROOT
import time
import numpy as np
from particle import Particle
from ROOT_analysis_functions import (
    caldeltaphi,
    particle_lineage_dfs,
    calpT,
    calphi,
    getParticleName,
    delta_phi_correlation,
    caleta,
    following_charm_dfs,
    fit_von_mises_region,
)

PI = np.pi

# Open the ROOT file
root_file = ROOT.TFile("/home/daniel/LibraFiles/CleanThesis/PythiaData/pythia_10Mil_Josey_events.root", "READ")
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
# dphi of c-cbar initial pair, pdg of last c descendant that still has charm, same for cbar, dphi of those 2
h1_all_ccbar_delta_phi = ROOT.TH1F("h1_all_ccbar_delta_phi", "#Delta#phi Correlation: c - cbar quarks; #Delta#phi (radians); Counts", 30, -PI/2, 3*PI/2)
h2_last_c_pdg = ROOT.TH1F("h2_last_c_pdg", "Last charm descendant PDG; PDG ID; Counts", 2001, -1000.5, 1000.5)
h3_last_cbar_pdg = ROOT.TH1F("h3_last_cbar_pdg", "Last anti-charm descendant PDG; PDG ID; Counts", 2001, -1000.5, 1000.5)
h4_last_ccbar_delta_phi = ROOT.TH1F("h4_last_ccbar_delta_phi", "#Delta#phi: Last charm vs last anti-charm descendant; #Delta#phi (radians); Counts", 30, -PI/2, 3*PI/2)
# now look at object that fills h2 and h3. get more specific and make pdg and pT and eta cuts
h5_last_c_descendant_D0_pT = ROOT.TH1F("h5", "Last charm descendant: pT of D0", 50, 0, 20)
h6_last_cbar_descendant_D0_pT = ROOT.TH1F("h6", "Last anti-charm descendant: pT of anti-D0", 50, 0, 20)
# low, med, high momentum D0 dphi
h7_low_pT_D0_descendants_dphi = ROOT.TH1F("h7", "#Delta#phi D0 (0 <= pT <= 2)", 30, -PI/2, 3*PI/2)
h8_med_pT_D0_descendants_dphi = ROOT.TH1F("h8", "#Delta#phi D0 (2 < pT <= 4)", 30, -PI/2, 3*PI/2)
h9_high_pT_D0_descendants_dphi = ROOT.TH1F("h9", "#Delta#phi D0 (4 < pT)", 30, -PI/2, 3*PI/2)

large_dict = {}
# Loop over all events
m_events = tree.GetEntries()
print(m_events, "events")
for i in range(500000): # m_entries, edit to look at single event
    ccbar_pair_count = 0
    ccbar_hadron_daughters = 0
    if i % 1000 == 0:
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

    # if len(c_quark_daughter_indices_list) != 0 or len(cbar_quark_daughter_indices_list) != 0:
    #     print(len(c_quark_daughter_indices_list), c_quark_daughter_indices_list)
    #     print(len(cbar_quark_daughter_indices_list), cbar_quark_daughter_indices_list)
    #     print()
    
    if c_quark_daughter_indices_list:
        last_c_descendant_index = c_quark_daughter_indices_list[-1]
        last_c_descendant_px = px[last_c_descendant_index]
        last_c_descendant_py = py[last_c_descendant_index]
        last_c_descendant_pT = calpT(last_c_descendant_px, last_c_descendant_py)
        last_c_descendant_pdg = pdg[last_c_descendant_index]
        if 0 <= last_c_descendant_index < n_particles:
            h2_last_c_pdg.Fill(last_c_descendant_pdg)
            if last_c_descendant_pdg == 421:
                h5_last_c_descendant_D0_pT.Fill(last_c_descendant_pT)

            # if pdg[last_c_index] == 4:
                # print(daughter_1[last_c_index], daughter_2[last_c_index])
                # print("weak?", pdg[daughter_1[last_c_index]])

    if cbar_quark_daughter_indices_list:
        last_cbar_descendant_index = cbar_quark_daughter_indices_list[-1]
        last_cbar_descendant_px = px[last_cbar_descendant_index]
        last_cbar_descendant_py = py[last_cbar_descendant_index]
        last_cbar_descendant_pT = calpT(last_cbar_descendant_px, last_cbar_descendant_py)
        last_cbar_descendant_pdg = pdg[last_cbar_descendant_index]
        if 0 <= last_cbar_descendant_index < n_particles:
            h3_last_cbar_pdg.Fill(last_cbar_descendant_pdg)
            if last_cbar_descendant_pdg == -421:
                h6_last_cbar_descendant_D0_pT.Fill(last_cbar_descendant_pT)

    if c_quark_daughter_indices_list and cbar_quark_daughter_indices_list:
        last_c_index = c_quark_daughter_indices_list[-1]
        last_cbar_index = cbar_quark_daughter_indices_list[-1]
        if 0 <= last_c_index < n_particles and 0 <= last_cbar_index < n_particles:
            last_c_phi = calphi(y_momentum=py[last_c_index], x_momentum=px[last_c_index])
            last_cbar_phi = calphi(y_momentum=py[last_cbar_index], x_momentum=px[last_cbar_index])
            last_descendants_dphi = caldeltaphi(last_c_phi, last_cbar_phi)
            h4_last_ccbar_delta_phi.Fill(last_descendants_dphi)
            pair_key = (pdg[last_c_index], pdg[last_cbar_index])
            large_dict.setdefault(pair_key, []).append(last_descendants_dphi)
            # if last_c_descendant_pdg == 421 and last_cbar_descendant_pdg == -421:
            #     if 0 <= last_c_descendant_pT <= 2 and 0 <= last_cbar_descendant_pT <= 2:
            #         h7_low_pT_D0_descendants_dphi.Fill(last_descendants_dphi)
            #     elif 2 < last_c_descendant_pT <= 4 and 2 < last_cbar_descendant_pT <= 4:
            #         h8_med_pT_D0_descendants_dphi.Fill(last_descendants_dphi)
            #     elif 4 < last_c_descendant_pT and 4 < last_cbar_descendant_pT:
            #         h9_high_pT_D0_descendants_dphi.Fill(last_descendants_dphi)



pair_histograms = []
pair_fit_functions = []
pair_fit_labels = []
for (pdg_from_c, pdg_from_cbar), delta_phi_values in large_dict.items():
    if not (pdg_from_c == 421 and pdg_from_cbar == -421): continue # Just D0 anti-D0
    hist_name = f"h_delta_phi_{pdg_from_c}_{pdg_from_cbar}"
    c_name = getParticleName(pdg_from_c)
    cbar_name = getParticleName(pdg_from_cbar)
    hist_title = f"#Delta#phi: D^{{0}}- #bar{{D}}^{{0}}; #Delta#phi (rad); D^{{0}}- #bar{{D}}^{{0}} pairs / bin"
    hist = ROOT.TH1F(hist_name, hist_title, 30, -PI/2, 3*PI/2)
    for value in delta_phi_values:
        hist.Fill(value)
    hist.SetLineWidth(2)
    pair_histograms.append(hist)

    fits_for_hist = []
    if hist.GetEntries() > 0:
        bin_margin = hist.GetBinWidth(1)
        near_window_min = -PI / 2 - bin_margin
        near_window_max = PI / 2 + bin_margin
        away_window_min = PI / 2 - bin_margin
        away_window_max = 3 * PI / 2 + bin_margin

        near_result = fit_von_mises_region(
            hist,
            f"{hist_name}_von_mises_near",
            near_window_min,
            near_window_max,
            0.0,
            ROOT.kAzure + 2,
            min_entries=15,
        )
        if near_result:
            near_fit, near_chi2, near_ndf = near_result
            pair_fit_functions.append(near_fit)
            fits_for_hist.append(("Near", near_fit, near_chi2, near_ndf))

        away_result = fit_von_mises_region(
            hist,
            f"{hist_name}_von_mises_away",
            away_window_min,
            away_window_max,
            PI,
            ROOT.kRed + 1,
            min_entries=15,
        )
        if away_result:
            away_fit, away_chi2, away_ndf = away_result
            pair_fit_functions.append(away_fit)
            fits_for_hist.append(("Away", away_fit, away_chi2, away_ndf))

    # for idx, (label, fit_obj, chi2, ndf) in enumerate(fits_for_hist):
    #     # label_text = f"{label}: #chi^{{2}}/ndf = {chi2:.1f}/{int(ndf) if ndf else 0}"
    #     text_y = 0.88 - 0.06 * idx
    #     latex = ROOT.TLatex(0.18, text_y, label_text)
    #     latex.SetNDC()
    #     latex.SetTextSize(0.035)
    #     latex.SetTextColor(fit_obj.GetLineColor())
    #     latex.SetTextFont(42)
    #     hist.GetListOfFunctions().Add(latex)
    #     pair_fit_labels.append(latex)

# Style for ROOT histogram
# h1_delta_phi.SetLineColor(ROOT.kBlue)
h1_all_ccbar_delta_phi.SetLineWidth(2)
h2_last_c_pdg.SetLineWidth(2)
h3_last_cbar_pdg.SetLineWidth(2)
h4_last_ccbar_delta_phi.SetLineWidth(2)
h5_last_c_descendant_D0_pT.SetLineWidth(2)
h6_last_cbar_descendant_D0_pT.SetLineWidth(2)

# Save the histogram as a ROOT file
output_file = ROOT.TFile("/home/daniel/LibraFiles/CleanThesis/RootOutputs/11_14_2025_grfp_500K.root", "RECREATE")

# Write out all histograms
h1_all_ccbar_delta_phi.Write()
h2_last_c_pdg.Write()
h3_last_cbar_pdg.Write()
h4_last_ccbar_delta_phi.Write()
h5_last_c_descendant_D0_pT.Write()
h6_last_cbar_descendant_D0_pT.Write()
h7_low_pT_D0_descendants_dphi.Write()
h8_med_pT_D0_descendants_dphi.Write()
h9_high_pT_D0_descendants_dphi.Write()
for hist in pair_histograms:
    hist.Write()
for fit_func in pair_fit_functions:
    fit_func.Write()


output_file.Close()
