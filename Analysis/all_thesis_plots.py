import ROOT
import time
import numpy as np
import ctypes
from particle import Particle
from datetime import date
from ROOT_analysis_functions import (
    caldeltaphi,
    # particle_lineage_dfs,
    calpT,
    calphi,
    getParticleName,
    # delta_phi_correlation,
    caleta,
    following_charm_dfs,
    fit_von_mises_region,
    is_charged_pdg,
    get_bin_label
)
start_time = time.perf_counter()
today = date.today()

PI = np.pi


def use_visible_errors(hist):
    hist.Sumw2()
    hist.SetOption("E1")
    hist.SetMarkerStyle(20)
    hist.SetMarkerSize(0.7)

MULT_CLASSES_MB = [7.9, 19.6, 300] # was 65.7 but upper bound should not matter
MULT_CLASSES_HARD = [29.2, 46.3, 300] # was 90.2 but upper bound should not matter
MULT_CLASSES_HI = [23.4, 51.1, 300] # was 128.6 but upper bound should not matter

# fpath = "/home/daniel/LibraFiles/CleanThesis/PythiaData/pythia_spring_MB_9M_events.root"
# fpath = "/home/daniel/LibraFiles/CleanThesis/PythiaData/pythia_spring_hard_9M_events.root"
# fpath = "/home/daniel/LibraFiles/CleanThesis/PythiaData/pythia_spring_minbias_1M_heavyion_events.root"

# file_list = ["/home/daniel/LibraFiles/CleanThesis/PythiaData/pythia_spring_MB_9M_events.root"]
file_list = ["/home/daniel/LibraFiles/CleanThesis/PythiaData/pythia_spring_hard_9M_events.root", "/home/daniel/LibraFiles/CleanThesis/PythiaData/pythia_spring_hard_9M_events_1.root", "/home/daniel/LibraFiles/CleanThesis/PythiaData/pythia_spring_hard_9M_events_2.root"]
# file_list = ["/home/daniel/LibraFiles/CleanThesis/PythiaData/pythia_spring_minbias_1M_heavyion_events.root"]

if "MB" in file_list[0]:
    file_type = "MB"
    mult_edges = MULT_CLASSES_MB
elif "hard" in file_list[0]:
    file_type = "hard"
    mult_edges = MULT_CLASSES_HARD
elif "heavyion" in file_list[0]:
    file_type = "HI"
    mult_edges = MULT_CLASSES_HI
else:
    raise ValueError(f"Could not determine file type from fpath: {file_list}")
# PT_BINS = [
#     ("low_pT", 0.0, 2.0),
#     ("mid_pT", 2.0, 4.0),
#     ("high_pT", 4.0, 8.0),
# ]
# make finer pT binning, also change pT_pretty
# PT_BINS = [
#     ("first_pT", 0.0, 1.0),
#     ("second_pT", 1.0, 2.0),
#     ("third_pT", 2.0, 3.0),
#     ("fourth_pT", 3.0, 4.0),
#     ("fifth_pT", 4.0, 5.0),
#     ("sixth_pT", 5.0, 6.0),
#     ("seventh_pT", 6.0, 7.0),
#     ("eighth_pT", 7.0, 8.0),
# ]
PT_BINS = [
    ("first_pT", 0.0, 0.5),
    ("second_pT", 0.5, 1.0),
    ("third_pT", 1.0, 1.5),
    ("fourth_pT", 1.5, 2.0),
    ("fifth_pT", 2.0, 2.5),
    ("sixth_pT", 2.5, 3.0),
    ("seventh_pT", 3.0, 3.5),
    ("eighth_pT", 3.5, 8.0),
]
MULT_BINS = [
    ("low_mult", 0.0, mult_edges[0]),
    ("mid_mult", mult_edges[0], mult_edges[1]),
    ("high_mult", mult_edges[1], mult_edges[2]),
]
h_dihadron_dphi = {}
for mult_label, mult_lo, mult_hi in MULT_BINS:
    for pt_label, pt_lo, pt_hi in PT_BINS:
        hname = f"h_dphi_D0antiD0_{pt_label}_{mult_label}"
        htitle = (
            f"#Delta#phi of D^{{0}}-#bar{{D}}^{{0}} pairs "
            f"({pt_label}, {mult_label});"
            f"#Delta#phi [rad];Counts"
        )
        hist = ROOT.TH1F(hname, htitle, 64, -ROOT.TMath.Pi()/2, 3*ROOT.TMath.Pi()/2)
        use_visible_errors(hist)
        h_dihadron_dphi[(pt_label, mult_label)] = hist
pt_pretty = {
    "first_pT": "0 < p_{T} #leq 0.5 GeV",
    "second_pT": "0.5 < p_{T} #leq 1 GeV",
    "third_pT": "1 < p_{T} #leq 1.5 GeV",
    "fourth_pT": "1.5 < p_{T} #leq 2 GeV",
    "fifth_pT": "2 < p_{T} #leq 2.5 GeV",
    "sixth_pT": "2.5 < p_{T} #leq 3.0 GeV",
    "seventh_pT": "3.0 < p_{T} #leq 3.5 GeV",
    "eighth_pT": "3.5 < p_{T} #leq 8 GeV",
}

mult_pretty = {
    "low_mult": "low multiplicity",
    "mid_mult": "mid multiplicity",
    "high_mult": "high multiplicity",
}

outpath = f"/home/daniel/LibraFiles/CleanThesis/RootOutputs/{today.month}_{today.day}_{today.year}_thesis_{file_type}.root"
width_outpath = (
    f"/home/daniel/LibraFiles/CleanThesis/RootOutputs/"
    f"{today.month}_{today.day}_{today.year}_thesis_{file_type}_von_mises_widths.txt"
)
# Open the ROOT file
# root_file = ROOT.TFile(fpath, "READ") # USER INPUT
#root_file = ROOT.TFile("pythia_100_events_pT_hard.root", "READ")
# tree = root_file.Get("events") # fetches the tree
tree = ROOT.TChain("events")

for fpath in file_list:
    tree.Add(fpath)


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
h1_all_ccbar_delta_phi = ROOT.TH1F("h1_all_ccbar_delta_phi", "Initial c-#bar{c} Pairs Correlation; #Delta#phi (radians); Counts", 30, -PI/2, 3*PI/2)
h2_last_c_pdg = ROOT.TH1F("h2_last_c_pdg", "PDG of Last Charmed Descendant of c; PDG ID; Counts", 2001, -1000.5, 1000.5)
h3_last_cbar_pdg = ROOT.TH1F("h3_last_cbar_pdg", "PDG of Last Charmed Descendant of #bar{c}; PDG ID; Counts", 2001, -1000.5, 1000.5)
h4_last_ccbar_delta_phi = ROOT.TH1F("h4_last_ccbar_delta_phi", "Last Charmed Descendant of c - Last Charmed Descendant of #bar{c} Correlation; #Delta#phi (radians); Counts", 30, -PI/2, 3*PI/2)

# now look at object that fills h2 and h3. get more specific and make pdg and pT and eta cuts
# h5_last_c_descendant_D0_pT = ROOT.TH1F("h5", "Last charm descendant: pT of D0", 50, 0, 20)
# h6_last_cbar_descendant_D0_pT = ROOT.TH1F("h6", "Last anti-charm descendant: pT of anti-D0", 50, 0, 20)
# low, med, high momentum D0 dphi
h5_D0_pT_dist = ROOT.TH1F("h5_D0_pT_dist", "p_{T} of D^{0};p_{T}(D^{0}) [GeV], Counts", 50, 0, 8)
h6_antiD0_pT_dist = ROOT.TH1F("h6_antiD0_pT_dist", "p_{T} of #bar{D}^{0};p_{T}(#bar{D}^{0}) [GeV], Counts", 50, 0, 8)

# h7_low_pT_D0_descendants_dphi = ROOT.TH1F("h7", "#Delta#phi D0 (0 <= pT <= 2)", 30, -PI/2, 3*PI/2)
# h8_med_pT_D0_descendants_dphi = ROOT.TH1F("h8", "#Delta#phi D0 (2 < pT <= 4)", 30, -PI/2, 3*PI/2)
# h9_high_pT_D0_descendants_dphi = ROOT.TH1F("h9", "#Delta#phi D0 (4 < pT <= 8)", 30, -PI/2, 3*PI/2)



h10_D0_vs_antiD0_pT = ROOT.TH2F(
    "h10_D0_vs_antiD0_pT",
    "p_{T}(D^{0}) vs p_{T}(#bar{D}^{0});p_{T}(D^{0}) [GeV];p_{T}(#bar{D}^{0}) [GeV];D^{0}-#bar{D}^{0} pairs / bin",
    50, 0, 8,          # X axis: D0 pT
    50, 0, 8           # Y axis: anti-D0 pT
)
h10_D0_vs_antiD0_pT.SetOption("COLZ")

h10_y_equals_x = ROOT.TLine(0, 0, 8, 8)
h10_y_equals_x.SetLineColor(ROOT.kRed + 1)
h10_y_equals_x.SetLineStyle(2)
h10_y_equals_x.SetLineWidth(2)
h10_D0_vs_antiD0_pT.GetListOfFunctions().Add(h10_y_equals_x)

h10_D0_vs_antiD0_pT_yx_asymmetry = ROOT.TH2F(
    "h10_D0_vs_antiD0_pT_yx_asymmetry",
    "Mirror-bin asymmetry for p_{T}(D^{0}) vs p_{T}(#bar{D}^{0});"
    "p_{T}(D^{0}) [GeV];p_{T}(#bar{D}^{0}) [GeV];"
    "(N_{ij} - N_{ji}) / (N_{ij} + N_{ji})",
    50, 0, 8,
    50, 0, 8
)
h10_D0_vs_antiD0_pT_yx_asymmetry.SetOption("COLZ")

large_dict = {}
# Loop over all events
m_events = tree.GetEntries()
print(m_events, "events in total")
print(file_type)
for i in range(1000000): # m_entries, edit to look at single event
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
    
    temp_mult = 0
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

        # filling the simple D0 and anti D0 pT histograms
        if particle_pdg == 421:
            h5_D0_pT_dist.Fill(particle_pT)
        elif particle_pdg == -421:
            h6_antiD0_pT_dist.Fill(particle_pT)

        if is_charged_pdg(particle_pdg) and particle_status > 0 and np.abs(particle_eta) < 1:
            temp_mult += 1

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
    temp_D0_antiD0 = [False, False]

    if c_quark_daughter_indices_list:
        last_c_descendant_index = c_quark_daughter_indices_list[-1]
        last_c_descendant_px = px[last_c_descendant_index]
        last_c_descendant_py = py[last_c_descendant_index]
        last_c_descendant_pT = calpT(last_c_descendant_px, last_c_descendant_py)
        last_c_descendant_pdg = pdg[last_c_descendant_index]
        if 0 <= last_c_descendant_index < n_particles:
            h2_last_c_pdg.Fill(last_c_descendant_pdg)
            if last_c_descendant_pdg == 421:
                # h5_last_c_descendant_D0_pT.Fill(last_c_descendant_pT)
                temp_D0_antiD0[0] = True

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
                # h6_last_cbar_descendant_D0_pT.Fill(last_cbar_descendant_pT)
                temp_D0_antiD0[1] = True
    
    # if no False, then we have D0 and antiD0. then just add pT's
    if False not in temp_D0_antiD0:
        h10_D0_vs_antiD0_pT.Fill(last_c_descendant_pT, last_cbar_descendant_pT)


    # if both lists exist, do a dihadron distribution i guess
    # if c_quark_daughter_indices_list and cbar_quark_daughter_indices_list:
    #     last_c_index = c_quark_daughter_indices_list[-1]
    #     last_cbar_index = cbar_quark_daughter_indices_list[-1]
    #     if 0 <= last_c_index < n_particles and 0 <= last_cbar_index < n_particles:
    #         last_c_eta = caleta(x_momentum=px[last_c_index], y_momentum=py[last_c_index], z_momentum=pz[last_c_index])
    #         last_cbar_eta = caleta(x_momentum=px[last_cbar_index], y_momentum=py[last_cbar_index], z_momentum=pz[last_cbar_index])
    #         if np.abs(last_c_eta) > 1 or np.abs(last_cbar_eta) > 1: 
    #             # print("not in detector visibility, cya!")
    #             continue
    #         if MULT_FLOOR < temp_mult <= MULT_CEILING:
    #             # only look at certain multiplicity range!
    #             continue
    #         last_c_pT = calpT(x_momentum=px[last_c_index], y_momentum=py[last_c_index])
    #         last_cbar_pT = calpT(x_momentum=px[last_cbar_index], y_momentum=py[last_cbar_index])

    #         last_c_phi = calphi(y_momentum=py[last_c_index], x_momentum=px[last_c_index])
    #         last_cbar_phi = calphi(y_momentum=py[last_cbar_index], x_momentum=px[last_cbar_index])
    #         last_descendants_dphi = caldeltaphi(last_c_phi, last_cbar_phi)
    #         h4_last_ccbar_delta_phi.Fill(last_descendants_dphi)
    #         pair_key = (pdg[last_c_index], pdg[last_cbar_index])
    #         large_dict.setdefault(pair_key, []).append(last_descendants_dphi)
    #         # if last_c_descendant_pdg == 421 and last_cbar_descendant_pdg == -421:
    #         #     if 0 <= last_c_descendant_pT <= 2 and 0 <= last_cbar_descendant_pT <= 2:
    #         #         h7_low_pT_D0_descendants_dphi.Fill(last_descendants_dphi)
    #         #     elif 2 < last_c_descendant_pT <= 4 and 2 < last_cbar_descendant_pT <= 4:
    #         #         h8_med_pT_D0_descendants_dphi.Fill(last_descendants_dphi)
    #         #     elif 4 < last_c_descendant_pT and 4 < last_cbar_descendant_pT:
    #         #         h9_high_pT_D0_descendants_dphi.Fill(last_descendants_dphi)
    if c_quark_daughter_indices_list and cbar_quark_daughter_indices_list:

        last_c_index = c_quark_daughter_indices_list[-1]
        last_cbar_index = cbar_quark_daughter_indices_list[-1]

        if not (0 <= last_c_index < n_particles and 0 <= last_cbar_index < n_particles):
            continue

        # Require exact D0 and anti-D0
        if not (pdg[last_c_index] == 421 and pdg[last_cbar_index] == -421):
            continue

        # eta acceptance
        last_c_eta = caleta(
            x_momentum=px[last_c_index],
            y_momentum=py[last_c_index],
            z_momentum=pz[last_c_index]
        )
        last_cbar_eta = caleta(
            x_momentum=px[last_cbar_index],
            y_momentum=py[last_cbar_index],
            z_momentum=pz[last_cbar_index]
        )

        # if np.abs(last_c_eta) > 1 or np.abs(last_cbar_eta) > 1:
        #     continue

        # multiplicity classification
        mult_label = get_bin_label(temp_mult, MULT_BINS)
        if mult_label is None:
            print("test?", temp_mult)
            continue

        # pT calculation
        last_c_pT = calpT(
            x_momentum=px[last_c_index],
            y_momentum=py[last_c_index]
        )
        last_cbar_pT = calpT(
            x_momentum=px[last_cbar_index],
            y_momentum=py[last_cbar_index]
        )

        # No longer require BOTH D0 and anti-D0 to be in the same pT class...
        pt_label_c = get_bin_label(last_c_pT, PT_BINS)
        pt_label_cbar = get_bin_label(last_cbar_pT, PT_BINS)

        if pt_label_c is None or pt_label_cbar is None:
            continue

        # ... because we commented this part out
        # if pt_label_c != pt_label_cbar:
        #     continue

        pt_label = (pt_label_c, pt_label_cbar)

        # delta phi
        last_c_phi = calphi(
            y_momentum=py[last_c_index],
            x_momentum=px[last_c_index]
        )
        last_cbar_phi = calphi(
            y_momentum=py[last_cbar_index],
            x_momentum=px[last_cbar_index]
        )
        last_descendants_dphi = caldeltaphi(last_c_phi, last_cbar_phi)

        # Fill only the relevant 1 of the 9 histograms
        # h_dihadron_dphi[(pt_label, mult_label)].Fill(last_descendants_dphi)

        pair_key = (
            pdg[last_c_index], 
            pdg[last_cbar_index], 
            pt_label_c, 
            pt_label_cbar,
            mult_label) # most likely this will error
        # pair_key = (pdg[last_c_index], pdg[last_cbar_index])
        large_dict.setdefault(pair_key, []).append(last_descendants_dphi)

symmetry_total_abs_diff = 0.0
symmetry_total_pairs = 0.0
symmetry_chi2 = 0.0
symmetry_ndf = 0
max_abs_asymmetry = 0.0
max_asymmetry_bin = None

for x_bin in range(1, h10_D0_vs_antiD0_pT.GetNbinsX() + 1):
    for y_bin in range(1, h10_D0_vs_antiD0_pT.GetNbinsY() + 1):
        content_xy = h10_D0_vs_antiD0_pT.GetBinContent(x_bin, y_bin)
        content_yx = h10_D0_vs_antiD0_pT.GetBinContent(y_bin, x_bin)
        pair_sum = content_xy + content_yx
        if pair_sum <= 0:
            continue

        asymmetry = (content_xy - content_yx) / pair_sum
        h10_D0_vs_antiD0_pT_yx_asymmetry.SetBinContent(x_bin, y_bin, asymmetry)

        if x_bin < y_bin:
            diff = content_xy - content_yx
            symmetry_total_abs_diff += abs(diff)
            symmetry_total_pairs += pair_sum
            symmetry_chi2 += diff * diff / pair_sum
            symmetry_ndf += 1

            abs_asymmetry = abs(asymmetry)
            if abs_asymmetry > max_abs_asymmetry:
                max_abs_asymmetry = abs_asymmetry
                max_asymmetry_bin = (x_bin, y_bin)

symmetry_norm_abs_diff = (
    symmetry_total_abs_diff / symmetry_total_pairs
    if symmetry_total_pairs > 0
    else 0.0
)
symmetry_chi2_per_ndf = symmetry_chi2 / symmetry_ndf if symmetry_ndf > 0 else 0.0

h10_D0_vs_antiD0_pT_yx_symmetry_check = ROOT.TH1F(
    "h10_D0_vs_antiD0_pT_yx_symmetry_check",
    "y=x symmetry check for h10_D0_vs_antiD0_pT;Metric;Value",
    4,
    0.5,
    4.5
)
h10_D0_vs_antiD0_pT_yx_symmetry_check.GetXaxis().SetBinLabel(1, "mirror bin pairs")
h10_D0_vs_antiD0_pT_yx_symmetry_check.GetXaxis().SetBinLabel(2, "norm abs diff")
h10_D0_vs_antiD0_pT_yx_symmetry_check.GetXaxis().SetBinLabel(3, "chi2/ndf")
h10_D0_vs_antiD0_pT_yx_symmetry_check.GetXaxis().SetBinLabel(4, "max abs asym")
h10_D0_vs_antiD0_pT_yx_symmetry_check.SetBinContent(1, symmetry_ndf)
h10_D0_vs_antiD0_pT_yx_symmetry_check.SetBinContent(2, symmetry_norm_abs_diff)
h10_D0_vs_antiD0_pT_yx_symmetry_check.SetBinContent(3, symmetry_chi2_per_ndf)
h10_D0_vs_antiD0_pT_yx_symmetry_check.SetBinContent(4, max_abs_asymmetry)
h10_D0_vs_antiD0_pT_yx_symmetry_check.SetOption("HIST TEXT0")

print("h10_D0_vs_antiD0_pT y=x symmetry check:")
print(f"  compared mirror bin pairs: {symmetry_ndf}")
print(f"  normalized absolute mirror-bin difference: {symmetry_norm_abs_diff:.6f}")
print(f"  chi2/ndf from mirror-bin differences: {symmetry_chi2_per_ndf:.6f}")
if max_asymmetry_bin is not None:
    x_bin, y_bin = max_asymmetry_bin
    x_center = h10_D0_vs_antiD0_pT.GetXaxis().GetBinCenter(x_bin)
    y_center = h10_D0_vs_antiD0_pT.GetYaxis().GetBinCenter(y_bin)
    print(
        "  largest fractional asymmetry: "
        f"{max_abs_asymmetry:.6f} at "
        f"pT(D0)={x_center:.3f} GeV, pT(anti-D0)={y_center:.3f} GeV"
    )

h10_D0_vs_antiD0_pT_surface = h10_D0_vs_antiD0_pT.Clone(
    "h10_D0_vs_antiD0_pT_surface"
)
h10_D0_vs_antiD0_pT_surface.SetTitle(
    "Surface plot of p_{T}(D^{0}) vs p_{T}(#bar{D}^{0});"
    "p_{T}(D^{0}) [GeV];p_{T}(#bar{D}^{0}) [GeV];"
    "D^{0}-#bar{D}^{0} pairs / bin"
)
h10_D0_vs_antiD0_pT_surface.SetOption("SURF2")



pair_histograms = []
pair_fit_functions = []
pair_fit_labels = []
pair_fit_annotations = []
width_rows = []
for (pdg_from_c, pdg_from_cbar, pt_label_c, pt_label_cbar, mult_label), delta_phi_values in large_dict.items():
    hist_name = f"h_delta_phi_{pdg_from_c}_{pdg_from_cbar}_{pt_label_c}_{pt_label_cbar}_{mult_label}"
    c_name = getParticleName(pdg_from_c)
    cbar_name = getParticleName(pdg_from_cbar)
    # hist_title = f"#Delta#phi: D^{{0}}- #bar{{D}}^{{0}}; #Delta#phi (rad); D^{{0}}- #bar{{D}}^{{0}} pairs / bin"
    # hist_title = f"#Delta#phi: {c_name} vs {cbar_name}; #Delta#phi (radians); Counts"
    hist_title = (
        f"#Delta#phi: D^{{0}} ({pt_pretty[pt_label_c]}) vs "
        f"#bar{{D}}^{{0}} ({pt_pretty[pt_label_cbar]}), "
        f"{mult_pretty[mult_label]};"
        f"#Delta#phi (radians);Counts"
        )
    hist = ROOT.TH1F(hist_name, hist_title, 30, -PI/2, 3*PI/2)
    use_visible_errors(hist)
    for value in delta_phi_values:
        hist.Fill(value)
    hist.SetLineWidth(2)
    pair_histograms.append(hist)

    fits_for_hist = []
    if hist.GetEntries() > 0:
        bin_margin = hist.GetBinWidth(1)
        # near_window_min = -PI / 2 - bin_margin
        # near_window_max = PI / 2 + bin_margin
        near_window_min = hist.GetXaxis().GetXmin()
        near_window_max = hist.GetXaxis().GetXmax()
        # away_window_min = PI / 2 - bin_margin
        # away_window_max = 3 * PI / 2 + bin_margin
        away_window_min = hist.GetXaxis().GetXmin()
        away_window_max = hist.GetXaxis().GetXmax()

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
    if fits_for_hist:
        axis = hist.GetXaxis()
        near_low_bin = axis.FindBin(-PI / 2 + 1e-6)
        near_high_bin = axis.FindBin(PI / 2 - 1e-6)
        away_low_bin = axis.FindBin(PI / 2 + 1e-6)
        away_high_bin = axis.FindBin(3 * PI / 2 - 1e-6)
        near_yield_err = ctypes.c_double(0.0)
        away_yield_err = ctypes.c_double(0.0)
        near_yield = hist.IntegralAndError(near_low_bin, near_high_bin, near_yield_err)
        away_yield = hist.IntegralAndError(away_low_bin, away_high_bin, away_yield_err)
        # x_center = 0.5 * (axis.GetXmin() + axis.GetXmax())
        y_max = hist.GetMaximum()
        if y_max <= 0:
            y_max = 1.0
        text_start_y = 0.8 * y_max
        text_step = 0.12 * y_max
        peak_colors = {
            "Near": ROOT.kAzure + 2,
            "Away": ROOT.kRed + 1,
        }
        for idx, (label, fit_obj, _, _) in enumerate(fits_for_hist):
            amplitude = fit_obj.GetParameter(0)
            kappa = fit_obj.GetParameter(2)
            kappa_err = fit_obj.GetParError(2)
            # Convert von Mises concentration into an approximate angular width.
            width = np.sqrt(1.0 / kappa) if kappa > 0 else 0.0
            width_err = 0.5 * kappa_err / (kappa ** 1.5) if kappa > 0 else 0.0
            peak_height = amplitude * np.exp(kappa)
            if label == "Near":
                peak_yield = near_yield
                peak_yield_err = near_yield_err.value
            else:
                peak_yield = away_yield
                peak_yield_err = away_yield_err.value
            width_rows.append({
                "hist_name": hist_name,
                "peak": label.lower(),
                "d0_pt_class": pt_label_c,
                "anti_d0_pt_class": pt_label_cbar,
                "multiplicity_class": mult_label,
                "yield_counts": peak_yield,
                "yield_err_counts": peak_yield_err,
                "width_rad": width,
                "width_err_rad": width_err,
                "kappa": kappa,
                "kappa_err": kappa_err,
                "amplitude": amplitude,
                "peak_height": peak_height,
            })
            text = f"{label}: Height={peak_height:.2f}, Width={width:.2f} rad"
            latex = ROOT.TLatex(3, text_start_y - idx * text_step, text)
            latex.SetTextAlign(21)  # center horizontally, align to top vertically
            latex.SetTextSize(0.03)
            latex.SetTextColor(peak_colors.get(label, ROOT.kBlack))
            latex.SetNDC(False)
            hist.GetListOfFunctions().Add(latex)
            pair_fit_annotations.append(latex)

# Style for ROOT histogram
# h1_delta_phi.SetLineColor(ROOT.kBlue)
h1_all_ccbar_delta_phi.SetLineWidth(2)
h2_last_c_pdg.SetLineWidth(2)
h3_last_cbar_pdg.SetLineWidth(2)
h4_last_ccbar_delta_phi.SetLineWidth(2)
h5_D0_pT_dist.SetLineWidth(2)
h6_antiD0_pT_dist.SetLineWidth(2)

# Save the histogram as a ROOT file
output_file = ROOT.TFile(outpath, "RECREATE")

# Write out all histograms
h1_all_ccbar_delta_phi.Write()
h2_last_c_pdg.Write()
h3_last_cbar_pdg.Write()
h4_last_ccbar_delta_phi.Write()
h5_D0_pT_dist.Write()
h6_antiD0_pT_dist.Write()
# h7_low_pT_D0_descendants_dphi.Write()
# h8_med_pT_D0_descendants_dphi.Write()
# h9_high_pT_D0_descendants_dphi.Write()
h10_D0_vs_antiD0_pT.Write()
h10_D0_vs_antiD0_pT_surface.Write()
h10_D0_vs_antiD0_pT_yx_asymmetry.Write()
h10_D0_vs_antiD0_pT_yx_symmetry_check.Write()
for hist in pair_histograms:
    hist.Write()
for fit_func in pair_fit_functions:
    fit_func.Write()

with open(width_outpath, "w", encoding="utf-8") as width_file:
    width_file.write(
        "hist_name\tpeak\td0_pt_class\tanti_d0_pt_class\tmultiplicity_class\t"
        "yield_counts\tyield_err_counts\twidth_rad\twidth_err_rad\t"
        "kappa\tkappa_err\tamplitude\tpeak_height\n"
    )
    for row in width_rows:
        width_file.write(
            f"{row['hist_name']}\t"
            f"{row['peak']}\t"
            f"{row['d0_pt_class']}\t"
            f"{row['anti_d0_pt_class']}\t"
            f"{row['multiplicity_class']}\t"
            f"{row['yield_counts']:.6f}\t"
            f"{row['yield_err_counts']:.6f}\t"
            f"{row['width_rad']:.6f}\t"
            f"{row['width_err_rad']:.6f}\t"
            f"{row['kappa']:.6f}\t"
            f"{row['kappa_err']:.6f}\t"
            f"{row['amplitude']:.6f}\t"
            f"{row['peak_height']:.6f}\n"
        )

end_time = time.perf_counter()
elapsed_time = end_time - start_time
# Print the result (formatted to 4 decimal places)
print(f"Program execution time: {elapsed_time:.4f} seconds")
print(f"Saved von Mises widths to {width_outpath}")

output_file.Close()
