import ROOT
import array
import numpy as np
from ROOT_analysis_functions import (
    caleta,
)

pdgDB = ROOT.TDatabasePDG.Instance()

def is_charged_pdg(pdgid: int) -> bool:
    part = pdgDB.GetParticle(pdgid)
    if not part:
        return False
    # Charge() returns charge in units of e/3 in ROOT's PDG DB
    return abs(part.Charge()) > 0.0

# Open the ROOT file
root_file = ROOT.TFile("/home/daniel/LibraFiles/CleanThesis/PythiaData/pythia_spring_MB_1M_heavyion_events.root", "READ")
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

# V0A estimating mult
hV0A = ROOT.TH1I("hV0A", "Final State Charged Particles", 100, 0, 200)

# mult_values = []  # optional: store values too

# Loop over all events
# to fill, we must make sure that at least one charged particle is found in |eta| < 1 for an event
# particles that fill the histogram will be in the V0A forward region, must be charged, and must be final
# final --> status code greater than 0
# charged --> 
m_events = tree.GetEntries()
print(m_events, "events")
for i in range(m_events): # m_entries, edit to look at single event
    ccbar_pair_count = 0
    ccbar_hadron_daughters = 0
    if i % 1000 == 0:
        print(str(i) + " events have completed")
    tree.GetEntry(i)
    

    temp_mult = 0
    n_particles = len(pdg)
    for j in range(n_particles): # particle loop, total number of particles in event "i"
        # Particle loop variables
        particle_pdg = pdg[j]
        particle_status = status[j]
        particle_px = px[j]
        particle_py = py[j]
        particle_pz = pz[j]
        particle_eta = caleta(x_momentum=particle_px, y_momentum=particle_py, z_momentum=particle_pz)
        if is_charged_pdg(particle_pdg) and particle_status > 0 and np.abs(particle_eta) < 1:
            temp_mult += 1

    hV0A.Fill(temp_mult)
        
        

# ---- inside your event loop ----
# multV0A = ...
# hV0A.Fill(multV0A)
# mult_values.append(multV0A)
# -------------------------------

# Example: edges for 0–10–30–50–100% (i.e. top 0-10% highest mult etc.)
# ROOT quantiles are "x such that fraction q is below x"
# qs = array.array('d', [0.01, 0.05, 0.10, 0.15, 0.20, 0.30, 0.40, 0.50, 0.70, 1.00])  # 50th, 70th, 90th percentiles
# xq = array.array('d', [0.0]*len(qs))

# hV0A.GetQuantiles(len(qs), xq, qs)

# print("Quantile edges:")
# for q, x in zip(qs, xq):
#     print(f"  q={q:.2f} -> Nch^V0A = {x:.1f}")

# lines = []
# for x in xq:
#     ln = ROOT.TLine(x, 0, x, hV0A.GetMaximum())
#     ln.SetLineStyle(2)
#     ln.SetLineWidth(2)
#     ln.Draw("same")
#     lines.append(ln)

# leg = ROOT.TLegend(0.55, 0.65, 0.88, 0.88)
# leg.AddEntry(hV0A, "V0A estimator", "l")

# leg.AddEntry(lines[0], "50%", "l")
# leg.AddEntry(lines[1], "70%", "l")
# leg.AddEntry(lines[2], "90%", "l")
# leg.Draw()



# # Save the histogram as a ROOT file
# output_file = ROOT.TFile("/home/daniel/LibraFiles/CleanThesis/RootOutputs/02_05_2026_multiplicity_intervals_1.root", "RECREATE")

# # Write out all histograms
# hV0A.Write()


# output_file.Close()

# # Interpretation:
# #  - events with multV0A >= xq[2] are in the top 10% (0–10% highest mult)
# #  - xq[1]..xq[2] are 10–30% highest mult, etc.

# --- quantiles ---
qs = array.array('d', [0.33, 0.66, 0.99])
xq = array.array('d', [0.0]*len(qs))

# Make sure histogram has entries!
if hV0A.GetEntries() == 0:
    raise RuntimeError("hV0A has 0 entries; cannot compute quantiles.")

hV0A.GetQuantiles(len(qs), xq, qs)

print("Quantile edges:")
for q, x in zip(qs, xq):
    print(f"  q={q:.2f} -> Nch^V0A = {x:.1f}")

# --- draw onto a canvas (this is what you'll save) ---
c = ROOT.TCanvas("cV0A", "Final state charged particles multiplicity with percentile cuts", 900, 700)
c.cd()

hV0A.SetStats(0)
hV0A.Draw("HIST")

# Use the drawn maximum (safer) and pad coords
ROOT.gPad.Update()
ymax = ROOT.gPad.GetUymax()

lines = []
for x in xq:
    ln = ROOT.TLine(x, 0.0, x, ymax)
    ln.SetLineStyle(2)
    ln.SetLineWidth(2)
    ln.Draw()
    lines.append(ln)

leg = ROOT.TLegend(0.55, 0.55, 0.88, 0.88)
leg.SetBorderSize(0)
leg.SetFillStyle(0)
leg.AddEntry(hV0A, "V0A estimator", "l")

# Label each line with the percentile corresponding to q
for ln, q in zip(lines, qs):
    leg.AddEntry(ln, f"{100*q:.0f}th pct", "l")

leg.Draw()
c.Modified()
c.Update()

# --- write everything to output file ---
output_file = ROOT.TFile("/home/daniel/LibraFiles/CleanThesis/RootOutputs/03_05_2026_multiplicity_intervals_1.root", "RECREATE")
output_file.cd()

hV0A.Write("hV0A")                 # histogram itself
c.Write("cV0A_with_quantiles")     # canvas with lines+legend embedded

# (Optional) also save an image for quick viewing
c.SaveAs("/home/daniel/LibraFiles/CleanThesis/RootOutputs/03_05_2026_multiplicity_intervals_HI.png")

output_file.Close()

