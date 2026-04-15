import ROOT
# # from ROOT_analysis_functions import fit_von_mises_region
# # import math

# # h = ROOT.TH1F("h_test", "Toy von Mises;#Delta#phi;Counts", 30, -math.pi/2, 3*math.pi/2)

# # # Fill with something vaguely von-Mises-shaped around 0 and pi:
# # for i in range(100000):
# #     # e.g. a Gaussian around 0, just to have a peak
# #     x = ROOT.gRandom.Gaus(0, 0.4)
# #     if x < -math.pi/2: x += 2*math.pi
# #     if x >= 3*math.pi/2: x -= 2*math.pi
# #     h.Fill(x)

# # near = fit_von_mises_region(
# #     h, "h_test_von_mises_near",
# #     -math.pi/2, math.pi/2, 0.0,
# #     ROOT.kBlue+1, min_entries=10,
# # )

# # c = ROOT.TCanvas("c", "test", 800, 600)
# # h.Draw("E")
# # if near:
# #     near[0].Draw("same")
# # c.SaveAs("test_von_mises.pdf")

# import ROOT

# pdgDB = ROOT.TDatabasePDG.Instance()

# def is_charged_pdg(pdgid: int) -> bool:
#     part = pdgDB.GetParticle(pdgid)
#     if not part:
#         return False
#     # Charge() returns charge in units of e/3 in ROOT's PDG DB
#     return abs(part.Charge()) > 0.0

# def charge_e(pdgid: int) -> float:
#     part = pdgDB.GetParticle(pdgid)
#     if not part:
#         return 0.0
#     return part.Charge() / 3.0

# print(is_charged_pdg(333))

f = ROOT.TFile.Open("/home/daniel/LibraFiles/CleanThesis/PythiaData/pythia_spring_MB_9M_events_1.root")
f.ls()
obj = f.Get("events")
print(obj)
