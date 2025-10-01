import pythia8
import ROOT

# Initialize Pythia
pythia = pythia8.Pythia()

# # Set the random seed (if desired)
# seed_value = 12345  # Example seed value
# pythia.readString(f"Random:seed = {seed_value}")

pythia.readString("Beams:eCM = 14000.")  # 14 TeV collisions
pythia.readString("HardQCD:all = on")    # Enable QCD processes
pythia.readString("PhaseSpace:pTHatMin = 20.")  # Minimum pT for hard scatterings

# Enable quarkonium production
# pythia.readString("Charmonium:all = on")
# pythia.readString("Bottomonium:all = on")

# from online:
# pythia.readString("Beams:eCM = 14000.")        # collision energy
# pythia.readString("HardQCD:hardccbar = on")    # enable direct heavy-quark production
# # alternatively: HardQCD:all = on  (this also gives ccbar but includes lots of light-QCD)
# pythia.readString("PhaseSpace:pTHatMin = 5.")  # tune to your pT-range of interest
# # keep parton shower and hadronization on (defaults). turn off for parton studies
# # optionally turn off quarkonium if you only want open charm:
# pythia.readString("Charmonium:all = off")     # don't generate c-cbar bound states (if undesired)
# # decays/hadronization default are on; you usually want that to see final D mesons

pythia.init()


# Create ROOT file and TTree
root_file = ROOT.TFile("pythia_100_events_pT_hard.root", "RECREATE")
tree = ROOT.TTree("events", "Pythia8 Full Event Record")

# Vectors to store event information
event_id = ROOT.std.vector('int')()
pid = ROOT.std.vector('int')()
status = ROOT.std.vector('int')()
px = ROOT.std.vector('float')()
py = ROOT.std.vector('float')()
pz = ROOT.std.vector('float')()
E = ROOT.std.vector('float')()
mother_1 = ROOT.std.vector('int')()
mother_2 = ROOT.std.vector('int')()
daughter_1 = ROOT.std.vector('int')()
daughter_2 = ROOT.std.vector('int')()

# Attach vectors to the ROOT tree
tree.Branch("event_id", event_id)
tree.Branch("pid", pid)  # Particle ID (PDG code)
tree.Branch("status", status)  # Pythia status code
tree.Branch("px", px)
tree.Branch("py", py)
tree.Branch("pz", pz)
tree.Branch("E", E)
tree.Branch("mother_1", mother_1)  # Mother indices
tree.Branch("mother_2", mother_2)
tree.Branch("daughter_1", daughter_1)  # Daughter indices
tree.Branch("daughter_2", daughter_2)

# Event loop
for iEvent in range(100):  # Generate N events
    if not pythia.next():
        continue

    # Clear vectors before filling them with new event data
    event_id.clear()
    pid.clear()
    status.clear()
    px.clear()
    py.clear()
    pz.clear()
    E.clear()
    mother_1.clear()
    mother_2.clear()
    daughter_1.clear()
    daughter_2.clear()

    # Loop over all particles
    for i in range(pythia.event.size()):
        p = pythia.event[i]

        event_id.push_back(iEvent)
        pid.push_back(p.id())
        status.push_back(p.status())
        px.push_back(p.px())
        py.push_back(p.py())
        pz.push_back(p.pz())
        E.push_back(p.e())

        # Store mother and daughter indices
        mother_1.push_back(p.mother1())  # First mother index
        mother_2.push_back(p.mother2())  # Second mother index
        daughter_1.push_back(p.daughter1())  # First daughter index
        daughter_2.push_back(p.daughter2())  # Second daughter index

    # Fill tree with the event data
    tree.Fill()

# Save and close the ROOT file
root_file.Write()
root_file.Close()

# Print Pythia statistics
pythia.stat()
