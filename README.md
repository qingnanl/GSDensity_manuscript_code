# GSDensity_manuscript_code

Simulation scRNA-seq data with SERGIO:
sergio_steady_gen.R: generate SERGIO input file (steady state mode)

run_sergio_code.py: run SERGIO (steady state mode)

sergio_dyn_gen.R: generate SERGIO input file (dynamic state mode)

run_sergio_dyn.1.py: run SERGIO(dynamic mode only)

normalize_sim_rand.R: combine random gene matrix with simulated gene matrix (only dynamic mode)

bMat_cID7.tab; bMat_cID10.mod.tab: cell state migration file for dynamic modes

Add_noise_to_data.ipynb: add noise to simulated data

sergio_analysis_functions.R; method.benchmarking.functions.R: a set of functions for simulated data benchmarking; 
data.sum <- benchmark.sergio(design = "5000.gene.3ct", version.n = 1, drop.out = "dp1", 
                             n.cell.types = 3, n.cells.per.type = 1500, n.genes = 5000)



Real world RNA-seq data:
curate_data.R: curation of real-world datasets and markers

Cell_marker_Mouse.xlsx; Cell_marker_Human.xlsx; st4.csv; PanglaoDB_markers_27_Mar_2020.tsv.gz: markers

method.benchmarking.functions.rd.R: a set of functions for real world data benchmarking; 
  run benchmark.real(fname = "pbmc3k") with "pbmc3k.rds" and "pbmc3k.markers.rds" in the dir.

Applications:

tnbc_preprocess.r; tnbc.subtype.cc.R: processing tnbc data and run cellchat

mb.r: mouse brain spatial

mouse_brain_dev.r: mouse brain development

panc.preprocess.R; panc.process.R: cancer visium data analysis
