'''
Shen et al. (2018) Science 362, 705-709

Individual design trajectories were filtered by the following criteria:
- the difference between the Rosetta energy of the bound (polymeric) and unbound (monomeric) states less than -15.0 Rosetta Energy Units
- interface surface area greater than 700 Å2
- shape complementarity greater than 0.62
- unsatisfied polar residues less than 5
'''


from biorazer_toolkit.apps.rosetta import RosettaApp

# Fill to the root directory of your Rosetta installation, like /mnt/c/WSLApplications/Rosetta/main
ROSETTA_DIR='/opt/rosetta.source.release-408/main' 

app = RosettaApp(
    app_dir=ROSETTA_DIR
)

score_jd2_bin = app.find_app("score_jd2")[0]
