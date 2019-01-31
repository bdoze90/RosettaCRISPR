"""This allows for a fully automated processing all the way from initial ensemble of 5 structures to the truncation and
RNA/DNA data output in raw txt data file for processing in either Analysis code or Excel."""

# Pseudocode is in comments, fill in the code as you go along refining the protocol from Ensemble_1 of 4UN3
# This code is at the PDB structure level.  Each new PDB will need to go through this process.
# Manual creation of the directories with the initial relaxed structure

# Run the SeqMutator and then the SeqMutatorOff function to create all the off-target sequences for each ensemble
# Run the