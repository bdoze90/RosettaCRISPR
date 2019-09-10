"""File for generating all the data for Cas12 structures.  Call this file in the command line for a single ensemble,
specifying which endonuclease you are using so that it pulls from the correct files in terms of rnaid/dnaid.
Calling multiple instances of this file is a bad idea, since the processes are already sufficiently parallelized
to maximize usage of a single node."""

from tempfile import mkstemp
import shutil
import os
from RosettaSub import RosettaSingleProcess, RosettaSubprocess
from PDBparse import PDB

class SeqMutatorCas12:

    def __init__(self,CasID,num_ensemble,base_dir,structure):
        self.CasID = CasID
        self.Ensemble = num_ensemble
        self.base_dir = base_dir + structure + "/"
        self.structure = structure

        self.rna_seqs = list()  # in order list of sequences from 1 to n
        self.dna_seqs = list()
        self.combinations = dict()  # relational database with keys as rna sequences and values as dna sequences
        self.on_target_combos = list()  # list of tuples containing the rna and dna combinations that comprise the on-target seqs
        # List of all the appropriate indexes for the mutated sequences and crystal structures
        self.cs_dict = {"5XUU": {"ChainB": ('r', 0, 20, '', 'AAUUUCUACUAAGUGUAGAU'),  # LbCas12
                                 "ChainA": ('protein', 'n', 'n', 'n', 'n'),
                                 "ChainC": ('d', 0, 20, 'rc', 'TGGAGGACG'),
                                 "ChainD": ('d', 'n', 'n', '', 'CGTCCTCCA')},

                        "5XUS": {"ChainB": ('r', 0, 20, '', 'AAUUUCUACUAAGUGUAGAU'),  # LbCas12
                                 "ChainA": ('protein', 'n', 'n', 'n', 'n'),
                                 "ChainC": ('d', 0, 20, 'rc', 'TAAAGGACG'),
                                 "ChainD": ('d', 'n', 'n', '', 'CGTCCTTTA')},

                        "5XUT": {"ChainB": ('r', 0, 20, '', 'AAUUUCUACUAAGUGUAGAU'),  # LbCas12
                                 "ChainA": ('protein', 'n', 'n', 'n', 'n'),
                                 "ChainC": ('d', 0, 20, 'rc', 'TAGAGGACG'),
                                 "ChainD": ('d', 'n', 'n', '', 'CGTCCTCTA')},

                        "5XH6": {"ChainB": ('r', 0, 'n', '', 'AAUUUCUACUCUUGUAGAU'),  # AsCas12
                                 "ChainA": ('protein', 'n', 'n', 'n', 'n'),
                                 "ChainC": ('d', 0, 'n', 'rc', 'TATAGGACTG'),
                                 "ChainD": ('d', 'n', 'n', '', 'CAGTCCTATA')},

                        "5XH7": {"ChainB": ('r', 0, 'n', '', 'AAUUUCUACUCUUGUAGAU'),  # AsCas12
                                 "ChainA": ('protein', 'n', 'n', 'n', 'n'),
                                 "ChainC": ('d', 0, 'n', 'rc', 'TGGAGGACTG'),
                                 "ChainD": ('d', 'n', 'n', '', 'CAGTCCTCCA')}
                        }

    # STEP 1: Collect the data from the raw files and place into objects.
    def collect_data(self):
        # Get the rna sequences
        f = open("/home/trinhlab/Desktop/RosettaCRISPR/rna_seqs_" + self.CasID + ".txt")
        for line in f:
            self.rna_seqs.append(line.split("\t")[1][-1])
        f.close()
        # Get the dna sequences
        y = open("/home/trinhlab/Desktop/RosettaCRISPR/dna_seqs_" + self.CasID + ".txt")
        for line in y:
            self.dna_seqs.append(line.split("\t")[1][-1])
        y.close()
        # Get the combinations
        z = open("/home/trinhlab/Desktop/RosettaCRISPR/combo_seqs_" + self.CasID + ".txt")
        for line in z:
            linelist = line.split("\t")
            if linelist[0] == "ON":
                self.on_target_combos.append((linelist[1],linelist[2]))
            else:
                if linelist[1] in self.combinations:
                    self.combinations[linelist[1]].append(linelist[2])
                else:
                    self.combinations[linelist[1]] = [linelist[2]]

    # STEP 2: Mutate the RNAs and DNAs from the base structure into the on targets sequences
    def mutate_to_on_target(self):
        # Run the Rosetta Subprocess for each sequence:
        rr = RosettaSingleProcess("rna_thread.default.linuxgccrelease")

        # Create the Base PDB object for extracting chains:
        base_pdb = PDB(self.base_dir, self.structure + ".pdb")

        # Get then mutate Chain D since it is always constant for all structures
        chainD = open(self.base_dir + "ChainD.pdb", 'w')
        chainD.write(base_pdb.return_chain("D"))
        rr.set_inputs(
            ["-s", self.base_dir + "ChainD.pdb", "-seq", self.cs_dict[self.structure]["ChainD"][4].lower(), "-o",
             self.base_dir + "ChainD_MUT.pdb"])
        rr.run_process()
        # Change the chain name for the file:
        self.change_chain_name(self.base_dir + "ChainD.pdb", "D")

        # Get the rest of the base chains extracted:
        chainB = open(self.base_dir + "ChainB.pdb",'w')
        chainB.write(base_pdb.return_chain("B"))
        chainB.close()
        chainC = open(self.base_dir + "ChainC.pdb","w")
        chainC.write(base_pdb.return_chain("C"))
        chainC.close()
        chainA = open(self.base_dir + "ChainA.pdb", 'w')
        chainA.write(base_pdb.return_chain("A"))
        chainA.close()

        # Mutate all the RNA sequences in the on-target-combos:
        for ids in self.on_target_combos:
            specs = self.cs_dict[self.structure]["ChainB"]
            if specs[2] == 'n':
                mysequence = self.rna_seqs[ids[0]][specs[1]:] + specs[4]
            else:
                mysequence = self.rna_seqs[ids[0]][specs[1]:specs[2]] + specs[4]
            # check to see if you need to run the sequence through revcom algorithm:
            if specs[3] == 'rc':
                mysequence = self.revcom(mysequence)
            rr.set_inputs(
                ["-s", self.base_dir + "ChainB.pdb", "-seq", mysequence.lower(), "-o",
                 self.base_dir + "ChainB_MUT/"+ ids[0] + ".pdb"])
            rr.run_process()

            # Change the chain name for the file:
            self.change_chain_name(self.base_dir + "ChainB_MUT/"+ ids[0] + ".pdb", "B")

        # Mutate the DNA sequences in the on-target-combos:
        for ids in self.on_target_combos:
            specs = self.cs_dict[self.structure]["ChainC"]
            if specs[2] == 'n':
                mysequence = self.dna_seqs[ids[1]][specs[1]:] + specs[4]
            else:
                mysequence = self.dna_seqs[ids[1]][specs[1]:specs[2]] + specs[4]
            # check to see if you need to run the sequence through revcom algorithm:
            if specs[3] == 'rc':
                mysequence = self.revcom(mysequence)
            rr.set_inputs(
                ["-s", self.base_dir + "ChainC.pdb", "-seq", mysequence.lower(), "-o",
                 self.base_dir + "ChainC_MUT/" + ids[1] + ".pdb"])
            rr.run_process()
            # Change the chain name for the file:
            self.change_chain_name(self.base_dir + "ChainC_MUT/" + ids[1] + ".pdb", "C")

        # Consolidating the mutated structures into full mutated PDB files:
        for combo in self.on_target_combos:
            a = self.base_dir + "ChainA.pdb"
            chains = [a]
            chains.append(self.base_dir + "ChainB_MUT/" + combo[0] + ".pdb")
            chains.append(self.base_dir + "ChainC_MUT/" + combo[1] + ".pdb")
            chains.append(self.base_dir + "ChainD_MUT.pdb")
            with open(self.base_dir + "FULL_MUT_PDBs/" + combo[0] + "-" + combo[1] + ".pdb", 'wb') as wfd:
                for f in chains:
                    with open(f, 'rb') as fd:
                        shutil.copyfileobj(fd, wfd)


    # STEP 3: Relax the on target structures, then add the _rel suffix.  Call this function again with the minimize
    def create_on_targets(self, process):
        runlist = list()
        os.chdir(self.base_dir + "FULL_MUT_PDBs/")
        for target in os.listdir(os.curdir):
            # Make sure not to open a .DS_Store
            if target.endswith(".pdb"):
                # Check for the minimization successive function:
                if process.startswith("minimize"):
                    if target.endswith("rel.pdb"):
                        runlist.append(target)
                else:
                    runlist.append(target)
        rsub = RosettaSubprocess(process, 16, runlist)
        if process.startswith("relax"):
            rsub.set_inputs(["-s", 'filler', "nstruct", "1", "relax:default_repeats", "5", "-out:suffix", "_rel"])
        else:
            rsub.set_inputs(["-s", 'filler', "-out:suffix", "min", "-run:min_tolerance", "0.0001"])
        rsub.run_batch()


    # STEP 4: Using the relaxed structure from the fully mutated structures, mutate all the dna sequences in the list,
    # while also creating the directories for all the off-targets in the OFF_TARGET folder.
    def off_target_structure_mutate(self):
        # Create the directories for storing the mutated files
        for rseq in self.rna_seqs:
            os.mkdir(self.base_dir + "OFF_TARGET/" + rseq)
        # Obtain the sequence from ChainC (The DNA that needs to be changed) for the mutation
        for targetfile in os.listdir(self.base_dir + "FULL_MUT_PDBs/"):
            rnaid = targetfile[0:8]
            target_pdb = PDB(self.base_dir + "FULL_MUT_PDBs/", targetfile)
            DNAchain = open(self.base_dir + "FULL_MUT_PDBs/" + "tempChainC.pdb")
            DNAchain.write(target_pdb.return_chain("C"))
            DNAchain.close()

            # Initialize the Rosetta mutator algorithm:
            rr = RosettaSingleProcess("rna_thread.default.linuxgccrelease")

            # Mutate the DNA to the respective off target sequences:
            for dnaid in self.combinations[rnaid]:
                specs = self.cs_dict[self.structure]["ChainC"]
                if specs[2] == 'n':
                    mysequence = self.dna_seqs[dnaid][specs[1]:] + specs[4]
                else:
                    mysequence = self.dna_seqs[dnaid][specs[1]:specs[2]] + specs[4]
                # check to see if you need to run the sequence through revcom algorithm:
                if specs[3] == 'rc':
                    mysequence = self.revcom(mysequence)
                rr.set_inputs(
                    ["-s", self.base_dir + "FULL_MUT_PDBs/tempChainC.pdb", "-seq", mysequence.lower(), "-o",
                     self.base_dir + "ChainC_MUT/" + dnaid + ".pdb"])
                rr.run_process()
                # Change the chain name for the file:
                self.change_chain_name(self.base_dir + "ChainC_MUT/" + dnaid + ".pdb", "C")

                # Reassemble the file
                f = open(self.base_dir + "OFF_TARGET/" + rnaid + "/" + rnaid + "-" + dnaid + ".pdb",'w')
                f.write(target_pdb.return_chain("A"))
                f.write(target_pdb.return_chain("B"))
                with open(self.base_dir + "FULL_MUT_PDBs/tempChainC.pdb", 'rb') as fd:
                    shutil.copyfileobj(fd, f)
                f.write(target_pdb.return_chain("D"))


    # STEP 5: Minimize all the off target structures across all the directories of off targets
    def minimize_off_target(self):
        runlist = list()
        os.chdir(self.base_dir + "OFF_TARGET/")
        # Iterate across all the off target groups
        for onbase in os.listdir(os.curdir):
            os.chdir(onbase)
            for offstruct in os.listdir(os.curdir):
                # Make sure not to open a .DS_Store
                if offstruct.endswith(".pdb"):
                    runlist.append(offstruct)
            rsub = RosettaSubprocess("minimize.default.linuxgccrelease", 16, runlist)
            rsub.set_inputs(["-s", 'filler', "-out:suffix", "_min", "-run:min_tolerance", "0.0001"])
            rsub.run_batch()


    # STEP 6a: Generate truncations by removing the nucleotides from the back end of the ChainB (RNA) sequence

    # STEP 6b: Generate truncations by removing the DNA nucleotides from the frond end of the ChainC (DNA) sequence


    #### START OF EXTRA FUNCTIONS NEEDED FOR THE PROGRAM TO RUN ####
    def revcom(self, sequence, complement=False):
        retseq = ""
        change = {'A': 'T',
                  'T': 'A',
                  'G': 'C',
                  'C': 'G'}
        for nt in sequence:
            rnt = change[nt.upper()]
            if complement:
                retseq += rnt
            else:
                retseq = rnt + retseq
        return retseq

    # Changes the chain name to be consistent with the new PDB file:
    def change_chain_name(self, change_file, new_chain_char):
        # Create temp file
        print("changing file: " + change_file)
        fh, abs_path = mkstemp()
        with os.fdopen(fh, 'w') as new_file:
            with open(change_file) as old_file:
                for line in old_file:
                    if line.find("TER") != -1 or line.find("HET") != -1:
                        new_file.write(line)
                    else:
                        new_line = line[:21] + new_chain_char + line[22:]  # chain char at position 21
                        new_file.write(new_line)
            # Remove original file
            os.remove(change_file)
            # Move new file
            shutil.move(abs_path, change_file)




SMC12 = SeqMutatorCas12("lbCas12","1","/home/trinhlab/Desktop/RosettaCRISPR/","5XUS")
SMC12.collect_data()  # step 1
SMC12.mutate_to_on_target() # step 2
SMC12.create_on_targets("relax.default.linuxgccrelease") # step 3a
SMC12.create_on_targets("minimize.default.linuxgccrelease") # step 3b
SMC12.off_target_structure_mutate() # step 4
SMC12.minimize_off_target() # step 5
# SMC12.generate_truncations("B", False)
# SMC12.generate_truncations("C", True)