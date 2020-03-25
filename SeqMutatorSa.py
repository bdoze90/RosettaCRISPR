"""This file runs through the complete mutations for SaCas9.  Most of this code is degenerate from the spCas9 Seq Mutator
and can eventually be combined into one file, but for now with such a large dataset it is easier to just keep the SaCas9
as a separate file even though the data structures as far as DNA/RNAIDs are the same."""

from tempfile import mkstemp
import shutil
import os
from RosettaSub import RosettaSingleProcess, RosettaSubprocess
from PDBparse import PDB

class SeqMutatorSa:

    def __init__(self,CasID,num_ensemble,base_dir,structure):
        self.CasID = CasID
        self.Ensemble = num_ensemble
        self.base_dir = base_dir + structure + "/Ensemble_" + str(num_ensemble) + "/"
        self.structure = structure

        self.rna_seqs = list()  # in order list of sequences from 1 to n
        self.dna_seqs = list()
        self.combinations = dict()  # relational database with keys as rna sequences and values as dna sequences
        self.on_target_combos = list()  # list of tuples containing the rna and dna combinations that comprise the on-target seqs
        # List of all the appropriate indexes for the mutated sequences and crystal structures
        self.cs_dict = {# Beginning of the saCas9 structures (25bp target sequence for baseline)
                        "5CZZ": {"ChainB": ('r', 0, -1, '', ''),
                                 "ChainA": ('protein', 'n', 'n', 'n', 'n'),
                                 "ChainC": ('d', 1, 25, 'rc', 'CTATTCAA'),
                                 "ChainD": ('d', 'n', 'n', '', 'TTGAATAG')},

                        "5AXW": {"ChainB": ('r', 0, -1, '', ''),
                                 "ChainA": ('protein', 'n', 'n', 'n', 'n'),
                                 "ChainC": ('d', 1, 26, 'rc', 'TTGGGTAG'),
                                 "ChainD": ('d', 0, 'n', '', 'TTGGGTAG')}
                        }

    # STEP 1: Collect the data from the raw files and place into objects.
    def collect_data(self):
        # Get the rna sequences
        f = open("/home/trinhlab/Desktop/RosettaCRISPR/rna_seqs_" + self.CasID + ".txt")
        for line in f:
            self.rna_seqs.append(line[:-1].split("\t")[1])
        f.close()
        # Get the dna sequences
        y = open("/home/trinhlab/Desktop/RosettaCRISPR/dna_seqs_" + self.CasID + ".txt")
        for line in y:
            self.dna_seqs.append(line[:-1].split("\t")[1])
        y.close()
        # Get the combinations
        z = open("/home/trinhlab/Desktop/RosettaCRISPR/combo_seqs_" + self.CasID + ".txt")
        for line in z:
            linelist = line[:-1].split("\t")
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
        self.change_chain_name(self.base_dir + "ChainD_MUT.pdb", "D")

        # Get the rest of the base chains extracted:
        chainB = open(self.base_dir + "ChainB.pdb",'w')
        chainB.write(base_pdb.return_chain("B"))
        chainB.close()
        chainA = open(self.base_dir + "ChainA.pdb", 'w')
        chainA.write(base_pdb.return_chain("A"))
        chainA.close()
        chainD = open(self.base_dir + "ChainC.pdb", 'w')
        chainD.write(base_pdb.return_chain("C"))
        chainD.close()

        # Mutate all the RNA sequences in the on-target-combos:
        for ids in self.on_target_combos:
            rnumid = int(ids[0][2:]) - 300001
            specs = self.cs_dict[self.structure]["ChainB"]
            if specs[2] == 'n':
                mysequence = self.rna_seqs[rnumid][specs[1]:] + specs[4]
            else:
                mysequence = self.rna_seqs[rnumid][specs[1]:specs[2]] + specs[4]
            # check to see if you need to run the sequence through revcom algorithm:
            if specs[3] == 'rc':
                mysequence = self.revcom(mysequence)
            print(len(mysequence))
            rr.set_inputs(
                ["-s", self.base_dir + "ChainB.pdb", "-seq", mysequence.lower(), "-o",
                 self.base_dir + "ChainB_MUT/"+ ids[0] + ".pdb"])
            rr.run_process()

            # Change the chain name for the file:
            self.change_chain_name(self.base_dir + "ChainB_MUT/"+ ids[0] + ".pdb", "A")

        # Mutate the DNA sequences in the on-target-combos:
        for ids in self.on_target_combos:
            dnumid = int(ids[1][2:]) -300001
            specs = self.cs_dict[self.structure]["ChainC"]
            if specs[2] == 'n':
                mysequence = self.dna_seqs[dnumid][specs[1]:] + specs[4]
            else:
                mysequence = self.dna_seqs[dnumid][specs[1]:specs[2]] + specs[4]
            # check to see if you need to run the sequence through revcom algorithm:
            if specs[3] == 'rc':
                mysequence = self.revcom(mysequence)
            print(len(mysequence))
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
        # remove the numbering:
        for target in os.listdir(os.curdir):
            if target.endswith("0001.pdb"):
                os.rename(os.getcwd() + "/" + target, os.getcwd() + "/" + target[:-9] + ".pdb")


    # STEP 4: Using the relaxed structure from the fully mutated structures, mutate all the dna sequences in the list,
    # while also creating the directories for all the off-targets in the OFF_TARGET folder.
    def off_target_structure_mutate(self):
        # Create the directories for storing the mutated files
        for i in range(len(self.rna_seqs)):
            if not os.path.isdir(self.base_dir + "OFF_TARGET/r_" + str(300001+ i)):
                os.mkdir(self.base_dir + "OFF_TARGET/r_" + str(300001 + i))
        # Obtain the sequence from ChainC (The DNA that needs to be changed) for the mutation
        for rnaid in self.combinations:
            #search for the corresponding DNAid
            cordna = str()
            for item in self.on_target_combos:
                if item[0] == rnaid:
                    cordna = item[1]
            targetfile = rnaid + "-" + cordna + "_rel.pdb"  # change to _min if not doing off-relaxation
            target_pdb = PDB(self.base_dir + "FULL_MUT_PDBs/", targetfile)
            DNAchain = open(self.base_dir + "FULL_MUT_PDBs/" + "tempChainC.pdb", "w")
            DNAchain.write(target_pdb.return_chain("C"))
            DNAchain.close()
            #DNAchain2 = open(self.base_dir + "FULL_MUT_PDBs/" + "tempChainD.pdb", "w")
            #DNAchain2.write(target_pdb.return_chain("D"))
            #DNAchain2.close()

            # Initialize the Rosetta mutator algorithm:
            rr = RosettaSingleProcess("rna_thread.default.linuxgccrelease")

            # Mutate the DNA to the respective off target sequences:
            for dnaid in self.combinations[rnaid]:
                dnumid = int(dnaid[2:]) - 300001
                specs = self.cs_dict[self.structure]["ChainC"]
                if specs[2] == 'n':
                    mysequence = self.dna_seqs[dnumid][specs[1]:] + specs[4]
                else:
                    mysequence = self.dna_seqs[dnumid][specs[1]:specs[2]] + specs[4]
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
                f = open(self.base_dir + "OFF_TARGET/" + rnaid + "/" + rnaid + "-" + dnaid + ".pdb", 'w')
                f.write(target_pdb.return_chain("A"))
                f.write(target_pdb.return_chain("B"))
                with open(self.base_dir + "FULL_MUT_PDBs/tempChainC.pdb", 'r') as fd:
                    content = fd.read()
                    f.write(content)
                f.write(target_pdb.return_chain("D"))


    # STEP 5: Minimize all the off target structures across all the directories of off targets
    def minimize_off_target(self):
        odir = self.base_dir + "OFF_TARGET/"
        # Iterate across all the off target groups
        for onbase in os.listdir(odir):
            print(onbase)
            runlist = []
            os.chdir(odir + onbase)
            for offstruct in os.listdir(os.curdir):
                # Make sure not to open a .DS_Store
                if offstruct.endswith(".pdb"):
                    runlist.append(offstruct)
            rsub = RosettaSubprocess("minimize.default.linuxgccrelease", 16, runlist)
            rsub.set_inputs(["-s", 'filler', "-out:suffix", "_min", "-run:min_tolerance", "0.0001"])
            rsub.run_batch()
            # Delete the non-minimized files
            for file in os.listdir(os.curdir):
                if file.endswith("min_0001.pdb"):
                    os.rename(file,file[:-9] + ".pdb")


    # STEP 6a: Generate truncations by removing the nucleotides from either the RNA or DNA sequence
    def generate_truncations(self, chain, direction):
        b = self.base_dir + "OFF_TARGET/"
        # Iterate through all off target groups:
        for onbase in os.listdir(b):
            if onbase.startswith("r"):
                os.chdir(b + onbase)
                for target in os.listdir(os.curdir):
                    if target.endswith("min.pdb"):
                        myp = PDB(os.getcwd() + "/", target)
                        myp.reassemble(target[:-7],24, chain, trunc_dir=direction)

    # STEP 6b: Generate truncations of the base files (relaxed not minimized) in the FullMut folder
    def generate_base_truncations(self, chain, direction):
        b = self.base_dir + "FULL_MUT_PDBs"
        # Iterate through full muts:
        for file in os.listdir(b):
            if file.endswith("rel.pdb"):
                myp = PDB(b + "/", file)
                myp.reassemble(file[:-7], 24, chain, trunc_dir=direction)

    # STEP 7: Score all the truncation files
    def score_truncations(self):
        # Iterate through all of the off target groups accessing the truncation folder:
        for onbase in os.listdir(self.base_dir + "OFF_TARGET/"):
            os.chdir(self.base_dir + "OFF_TARGET/")
            if onbase.startswith("r"):  # Making sure to only pull in the actual directories
                mylist = list()
                if not os.path.isdir(onbase + "/truncs_from_min"):
                    continue
                os.chdir(onbase + "/truncs_from_min")
                for item in os.listdir(os.curdir):
                    if item.endswith(".pdb"):
                        print(item)
                        mylist.append(os.getcwd() + "/" + item)
                R = RosettaSubprocess("score_jd2.default.linuxgccrelease", 16, mylist)
                R.set_inputs(["-in:file:s",'filler', "-out:pdb", "-out:suffix", "_scr"])
                R.run_batch(sleeptime=15)
            # Delete the non-scored truncation files so they don't cloud memory:
            for file in os.listdir(os.curdir):
                if file.endswith("scr_0001.pdb"):
                    os.rename(file,file[:-9] + ".pdb")
            for file in os.listdir(os.curdir):
                if file.endswith(".pdb") and not file.endswith("scr.pdb"):
                    os.remove(file)

    def score_truncations_base(self):
        os.chdir(self.base_dir + "FULL_MUT_PDBs/truncs_from_min")
        mylist = list()
        for item in os.listdir(os.curdir):
            if item.endswith(".pdb"):
                print(item)
                mylist.append(os.getcwd() + "/" + item)
        R = RosettaSubprocess("score_jd2.default.linuxgccrelease", 16, mylist)
        R.set_inputs(["-in:file:s", 'filler', "-out:pdb", "-out:suffix", "_scr"])
        R.run_batch(sleeptime=15)
        # Delete the non-scored truncation files so they don't cloud memory:
        for file in os.listdir(os.curdir):
            if file.endswith("scr_0001.pdb"):
                os.rename(file, file[:-9] + ".pdb")
        for file in os.listdir(os.curdir):
            if file.endswith(".pdb") and not file.endswith("scr.pdb"):
                os.remove(file)


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



SMCa9 = SeqMutatorSa("saCas9","1","/home/trinhlab/Desktop/RosettaCRISPR/","5CZZ")
SMCa9.collect_data()  # step 1
#SMCa9.mutate_to_on_target() # step 2
#SMCa9.create_on_targets("relax.default.linuxgccrelease") # step 3a
#SMCa9.create_on_targets("minimize.default.linuxgccrelease") # step 3b
SMCa9.off_target_structure_mutate() # step 4
SMCa9.minimize_off_target() # step 5
SMCa9.generate_truncations("B", False) # step 6a
#SMC12.generate_truncations("C", True)  # step 6b
SMCa9.score_truncations()  # step 7