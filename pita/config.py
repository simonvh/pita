import logging
import yaml
import os
import sys
import pysam
from pita.io import _create_tabix

SAMTOOLS = "samtools"
TSS_FOUND = "v"
TSS_UPSTREAM = "u"
TSS_DOWNSTREAM = "a"
TSS_NOTFOUND = "x"
SEP = ":::"
VALID_TYPES = ["bed", "gff", "gff3", "gtf"]
DEBUG_LEVELS = ["CRITICAL", "ERROR", "WARNING", "INFO", "DEBUG"]

class PitaConfig(object):
    def __init__(self):
        """ fname: name of yaml configuration file
        """
        self.logger = logging.getLogger("pita")
        
        self.maxentpath = ""

    def load(self, fname,  reannotate=False):
        # Parse YAML config file
        f = open(fname, "r")
        self.config = yaml.load(f)
        f.close()

        self.db_conn = "sqlite:///pita_database.db"
        if "database" in self.config:
            self.db_conn = self.config["database"]
        
        # Data directory
        self.base = "."
        if "data_path" in self.config:
            self.base = self.config["data_path"]

        # Prune overlaps
        self.prune = None
        if "prune" in self.config:
            self.prune = self.config["prune"]
        
        self.keep = []
        self.filter = []
        self.experimental = []
    
        # MaxEnt directory
        self.maxentpath = "" 
        if "maxent" in self.config:
            self.maxentpath = self.config["maxent"]
            
        # Pita UTR
        self.pitaUTR = False
        if "extendUtr" in self.config:
            self.pitaUTR = True

        # Scoring weight
        self.weight = {}
        if "scoring" in self.config:
            self.weight = self.config["scoring"]

        self._parse_repeats()
       
        # load annotation files
        self._parse_annotation(reannotate)
  
        # only use chromosome specified in config file
        self.chroms = self.chroms.keys()
        if "chromosomes" in self.config and self.config["chromosomes"]:
            if type(self.config["chromosomes"]) == type([]):
                self.chroms = self.config["chromosomes"]
            else:
                self.chroms = [self.config["chromosomes"]]

        # check the data files
        self._check_data_files()
        
        # output option
        self.min_protein_size = 20 
    
    def _parse_repeats(self):
        self.repeats = []
        if "repeats"in self.config:
            for d in self.config["repeats"]:
                fname = os.path.join(self.base, d["path"])
                tabix_fname = _create_tabix(fname, "bed")
                d["path"] = fname
                d["tabix"] = tabix_fname
 
                self.repeats.append(d)

    def _parse_annotation(self, reannotate=False):

        if not self.config.has_key("annotation") or len(self.config["annotation"]) == 0:
            self.logger.error("No annotation files specified.")
            sys.exit(1)
        
        self.anno_files = []
        self.chroms = {}
        for d in self.config["annotation"]:
            self.logger.debug("annotation: %s", d)
            fname = os.path.join(self.base, d["path"])
            t = d["type"].lower()
            min_exons = 2
            if d.has_key("min_exons"):
                min_exons = d["min_exons"]
            
            if d.has_key("keep") and d["keep"]:
                self.keep.append(fname)
           
            if d.has_key("filter") and d["filter"]:
                self.filter.append(fname)
            
            if d.has_key("experimental") and d["experimental"]:
                self.experimental.append(fname)

            if not t in VALID_TYPES:
                self.logger.error("Invalid type: %s", t)
                sys.exit(1)
            if not os.path.exists(fname):
                self.logger.error("File does not exist: %s", fname)
                sys.exit(1)
            else:
                tabix_file = ""
                if not reannotate:
                    tabix_file = _create_tabix(fname, t)
                    
                    # Save chromosome names
                    for chrom in pysam.Tabixfile(tabix_file).contigs:
                        self.chroms[chrom] = 1

                # Add file info
                self.anno_files.append([d["name"], fname, tabix_file, t, min_exons])
       
    def _check_data_files(self):
        # data  config
        self.logger.info("Checking data files")
        self.data = []
        if self.config.has_key("data") and self.config["data"]:
            for d in self.config["data"]:
                self.logger.debug("data: %s", d)
                d.setdefault("up", 0)
                d.setdefault("down", 0)
                if type("") == type(d["path"]):
                    d["path"] = [d["path"]]
       
                d.setdefault("feature", "all")
                if d["feature"] not in ["all", "start", "end", "splice"]:
                    self.logger.error("Incorrect span: %s", ["feature"])
                    sys.exit(1)

                fnames = [os.path.join(self.base, x) for x in d["path"]]
                for fname in fnames:
                    if not os.path.exists(fname):
                        self.logger.error("File does not exist: %s", fname)
                        sys.exit(1)
                  
                    if fname.endswith("bam") and not os.path.exists(fname + ".bai"):
                        self.logger.error("BAM file %s needs to be indexed!", fname)
                        sys.exit(1)

                    #if fname.endswith("bam"):
                    #    names_and_stats.append((fname, read_statistics(fname)))
                    #else:
                     #   names_and_stats.append((fname, None))
                row = [d["name"], fnames, d["feature"], (int(d["up"]), int(d["down"]))]
                self.data.append(row)

config = PitaConfig()
