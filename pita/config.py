import logging
import yaml
import os
import sys
import pysam
import subprocess
from tempfile import NamedTemporaryFile

SAMTOOLS = "samtools"
TSS_FOUND = "v"
TSS_UPSTREAM = "u"
TSS_DOWNSTREAM = "a"
TSS_NOTFOUND = "x"
SEP = ":::"
VALID_TYPES = ["bed", "gff", "gff3", "gtf"]
DEBUG_LEVELS = ["CRITICAL", "ERROR", "WARNING", "INFO", "DEBUG"]

class PitaConfig:
    def __init__(self, fname):
        """ fname: name of yaml configuration file
        """
        
        self.logger = logging.getLogger("pita")

        # Parse YAML config file
        f = open(fname, "r")
        self.config = yaml.load(f)
        f.close()

        self.db_conn = "sqlite:///pita_database.db"
        if self.config.has_key("database"):
            self.db_conn = self.config["database"]
        
        # Data directory
        self.base = "."
        if self.config.has_key("data_path"):
            self.base = self.config["data_path"]

        # Prune overlaps
        self.prune = None
        if self.config.has_key("prune_overlap"):
            self.prune = self.config["prune_overlap"]
        
        # Scoring weight
        self.weight = {}
        if self.config.has_key("scoring"):
            self.weight = self.config["scoring"]

        # load annotation files
        self._parse_annotation()
  
        # only use chromosome specified in config file
        self.chroms = self.chroms.keys()
        if self.config.has_key("chromosomes") and self.config["chromosomes"]:
            if type(self.config["chromosomes"]) == type([]):
                self.chroms = self.config["chromosomes"]
            else:
                self.chroms = [self.config["chromosomes"]]

        # check the data files
        self._check_data_files()

    def _parse_annotation(self):

        if not self.config.has_key("annotation") or len(self.config["annotation"]) == 0:
            self.logger.error("No annotation files specified.")
            sys.exit(1)
        
        self.anno_files = []
        self.chroms = {}
        for d in self.config["annotation"]:
            self.logger.debug("annotation: {0}".format(d))
            fname = os.path.join(self.base, d["path"])
            t = d["type"].lower()
            min_exons = 2
            if d.has_key("min_exons"):
                min_exons = d["min_exons"]
            if not t in VALID_TYPES:
                self.logger.error("Invalid type: {0}".format(t))
                sys.exit(1)
            if not os.path.exists(fname):
                self.logger.error("File does not exist: {0}".format(fname))
                sys.exit(1)
            else:
                self.logger.info("Creating tabix index for {0}".format(os.path.basename(fname)))
                self.logger.debug("Preparing {0} for tabix".format(fname))
                tmp = NamedTemporaryFile(prefix="pita")
                preset = "gff"
                if t == "bed":
                    cmd = "sort -k1,1 -k2g,2 {0} | grep -v track | grep -v \"^#\" > {1}"
                    preset = "bed"
                elif t in ["gff", "gff3", "gtf3"]:
                    cmd = "sort -k1,1 -k4g,4 {0} | grep -v \"^#\" > {1}"
                
                # Sort the input file
                self.logger.debug(cmd.format(fname, tmp.name))
                subprocess.call(cmd.format(fname, tmp.name), shell=True)
                # Compress using bgzip
                self.logger.debug("compressing {0}".format(tmp.name))
                tabix_file = tmp.name + ".gz"
                pysam.tabix_compress(tmp.name, tabix_file)
                tmp.close()
                # Index (using tabix command line, as pysam.index results in a Segmentation fault
                self.logger.debug("indexing {0}".format(tabix_file))
                subprocess.call("tabix {0} -p {1}".format(tabix_file, preset), shell=True)
                
                #fobj = pysam.Tabixfile(tabix_file)
                # Add file info
                self.anno_files.append([d["name"], fname, tabix_file, t, min_exons])
                # Save chromosome names
                for chrom in pysam.Tabixfile(tabix_file).contigs:
                    self.chroms[chrom] = 1

    def _check_data_files(self):
        # data  config
        self.logger.info("Checking data files")
        self.data = []
        if self.config.has_key("data") and self.config["data"]:
            for d in self.config["data"]:
                self.logger.debug("data: {0}".format(d))
                d.setdefault("up", 0)
                d.setdefault("down", 0)
                if type("") == type(d["path"]):
                    d["path"] = [d["path"]]
       
                d.setdefault("feature", "all")
                if d["feature"] not in ["all", "start", "end", "splice"]:
                    self.logger.error("Incorrect span: {}".format(d["feature"]))
                    sys.exit(1)

                names_and_stats = []
                fnames = [os.path.join(self.base, x) for x in d["path"]]
                for fname in fnames:
                    if not os.path.exists(fname):
                        self.logger.error("File does not exist: {0}".format(fname))
                        sys.exit(1)
                  
                    if fname.endswith("bam") and not os.path.exists(fname + ".bai"):
                        self.logger.error("BAM file {0} needs to be indexed!".format(fname))
                        sys.exit(1)

                    #if fname.endswith("bam"):
                    #    names_and_stats.append((fname, read_statistics(fname)))
                    #else:
                     #   names_and_stats.append((fname, None))
                row = [d["name"], fnames, d["feature"], (int(d["up"]), int(d["down"]))]
                self.data.append(row)
