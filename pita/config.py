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

        # Parse YAML config file
        f = open(configfile, "r")
        self.config = yaml.load(f)

        # Data directory
        base = "."
        if config.has_key("data_path"):
            self.base = config["data_path"]

        # Prune overlaps
        self.prune = None
        if config.has_key("prune_overlap"):
            self.prune = config["prune_overlap"]

if not config.has_key("annotation") or len(config["annotation"]) == 0:
    logger.error("No annotation files specified.")
    sys.exit(1)

anno_files = []
chroms = {}
for d in config["annotation"]:
    logger.debug("annotation: {0}".format(d))
    fname = os.path.join(base, d["path"])
    t = d["type"].lower()
    min_exons = 2
    if d.has_key("min_exons"):
        min_exons = d["min_exons"]
    if not t in VALID_TYPES:
        logger.error("Invalid type: {0}".format(t))
        sys.exit(1)
    if not os.path.exists(fname):
        logger.error("File does not exist: {0}".format(fname))
        sys.exit(1)
    else:
        logger.info("Creating tabix index for {0}".format(os.path.basename(fname)))
        logger.debug("Preparing {0} for tabix".format(fname))
        tmp = NamedTemporaryFile(prefix="pita")
        preset = "gff"
        if t == "bed":
            cmd = "sort -k1,1 -k2g,2 {0} | grep -v track | grep -v \"^#\" > {1}"
            preset = "bed"
        elif t in ["gff", "gff3", "gtf3"]:
            cmd = "sort -k1,1 -k4g,4 {0} | grep -v \"^#\" > {1}"
        
        # Sort the input file
        logger.debug(cmd.format(fname, tmp.name))
        subprocess.call(cmd.format(fname, tmp.name), shell=True)
        # Compress using bgzip
        logger.debug("compressing {0}".format(tmp.name))
        tabix_file = tmp.name + ".gz"
        pysam.tabix_compress(tmp.name, tabix_file)
        tmp.close()
        # Index (using tabix command line, as pysam.index results in a Segmentation fault
        logger.debug("indexing {0}".format(tabix_file))
        subprocess.call("tabix {0} -p {1}".format(tabix_file, preset), shell=True)
        
        #fobj = pysam.Tabixfile(tabix_file)
        # Add file info
        anno_files.append([d["name"], tabix_file, t, min_exons])
        # Save chromosome names
        for chrom in pysam.Tabixfile(tabix_file).contigs:
            chroms[chrom] = 1

# data  config
logger.info("Checking data files")
data = []
if config.has_key("data") and config["data"]:
    for d in config["data"]:
        logger.debug("data: {0}".format(d))
        d.setdefault("up", 0)
        d.setdefault("down", 0)
        if type("") == type(d["path"]):
            d["path"] = [d["path"]]
       

        names_and_stats = []
        fnames = [os.path.join(base, x) for x in d["path"]]
        for fname in fnames:
            if not os.path.exists(fname):
                logger.error("File does not exist: {0}".format(fname))
                sys.exit(1)
          
            if fname.endswith("bam") and not os.path.exists(fname + ".bai"):
                logger.error("BAM file {0} needs to be indexed!".format(fname))
                sys.exit(1)

            #if fname.endswith("bam"):
            #    names_and_stats.append((fname, read_statistics(fname)))
            #else:
             #   names_and_stats.append((fname, None))
        row = [d["name"], fnames, d["feature"], (d["up"], d["down"])]
        data.append(row)

weight = {}
if config.has_key("scoring"):
    weight = config["scoring"]

line_format = "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}"
print 'track name="{0}"'.format(basename)
chroms = chroms.keys()

if config.has_key("chromosomes") and config["chromosomes"]:
    if type(config["chromosomes"]) == type([]):
        chroms = config["chromosomes"]
    else:
        chroms = [config["chromosomes"]]

def init_worker():
    signal.signal(signal.SIGINT, signal.SIG_IGN)

def print_output(genename, exons, lock=None):
    print model_to_bed(exons, genename)
    #alog.log_to_file(genename, exons) 

    if lock:
        lock.acquire()
    # Print sequences
    cdna_fh.write(">{0}\n{1}\n".format(genename, exons_to_seq(exons)))
    protein_fh.write(">{0}\n{1}\n".format(genename, longest_orf(exons, do_prot=True)))
    cdna_fh.flush()
    protein_fh.flush()
    if lock:
        lock.release()

def listener(q, lock=None):
    '''listens for messages on the q, writes to file. '''
    
    while 1:
        m = q.get()
        if m == 'kill':
            break
         
        genename, exons = m
        logger.debug("calling print_output for {0}".format(genename))
        print_output(genename, exons, lock)

def annotate_chrom(chrom, q, anno_files, data, weight, prune, index):
    logger.info("Chromosome {0} started".format(chrom))
    for genename, best_exons in get_chrom_models(chrom, anno_files, data, weight, prune, index):
        #results.append([genename, best_exons])
        logger.debug("Putting {0} in print queue".format(genename))
        q.put([genename, best_exons])

    logger.info("Chromosome {0} finished".format(chrom))

if threads > 1:
    logger.info("Starting threaded work")
    manager = mp.Manager()
    lock = manager.Lock()
    q = manager.Queue()
    pool = mp.Pool(threads, init_worker) 
    
    try:  
        
        #put listener to work first
        watcher = pool.apply_async(listener, args=(q, lock) )
        
        # do the main work 
        partialAnnotate = partial(annotate_chrom, q=q, anno_files=anno_files, data=data, weight=weight, prune=prune, index=index)
        pool.map(partialAnnotate, chroms) 
        
        # kill the queue!
        q.put('kill')
        pool.close()
        pool.join()

    except KeyboardInterrupt:
        logger.exception("Caught KeyboardInterrupt, terminating workers")
        pool.terminate()
        pool.join()
else:
    for chrom in chroms:
        for genename,best_exons in get_chrom_models(chrom, anno_files, data, weight, prune, index):
            print_output(genename,best_exons)

cdna_fh.close()
protein_fh.close()
