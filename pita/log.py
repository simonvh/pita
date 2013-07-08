from pita.collection import get_updated_exons
#import itertools
import logging

class AnnotationLog:
    header = "Model\tNr. Exons\tExons in best model\tExons in other models\tUpdated 5'\tUpdated 3'\tOriginal models\n"
    log_str = "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\n"
    def __init__(self, append):
        self.files = {}
        self.logger = logging.getLogger('pita')
        self.append = append 

    def add(self, name):
        mode = "w"
        if self.append:
            mode = "a"
        self.files[name] = open("pita.{0}.log".format(name), mode)
        self.files[name].write(self.header)
         
    def log_to_file(self, genename, model, ev, best_ev, other_ev):
#        best_exons = [e for e in model]
#        best_ev = {}
#        for e in best_exons:
#            for ev in set([x.split(":")[0] for x in e.evidence]):
#                best_ev[ev] = best_ev.setdefault(ev, 0) + 1
#        other_exons = []
#        other_ev = {}
#        #self.logger.debug("cluster {0}".format(len(cluster)))
#        
#        # Fast way to collapse
#        other_exons = [e for e in set(itertools.chain.from_iterable(cluster)) if not e in best_exons]
#        for e in other_exons:
#            for ev in set([x.split(":")[0] for x in e.evidence]):
#                other_ev[ev] = other_ev.setdefault(ev, 0) + 1
#    
#        ev = []
#        for e in best_exons + other_exons:
#            for evidence in e.evidence:
#                ev.append(evidence.split(":"))
       
        # ev, model, best_exons
        for name, f in self.files.items():
            orig_models = {}
            for (origin,orig_name) in ev:
                if origin == name:
                    orig_models[orig_name] = orig_models.setdefault(orig_name, 0) + 1
    
            u5, u3 = get_updated_exons(model, name)
            f.write(self.log_str.format(
                                       genename,
                                       len(model),
                                       best_ev.setdefault(name, 0),
                                       other_ev.setdefault(name, 0),
                                       u5,
                                       u3,
                                       ",".join(orig_models.keys())
                                      )
                        )
            f.flush()       

    def __del__(self):
        for f in self.files.values():
            f.close()


