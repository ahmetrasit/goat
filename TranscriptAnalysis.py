from collections import Counter
import os
import json
import re

class TranscriptAnalysis:
    #'s' for sense, 'a' for antisense
    gene2element = {}
    gene2name = {}
    gene2alias = {}
    trans2pos = {}
    exon2pos = {}
    rest2pos = {}
    
    sam_file = None
    type2parse = {}
    calculateWith = {}
    normalizeBy = {}
    
    def __init__(self, sam_file):
        self.sam_file = sam_file
        self.type2parse = {
            'PC_EXON':self.parseExon,
            'TRANSPO':self.parseTransposon,
            'TRANSCR':self.parseTranscript,
            'PSEUDOG':self.parsePseudogene}
        self.calculateWith = {'everything':self.calcTotal, 'structural':self.calcTotalExcludeStructural}
        self.normalizeBy = {'everything':self.normAll, 'tRNA':self.normtRNA, 'miRNA':self.normmiRNA}
        
        
    def run(self):
        seq2genes, seq2count = self.parseSamFile(self.sam_file)

        unique, multi = self.findMappers(seq2genes)
        gene2seq = self.getGene2Seq(seq2genes, unique, multi)
        #shorts = [curr for curr in gene2seq if len(curr) == 1]

        return seq2genes, gene2seq, seq2count, unique, multi
    
    
                        

    ########################################
    ### Pre-processing and Normalization ###
    
    def parseSamFile(self, sam_file):
        seq2genes = {}
        seq2count = {}
        with open(sam_file) as f:
            for line in f:
                if not line.startswith('@'):
                    seq_count, flag, target, pos, _, cigar, _ = line.split('\t', 6)
                    if not target.startswith('*'):
                        seq, count = seq_count.split(':')
                        strand = self.flag2strand(int(flag))
                        seq2count[seq] = int(count)
                        gene = self.parseMain(target, pos, seq, strand)
                        if len(gene) < 2:
                            print(gene, target)
                        if seq in seq2genes:
                            try:
                                seq2genes[seq][strand].add(gene)
                            except:
                                seq2genes[seq][strand] = set([gene])
                        else:
                            try:
                                seq2genes[seq] = {strand: set([gene])}
                            except Exception as e:
                                print(f'Error parsing a line:{target, pos, seq, strand}, {e}')
        #shorts = set([gene for seq in seq2genes for strand in seq2genes[seq] for gene in seq2genes[seq][strand] if len(gene) == 1])
        #print('!!>>S', shorts)
        return seq2genes, seq2count
    
    
    def calcTotal(self, gene2seq, seq2genes, seq2count, multi):
        gene2total_ppm = {}
        gene2total = self.getGene2Total(gene2seq, seq2genes, seq2count, multi)
        total_read_count = 0
        for gene in gene2total:
            total_read_count += sum([gene2total[gene][strand][mapper] for strand in gene2total[gene] for mapper in gene2total[gene][strand]])
        ppm_ratio = total_read_count/1000000
        for gene in gene2total:
            curr_values = gene2total[gene]
            gene2total_ppm[gene] = {'a':{'u':curr_values['a']['u']/ppm_ratio, 'm':curr_values['a']['m']/ppm_ratio}, 
                                    's':{'u':curr_values['s']['u']/ppm_ratio, 'm':curr_values['s']['m']/ppm_ratio}}

        seq2ppm = {}
        for seq in seq2count:
            seq2ppm[seq] = seq2count[seq]/ppm_ratio
        
        return gene2total_ppm, seq2ppm, total_read_count
    
    
    def calcTotalExcludeStructural(self, gene2seq, seq2genes, seq2count, _):
        exclude_set = set(['ASRNA', 'LINCRNA', 'RRNA', 'SCRNA', 'SNORNA', 'SNRNA'])

        def filterStructural(exclude_set, seq2genes):
            exclude_gene_set = set([gene for gene in gene2element if gene2element[gene] in exclude_set])
            filtered_seq2genes = {}
            for seq in seq2genes:
                targets = set([gene for strand in seq2genes[seq] for gene in seq2genes[seq][strand]])
                if len(targets & exclude_gene_set) == 0:
                    filtered_seq2genes[seq] = seq2genes[seq]
            filtered_unique, filtered_multi = self.findMappers(filtered_seq2genes)
            filtered_gene2seq = self.getGene2Seq(filtered_seq2genes, filtered_unique, filtered_multi)
            return filtered_gene2seq, filtered_seq2genes, filtered_multi

        filtered_gene2seq, filtered_seq2genes, filtered_multi = filterStructural(exclude_set, seq2genes)

        return self.calcTotal(filtered_gene2seq, filtered_seq2genes, seq2count, filtered_multi)   

    
    ### Normalization ###

    def sumOfType(self, gene2total_ppm, _type):
        total = 0
        for gene in gene2total_ppm:
            if self.gene2element[gene] == _type:
                total += sum([gene2total_ppm[gene][strand][mapper] for strand in gene2total_ppm[gene] for mapper in gene2total_ppm[gene][strand]])

        return total


    def normByType(self, gene2total_ppm, seq2ppm, _type):
        norm_gene2total_ppm = {}
        norm_seq2ppm = {}
        ref = self.sumOfType(gene2total_ppm, _type.upper())
        ratio = 1000000/ref
        print(ref)
        for gene in gene2total_ppm:
            curr = gene2total_ppm[gene]
            norm_gene2total_ppm[gene] = {'a': {'u':curr['a']['u']*ratio, 'm':curr['a']['m']*ratio},
                                         's': {'u':curr['s']['u']*ratio, 'm':curr['s']['m']*ratio}}
        for seq in seq2ppm:
            norm_seq2ppm[seq] = seq2ppm[seq] * ratio

        return norm_gene2total_ppm, norm_seq2ppm


    def normAll(self, gene2total_ppm, seq2ppm):
        return gene2total_ppm, seq2ppm


    def normtRNA(self, gene2total_ppm, seq2ppm):
        return self.normByType(gene2total_ppm, seq2ppm, 'TRNA')


    def normmiRNA(self, gene2total_ppm, seq2ppm):
        return self.normByType(gene2total_ppm, seq2ppm, 'MIRNA')

    
        
    
    
    ############################################
    ### Identifying Unique and Multi Mappers ###

    def findMappers(self, seq2genes):
        #find unique and multi mappers
        unique = set([])
        multi = set([])
        for seq in seq2genes:
            if len(seq) >= 10:
                curr = seq2genes[seq]
                if len(list(curr.values())[0]) == 1:
                    unique.add(seq)
                else:
                    multi.add(seq)
        return unique, multi

    def getMultiCategory(self, seq, seq2genes, gene2seq):
        count_existing = 0
        genes = [gene for strand in seq2genes[seq] for gene in seq2genes[seq][strand]]
        for strand in seq2genes[seq]:
            for gene in seq2genes[seq][strand]:
                try:
                    #gene2seq should only have unique mappers
                    count_existing += 1 if len(gene2seq[gene][strand]['u']) > 0 else 0
                except: #the gene has not been met, because it does not have any unique mappers
                    pass

        if count_existing == 0:
            return 'none'
        elif count_existing == len(genes):
            return 'all'
        else:
            return 'some'


    def getGene2Seq(self, seq2genes, unique, multi):
        #creating gene based sequence information
        gene2seq = {}
        for seq in unique:
            for strand in seq2genes[seq]:
                for gene in seq2genes[seq][strand]:
                    if gene not in gene2seq:
                        gene2seq[gene] = {
                            'a':{'u':set([]), 'm':{'all':set([]), 'some':set([]), 'none':set([])}},
                            's':{'u':set([]), 'm':{'all':set([]), 'some':set([]), 'none':set([])}}}
                    gene2seq[gene][strand]['u'].add(seq)

        for index, seq in enumerate(multi):
            for strand in seq2genes[seq]:
                for gene in seq2genes[seq][strand]:                
                    if gene not in gene2seq:
                        gene2seq[gene] = {
                            'a':{'u':set([]), 'm':{'all':set([]), 'some':set([]), 'none':set([])}},
                            's':{'u':set([]), 'm':{'all':set([]), 'some':set([]), 'none':set([])}}}
                    category = self.getMultiCategory(seq, seq2genes, gene2seq)
                    gene2seq[gene][strand]['m'][category].add(seq)

        return gene2seq

    
    def getGene2Total(self, gene2seq, seq2genes, seq2count, multi, filterGene2Seq=None, filterGene2Seq_args=None):
        #calculating total for unique mappers
        gene2total = {gene:{'a':{'u':0, 'm':0}, 's':{'u':0, 'm':0}} for gene in gene2seq}
        filtered_gene2seq = gene2seq if not filterGene2Seq else filterGene2Seq(gene2seq, filterGene2Seq_args)
        for gene in filtered_gene2seq:
            for strand in filtered_gene2seq[gene]:
                for seq in filtered_gene2seq[gene][strand]['u']:
                    gene2total[gene][strand]['u'] += seq2count[seq]

        #calculating total of multimappers, in regards to the ratio of unique mappers of a given gene
        multi2total = {}
        for seq in multi:
            multi2total[seq] = 0
            for strand in seq2genes[seq]:
                if 'a' in seq2genes[seq]:
                    multi2total[seq] += sum([gene2total[gene]['a']['u']+gene2total[gene]['s']['u'] for gene in seq2genes[seq]['a']])
                if 's' in seq2genes[seq]:
                    multi2total[seq] += sum([gene2total[gene]['a']['u']+gene2total[gene]['s']['u'] for gene in seq2genes[seq]['s']])

        for gene in gene2seq:
            for strand in gene2seq[gene]:
                for category in gene2seq[gene][strand]['m']:
                    for seq in gene2seq[gene][strand]['m'][category]:
                        try:
                            gene2total[gene][strand]['m'] += self.addByCategory(category, seq2count[seq], gene2total[gene][strand]['u'], multi2total[seq], seq2genes[seq]) 
                        except:
                            pass

        return gene2total


    def addByCategory(self, category, seq_count, gene_strand_u_count, genes_total_count, genelist):
        def countAll(seq_count, gene_strand_u_count, genes_total_count, genelist):
            return seq_count * (gene_strand_u_count / genes_total_count)

        def countNone(seq_count, gene_strand_u_count, genes_total_count, genelist):
            return seq_count * (1 / len([gene for strand in genelist for gene in genelist[strand]]))

        def countSome(seq_count, gene_strand_u_count, genes_total_count, genelist):
            return seq_count * ((gene_strand_u_count if gene_strand_u_count > 0 else 1 ) / genes_total_count)

        return {'all':countAll, 'some':countSome, 'none':countNone}[category](seq_count, gene_strand_u_count, genes_total_count, genelist)


    def filterGene2Seq(self, gene2seq, nt, len_min, len_max):
        def accepted(seq):
            return True if seq[0] in nt and len_max >= len(seq) >= len_min else False

        filtered = {}
        for gene in gene2seq:
            filtered[gene] = {
                'a':{'u':set([]), 'm':{'all':set([]), 'some':set([]), 'none':set([])}},
                's':{'u':set([]), 'm':{'all':set([]), 'some':set([]), 'none':set([])}}}
            for strand in gene2seq[gene]:
                filtered[gene][strand]['u'] = set([seq for seq in gene2seq[gene][strand]['u'] if accepted(seq)])
                for category in filtered[gene][strand]['m']:
                    filtered[gene][strand]['m'][category] = set([seq for seq in gene2seq[gene][strand]['m'][category] if accepted(seq)])

        return filtered

    
    def updatePosDicts(self, seq2ppm):
        gene2pos = {}
        def updateSelected(curr_dict, seq2ppm):
            selected = {}
            for curr in curr_dict:
                selected[curr] = {'a': {}, 's': {}}
                for strand in 'as':
                    for pos in curr_dict[curr][strand]:
                        for seq in curr_dict[curr][strand][pos]:
                            try:
                                selected[curr][strand][int(pos)].append({seq: round(seq2ppm[seq], 4)})
                            except:
                                selected[curr][strand][int(pos)] = {seq: round(seq2ppm[seq], 4)}
            return selected
        pairs = {'trans': self.trans2pos, 'exon': self.exon2pos, 'rest': self.exon2pos}
        gene2pos = {'trans': {}, 'exon': {}, 'rest': {}}

        for category in pairs:
            gene2pos[category] = updateSelected(pairs[category], seq2ppm)

        return gene2pos


    
    ########################
    ### Parsing SAM File ###
    
    def getTranscriptGeneName(self, name):
        match = re.match(r'(\w+\.t?\d+)((\.\d+)|([a-z](\.\d+)*))*', name.split(":")[-1])
        if match:
            return match.group(1)
        return '-'

    def flag2strand(self, flag):
        if flag % 256 == 0:
            return 's'
        elif flag % 256 == 16:
            return 'a'
        else:
            return 'error'

    def parseTransposon(self, meta, pos, seq, strand):
        gene = meta.split(':')[2]
        self.gene2element[gene] = 'TRANSPOSON'
        self.gene2name[gene] = gene
        self.gene2alias[gene] = gene
        try:
            self.rest2pos[gene][strand][int(pos)].append(seq)
        except:
            if gene not in self.rest2pos:
                self.rest2pos[gene] = {'a': {}, 's': {}}
            self.rest2pos[gene][strand] = {int(pos): [seq]}
        if len(gene) < 2:
            print(gene, 'transposon')
        return gene

    def parseExon(self, meta, pos, seq, strand):
        _, exon_id, genes = meta.split(':', 2)
        for pair in genes.split("|"):
            gene, alias_name = pair.split(":")
            self.gene2element[gene] = 'PC'
            self.gene2alias[gene] = alias_name.split("_")[0]
            self.gene2name[gene] = alias_name.split("_")[-1]
        try:
            self.exon2pos[exon_id][strand][int(pos)].append(seq)
        except:
            if exon_id not in self.exon2pos:
                self.exon2pos[exon_id] = {'a':{}, 's':{}}
            self.exon2pos[exon_id][strand] = {int(pos):[seq]}
        if len(gene) < 2:
            print(gene, 'exon')
        return gene

    def parseTranscript(self, meta, pos, seq, strand):
        _, trans, gene = meta.split(':')
        self.gene2element[gene] = 'PC'
        name = self.getTranscriptGeneName(trans)
        self.gene2name[gene] = name
        try:
            self.trans2pos[trans][strand][int(pos)].append(seq)
        except:
            if trans not in self.trans2pos:
                self.trans2pos[trans] = {'a': {}, 's': {}}
            self.trans2pos[trans][strand] = {int(pos): [seq]}
        if len(gene) < 2:
            print(gene, 'transcript')
        return gene

    def parsePseudogene(self, meta, pos, seq, strand):
        _, name, gene = meta.split(':')
        self.gene2element[gene] = 'PSEUODOGENE'
        self.gene2name[gene] = name
        self.gene2alias[gene] = name
        try:
            self.rest2pos[gene][strand][int(pos)].append(seq)
        except:
            if gene not in self.rest2pos:
                self.rest2pos[gene] = {'a': {}, 's': {}}
            self.rest2pos[gene][strand] = {int(pos): [seq]}
        if len(gene) < 2:
            print(gene, 'pseudogene')
        return gene

    def parseRest(self, meta, pos, seq, strand):
        #print(meta)
        try:
            element, _, gene, alias_name = meta.split(':')
            self.gene2alias[gene] = alias_name.split("_")[0]
            self.gene2name[gene] = alias_name.split("_")[-1]
            self.gene2element[gene] = element
            try:
                self.rest2pos[gene][strand][int(pos)].append(seq)
            except:
                if gene not in self.rest2pos:
                    self.rest2pos[gene] = {'a': {}, 's': {}}
                self.rest2pos[gene][strand] = {int(pos): [seq]}
        except:
            print(meta)
        if len(gene) < 2:
            print(gene, 'rest')
        return gene

    def parseMain(self, meta, pos, seq, strand):
        tag = meta[:7]
        if tag in self.type2parse:
            return self.type2parse[tag](meta, pos, seq, strand)
        else:
            return self.parseRest(meta, pos, seq, strand)
        
        
    ###################
    ### Accessories ###  
    def createRatioMatrix(gene2total_ppm):
        ref_order = [
            'PC', 'MIRNA', 'NCRNA', 'PSEUODOGENE', 'TRANSPOSON', 'RRNA',
            'LINCRNA', 'PIRNA', 'TRNA', 'SNORNA', 'ASRNA', 'SNRNA', 'SCRNA']
        ratios = sorted([(int(self.sumOfType(gene2total_ppm, _type)), _type) for _type in set(gene2element.values())], reverse=True)
        curr_order_elements = set([y for x,y in ratios])
        order = []
        for curr in ref_order:
            if curr in curr_order_elements:
                order.append(curr)

        melted = []
        for val1, _type1 in ratios:
            for val2, _type2 in ratios:
                melted.append([_type1, _type2, round(val1/val2, 2)])
        df = pd.DataFrame(melted)
        df.columns = 'type1', 'type2', 'ratio'
        matrix = df.pivot_table(index='type1', columns='type2', values='ratio').reindex(order, axis=0).reindex(order, axis=1)
        ratio_matrix = matrix.where(np.triu(np.ones(matrix.shape)).astype(np.bool))

        return ratio_matrix
    
    ###########################################
    ### Split data into smaller json chunks ###
    def splitNormalizedSample(self, gene2seq, norm_seq2ppm, sample_folder):
        self.createFilteredJson(gene2seq, norm_seq2ppm, sample_folder)


    def createFilteredJson(self, gene2seq, norm_seq2ppm, folder):
        for nt in 'ATGC':
            for length in range(10, 36):
                species = f'{length}{nt}'
                len_start = length
                len_end = length
                filtered = self.filterGenesBySpecies(gene2seq, len_start, len_end, nt, norm_seq2ppm)
                if filtered:
                    with open(os.path.join(folder, f'{species}.json'), 'w') as f:
                        json.dump(filtered, f)


    def filterSeqBySpecies(self, seq_set, len_start, len_end, nt, norm_seq2ppm):
            filtered = {seq:round(norm_seq2ppm[seq], 4) for seq in seq_set if seq.startswith(nt) and len_start >= len(seq) >= len_end}
            return filtered, bool(filtered)

    def filterGenesBySpecies(self, gene2seq, len_start, len_end, nt, norm_seq2ppm):
        filtered_gene2seq = {}
        for gene in gene2seq:
            found_set = set([])
            filtered_gene2seq[gene] = { 'a':{'u':set([]), 'm':{'all':set([]), 'some':set([]), 'none':set([])}},
                                        's':{'u':set([]), 'm':{'all':set([]), 'some':set([]), 'none':set([])}}}
            for strand in 'as':
                filtered_gene2seq[gene][strand]['u'], curr_found = self.filterSeqBySpecies(gene2seq[gene][strand]['u'], len_start, len_end, nt, norm_seq2ppm)
                found_set.add(curr_found)
                for multi_type in ['all', 'some', 'none']:
                    filtered_gene2seq[gene][strand]['m'][multi_type], curr_found = self.filterSeqBySpecies(gene2seq[gene][strand]['m'][multi_type], len_start, len_end, nt, norm_seq2ppm)
                    found_set.add(curr_found)

            #remove empty genes
            if not any(found_set):
                del filtered_gene2seq[gene]

        return filtered_gene2seq    


    def createDir(self, folder):
        try:
            os.mkdir(folder)
        except Exception as e:
            print(f'folder already exists: {e}, {folder}')
        return folder