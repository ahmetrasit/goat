import json
import os
from markupsafe import escape
import html
import re
import logging


class Process:
    def __init__(self):
        print('ready')

    def setops(self, formdata):
        groupA = formdata['groupA']
        groupB = formdata['groupB']
        groupC = formdata['groupC']
        operation = formdata['hidden_operation']
        save_filename = escape(formdata['save']).replace('\s', '_')
        id_sets = self.prepIdSets()
        set_a = self.getGeneListInLocusName(groupA, id_sets)
        set_b = self.getGeneListInLocusName(groupB, id_sets)
        set_c = self.getGeneListInLocusName(groupC, id_sets) if groupC else []
        print('>>', operation)
        output = self.applyOperation(set_a, set_b, set_c, operation)
        if output:
            saved_filename = self.saveFile(f'data/genelist/{save_filename}.json', list(output))
            return '', f'Genelist saved as "{saved_filename}", having {len(output)} genes.'
        else:
            return '', '!No genes found after set operation, file is not saved.'

    def filterData(self, formdata):
        sample = formdata['dataset'].split('/')[-1]
        len_start = int(formdata['length_start'])
        len_end = int(formdata['length_end'])
        nts = ''.join([nt for nt in formdata.getlist('nucleotides')])
        strand = formdata['strand']

        rules = {}
        for key in formdata.keys():
            if key.startswith('hidden_rule'):
                group_id, rule_id = key.replace('hidden_rule_', '').split('_')
                rule = formdata[key].split(';')
                try:
                    rules[group_id].append(rule)
                except:
                    rules[group_id] = [rule]

        message = ''
        if sample.endswith('json'):
            options = {'nts': nts, 'len_start': len_start, 'len_end': len_end}
            data, message = self.filterJsonFile(f'data/{formdata["dataset"]}', {}, options)
        else:
            options = {'rules': rules}
            data = {}
            for nt in nts:
                for length in range(len_start, len_end + 1):
                    file = f'data/{formdata["dataset"]}/{length}{nt}.json'
                    data, message = self.filterFileBundle(file, data, strand)

        data = self.furtherFilter(data, rules)
        data = self.removeEmptyGenes(data)

        if data:
            len_interval = len_start if len_start == len_end else f'{len_start}-{len_end}'
            strand_long = {'a': 'asense', 's': 'sense'}[strand]

            save_filename = escape(formdata['save']).replace('\s', '_')
            filename = f'data/filtered/{save_filename}_{sample.replace(".json", "")}_{len_interval}{nts}_{strand_long}.json'
            saved_filename = self.saveFile(filename, data)

            return data, f'Saved as {saved_filename}'
        else:
            return [], 'Nothing found, file is not created!'

    def speciesFilter(self, seq_dict, options):
        filtered = {}
        for seq in seq_dict:
            if seq[0] in options['nts'] and options['len_end'] >= len(seq) >= options['len_start']:
                filtered[seq] = seq_dict[seq]
        return filtered

    def filterJsonFile(self, file, data, options):
        with open(file) as f:
            curr_data = json.load(f)
            for gene in curr_data:
                curr_gene_data = self.speciesFilter(curr_data[gene], options)

                if len(curr_gene_data) > 0:
                    if gene in data:
                            data[gene] = {**data[gene], **curr_gene_data}
                    else:
                        data[gene] = curr_gene_data

        return data, 'Success'


    def filterFileBundle(self, file, data, strand):
        with open(file) as f:
            curr_data = json.load(f)
            for gene in curr_data:
                u_data = curr_data[gene][strand]['u']
                m_data = {}
                for multi_type in curr_data[gene][strand]['m']:
                    m_data = {**m_data, **curr_data[gene][strand]['m'][multi_type]}

                curr_gene_data = {'u': u_data, 'm': m_data}

                if len(u_data) > 0 or len(m_data) > 0:
                    if gene in data:
                        for mapper in 'um':
                            data[gene][mapper] = {**data[gene][mapper], **curr_gene_data[mapper]}
                    else:
                        data[gene] = curr_gene_data

        return data, 'Success'


    #####################
    ### FurtherFilter ###
    def furtherFilter(self, data, rules):
        filtered = {}
        initial_gene_len = len(data)
        seq_rules, gene_rules = self.parseRules(rules)

        for gene in data:
            passed_seq_rules = {}
            seq_dict = {**data[gene]['u'], **data[gene]['m']} if 'u' in data[gene] else data[gene]
            for seq in seq_dict:
                count = seq_dict[seq]
                if self.passedSeqRules(seq, count, seq_rules):
                    passed_seq_rules[seq] = count
            passed_gene_rules = self.passedGeneRules(passed_seq_rules, gene_rules)
            filtered[gene] = passed_gene_rules
        print(f'{len(filtered)}/{initial_gene_len} genes filtered using rules')
        return filtered

    def passedSeqRules(self, seq, count, rules):
        if not rules:
            return seq
        for rule in rules:
            ref_value = float(rule[1])
            operator = {'>': float(count) > ref_value,
                        '<': float(count) < ref_value,
                        '=': float(count) == ref_value,
                        '≥': float(count) >= ref_value,
                        '≤': float(count) <= ref_value
                        }
            return operator[rule[0]]

    def passedGeneRules(self, seq_dict, rules):
        if not rules:
            return seq_dict
        passed = set([])
        run_count = 0
        for rule in rules:
            sorted_seq_list = sorted(seq_dict, key=seq_dict.get)
            ratio = float(rule[1]) / 100
            levels = {'Top': sorted_seq_list[::-1][:int(len(sorted_seq_list) * ratio) + 1],
                      'Bottom': sorted_seq_list[:int(len(sorted_seq_list) * ratio) + 1]
                      }
            curr_passed = set(levels[rule[0]])
            if run_count == 0:
                passed = curr_passed
            else:
                passed = passed & curr_passed
            run_count += 1
        return {curr_seq: seq_dict[curr_seq] for curr_seq in passed}

    def parseRules(self, rules):
        gene_rules = []
        seq_rules = []
        for group_id in rules:
            for rule in rules[group_id]:
                if 'Top' in rule[0] or 'Bottom' in rule[0]:
                    gene_rules.append(rule)
                else:
                    seq_rules.append(rule)
        return seq_rules, gene_rules

    def removeEmptyGenes(self, data):
        cleaned_data = {}
        for gene in data:
            if data[gene]:
                cleaned_data[gene] = data[gene]
        return cleaned_data

    def binData(self, formdata):
        mappers = ''
        mappers += 'u' if 'unique_mappers' in formdata else ''
        mappers += 'm' if 'multi_mappers' in formdata else ''
        mappers = 'um' if mappers=='' else mappers
        gene_sets =  formdata.getlist('gene_sets')
        gene_types = formdata.getlist('gene_types')

        has_sequence = True if 'hasSequence' in formdata else False
        is_transcript = True if 'isTranscript' in formdata else False
        save_filename = escape(formdata['save']).replace('\s', '_')
        dataset = formdata['dataset']
        sample = dataset.split('/')[-1]

        binned = {}
        with open(f'data/{dataset}') as f:
            curr_data = json.load(f)
            for gene in curr_data:
                genename = gene if not is_transcript else self.getTranscriptGeneName(gene)
                if not has_sequence:    #assuming data does not have strand information
                    try:
                        binned[genename] += float(curr_data[gene])
                    except:
                        binned[genename] = float(curr_data[gene])
                else:
                    if 'u' in curr_data[gene]:   #means there's mapper information in the file
                        for mapper in mappers:
                            for seq in curr_data[gene][mapper]:
                                try:
                                    binned[genename] += float(curr_data[gene][mapper][seq])
                                except:
                                    binned[genename] = float(curr_data[gene][mapper][seq])
                    else:
                        for seq in curr_data[gene]:
                            try:
                                binned[genename] += float(curr_data[gene][seq])
                            except:
                                binned[genename] = float(curr_data[gene][seq])


        binned_gene_set = set(binned.keys())
        id_sets = self.prepIdSets()
        binned_id_type, perc_common = self.getIdType(id_sets, binned_gene_set)

        if gene_types:
            passing_gene_type_filter = self.filterByGeneTypes(binned_gene_set, binned_id_type, gene_types)
        else:
            passing_gene_type_filter = binned_gene_set

        #print('@passing_gene_type_filter', len(passing_gene_type_filter), list(passing_gene_type_filter)[:10])

        if gene_sets:
            passing_gene_set_filter = self.filterByGeneSets(passing_gene_type_filter, binned_id_type, gene_sets, id_sets)
        else:
            passing_gene_set_filter = passing_gene_type_filter

        stats = f'Started with {len(binned_gene_set)} genes, {perc_common}% of them can be mapped. {len(passing_gene_type_filter)} genes passed gene type filter. {len(passing_gene_set_filter)} genes passed gene set filter.'
        if passing_gene_set_filter:
            filename = f'data/binned/{save_filename}_{sample.replace(".json", "")}_{mappers}_binned.json'
            data = {gene:binned[gene] for gene in passing_gene_set_filter}
            saved_filename = self.saveFile(filename, data)

            return passing_gene_set_filter, stats + f' Saved as {saved_filename}'
        else:
            return [], stats + ' No file created.'


    def compareDataset(self, formdata):
        def applyRule(dataA, dataB, element, rule):
            def getDatumFromRule(rule_first, currA, currB):
                datum = {'A': currA, 'B': currB, 'A/B': currA / currB if currB > 0 else 0,
                         'B/A': currB / currA if currA > 0 else 0}
                return datum[rule_first]

            def returnRule(datum, rule):
                ref_value = float(rule[2])
                operator = {'>': datum > ref_value,
                            '<': datum < ref_value,
                            '=': datum == ref_value,
                            '≥': datum >= ref_value,
                            '≤': datum <= ref_value
                            }
                return operator[rule[1]]

            currA = float(dataA[element]) if element in dataA else 0
            currB = float(dataB[element]) if element in dataB else 0
            datum = getDatumFromRule(rule[0], currA, currB)
            return returnRule(datum, rule)

        def getResult(dataA, dataB, rules):
            output = set([])
            elements = set(dataA.keys()) | set(dataB.keys())
            for element in elements:
                if all([applyRule(dataA, dataB, element, rule) for rule in rules]):
                    output.add(element)
            return output

        rules = {}
        for key in formdata.keys():
            if key.startswith('hidden_rule'):
                group_id, rule_id = key.replace('hidden_rule_', '').split('_')
                rule = formdata[key].split(';')
                try:
                    rules[group_id].append(rule)
                except:
                    rules[group_id] = [rule]

        with open(f'data/{formdata["groupA"]}') as f:
            dataA = json.load(f)
        with open(f'data/{formdata["groupB"]}') as f:
            dataB = json.load(f)

        output = set([])
        for group in rules:
            output |= getResult(dataA, dataB, rules[group])

        save_filename = escape(formdata['save']).replace('\s', '_')
        saved_filename = ''
        if len(output) > 0:
            filename = f'data/genelist/{save_filename}.json'
            saved_filename = self.saveFile(filename, list(output))

        return json.dumps(list(output)), f'{len(output)} genes passing the filter.' + f'Saved as {saved_filename}' if saved_filename else 'No genes passed the filter, data not saved.'


    def getTranscriptGeneName(self, name):
        match = re.match(r'(\w+\.t?\d+)((\.\d+)|([a-z](\.\d+)*))*', name.split(":")[-1])
        if match:
            return match.group(1)
        return '-'


    def prepIdSets(self):
        id_sets = {}
        for curr in ['alias2name', 'name2gene', 'gene2name']:
            with open(f'mappers/{curr}.json') as f:
                id_sets[curr] = set(json.load(f).keys())
        return id_sets

    def prepIdSetsWithConverters(self):
        id_sets = {}
        converters = {}
        for curr in ['alias2name', 'name2gene', 'gene2name']:
            with open(f'mappers/{curr}.json') as f:
                curr_data = json.load(f)
                id_sets[curr] = set(curr_data.keys())
                converters[curr] = curr_data
        return id_sets, converters



    def getIdType(self, id_sets, gene_set):
        max = 0
        found_type = ''
        for curr_type in id_sets:
            common = len(id_sets[curr_type] & set(gene_set))
            if common > max:
                max = common
                found_type = curr_type
        return found_type.split('2')[0], 100*max//len(gene_set)


    def filterByGeneTypes(self, gene_set, id_type, selected_gene_types):
        with open(f'mappers/{id_type}2type.json') as f:
            id2type = json.load(f)
        return set([gene for gene in gene_set if gene in id2type and id2type[gene] in selected_gene_types])

    def getTypeDict(self, id_type):
        with open(f'mappers/{id_type}2type.json') as f:
            converter = json.load(f)
            return converter


    def filterByGeneSets(self, gene_set, id_type, selected_gene_sets, id_sets):
        genes_in_selected_sets = set([])
        for selected in selected_gene_sets:
            with open(f'data/genelist/{selected}.json') as f:
                curr_data = json.load(f)
                curr_id_type, perc_common = self.getIdType(id_sets, curr_data)
                converted_ids = self.convert(curr_data, curr_id_type, id_type) if curr_id_type != id_type else set(curr_data)
                genes_in_selected_sets |= converted_ids

        return set(gene_set) & genes_in_selected_sets


    def convert(self, genes, origin, target):
        with open(f'mappers/{origin}2{target}.json') as f:
            converter = json.load(f)
        return set([converter[gene] for gene in genes if gene in converter])


    def getConverter(self, origin, target):
        with open(f'mappers/{origin}2{target}.json') as f:
            converter = json.load(f)
        return converter


    def saveFile(self, filename, data):
        print(filename)
        filename = re.sub('.json$', '', filename)
        *folder_fields, file = filename.split('/')
        folder = os.path.join(*folder_fields)
        file = escape(file)
        file = re.sub(r'[^\w\s-]', '', file.lower())
        file = re.sub(r'[-\s]+', '-', file).strip('-_')
        print(file)
        file_path = os.path.join(folder, file) + '.json'
        print('>FP', file_path)
        while os.path.isfile(file_path):
            file = re.sub('.json$', '', file)
            file += '_1'
            file += '.json'
            file_path = os.path.join(folder, file)
        print('>FP2', file_path)
        with open(file_path, 'w') as f:
            json.dump(data, f)
        return f'{file}'


    def getGeneListInLocusName(self, selected_gene_set_name, id_sets):
        try:
            with open(f'data/genelist/{selected_gene_set_name}.json') as f:
                curr_data = json.load(f)
                curr_id_type, perc_common = self.getIdType(id_sets, curr_data)
                converted_ids = self.convert(curr_data, curr_id_type, 'name') if curr_id_type != 'name' else set(curr_data)
                return converted_ids
        except Exception as e:
            print(f'>>>>>>>Error with genelist {selected_gene_set_name}:{e}')
            return set([])

    def applyOperation(self, set_a, set_b, set_c, operation):
        operations = {'and': set_a & set_b, 'or': set_a | set_b, 'a-b': set_a - set_b, 'b-a': set_b - set_a}
        return operations[operation]

