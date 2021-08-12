import json
import os
from markupsafe import escape
import html
import subprocess
import multiprocessing
import time
import gzip
import shutil
from TranscriptAnalysis import TranscriptAnalysis
import random
import string

class Preprocess:
    def __init__(self):
        pass


    def main(self, formdata):
        if 'filename_0' in formdata:
            path = formdata['path'].strip()
            exp_name = formdata['experiment_name'].strip()
            if len(exp_name) == 0:
                letters = string.ascii_lowercase
                exp_name = ''.join(random.choice(letters) for i in range(10))
            file2alias = {}
            alias_list = [formdata[key].strip() for key in formdata.keys() if key.startswith('alias_')]
            alias_list_ok = True if len(alias_list) == len(set(alias_list)) else False
            for key in formdata.keys():
                if key.startswith('filename_'):
                    curr_filename = formdata[key]
                    if alias_list_ok:
                        file_id = key.split('_')[1]
                        curr_alias = formdata['alias_'+file_id].strip()
                        file2alias[curr_filename] = curr_alias if len(curr_alias) > 0 else formdata[key]
                    else:
                        file2alias[curr_filename] = curr_filename
            exp_folder = self.createDir(os.path.join(path, exp_name))
            with open(os.path.join(exp_folder, 'file2alias.json'), 'w') as f:
                json.dump(file2alias, f)

            if 'use_counts_fa' in formdata:
                print('>>>>counts')
                self.align(path, exp_folder, file2alias, exp_name)
            else:
                adapter_seq = formdata['adapter_seq']
                min_seq_len = formdata['min_seq_len']
            return formdata, ''
        else:
            return '', 'No files found in the path'



    def align(self, data_folder, exp_folder, file2alias, exp_name):
        process_dict = {}
        print('>>>Bowtie', data_folder)
        sam_folder = os.path.join(exp_folder, 'sam')
        split_folder = os.path.join(exp_folder, 'split')
        self.createDir(sam_folder)
        self.createDir(split_folder)

        for file in os.listdir(data_folder):
            if file.endswith('.fa'):
                file_path = os.path.join(data_folder, file)
                sam_path = os.path.join(sam_folder, file2alias[file])
                process_dict[file] = ['bowtie2', f'-a -f -x mappers/merged.ws274 -U {file_path} | grep', "-v '^@'",
                                      f'> {sam_path}.sam']

        p = multiprocessing.Process(target=self.processController, args=(process_dict, sam_folder, split_folder, file2alias, exp_name))
        p.start()

    def processController(self, process_dict, output_folder, split_folder, file2alias, exp_name):
        print('inside process controller')
        q = multiprocessing.Queue()
        global_start = time.time()
        alignment_report = {}
        proc_list = []
        for fi, file in enumerate(process_dict):
            proc_list.append([f'{fi}/{len(process_dict)}', file, process_dict[file], output_folder, split_folder, file2alias])
        pool = multiprocessing.Pool(os.cpu_count() // 2)
        report = pool.map(func=self.handlePreProcessing, iterable=proc_list, chunksize=1)
        pool.close()
        pool.join()
        for file, output in report:
            alignment_report[file] = output

        with open(os.path.join(output_folder, 'alignment_report.json'), 'w') as f:
            json.dump(alignment_report, f)

        self.moveFolder(split_folder, f'data/original/{exp_name}')
        print('finished with alignment and files are moved', f'{int(time.time() - global_start)} sec elapsed')

    def handlePreProcessing(self, args):
        fi, file, cmd, output_folder, split_folder, file2alias = args
        start = time.time()
        print(f'{file} ({fi})', 'starting..')
        sp = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        out, err = sp.communicate()
        print(f'{file} ({fi})', 'bowtie2', f'{int(time.time() - start)} sec elapsed')
        sam_file = os.path.join(output_folder, f'{file2alias[file]}.sam')
        ta = TranscriptAnalysis(sam_file)
        seq2genes, gene2seq, seq2count, unique, multi = ta.run()
        print(f'{file} ({fi})', 'ta.run()', f'{int(time.time() - start)} sec elapsed')
        gene2total_ppm, seq2ppm, total_read_count = ta.calculateWith['everything'](gene2seq, seq2genes, seq2count,
                                                                                   multi)
        norm_gene2total_ppm, norm_seq2ppm = ta.normalizeBy['everything'](gene2total_ppm, seq2ppm)
        ta.splitNormalizedSample(gene2seq, norm_seq2ppm, os.path.join(split_folder, file2alias[file]))
        print(f'{file} ({fi})', 'finished', f'{int(time.time() - start)} sec elapsed')
        self.compressFile(sam_file)
        return file, err.decode("utf-8").split('\n')


    def moveFolder(self, source_folder, target_folder):
        try:
            subprocess.Popen(['mv', f'{source_folder}', f'{target_folder}'])
        except Exception as e:
            print(f'Error moving folder {source_folder} to {target_folder} :{e}')

    def compressFile(self, file):
        try:
            subprocess.Popen(['gzip', f'{file}'])
        except Exception as e:
            print(f'Error gzipping {file}: {e}')


    def createDir(self, folder):
        try:
            os.mkdir(folder)
        except Exception as e:
            print(f'folder already exists: {e}, {folder}')
        return folder