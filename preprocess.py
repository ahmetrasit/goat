import json
import os
import shlex
from markupsafe import escape
import html
from subprocess import Popen, PIPE
import multiprocessing
import time
import shutil
from TranscriptAnalysis import TranscriptAnalysis
import random
import string
import gzip
import logging



class Preprocess:

    def __init__(self):
        pass


    def main(self, formdata):
        if 'filename_0' in formdata:
            path = formdata['path'].strip()
            username = formdata['username'].strip()
            exp_name = username + '_' + formdata['experiment_name'].strip()
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
            self.createDir(f'data/original/{exp_name}')
            with open(os.path.join(exp_folder, 'file2alias.json'), 'w') as f:
                json.dump(file2alias, f)

            if 'use_counts_fa' in formdata:
                print('>>>>counts')
                self.align(path, exp_folder, file2alias, exp_name)
            else:
                print('time for JBrowse')
                adapter_seq = formdata['adapter_seq']
                min_seq_len = formdata['min_seq_len']
                first_n = int(formdata['first_n']) + 1
                self.jbrowse_prep(path, exp_folder, file2alias, exp_name, adapter_seq, min_seq_len, first_n)
            return formdata, ''
        else:
            return '', 'No files found in the path'


    def jbrowse_prep(self, data_folder, exp_folder, file2alias, exp_name, adapter_seq, min_seq_len, first_n):
        p = multiprocessing.Process(target=self.doJBrowseProcessing,
                                    args=(data_folder, exp_folder, file2alias, exp_name, adapter_seq, min_seq_len, first_n))
        p.start()


    def doJBrowseProcessing(self, data_folder, exp_folder, file2alias, exp_name, adapter_seq, min_seq_len, first_n):
        #a separate main function for jbrowse processing
        folders = {'data_folder':data_folder, 'exp_folder':exp_folder}
        folders['trimmed_folder'] = self.createDir(os.path.join(exp_folder, 'trimmed'))
        folders['jbrowse_folder'] = self.createDir(os.path.join(exp_folder, 'jbrowse_data'))
        folders['counts_folder'] = self.createDir(os.path.join(exp_folder, 'counts'))
        folders['exp_folder'] = exp_folder

        proc_list = []
        fastq_gz_files = [file for file in os.listdir(data_folder) if file.endswith('.fastq.gz')]
        for fi, file in enumerate(fastq_gz_files):
            if file.endswith('.fastq.gz'):
                proc_list.append([f'{fi+1}/{len(fastq_gz_files)}', file, folders, file2alias, exp_name, adapter_seq, min_seq_len, first_n])
        pool = multiprocessing.Pool(os.cpu_count() // 2)
        report = pool.map(func=self.singleFileJBrowse, iterable=proc_list, chunksize=1)
        pool.close()
        pool.join()

        jbrowse_report = {}
        for file, *output in report:
            jbrowse_report[file] = output

        with open(os.path.join(folders['exp_folder'], 'alignment_report.json'), 'w') as f:
            json.dump(jbrowse_report, f)


    def singleFileJBrowse(self, args):
        fi, file, folders, file2alias, exp_name, adapter_seq, min_seq_len, first_n = args
        print('>before cutadapt')
        ca_out, ca_err = self.cutadaptSamples(fi, file, folders, file2alias, exp_name, adapter_seq, min_seq_len)

        cc_out, cc_err = self.createCountsFa(fi, file, folders, file2alias, first_n)
        _, a_err = self.singleSampleAlignProcess(fi, file, folders, file2alias, exp_name)
        cbbw_out, cbbw_err = self.createBamBigWig(fi, file, folders, file2alias)

        return file, [ca_out, ca_err], [cc_out, cc_err], [cbbw_out, cbbw_err]

    def singleSampleAlignProcess(self, fi, file, folders, file2alias, exp_name):
        sam_folder = self.createDir(os.path.join(folders['exp_folder'], 'sam'))
        split_folder = self.createDir(os.path.join(folders['exp_folder'], 'split'))
        file_path = os.path.join(folders['counts_folder'], file2alias[file]) + '.counts.fa'
        sam_path = os.path.join(sam_folder, file2alias[file])

        process = ["/usr/bin/bowtie2", f'-a -f -p 10 -x mappers/merged.ws274 -U {file_path}', f'> {sam_path}.sam']
        file, error = self.handlePreProcessing([f'{fi}', file, process, sam_folder, split_folder, file2alias])
        orig_split_folder = os.path.join(split_folder, file2alias[file])
        self.moveFolder(orig_split_folder, f'data/original/{exp_name}/{file2alias[file]}')

        return file, error

    def createCountsFa(self, fi, file, folders, file2alias, first_n):
        trimmed_file_path = os.path.join(folders['trimmed_folder'], file2alias[file]) + '.trimmed.fastq.gz'
        counts_file_prefix = os.path.join(folders['counts_folder'], file2alias[file])

        counts_fa_prep_cmd = ["scripts/create_counts_fa.sh", counts_file_prefix, trimmed_file_path, str(first_n)]
        print(counts_fa_prep_cmd)
        sp = Popen(counts_fa_prep_cmd, stdout=PIPE, stderr=PIPE)
        out, err = sp.communicate()
        for curr in [out, err]:
            for currline in str(curr).split('\n'):
                print(curr)
            print()
        return str(out), str(err)

    def createBamBigWig(self, fi, file, folders, file2alias):
        trimmed_file_path = os.path.join(folders['trimmed_folder'], file2alias[file]) + '.trimmed.fastq.gz'
        jbrowse_file_prefix = os.path.join(folders['jbrowse_folder'], file2alias[file])
        jupyter_prep_cmd = ["scripts/jupyter_prep.sh", jbrowse_file_prefix, trimmed_file_path, 'mappers/genome', '/users/ahmetrasit/ucsc-tools']
        sp = Popen(jupyter_prep_cmd, stdout=PIPE, stderr=PIPE, shell=False)
        out, err = sp.communicate()
        for curr in [out, err]:
            for currline in str(curr).split('\n'):
                for x in currline.split('\\n'):
                    print(x)
            print()
        return str(out).split('\n'), str(err).split('\n')

    def cutadaptSamples(self, fi, file, folders, file2alias, exp_name, adapter_seq, min_seq_len):
        #do not use multithreading, it may clog when multiple users submit jobs, and I'm already using cpu/2
        trimmed_file_path = os.path.join(folders['trimmed_folder'], file2alias[file]) + '.trimmed.fastq.gz'
        fastq_file_path = os.path.join(folders['data_folder'], file)
        cutadapt_cmd = ["cutadapt", "-a", adapter_seq, "-j", "10", "-m", min_seq_len, "-o", trimmed_file_path, fastq_file_path]
        print(' '.join(cutadapt_cmd))
        sp = Popen(cutadapt_cmd, stdout=PIPE, stderr=PIPE, shell=False)
        out, err = sp.communicate()
        return str(out).split('\n'), str(err).split('\n')

    def align(self, data_folder, exp_folder, file2alias, exp_name):
        process_dict = {}
        sam_folder = self.createDir(os.path.join(exp_folder, 'sam'))
        split_folder = self.createDir(os.path.join(exp_folder, 'split'))

        for file in os.listdir(data_folder):
            if file.endswith('.fa'):
                file_path = os.path.join(data_folder, file)
                sam_path = os.path.join(sam_folder, file2alias[file])
                process_dict[file] = ['/usr/bin/bowtie2', '-a',  '-f', 'p', '10', '-x mappers/merged.ws274',  f'-U {file_path}']

        p = multiprocessing.Process(target=self.processController, args=(process_dict, sam_folder, split_folder, file2alias, exp_name))
        p.start()

    def processController(self, process_dict, output_folder, split_folder, file2alias, exp_name):
        #print('inside process controller')
        q = multiprocessing.Queue()
        global_start = time.time()
        alignment_report = {}
        proc_list = []
        for fi, file in enumerate(process_dict):
            proc_list.append([f'{fi+1}/{len(process_dict)}', file, process_dict[file], output_folder, split_folder, file2alias])
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
        sp = Popen(cmd, stdout=PIPE, stderr=PIPE)
        out, err = sp.communicate()
        for curr in [out, err]:
            for currline in str(curr).split('\n'):
                for x in currline.split('\\n'):
                    print(x)
            print()
        print(f'{file} ({fi})', 'bowtie2', f'{int(time.time() - start)} sec elapsed')
        sam_file = os.path.join(output_folder, f'{file2alias[file]}.sam')
        ta = TranscriptAnalysis(sam_file)
        seq2genes, gene2seq, seq2count, unique, multi = ta.run()
        print(f'{file} ({fi})', 'ta.run()', f'{int(time.time() - start)} sec elapsed')
        gene2total_ppm, seq2ppm, total_read_count = ta.calculateWith['everything'](gene2seq, seq2genes, seq2count,                                                                                   multi)
        norm_gene2total_ppm, norm_seq2ppm = ta.normalizeBy['everything'](gene2total_ppm, seq2ppm)
        gene2pos = ta.updatePosDicts(norm_seq2ppm)
        sample_split_folder = self.createDir(os.path.join(split_folder, file2alias[file]))
        self.saveFile(sample_split_folder + f'/{file2alias[file]}.gene2pos.json', gene2pos)
        ta.splitNormalizedSample(gene2seq, norm_seq2ppm, sample_split_folder)
        print(f'{file} ({fi})', 'finished', f'{int(time.time() - start)} sec elapsed')
        self.compressFile(sam_file)
        return file, err.decode("utf-8").split('\n')


    def saveFile(self, path, data):
        with open(path, 'w') as f:
            json.dump(data, f)

    def moveFolder(self, source_folder, target_folder):
        try:
            Popen(['mv', f'{source_folder}', f'{target_folder}'])
        except Exception as e:
            print(f'Error moving folder {source_folder} to {target_folder} :{e}')

    def compressFile(self, file):
        try:
            Popen(['gzip', f'{file}'])
        except Exception as e:
            print(f'Error gzipping {file}: {e}')


    def createDir(self, folder):
        try:
            os.mkdir(folder)
        except Exception as e:
            print(f'folder already exists: {e}, {folder}')
        return folder