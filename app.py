from flask import Flask, render_template, request, flash
from process import Process
import subprocess
import multiprocessing
import time
import gzip
import shutil
from TranscriptAnalysis import TranscriptAnalysis




import os
import json

app = Flask(__name__)
app.secret_key = b'_5#y2L"F4Q8z\n\xec]/kkk'





def getFileList(include_list=[]):
    try:
        parents = [folder for folder in os.listdir('./data/') if os.path.isdir(f'./data/{folder}')]
        parents = set(parents) & set(include_list) if include_list else parents
        parent2child = {parent:[folder for folder in os.listdir(f'./data/{parent}')] for parent in parents}
    except Exception as e:
        return f'<p>Error somewhere</p>{e}'
    return [f'{parent}/{child}' for parent in parent2child for child in parent2child[parent]]


def getGeneSets():
    try:
        gene_set_list = [file.replace('.json', '') for file in os.listdir('./data/genelist') if file.endswith('json')]
        return gene_set_list
    except Exception as e:
        print(f'Error getting gene sets: {e}')
        return []



def organizeFileList(file_list):
    organized_list = []
    for parent in sorted(set([path.split('/')[0] for path in file_list])):
        organized_list.append((f'<h5 class="mt-3"> {parent} </h5>', ''))
        for file in [file for file in file_list if file.startswith(f'{parent}/')]:
            organized_list.append((file.split('/')[-1], file))
        organized_list.append(('', ''))
    return organized_list


process_instance = Process()
process_func = {'filter':process_instance.filterData, 'bin':process_instance.binData, 'compare':process_instance.compareDataset}


@app.route("/", methods=['GET', 'POST'])
def home():
    if request.method == 'POST':
        dummy, message = process_func[request.form['type']](request.form)
        #return request.form
        #return dummy
        flash(message)
    file_list = getFileList()
    return render_template('home.html', data=organizeFileList(file_list))


@app.route("/filter")
def filter():
    # add ppm filter
    # add gene type filter
    file_list = getFileList(['original', 'filtered'])
    return render_template('filter.html', data=file_list)


@app.route("/bin")
def bin():
    # add ppm filter
    # add gene type filter
    file_list = getFileList(['filtered'])
    gene_sets = getGeneSets()
    return render_template('bin.html', data=file_list, gene_sets=gene_sets)


@app.route("/compare")
def compare():
    file_list = getFileList(['binned'])
    return render_template('compare.html', data=file_list)


@app.route("/annotate")
def annotate():
    file_list = getFileList()
    return render_template('annotate.html', data=file_list)


@app.route("/align", methods=['GET', 'POST'])
def align():
    if request.method == 'POST':
        print(request.form['path'])
        bowtie2(request.form['path'])
        flash('Alignment process started!')
    return render_template('preprocess.html')


def createDir(folder):
    try:
        os.mkdir(folder)
    except Exception as e:
        print(f'folder already exists: {e}, {folder}')

def bowtie2(folder):
    process_dict = {}
    print('>>>Bowtie', folder)
    sam_folder = os.path.join(folder,'sam')
    split_folder = os.path.join(folder, 'split')
    createDir(sam_folder)
    createDir(split_folder)

    for file in os.listdir(folder):
        if file.endswith('.fa'):
            file_path = os.path.join(folder, file)
            sam_path = os.path.join(sam_folder, file)
            process_dict[file] = ['bowtie2', f'-a -f -x mappers/merged.ws274 -U {file_path} | grep', "-v '^@'",
                                  f'> {sam_path}.sam']

    p = multiprocessing.Process(target=processController, args=(process_dict, sam_folder, split_folder))
    p.start()

def processController(process_dict, output_folder, split_folder):
    print('inside process controller')
    q = multiprocessing.Queue()
    global_start = time.time()
    alignment_report = {}
    proc_list = []
    for fi, file in enumerate(process_dict):
        proc_list.append([f'{fi}/{len(process_dict)}', file, process_dict[file], output_folder, split_folder])
        #alignment_report[file] = handlePreProcessing(file, process_dict[file], output_folder, split_folder)
    pool = multiprocessing.Pool(os.cpu_count()//2)
    report = pool.map(func=handlePreProcessing, iterable=proc_list, chunksize=1)
    pool.close()
    pool.join()
    for file, output in report:
        alignment_report[file] = output

    with open(os.path.join(output_folder, 'alignment_report.json'), 'w') as f:
        json.dump(alignment_report, f)
    print('finished with alignment', f'{int(time.time() - global_start)} sec elapsed')



def handlePreProcessing(args):
    fi, file, cmd, output_folder, split_folder = args
    start = time.time()
    print(f'{file} ({fi})', 'starting..')
    sp = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = sp.communicate()
    print(f'{file} ({fi})', 'bowtie2', f'{int(time.time() - start)} sec elapsed')
    sam_file = os.path.join(output_folder, f'{file}.sam')
    ta = TranscriptAnalysis(sam_file)
    seq2genes, gene2seq, seq2count, unique, multi = ta.run()
    print(f'{file} ({fi})', 'ta.run()', f'{int(time.time() - start)} sec elapsed')
    gene2total_ppm, seq2ppm, total_read_count = ta.calculateWith['everything'](gene2seq, seq2genes, seq2count, multi)
    norm_gene2total_ppm, norm_seq2ppm = ta.normalizeBy['everything'](gene2total_ppm, seq2ppm)
    ta.splitNormalizedSample(gene2seq, norm_seq2ppm, os.path.join(split_folder, file.replace('.fa', '')))
    print(f'{file} ({fi})', 'finished', f'{int(time.time() - start)} sec elapsed')
    compressFile(sam_file)
    print(f'{sam_file} compressed')
    return file, err.decode("utf-8").split('\n')


def compressFile(file):
    subprocess.Popen(['gzip', f'{file}'])
    return True
    with open(file, 'rb') as f_in:
        with gzip.open(f'{file}.gz', 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)

