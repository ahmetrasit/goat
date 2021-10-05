from flask import Flask, render_template, request, flash, jsonify
from process import Process
from preprocess import Preprocess
from plot import Plot
import logging
import re

import subprocess
import multiprocessing
import time
import gzip
import shutil
from TranscriptAnalysis import TranscriptAnalysis





import os
import json

app = Flask(__name__)
app.config['SEND_FILE_MAX_AGE_DEFAULT'] = 0

app.secret_key = b'_5#y2L"F4Q8z\n\xec]/kkk'








process_instance = Process()
preprocess_instance = Preprocess()
plot_instance = Plot()
process_func = {'filter':process_instance.filterData,
                'bin':process_instance.binData,
                'compare':process_instance.compareDataset,
                'preprocess':preprocess_instance.main,
                'setops':process_instance.setops,
                }


@app.route("/", methods=['GET', 'POST'])
def home():
    if request.method == 'POST':
        print(request.form)
        dummy, message = process_func[request.form['type']](request.form)
        flash(message)
    file_list = getOrigFolderList() + getFileList()
    return render_template('home.html', data=organizeFileList(file_list))


@app.route("/filter")
def filter():
    # add ppm filter
    # add gene type filter
    orig_file_list = sorted(getOrigFolderList())
    file_list = getFileList(['original', 'filtered'])
    return render_template('filter.html', data=file_list+orig_file_list)


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
    return render_template('compare.html', main_file_list=file_list)


@app.route("/setops")
def setops():
    file_list = getFileList(['genelist'])
    return render_template('setops.html', data=file_list)


@app.route("/preview/<groupA>/<groupB>/<groupC>/<operation>")
def preview(groupA, groupB, groupC, operation):
    id_sets = process_instance.prepIdSets()
    set_a = process_instance.getGeneListInLocusName(groupA.replace('|', ''), id_sets) if groupA != '|' else []
    set_b = process_instance.getGeneListInLocusName(groupB.replace('|', ''), id_sets) if groupB != '|' else []
    set_c = process_instance.getGeneListInLocusName(groupC.replace('|', ''), id_sets) if groupC != '|' else []
    output = process_instance.applyOperation(set_a, set_b, set_c, operation)
    return jsonify({'output':list(output), 'list_a':list(set_a), 'list_b':list(set_b), 'list_c':list(set_c)})


@app.route("/plot")
def plot():
    template_list = getPlotTemplateList()
    return render_template('plot.html', plot_templates=template_list)


@app.route("/getDataPair/<fileA>/<fileB>/<folder>")
def getDataPair(fileA, fileB, folder):
    def getConverter(data, id_type, converters):
        if id_type == 'name':
            return {gene:gene for gene in data}
        curr_cunverter = converters[f'{id_type}2name']
        return {gene:curr_cunverter[gene] if gene in curr_cunverter else gene for gene in data}

    pair = []
    data_a = getData(fileA, folder)
    data_b = getData(fileB, folder)
    id_sets, converters = process_instance.prepIdSetsWithConverters()
    id_type_a, perc_common_a = process_instance.getIdType(id_sets, data_a)
    id_type_b, perc_common_b = process_instance.getIdType(id_sets, data_b)
    converter_a = getConverter(data_a, id_type_a, converters)
    converter_b = getConverter(data_b, id_type_b, converters)
    type_converter = process_instance.getTypeDict('name')

    for gene in set(data_a.keys()) | set(data_b.keys()):
        converted_gene_name = converter_a[gene] if gene in converter_a else converter_b[gene]
        pair.append({'gene':converted_gene_name,
                     'type':type_converter[converted_gene_name] if converted_gene_name in type_converter else 'NA',
                     'x':round(data_a[gene], 1)+0.1 if gene in data_a else 0.1,
                     'y':round(data_b[gene], 1)+0.1 if gene in data_b else 0.1})

    return jsonify(pair)


@app.route("/plot/<plot_name>")
def plotSelected(plot_name):
    template_list = getPlotTemplateList()
    binned_file_list = getFileList(['binned'])
    gene_file_list = getFileList(['genelist'])
    return render_template('plot_custom.html', selected_plot=plot_name, plot_templates=template_list, main_file_list=binned_file_list, gene_file_list=gene_file_list)


@app.route("/discover")
def discover():
    file_list = getFileList(['genelist'])
    genelist2data = prepGeneLists(path.split('/')[-1].replace('.json', '') for path in file_list)
    return render_template('explain.html', data=genelist2data)

@app.route("/venn")
def venn():
    return render_template('venn_dev.html')

def prepGeneLists(genelists):
    list2genes = {}
    id_sets = process_instance.prepIdSets()
    for genelist in genelists:
        list2genes[genelist] = list(process_instance.getGeneListInLocusName(genelist, id_sets))
    return list2genes


@app.route("/annotate")
def annotate():
    file_list = getFileList()
    return render_template('annotate.html', data=file_list)


@app.route("/getgenelist/<genelist>")
def getGeneList(genelist):
    id_sets = process_instance.prepIdSets()
    gene_set = process_instance.getGeneListInLocusName(genelist, id_sets)
    return jsonify(list(gene_set))


@app.route("/getlist")
def getlist():
    path = request.args.get('path')
    if os.path.isdir(path):
        file_type = request.args.get('file_type')
        if file_type in ['.fastq.gz', '.fa']:
            files = [file for file in os.listdir(path) if file.endswith(file_type)]
            #print(files)
            return jsonify(files)
    else:
        return jsonify([])


@app.route("/preprocess", methods=['GET', 'POST'])
def preprocess():
    if request.method == 'POST':
        print(request.form['path'])

    return render_template('preprocess.html')


if __name__ == "__main__":
    gunicorn_logger = logging.getLogger('gunicorn.error')
    app.logger.handlers = gunicorn_logger.handlers
    app.logger.setLevel(gunicorn_logger.level)
    app.run(host='0.0.0.0', debug=True)

####################################
####################################

def getData(file, folder):
    try:
        with open(f'{folder.replace("_", "/")}/{file}') as f:
            data = json.load(f)
            return data
    except:
        return []

def getFileList(include_list=[]):
    try:
        parents = [folder for folder in os.listdir('./data/') if os.path.isdir(f'./data/{folder}') and folder != 'original']
        parents = set(parents) & set(include_list) if include_list else parents
        parent2child = {parent:[folder for folder in os.listdir(f'./data/{parent}')] for parent in parents}
    except Exception as e:
        return f'<p>Error somewhere</p>{e}'
    return [f'{parent}/{child}' for parent in parent2child for child in sorted(parent2child[parent])]

def getOrigFolderList():
    try:
        parents = [folder for folder in os.listdir('./data/original')]
        parent2child = {parent:[folder for folder in os.listdir(f'./data/original/{parent}')] for parent in parents}
    except Exception as e:
        return f'<p>Error somewhere</p>{e}'
    return [f'original/{parent}/{child}' for parent in parent2child for child in parent2child[parent]]

def getGeneSets():
    try:
        gene_set_list = [file.replace('.json', '') for file in os.listdir('./data/genelist') if file.endswith('json')]
        return gene_set_list
    except Exception as e:
        #print(f'Error getting gene sets: {e}')
        return []

def organizeFileList(file_list):
    organized_list = []
    for parent in sorted(set([path.split('/')[0] for path in file_list])):
        organized_list.append((f'<h5 class="mt-3"> {parent} </h5>', ''))
        for file in [file for file in file_list if file.startswith(f'{parent}/')]:
            #organized_list.append((file.split('/')[-1], file))
            if file.startswith('original'):
                organized_list.append(('/'.join(file.split('/')[-2:]), file))
            else:
                organized_list.append((file.split('/')[-1].replace('.json', ''), file))
        organized_list.append(('', ''))
    return organized_list


def getPlotTemplateList():
    return [(re.sub(r'\.html$', '', file).replace('_', ' '), re.sub(r'\.html$', '', file)) for file in os.listdir('templates/plot_templates')]