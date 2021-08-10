from flask import Flask, render_template, request, flash
from process import Process


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