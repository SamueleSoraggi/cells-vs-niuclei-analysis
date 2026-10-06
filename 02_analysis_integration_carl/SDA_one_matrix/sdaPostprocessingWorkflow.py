from gwf import Workflow
import os, sys

gwf = Workflow()




def SDAanalysis(COMPONENT, notebook, outputMap, cores='1', memory='16g', walltime='00:10:00', account='Carl'):
    inputs=[notebook+'.ipynb']
    outputs=[f'{outputMap}/Genes_component_{COMPONENT}.log'] 
    options = {
        'cores': cores,
        'memory': memory,
        'walltime':walltime,
        'account': account
    }

    spec = f'''
    eval "$(conda shell.bash hook)"
    conda activate /faststorage/project/testis_singlecell/Workspaces/samuele/Environments/Carl
    mkdir -p {outputMap}
    papermill {notebook}.ipynb ./SDA_components_notebooks/SDA_C{COMPONENT}.ipynb --parameters comp {COMPONENT}
    echo "hello" > {outputMap}/Genes_component_{COMPONENT}.log
    '''
    return inputs, outputs, options, spec


##Run SDA plots for each component for data on one concatenated matrix. Output notebooks
##component by component in the chosen output folder
inputFile = ['./Part02_SDA_component_plots']
outputMap = ['./SDA_components_notebooks']
for filename, outmap in zip(inputFile,outputMap):
    for component in range(100):
        name = filename.split('/')[-1]
        gwf.target_from_template( f'sda_{name}_C{component}' , SDAanalysis(notebook=filename,
                                                outputMap=outmap,
                                                COMPONENT=component) )

