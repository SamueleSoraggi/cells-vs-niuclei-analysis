from gwf import Workflow
import os, sys

gwf = Workflow()




def SDAanalysis(COMPONENT, notebook='MAST_TEMPLATE', outputMap='MAST', cores='1', memory='50g', walltime='1:00:00', account='primatescrna'):
    if COMPONENT==0:
        inputs=[notebook+'.ipynb']
    else:
        OLDCOMP=COMPONENT-1
        inputs=[f'{outputMap}/Genes_component_{OLDCOMP}.log']
    outputs=[f'{outputMap}/Genes_component_{COMPONENT}.log']
    options = {
        'cores': cores,
        'memory': memory,
        'walltime':walltime,
        'account': account
    }

    spec = f'''
    source activate SCdge
    mkdir -p {outputMap}
    COMPONENT={COMPONENT} runipy {notebook}.ipynb {outputMap}/Genes_component_{COMPONENT}.ipynb
    echo "hello" > {outputMap}/Genes_component_{COMPONENT}.log
    '''
    return inputs, outputs, options, spec


##Run SDA plots for each component for data on one concatenated matrix. Output notebooks
##component by component in the chosen output folder
inputFile = ['./Part02_SDA_component_plots']
outputMap = ['./SDA_component_plots']
for filename, outmap in zip(inputFile,outputMap):
    for component in range(15):
        name = filename.split('/')[-1]
        gwf.target_from_template( f'sda_{name}_C{component}' , SDAanalysis(notebook=filename,
                                                outputMap=outmap,
                                                COMPONENT=component) )

