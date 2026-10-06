from gwf import Workflow
import os, sys

gwf = Workflow()

########
#run sda on a tensor. This can be produced by exporting each matrix (same nr cells and genes) and concatenating it. File must be tab separated without any
#header or row name. However keep track of cell names to see their properties in the analysis. Example: save each matrix using pandas in a space separated
#csv file with row names as cell names. Then concatenate them excluding the row names - so cells can be tracked back again.
def sda(inputFile, nInd, components=10, outFile=None, cores=8, memory='100g', walltime='48:00:00', account='Carl'):
    inFile = inputFile+'.csv'
    if outFile==None:
        outFile = inputFile
    inputs = [inFile]
    outputs =  [outFile + '/it2000/'+ i for i in ['A','S1','X1','free_energy'] ]
    options = {
        'cores': cores,
        'memory': memory,
        'walltime': walltime,
        'account': account
    }

    spec ='''
    mkdir -p {outFile}
    ./SDA/sda_static_linux --data {inFile} --out {outFile} \
    --N {nInd} \
    --num_comps {components} \
    --eigen_parallel true \
    --num_blocks {cores_10} \
    num_openmp_threads {cores} \
    --remove_zero_comps true \
    --set_seed 21090 41205
'''.format(inputFile=inputFile, inFile=inFile, outFile=outFile, nInd=nInd, cores_10=cores*10, cores=cores, components=components)

    return inputs, outputs, options, spec

#SDA with scTransform-ed and gaussian-normalized data. All data in a single matrix as they do in Jung et al.
#This should remove the extra matrix with tissue scores - but the principle is the same, just one
#score less to multiply with its relative sign. Maybe using only one matrix will better highlight the
#presence of batch components that isolate technical artifacts in the data. I found none in my analysys.
#Using scTransformed data might have removed such effects? Or you need more components to isolate it?
inputFile = ['./CELLS_scaled_data']
for filename in inputFile:
    name = filename.split('/')[-1]
    gwf.target_from_template( f'sda_{name}_50C' , sda(inputFile=filename,
                                                    outFile=f'{filename}_50comps',
                                                    nInd=2364, ########### WRITE HERE THE NUMBER OF CELLS OF YOUR DATASET
                                                    components=50,
                                                    cores=10) )

