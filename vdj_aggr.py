import re, sys, getopt

usage_prompt = """
python vdj_aggr.py [-h] [-v] -i <input paths> -f <FASTA output> -a <annotation output>
	    <input paths>
                input paths file
        <fasta output>
                output FASTA file path
        <annotation output>
                output annotation CSV file path
options:
	-h, --help
         	prints out the above usage statement
    -v
            verbose, prints out additional debugging statements during runtime
"""
VERBOSE = False

def process_vdj_fasta(fastaFilePath, keptTranscriptsFile, outFile, sampleLabel):
    transcriptsFile = open(keptTranscriptsFile, 'r')
    fastaFile = open(fastaFilePath, 'r')

    transcriptList = []
    
    for line in transcriptsFile:
        transcriptList.append('>' + line.rstrip())
    if VERBOSE:
        print('process_vdj_fasta:', transcriptList)
    transcriptsFile.close()

    keep_transcript = False
    for line in fastaFile:
        if '>' in line:
            splitLine = re.split(r'-\d', line.rstrip())
            transcript = splitLine[0]
            contig_number = splitLine[1]
            if VERBOSE:
                print('process_vdj_fasta:', transcript)
            if transcript in transcriptList:
                reconstructedLine = splitLine[0] + '-' + sampleLabel + contig_number + '\n'
                outFile.write(reconstructedLine)
                keep_transcript = True
        elif keep_transcript:
            outFile.write(line)
            keep_transcript = False
    
    fastaFile.close()

def process_vdj_annotation(annotationFilePath, keptTranscriptsFile, outFile, sampleLabel, is_first):
    transcriptsFile = open(keptTranscriptsFile, 'r')
    annotationFile = open(annotationFilePath, 'r')

    transcriptList = []
    
    for line in transcriptsFile:
        transcriptList.append(line.rstrip())
    if VERBOSE:
        print('process_vdj_annotation:', transcriptList)
    transcriptsFile.close()
    first_row = is_first
    for line in annotationFile:
        if first_row:
                outFile.write(line)
                first_row = False
        else:
            splitLine = line.rstrip().split(',', 4)
            
            transcript = splitLine[0].split('-')[0]
            
            is_cell = splitLine[1]
            
            contig_id = splitLine[2]
            contig_id_split = re.split(r'-\w+_contig_', contig_id)
            print(contig_id_split)
            contig_id_transcript = contig_id_split[0]
            contig_id_contig_num = contig_id_split[1]

            high_confidence = splitLine[3]
            
            remainder = splitLine[4]

            if VERBOSE:
                print('process_vdj_annotation:', transcript)
            if transcript in transcriptList:
                reconstructedLine = transcript + '-' + sampleLabel + ',' + is_cell + ',' + contig_id_transcript + '-' + sampleLabel + '_contig_' +contig_id_contig_num + ',' + high_confidence + ',' + remainder + '\n'
                if VERBOSE:
                    print('process_vdj_annotation:', reconstructedLine)
                outFile.write(reconstructedLine)

    annotationFile.close()

def process_vdj(inRefPaths, inFastaPaths, inAnnoPaths, outFastaPath, outputAnnotationPath, sampleLabels):
    outFastaFile = open(outFastaPath, 'w')
    outAnnotationFile = open(outputAnnotationPath, 'w')

    numSamples = len(sampleLabels)
    for i in range(numSamples):
        process_vdj_fasta(inFastaPaths[i], inRefPaths[i], outFastaFile, sampleLabels[i])
        process_vdj_annotation(inAnnoPaths[i], inRefPaths[i], outAnnotationFile, sampleLabels[i], i == 0)
    outFastaFile.close()
    outAnnotationFile.close()

def main():
    
    inputFastaFiles = []
    inputAnnotationFiles = []
    inputReferenceFiles = []
    sampleLabels = []
    outputFastaFile = 'contigs.fasta'
    outputAnnotationFile = 'contigs_annotation.csv'
    inputsGiven = False

    try:
        opts, args = getopt.getopt(sys.argv[1:], "hi:f:a:v", ["help"])
    except getopt.GetoptError as err:
        print(usage_prompt)
        print(str(err))
        return
    for opt, arg in opts:
        if opt in ('-h', '--help'):
            print(usage_prompt)
            return
        elif opt == '-i':
            inputsGiven = True
            inputConfigPath = arg
            inputConfigFile = open(inputConfigPath)
            for line in inputConfigFile:
                splitLine = line.rstrip().split(',')
                sampleLabels.append(splitLine[0])
                inputReferenceFiles.append(splitLine[1])
                inputFastaFiles.append(splitLine[2])
                inputAnnotationFiles.append(splitLine[3])
            inputConfigFile.close()
        elif opt == '-f':
            outputFastaFile = arg
        elif opt == '-a':
            outputAnnotationFile = arg
        elif opt == '-v':
           global VERBOSE
           VERBOSE = True
        else:
            print(usage_prompt)
            return
    if inputsGiven == False:
        print(usage_prompt)
        print('ERROR: Not all required options were added')
        return
    process_vdj(inputReferenceFiles, inputFastaFiles, inputAnnotationFiles, outputFastaFile, outputAnnotationFile, sampleLabels)

if __name__ == "__main__":
        main()