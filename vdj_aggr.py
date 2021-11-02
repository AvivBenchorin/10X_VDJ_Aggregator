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

VERBOSE = True

def process_vdj_fasta(fastaFilePath, keptTranscriptsFile, outputPath, sampleNumber):
    transcriptsFile = open(keptTranscriptsFile, 'r')
    fastaFile = open(fastaFilePath, 'r')
    outFile = open(outputPath, 'a')

    transcriptList = []
    
    for line in transcriptsFile:
        transcriptList.append('>' + line)
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
                reconstructedLine = splitLine[0] + '-' + str(sampleNumber) + contig_number
                outFile.write(reconstructedLine)
                keep_transcript = True
        elif keep_transcript:
            outFile.write(line)
            keep_transcript = False
    
    outFile.close()
    fastaFile.close()

def process_vdj_annotation(annotationFilePath, keptTranscriptsFile, outputPath, sampleNumber):
    transcriptsFile = open(keptTranscriptsFile, 'r')
    annotationFile = open(annotationFilePath, 'r')
    outFile = open(outputPath, 'a')

    transcriptList = []
    
    for line in transcriptsFile:
        transcriptList.append('>' + line)
    if VERBOSE:
        print('process_vdj_annotation:', transcriptList)
    transcriptsFile.close()
    
    first_row = True
    for line in annotationFile:
        if first_row:
            first_row = False
            pass
        else:
            splitLine = line.rstrip().split(',', maxsplit=1)
            transcript = splitLine[0].split('-')[0]
            remainder = splitLine[1]
            if VERBOSE:
                print('process_vdj_annotation:', transcript)
            reconstructedLine = transcript + '-' + str(sampleNumber) + remainder
            outFile.write(reconstructedLine)

    outFile.close()
    annotationFile.close()

def process_vdj(inRefPaths, inFastaPaths, inAnnoPaths, outFastaPath, outputAnnotationPath):
    numSamples = len(inRefPaths)
    for i in range(numSamples):
        process_vdj_fasta(inFastaPaths[i], inRefPaths[i], outFastaPath, i)
        process_vdj_annotation(inAnnoPaths[i], inRefPaths[i], outputAnnotationPath, i)

def main():
    
    inputFastaFiles = []
    inputAnnotationFiles = []
    inputReferenceFiles = []
    outputFastaFile = 'contigs.fasta'
    outputAnnotationFile = 'contigs_annotation.csv'
    
    try:
        opts, args = getopt.getopt(sys.argv[1:], "hi:f:a:v", ["help"])
        print(opts, args)
    except getopt.GetoptError as err:
        print(usage_prompt)
        print(str(err))
        return
    for opt, arg in opts:
        if opt in ('-h', '--help'):
            print(usage_prompt)
            return
        elif opt == '-i':
            inputConfigPath = arg
            inputConfigFile = inputConfigPath.open()
            for line in inputConfigFile:
                splitLine = line.split('->')
                referencePath = splitLine[0]
                inputReferenceFiles.append(referencePath)
                paths = splitLine[1].split(',')

                inputFastaFiles.append(paths[0])
                inputAnnotationFiles.append(paths[1])
            inputConfigFile.close()
        elif opt == '-f':
            outputFastaFile = arg
        elif opt == '-a':
            outputAnnotationFile = arg
        elif opt == 'v':
            VERBOSE = True
        else:
            print(usage_prompt)
            return
    if len(opts) != 3:
        print(usage_prompt)
        print('ERROR: Not all required options were added')
    process_vdj(inputReferenceFiles, inputFastaFiles, inputAnnotationFiles, outputFastaFile, outputAnnotationFile)

if __name__ == "__main__":
        main()