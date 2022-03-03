import re, sys, getopt, argparse, os
VERBOSE = False

def process_vdj_fasta(fastaFilePath, keptTranscriptsFile, outFile, sampleLabel):
    # build list of transcripts to keep
    transcriptsFile = open(keptTranscriptsFile, 'r')
    transcriptList = []
    for line in transcriptsFile:
        transcriptList.append('>' + line.rstrip().split(',')[0])
    if VERBOSE:
        print('process_vdj_fasta, building transcript list:', transcriptList)
    transcriptsFile.close()

    fastaFile = open(fastaFilePath, 'r')
    keep_transcript = False
    for line in fastaFile:
        if '>' in line:
            splitLine = re.split(r'-\d', line.rstrip())
            transcript = splitLine[0]
            contig_number = splitLine[1]
            if VERBOSE:
                print('process_vdj_fasta, found in fasta:', transcript)
            if transcript in transcriptList:
                reconstructedLine = splitLine[0] + '-' + str(sampleLabel) + contig_number + '\n'
                outFile.write(reconstructedLine)
                keep_transcript = True
        elif keep_transcript:
            outFile.write(line)
            keep_transcript = False
    
    fastaFile.close()

def process_vdj_annotation(annotationFilePath, keptTranscriptsFile, outFile, sampleLabel, is_first_sample, metadataLabels = None):
    metadataDict = None
    if metadataLabels != None:
        metadataDict = dict()
        
    transcriptsFile = open(keptTranscriptsFile, 'r')
    transcriptList = []
    for line in transcriptsFile:
        if metadataDict == None:
            assert len(line.rstrip().split(',')) == 1
            transcriptList.append(line.rstrip())
        else:
            splitTranscriptLine = line.rstrip().split(',')
            transcript, metadataValueList = splitTranscriptLine[0], splitTranscriptLine[1:]
            transcriptList.append(transcript)
            assert len(metadataValueList) == len(metadataLabels)
            metadataDict[transcript] = metadataValueList
    if VERBOSE:
        print('process_vdj_annotation:', transcriptList)
    transcriptsFile.close()
    
    annotationFile = open(annotationFilePath, 'r')
    first_row = True
    contig_id_regex = re.compile(r'-\w+_contig_')
    for line in annotationFile:
        # Only write the header row to the outfile once, during the first sample file read
        if first_row:
            if is_first_sample:
                if metadataLabels != None:
                    line = ','.join([line.rstrip(), ','.join(metadataLabels)]) + '\n'
                outFile.write(line)
            first_row = False
        else:
            splitLine = line.rstrip().split(',', 4)
            
            transcript_and_sample = splitLine[0].split('-')
            transcript = transcript_and_sample[0]
            original_sample = transcript_and_sample[1]
            
            is_cell = splitLine[1]
            
            contig_id = splitLine[2]
            contig_id_split = contig_id_regex.split(contig_id)
            contig_id_transcript = contig_id_split[0]
            contig_id_contig_num = contig_id_split[1]

            high_confidence = splitLine[3]
            
            remainder = splitLine[4]

            if VERBOSE:
                print('process_vdj_annotation:', transcript)
            if transcript in transcriptList:
                if sampleLabel == '':
                    sampleLabel = original_sample
                reconstructedLine = transcript + '-' + sampleLabel + ',' + is_cell + ',' + contig_id_transcript + '-' + sampleLabel + '_contig_' +contig_id_contig_num + ',' + high_confidence + ',' + remainder
                if metadataDict != None:
                    metadataToAdd = ','.join(metadataDict[transcript])
                    if VERBOSE:
                        print('process_vdj_annotation:', metadataToAdd)
                    reconstructedLine = ','.join([reconstructedLine, metadataToAdd])
                if VERBOSE:
                    print('process_vdj_annotation:', reconstructedLine)
                outFile.write(reconstructedLine + '\n')

    annotationFile.close()

def process_vdj(inRefPaths, metadataList, inFastaPaths, inAnnoPaths, outFastaPath, outAnnotationPath, sampleLabels):
    numSamples = len(sampleLabels)
    if outAnnotationPath != None:
        outAnnotationFile = open(outAnnotationPath, 'w')
        for i in range(numSamples):
            process_vdj_annotation(inAnnoPaths[i], inRefPaths[i], outAnnotationFile, sampleLabels[i], i == 0, metadataList)
        outAnnotationFile.close()
    if outFastaPath != None:
        outFastaFile = open(outFastaPath, 'w')
        for i in range(numSamples):
            process_vdj_fasta(inFastaPaths[i], inRefPaths[i], outFastaFile, sampleLabels[i])
        outFastaFile.close()

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('config', help='CSV with the input paths and metadata label')
    parser.add_argument('-a', '--annotation_output', help='path to write annotation CSV file to', default='contigs_annotation.csv')
    parser.add_argument('-f', '--fasta_output', help='path to write contigs fasta file to', default='contigs.fasta')
    parser.add_argument("-v", "--verbose", help="increase logging verbosity",
                    action="store_true", default=False)
    args = parser.parse_args()
    
    compareAnnotation = False
    compareFasta = False
    
    metadataList = None
    
    inputAnnotationFiles = None
    inputFastaFiles = None
    
    outputFastaFile = None
    outputAnnotationFile = None
    
    sampleLabels = []
    inputReferenceFiles = []
    
    if os.path.exists(args.config):
        with open(args.config) as configFile:
            generalConfig = configFile.readline().rstrip().split(',')
            compareAnnotation = generalConfig[0] == 'True'
            compareFasta = generalConfig[1] == 'True'
            
            if not compareAnnotation and not compareFasta:
                print('ERROR: must aggregate annotations, fastas, or both')
                return
            else:
                if compareAnnotation:
                    inputAnnotationFiles = []
                if compareFasta:
                    inputFastaFiles = []
            
            if len(generalConfig) > 2:
                metadataList = generalConfig[2:]
            
            for line in configFile:
                splitLine = line.rstrip().split(',')
                numConfigs = len(splitLine)
                
                sampleLabels.append(splitLine[0])
                inputReferenceFiles.append(splitLine[1])
                if numConfigs == 3:
                    if compareAnnotation and not compareFasta:
                        inputAnnotationFiles.append(splitLine[2])
                    elif compareFasta and not compareAnnotation:
                        inputFastaFiles.append(splitLine[2])
                    else:
                        print('ERROR: three configuration values given, when four expected')
                        return
                elif numConfigs == 4:
                    if compareAnnotation and compareFasta:
                        inputAnnotationFiles.append(splitLine[2])
                        inputFastaFiles.append(splitLine[3])
                    else:
                        print('ERROR: four configuration values given, when three expected')
                        return
                else:
                    print('ERROR: incorrect number of configuration values given', numConfigs)
                    return
    else:
        print('ERROR: the configuration file could not be found')
        return
    if compareAnnotation:
        outputAnnotationFile = args.annotation_output
    if compareFasta:
        outputFastaFile = args.fasta_output 
    if args.verbose:
        global VERBOSE
        VERBOSE = True 

    process_vdj(inputReferenceFiles, metadataList, inputFastaFiles, inputAnnotationFiles, outputFastaFile, outputAnnotationFile, sampleLabels)

if __name__ == "__main__":
        main()